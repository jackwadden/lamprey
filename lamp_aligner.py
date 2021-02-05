#!/usr/bin/python

import sys, os, subprocess, argparse, collections

import random as rng

# use local swalign?
sys.path.insert(0, '/home/dna/software/cswalign')
import swalign

import mappy
import pysam
import vcfpy
import numpy as np
import statsmodels.api as sm

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


#########################
class bcolors:
    HEADER = '\033[95m'
    BLUE = '\033[94m'
    CYAN = '\033[96m'
    GREEN = '\033[32;1m'
    BGREEN = '\033[32;7m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def print_blue(text):
    sys.stdout.flush()
    sys.stdout.write(bcolors.BLUE + text + bcolors.ENDC)

def print_cyan(text):
    sys.stdout.flush()
    sys.stdout.write(bcolors.CYAN + text + bcolors.ENDC)

def print_red(text):
    sys.stdout.flush()
    sys.stdout.write(bcolors.FAIL + text + bcolors.ENDC)

def print_green(text):
    sys.stdout.flush()
    sys.stdout.write(bcolors.GREEN + text + bcolors.ENDC + "\n")
    
#########################
class PrimerSet:
    '''
    Stores bases in a primer set.
    '''
    
    def __init__(self, primer_set_filename):

        self.primer_dict = dict()
        self.fwd_primer_order = ['F3', 'F2', 'F1', 'T', 'B1c', 'B2c', 'B3c']
        self.rev_primer_order = ['B3', 'B2', 'B1', 'Tc', 'F1c', 'F2c', 'F3c']
        
        # initialize primer set given file
        with open(primer_set_filename) as fp:
            for line in fp:

                # read lines in format <primer_name> <primer_seq>
                line = line.strip()
                if not line: continue
                if line.startswith("#"): continue
                fields = line.split(' ')
                primer_name = fields[0]

                # parse fwd_primer_order line
                if primer_name == 'fwd_primer_order':

                    # clear default
                    self.fwd_primer_order.clear()

                    for field in fields[1:]:
                        self.fwd_primer_order.append(field)

                    continue

                # parse rev_primer_order line
                if primer_name == 'rev_primer_order':

                    # clear default
                    self.rev_primer_order.clear()
                    
                    for field in fields[1:]:
                        self.rev_primer_order.append(field)

                    continue

                # parse normal primer
                primer_string = fields[1]

                # also store revcomp version
                if primer_name.endswith('c'):
                    self.primer_dict[primer_name[:-1]] = rc(primer_string)
                else:
                    self.primer_dict[primer_name] = primer_string
                        
#########################
class Alignment:
    '''
    Stores alignment info.
    '''

    def __init__(self, start=0, end=0, identity=0.0, primer_name=''):
        self.start = start
        self.end = end
        self.identity = identity
        self.primer_name = primer_name

    def __repr__(self):
        return "{} at ({}, {}) identity: {}".format(
				self.primer_name, self.start, self.end, self.identity)
        
    def __str__(self):
        return "{} at ({}, {}) identity: {}".format(
				self.primer_name, self.start, self.end, self.identity)

    def __lt__(self, other):
        return self.start < other.start
    
    def __eq__(self, other):
        self.start == other.start
    
#######
def alignPrimer(aligner, seq, primer, primer_name):
    '''
    Aligns a primer to a DNA sequence.
    '''

    alignment_tmp = aligner.align(str(seq), primer)
    alignment = Alignment(alignment_tmp.r_pos,
                          alignment_tmp.r_end,
                          alignment_tmp.identity,
                          primer_name)

    return alignment

#######
def findPrimerAlignments(aligner, seq, primer, primer_name, identity_threshold):
    ''' 
    Greedy approach which finds the best primer alignment for a long sequence, 
    then uses left/right recursion to find all other possible (non-overlapping) 
    alignments. 
    '''

    # find optimal alignment
    alignment = alignPrimer(aligner, seq, primer, primer_name)
    alignment_len = alignment.end - alignment.start
    
    # TODO: improve heuristic (50% threshold -> 25% match passes worst case)
    alignment_list = list()
    if alignment.identity >= identity_threshold and \
			alignment_len >= len(primer) * identity_threshold:
        
        # add the optimal alignment to list of valid alignments
        alignment_list.append(alignment)
        
        # split the sequence into left/right sections based on the alignment
        left_seq = seq[:alignment.start]
        right_seq = seq[alignment.end:]
        
        # recurse left
        if len(left_seq) >= len(primer):
            left_alignments = findPrimerAlignments(
                    aligner, left_seq, primer, primer_name, identity_threshold) 
            alignment_list.extend(left_alignments)
        
        # recurse right
        if len(right_seq) >= len(primer):
            right_alignments = findPrimerAlignments(
                    aligner, right_seq, primer, primer_name, identity_threshold)

            # adjust right alignment positioning (since we cropped out seq start)
            for right_alignment in right_alignments:
                right_alignment.start += alignment.end
                right_alignment.end += alignment.end

            # add all right alignments
            alignment_list.extend(right_alignments)

    return alignment_list

#########################
# from stack overflow: https://stackoverflow.com/questions/25188968/reverse-complement-of-dna-strand-using-python
def rc(seq):
    ''' 
    Reverse-complement DNA sequence.
    '''
    
    # replace bases with complement, defaulting to same identifier if not found
    bases = list(seq)
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    return bases

############################
def pruneRapidAdapters(aligner, seq):
    '''
    Fancy printing of suspected rapid adapters still in sequence.
    '''

    # TODO: make adapter seq and thresholds variable
    adapter_ytop='GGCGTCTGCTTGGGTGTTTAACCTTTTTTTTTTAATGTACTTCGTTCAGTTACGTATTGCT'
    adapter_ybot='GCAATACGTAACTGAACGAAGT'

    # check ytop
    alignment = aligner.align(adapter_ytop, read_string)
    print(alignment.identity)
    print(alignment.q_pos)
    print(alignment.q_end)
    if alignment.q_end < 100 and alignment.identity > 0.55:
        print("Found possible adapter top of fwd:")
        alignment.dump()

    # check rc ytop
    alignment = aligner.align(adapter_ytop, rc(read_string))
    if alignment.q_end < 100 and alignment.identity > 0.55:
        print("Found possible adapter top on rc:")
        alignment.dump()

    # check ybot
    alignment = aligner.align(adapter_ybot, read_string)
    if alignment.q_pos > len(read_string) - (len(adapter_ybot) * 2) and alignment.identity > 0.55:
        print("Found possible adapter bot on fwd:")
        alignment.dump()

    # check rc ybot
    alignment = sw.align(adapter_ybot, rc(read_string))
    if alignment.q_pos > len(read_string) - (len(adapter_ybot) * 2) and alignment.identity > 0.55:
        print("Found possible adapter bot on rc:")
        alignment.dump()

##################
def findAllPrimerAlignments(aligner, seq, primers, identity_threshold, args):

    alignment_list = list()

    alignment_list.extend(findPrimerAlignments(aligner, seq, primers.primer_dict["T"], "T", identity_threshold))
    alignment_list.extend(findPrimerAlignments(aligner, seq, rc(primers.primer_dict["T"]), "Tc", identity_threshold))

    # bail if we didn't find a target in this read
    # this is a lazy shortcut optimization
    #if len(alignment_list) == 0:
    #    return(sorted(alignment_list))

    if args.high_confidence and len(alignment_list) < 2:
        return(sorted(alignment_list), 0.0)
    
    for primer_name, primer_seq in primers.primer_dict.items():
        if primer_name == 'T':
            continue
        else:
            # fwd
            alignment_list.extend(findPrimerAlignments(aligner, seq, primer_seq, primer_name, identity_threshold))
            # rc
            alignment_list.extend(findPrimerAlignments(aligner, seq, rc(primer_seq), primer_name + "c", identity_threshold))


    alignment_list = sorted(alignment_list)
            
    # get alignment coverage
    seq_length = len(seq)
    primer_coverage = 0
    for alignment in alignment_list:
        primer_coverage = primer_coverage + (alignment.end - alignment.start)
        #print("Primer coverage: " + str(primer_coverage))

    alignment_coverage = float(primer_coverage)/float(seq_length)
     
    return(alignment_list, alignment_coverage)
    
###################
def clearColor():
    sys.stdout.write(bcolors.ENDC)
    
def setPrimerColor(primer_string):

    color = ''
    
    if primer_string == 'T':
        color = bcolors.BGREEN
    elif primer_string == 'Tc':
        color = bcolors.GREEN
    elif primer_string == 'F2':
        color = '\u001b[38;5;4;7m'
    elif primer_string == 'F2c':
        color = '\u001b[38;5;4m'
    elif primer_string == 'F1':
        color = '\u001b[38;5;75;7m'
    elif primer_string == 'F1c':
        color = '\u001b[38;5;75m'
    elif primer_string == 'B2':
        color = '\u001b[38;5;1;7m'
    elif primer_string == 'B2c':
        color = '\u001b[38;5;1m'
    elif primer_string == 'B1':
        color = '\u001b[38;5;204;7m'
    elif primer_string == 'B1c':
        color = '\u001b[38;5;204m'
    sys.stdout.write(color)

def printPrimerAlignments(seq, alignments):
    '''
    Fancy printing of all primers aligned to sequence.
    '''

    print("", flush=True)
    print(seq)
    print(alignments)

    wrap = 100
    wrap_counter = 0
    
    found_primer = False
    primer_string = ''
    primer_string_counter = 0

    primer_queue = []
    
    for pos in range(0, len(seq)):
        sys.stdout.write(seq[pos])

        # print the alignment if we're at the wrap factor OR end of the string
        if (pos > 0 and (pos+1) % wrap == 0) or pos == len(seq)-1:
            sys.stdout.write('\n')

            # print corresponding alignment
            for aln_pos in range(wrap_counter * wrap, pos+1):

                found_start = False
                found_end = False

                    
                # do we need to emit a spec char?
                for alignment in alignments:
                    if alignment.start == aln_pos:                    
                        found_start = True
                        primer_string = alignment.primer_name
                        primer_string_counter = 0
                        primer_queue.append(alignment)
                    elif alignment.end == aln_pos:
                        found_end = True
                        primer_queue.pop()

                # reset color
                if found_primer:
                    setPrimerColor(primer_string)
                else:
                    sys.stdout.write(bcolors.ENDC)
           
                        
                # emit special char
                if found_start and found_end:
                    sys.stdout.write('X')
                    found_primer = True
                elif found_start:
                    # turn on color
                    setPrimerColor(primer_string)
                    sys.stdout.write('$')
                    found_primer = True
                elif found_end:
                    sys.stdout.write('^')
                    if len(primer_queue) == 0:
                        found_primer = False
                else:

                    # if no special char, print either a primer string, -, 
                    # or space if not currently in a primer
                    if found_primer:
                        if primer_string_counter < len(primer_string):
                            sys.stdout.write(primer_string[primer_string_counter])
                            primer_string_counter = primer_string_counter + 1
                        else:
                            sys.stdout.write('-')
                    else:
                        sys.stdout.write(' ')
                        
            # end alignment line
            sys.stdout.write(bcolors.ENDC)
            sys.stdout.write('\n\n')

            # increment wrap counter
            wrap_counter = wrap_counter + 1
            
    sys.stdout.write('\n')

################################################
def removeOverlappingPrimerAlignments(alignments, allowed_overlap):

    print("* Removing overlapping primer alignments...")
    
    new_alignments = list()
    overlapping_alignments = list()

    overlapping_alignments = 0

    last_start_pos = 0
    last_end_pos = 0
    last_alignment = alignments[0]
    
    for idx,alignment in enumerate(alignments):

        # skip first alignment
        if idx == 0:
            continue

        # if this alignment overlaps with the last alignment, pick a winner
        if alignment.start < last_alignment.end - allowed_overlap:

            print("  - Found overlapping primer alignment: {}".format(str(alignment.start) + " : " + str(alignment.end)), flush=True)

            if alignment.identity > last_alignment.identity:
                # erase last alignment without adding it
                last_alignment = alignment
            else:
                # skip this alignment altogether
                last_alignment = last_alignment
        else:
            # add last alignment... it's cleared
            new_alignments.append(last_alignment)
            last_alignment = alignment

    # append last alignment
    new_alignments.append(last_alignment)

    return new_alignments

    
###################
def extractAmpliconAroundTarget(primer_set, alignments, target):
    '''
    Target sequence was found in lamplicon. Extend this sequence as far as 
    possible both left and right (including primers) to increase mappability.
    This is done using expected next primer, allowing one mismatch.
    '''

    amplicon_start = 0
    amplicon_end = 0
    
    # is target forward or reverse?
    fwd_strand = target.primer_name == "T"

    # find target alignment index
    target_index = alignments.index(target)

    #####
    # look to the right 
    #####
    allowed_mismatches = 0
    primer_counter = primer_set.fwd_primer_order.index('T') if fwd_strand else primer_set.rev_primer_order.index('Tc')
    primer_counter = primer_counter + 1

    all_matched = True
    for i in range(target_index + 1, len(alignments)):

        # primer name
        primer_name = alignments[i].primer_name

        # expected primer name
        if fwd_strand:
            expected_primer_name = primer_set.fwd_primer_order[primer_counter]
        else:
            expected_primer_name = primer_set.rev_primer_order[primer_counter]

        # on a mismatch, end the search
        if primer_name != expected_primer_name:

            if allowed_mismatches > 0:
                allowed_mismatches = allowed_mismatches - 1
            else:
                # we've found our end, point, so make our end the last match end
                amplicon_end = alignments[i-1].end
                all_matched = False
                break
        else:
            primer_counter = primer_counter + 1

            # if we matched all primers, we have to quit no matter what
            if primer_counter == len(primer_set.fwd_primer_order):
                break

    # if we matched all in the primer sequence, extend to the end of the entire read
    if all_matched:
        amplicon_end = -1

    #####
    # look to the left
    #####
    allowed_mismatches = 0
    primer_counter = primer_set.fwd_primer_order.index('T') if fwd_strand else primer_set.rev_primer_order.index('Tc')
    primer_counter = primer_counter - 1

    all_matched = True
    for i in range(target_index - 1, 0, -1):

        # primer name
        primer_name = alignments[i].primer_name

        # expected primer name
        if fwd_strand:
            expected_primer_name = primer_set.fwd_primer_order[primer_counter]
        else:
            expected_primer_name = primer_set.rev_primer_order[primer_counter]

        # on a mismatch, end the search
        if primer_name != expected_primer_name:
            if allowed_mismatches > 0:
                allowed_mismatches = allowed_mismatches - 1
                
            else:
                # we've found our end, point, so make our end the last match end
                amplicon_start = alignments[i+1].start
                all_matched = False
                break
        else:
            primer_counter = primer_counter - 1

    if all_matched:
        amplicon_start = 0

    return amplicon_start, amplicon_end


################################################
def getI(read):
    return read.get_cigar_stats()[0][1]

def getD(read):
    return read.get_cigar_stats()[0][2]

def getEq(read):
    return read.get_cigar_stats()[0][7]

def getX(read):
    return read.get_cigar_stats()[0][8]

def getAccuracy(read, alt_pos=None, alt_base=None):
    '''
    Calculate accuracy for a single read, ignoring suspected alternate bases.
    '''

    # get CIGAR stats
    Eq = getEq(read)
    X = getX(read)
    I = getI(read)
    D = getD(read)

    # return easy answer if we don't care about alternate reads
    if alt_pos is None or alt_base is None: 
        return (float(Eq) / float(Eq + I + X + D))

    # adjust Eq and X if this is a suspected alt read
    # find the alignment pair for the ref position
    # query the query pos that corresponds, and check if the base is the alt base
    for pair in read.get_aligned_pairs():
        if pair[0] == None or pair[1] == None:
            continue
        if pair[1] == (alt_pos - 1):
            if read.query_alignment_sequence[pair[0]-read.query_alignment_start] == alt_base:
                Eq = Eq + 1
                X = X - 1
                
    return (float(Eq) / float(Eq + I + X + D))

#### parses CIGAR string and calc total errors.
# ignores leading/trailing errors
# can accept one alt base possibility/position pair to account for tumor reads
# if an alt base is provided, a consensus vote for the alt at that position will
# not be counted as an error
def calcErrorRate(ref_fn, sam_fn, alt_pos=None, alt_base=None):
    '''
    Calculate overall error rate in SAM file.
    '''

    # init
    samfile = pysam.AlignmentFile(sam_fn, 'r')
    cumulative_err_rate = 0.0
    mapped_read_count = 0

    # accumulate errors
    for read in samfile:
        if not read.is_unmapped:
            mapped_read_count += 1
            cumulative_err_rate += getAccuracy(read, alt_pos, alt_base)

    # if no mapped reads, cumulative accuracy is 0
    if mapped_read_count == 0: return 0

    # return average error
    avg_err_rate = cumulative_err_rate / mapped_read_count
    return avg_err_rate


################################################
def parsePileupCodeString(code_string, variant):

    variant_count = 0
    ref_count = 0
    total_count = 0
    skip_counter = 0
    set_skip = False

    for c in code_string:

        #print(c)                                                                                                                                                            

        if(skip_counter > 0):
            skip_counter = skip_counter - 1
            continue

        if(set_skip):
            skip_counter = int(c)
            set_skip = False
            continue

        total_count = total_count + 1

        if c == ',' or c == '.':
            ref_count = ref_count + 1

        elif c == variant or c == variant.lower():
            variant_count = variant_count + 1

        elif c == '-' or c == '+':
            set_skip = True

        elif c == '^':
            skip_counter = 1

    return ref_count, variant_count
                          

################################################
def extractPileupInfo(pileup_fn, contig, target_locus, limit):
    # parse lines in the pileup, extract the depth at the locus

    ref = ''
    depth = 0
    code_string = ''
    covers_roi = False
    within_roi = False
    with open(pileup_fn) as fp:
        for line in fp:
            fields = line.split('\t')
            if fields[0] == contig:

                pileup_locus = int(fields[1])
                
                if pileup_locus == target_locus:
                    ref = fields[2]
                    depth = int(fields[3])
                    code_string = fields[4]
                    covers_roi = True
                elif pileup_locus < (target_locus + limit) and pileup_locus > (target_locus - limit):
                    within_roi = True
                
    return ref, depth, code_string, covers_roi, within_roi

################################################
def generatePileup(bam_to_pileup_sh, ref_fn, bam_fn, out_fn):

    subprocess.run([bam_to_pileup_sh, ref_fn, bam_fn], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

################################################
def proportionConfidenceInterval(mut, wt, conf):


    conf_int = sm.stats.proportion_confint(count=mut,
                                           nobs=mut+wt,
                                           alpha=(1 - conf))

    samp_mean = float(mut)/float(mut + wt) * 100.0
    if conf_int[0] == 0.0:
        low = 0.0
    else:
        low = samp_mean - conf_int[0] * 100.0
    high = conf_int[1] * 100.0 - samp_mean

    #print("conf:", conf, "%.2f" % low, " -- ", "%.2f" % samp_mean, " -- ", "%.2f" % high, ")")
    print("conf:", conf, "%.4f" % (conf_int[0]*100.0), " -- ", "%.2f" % samp_mean, " -- ", "%.4f" % (conf_int[1]*100.0), ")")
    
################################################
def parseVCF(vcf_fn):

    reader = vcfpy.Reader.from_path(vcf_fn)

    # extract variant roi
    for record in reader:
        return record

    return
    
################################################
def is_hit_near_target(alignment, vcf_record, limit):
    contig = alignment.ctg
    ref_start = alignment.r_st
    ref_end = alignment.r_en

    target_contig = vcf_record.CHROM
    target_locus = vcf_record.POS

    if target_contig == contig:
        # fully contained to be considered a fragment
        if target_locus + limit > ref_end and target_locus - limit < ref_start:
            return True

    return False


################################################
def processLamplicon(sw, process_candidates_path, generate_consensus_path, minimap2, lamp_idx, lamplicon, primers, ref_vcf_record, target_vcf_record, args):


    # fragment limit
    # the distance between the locus and th start or end of the alignment to be considered a fragment
    fragment_limit = 1000

    # track how many votes each has
    mut = 0
    wt = 0
    
    alignments, alignment_coverage = findAllPrimerAlignments(sw, lamplicon.seq, primers, args.threshold, args)
    if len(alignments) > 0:
        # remove alignments that overlap by more than 4bp
        alignments = removeOverlappingPrimerAlignments(alignments, 4)

    # optionally print alignments of all primers
    if args.debug_print:
        printPrimerAlignments(lamplicon.seq, alignments)

    # does this lamplicon possibly cover the target?
    num_targets = 0
    num_ont_seqs = 0
    for alignment in alignments:
        if alignment.primer_name in ['T', 'Tc']:
            num_targets = num_targets + 1

        # track high quality ont sequences
        if alignment.primer_name.startswith('ont') and alignment.identity >= args.threshold:
            num_ont_seqs = num_ont_seqs + 1

    # if high confidence mode, only proceed if there are 2+ targets
    if args.high_confidence and num_targets < 2:
        return 0, 'unknown', 0, 0

    # low alignment coverage indicates this may be a genomic read
    if alignment_coverage < 0.15:
        print("* Low alignment coverage. Suspected genomic read. Aligning to human ref...")
        for hit in minimap2.map(lamplicon.seq):
            print(hit)
            
            # skip non-primary alignments
            if not hit.is_primary:
                continue
            
            # if the hit is within fragment_limit of the target, mark it as a fragment
            if not is_hit_near_target(hit, ref_vcf_record, fragment_limit):
                print("  - Alignment not near reference target. Assuming background genomic read.")
                classification = 'background'                
                return num_targets, classification, mut, wt

        print("  - Could not confirm genomic read, diagnosing further...")
                
    # if there were no targets identified over the threshold
    # classify the read
    if num_targets == 0:        
        
        print("* Didn't find any targets.... diagnosing...")

        # this is maybe a background genomic read or a lamplicon "fragment"
        # fragments occur due to fragmentation-based library prep or shearing
        # try to align the read to see if it's a background genomic read
        # if we have any hits, classify as background
        for hit in minimap2.map(lamplicon.seq):

            print(hit)
            
            # skip non-primary alignments
            if not hit.is_primary:
                continue
            
            # if the hit is within fragment_limit of the target, mark it as a fragment
            if is_hit_near_target(hit, ref_vcf_record, fragment_limit):
                classification = 'fragment'
                return num_targets, classification, mut, wt
            
            # else, just mark it as background
            classification = 'background'
            return num_targets, classification, mut, wt

        
        # classify read as "ont" if it contains > 1 ONT sequence and fewer than 2 other suspected primer seqs and greater than 50% coverage
        if num_ont_seqs > 0 and (len(alignments) - 2) <= num_ont_seqs :

            # get coverage
            if alignment_coverage >= 0.15:
                classification = 'ont'
            else:
                # assume unknown
                classification = 'unknown'

                # try to align the read to see if it's a background genomic read
                for hit in minimap2.map(lamplicon.seq):
                    classification = 'background'

            return num_targets, classification, mut, wt
                    
            
        # classify the read as spurious (bad/erroneous primer amplification) if it has at least 3 primer sequences, didn't align at the locus, and no target
        if len(alignments) >= 3:
            classification = 'spurious'
            return num_targets, classification, mut, wt

        
        # if there aren't any ont sequences or target/lamp sequences, then let's try to align this read to the human genome
        # write separated targets to new FASTA file
        #fastq_fn = "{}/tmp.fastq".format(args.output_dir)
        #with open(fastq_fn, "w") as output_handle:
        #    SeqIO.write(lamplicon, output_handle, "fastq")
        
        # map each target sequence; ignore secondary mappings 
        #out_sam_tmp = "{}/tmp.sam".format(args.output_dir)
        #subprocess.run(["minimap2","-a","-xmap-ont","--eqx","-t 1",
        #                "-n 1", "-N 5", "--secondary=no", "-m 0", "-o{}".
        #                format(out_sam_tmp),args.human_ref_fn,fastq_fn])

        
        #               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        # parse sam file and see if there was a mapping to the human genome

        #samfile = pysam.AlignmentFile(out_sam_tmp, 'r')
        #for read in samfile:
        #    if not read.is_unmapped:
        #        classification = 'off_target'

        # suspiciously short sequences might just be a fragment that can't align with enough identity
        if len(lamplicon.seq) < 60:
            classification = 'short'
        else:
            # final default to unknown
            classification = 'unknown'

        # Return the final classification
        return num_targets, classification, mut, wt

    ####################################
    # if we found a target....

    # extend each target into a single amplicon and emit as fasta
    target_idx = 0
    targets = []
    alignment_intervals = list()
    for alignment in alignments:

        if alignment.primer_name == 'T' or alignment.primer_name == 'Tc':
        
            # extend sequence around target for better mapping
            seq_start, seq_end = extractAmpliconAroundTarget(primers, alignments, alignment)
        
            # if we get a -1 for seq end, we never matched to the right, 
            # so just grab the whole read
            if seq_end == -1:
                seq_end = len(lamplicon.seq) - 1

            
            alignment_intervals.append(str(seq_start) + " : " + str(seq_end))
            
            # generate new read from lamplicon.seq substring
            record = SeqRecord(
                lamplicon.seq[seq_start:seq_end],
                "{}_{}".format(lamp_idx, target_idx),
                "{}_{}".format(lamp_idx, target_idx),
                "lamplicon {} target {}".format(lamp_idx, target_idx)
            )

            targets.append(record)
            target_idx += 1

    # Print intervals
    print("* Found {} candidate sub-reads:".format(len(alignment_intervals)))
    for interval in alignment_intervals:
        print("  " + interval)
    print("", flush=True)
          
    # write separated targets to new FASTA file
    fasta_fn = "{}/{}.fasta".format(args.output_dir, lamp_idx)
    with open(fasta_fn, "w") as output_handle:
        SeqIO.write(targets, output_handle, "fasta")


    ##################################################################
    #  STAGE 2: Map each sub-read against a smaller target reference #
    ##################################################################
    # map each sub read to the target sequence; ignore secondary mappings 
    print("* Aligning candidate sub-reads...")
    out_sam_all = "{}/{}_all.sam".format(args.output_dir, lamp_idx)
    out_bam = "{}/{}.bam".format(args.output_dir, lamp_idx)
    subprocess.run([process_candidates_path, fasta_fn, out_sam_all, args.target_ref_fn], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    # output pileup generated by process_candidates.sh
    pileup_fn = "{}/{}.pileup".format(args.output_dir, lamp_idx)
        
    # extract pileup info
    limit = 1000
    ref, depth, code_string, covers_roi, within_roi = extractPileupInfo(pileup_fn, target_vcf_record.CHROM, target_vcf_record.POS, limit)

    for alt in target_vcf_record.ALT:
        variant = alt.value.upper()

    print("  - Found {} reads that cover the target".format(depth))
          
    print("  - Pileup result: ",ref, depth, code_string, covers_roi)

    # in high confidence mode, ignore pileups that have less than 2 entries
    if args.high_confidence and depth < 2:
        return 0, 'unknown', 0, 0
    
    # if we don't cover the roi, then this was probably a fragment, off target, or spurious
    if not covers_roi:
        if within_roi:
            classification = 'fragment'
        else:
            classification = 'spurious'
        return num_targets, classification, mut, wt
        
    # parse code string and set mut/wt counts
    wt, mut = parsePileupCodeString(code_string, variant)

    print("  - Found {} mut and {} wt calls".format(mut,wt))
            
    # num targets is now the actual depth from the pileup
    num_targets = depth

    ## TODO ## re-implement consensus polishing here. External call to Medaka?

    # if we don't have any targets after splitting, and aligning
    if num_targets == 0:

        print("WHY!?!?!")
        exit(1)

        if len(alignments) >= 2:
            classification = 'spurious'
        elif len(lamplicon.seq) < 60:
            classification = 'short'
        else:
            classification = 'unknown'

    # if we do cover the roi and have at least one target, we are  target read
    else:
        classification = 'target'

    return num_targets, classification, mut, wt
        
        
def main(args):

    # get file directory for relative script paths
    dirname = os.path.dirname(__file__)
    process_candidates_path = os.path.join(dirname, 'scripts/process_candidates.sh')
    generate_consensus_path = os.path.join(dirname, 'scripts/generate_consensus.sh')
    os.makedirs(args.output_dir, exist_ok=True)

    # parse all lamplicons from FASTQ file
    print("> parsing lamplicons")
    #lamplicon_seqs = [l.seq for l in SeqIO.parse(args.lamplicons_fn, "fastq")]
    lamplicons = SeqIO.parse(args.lamplicons_fn, "fastq")
    
    # parse primers from config file
    print("> parsing primers")
    primers = PrimerSet(args.primer_set_fn)

    # parse ref VCF file if there is one
    if args.ref_vcf_fn != '':
        ref_vcf_record = parseVCF(args.ref_vcf_fn)

    # parse VCF file if there is one
    if args.target_vcf_fn != '':
        target_vcf_record = parseVCF(args.target_vcf_fn)
    
    # initialize aligner
    print("> initializing aligners")
    match = 2
    mismatch = -1
    scoring = swalign.NucleotideScoringMatrix(match, mismatch)
    sw = swalign.LocalAlignment(scoring)  # you can also choose gap penalties, etc...

    # set up a mappy instance to align reads
    minimap2 = mappy.Aligner(args.human_ref_fn)  # load or build index
    if not minimap2: raise Exception("ERROR: failed to load/build index")
    
    # iterate over all lamplicon sequences
    print("> splitting and mapping lamplicons")
    target_records = list()
    fragment_records = list()
    spurious_records = list()
    ont_records = list()
    background_records = list()
    short_records = list()
    unknown_records = list()

    # histogram data structures
    target_sequence_map = dict()
    
    # process each record
    class_counters = dict()
    class_counters['target'] = 0
    class_counters['fragment'] = 0
    class_counters['spurious'] = 0
    class_counters['ont'] = 0
    class_counters['background'] = 0
    class_counters['short'] = 0
    class_counters['unknown'] = 0


    # seed PRNG for VAF calcs
    simulate_vaf = False
    if simulate_vaf:
        rng.seed(7)
    
    mut_total = 0
    wt_total = 0

    mut_reads = 0
    wt_reads = 0

    target_lamplicon_counter = 0

    for lamp_idx,lamplicon in enumerate(lamplicons):
        
        # print status, stop early if we've found enough
        if args.max_lamplicons > 0 and lamp_idx >= args.max_lamplicons: break
        
        print_green("Processing Lamplicon {} \r".format(lamp_idx + 1))
        #print("Processing Lamplicon {} \r".format(lamp_idx + 1))
        sys.stdout.flush()
        
        num_targets,classification,mut,wt = processLamplicon(sw,
                                                             process_candidates_path,
                                                             generate_consensus_path,
                                                             minimap2,
                                                             lamp_idx,
                                                             lamplicon,
                                                             primers,
                                                             ref_vcf_record,
                                                             target_vcf_record,
                                                             args);


        # here, we can hijack the classification and inject our own VAF
        if simulate_vaf:
            target_vaf = 0.1
            if classification == "target":
                # roll the dice based on the VAF
                # adjust mut/wt calls depending on the result of the dice
                if rng.random() < target_vaf:
                    mut = mut + wt
                    wt = 0
                else:
                    # IMPORTANT: don't do anything to preserve error in wt dataset
                    mut = mut
                    wt = wt
            
        if classification == "target":

            target_lamplicon_counter = target_lamplicon_counter + 1

            # disaggregated totals (each target gets a vote)
            mut_total = mut_total + mut
            wt_total = wt_total + wt

            # aggregated totals (each lamplicon gets a vote)
            if mut > 0 and mut > wt:
                mut_reads = mut_reads + 1
            elif wt > 0 and wt >= mut:
                wt_reads = wt_reads + 1

        print("  - Classified as: " + classification)
        class_counters[classification] = class_counters[classification] + 1
        sys.stdout.write("  ")
        print(class_counters)

        
        if (mut_total + wt_total) > 0:
            print("  - Disagg VAF: {}/{} {}".format(mut_total, mut_total + wt_total, float(mut_total)/(float(mut_total)+float(wt_total))))
            print("  - Aggreg VAF: {}/{} {}".format(mut_reads, mut_reads + wt_reads, float(mut_reads)/(float(mut_reads)+float(wt_reads))))


            print("* VAF Confidence Intervals: ")
            sys.stdout.write("  ")
            proportionConfidenceInterval(mut_reads, wt_reads, 0.95)
            sys.stdout.write("  ")
            proportionConfidenceInterval(mut_reads, wt_reads, 0.99)
            sys.stdout.write("  ")
            proportionConfidenceInterval(mut_reads, wt_reads, 0.999)
            sys.stdout.write("  ")
            proportionConfidenceInterval(mut_reads, wt_reads, 0.9999)
            sys.stdout.flush()

        # Add record to list
        if classification == 'target':
            target_records.append(lamplicon)
        elif classification == 'fragment':
            fragment_records.append(lamplicon)
        elif classification == 'spurious':
            spurious_records.append(lamplicon)
        elif classification == 'ont':
            ont_records.append(lamplicon)
        elif classification == 'background':
            background_records.append(lamplicon)
        elif classification == 'short':
            short_records.append(lamplicon)
        elif classification == 'unknown':
            unknown_records.append(lamplicon)


        # keep track of each sequence based on how many targets they have
        if classification == "target":
            if num_targets not in target_sequence_map.keys():
                target_sequence_map[num_targets] = list()

            target_sequence_map[num_targets].append(lamplicon)
        
        # sort the target count histogram
        target_sequence_map = collections.OrderedDict(sorted(target_sequence_map.items()))
        print("")
        
    print("")

    
    if(args.save_target_reads):
        # write out on-target/off-target splits
        with open("{}/target_records.fastq".format(args.output_dir), "w") as output_handle:
            SeqIO.write(target_records, output_handle, "fastq")

        with open("{}/ont_records.fastq".format(args.output_dir), "w") as output_handle:
            SeqIO.write(ont_records, output_handle, "fastq")

        with open("{}/unknown_records.fastq".format(args.output_dir), "w") as output_handle:
            SeqIO.write(unknown_records, output_handle, "fastq")

    
    # print summary stats
    if args.print_summary_stats:
        total_records = len(target_records) + len(fragment_records) + len(spurious_records) + len(ont_records) + len(background_records) + len(short_records) + len(unknown_records)
        results_str = ''
        results_str = results_str + "> SUMMARY Statistics:\n"
        results_str = results_str + "> Examined {} lamplicons\n".format(total_records)
        results_str = results_str + ">  Found {}/{} suspected target lamplicons\n".format(len(target_records), total_records)
        results_str = results_str + ">  Found {}/{} suspected ont sequences\n".format(len(ont_records), total_records)
        results_str = results_str + ">  Found {}/{} sequences of unknown origin\n".format(len(unknown_records), total_records)
        results_str = results_str + ">  Disaggregate VAF (each target gets a vote): {}/{}\n".format(mut_total,(mut_total+wt_total))
        results_str = results_str + ">  Aggregate VAF (each lamplicon gets a vote): {}/{}\n".format(mut_reads,(mut_reads+wt_reads))
        results_str = results_str + ">  Concatemer Length Histogram:\n"
        for num_targets, target_read_list in target_sequence_map.items():
            results_str = results_str + ">    {} : {}\n".format(num_targets, len(target_read_list))

        print(results_str)

    # write stats to output file
    with open("{}/summary_stats.txt".format(args.output_dir), "w") as output_handle:
        output_handle.write(results_str)
            
    if args.save_concatemers:
        for num_targets, target_read_list in target_sequence_map.items():
            with open("{}/{}_targets.fastq".format(args.output_dir, num_targets), "w") as output_handle:
                SeqIO.write(target_read_list, output_handle, "fastq")

    

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("target_ref_fn")
    parser.add_argument("human_ref_fn")
    parser.add_argument("primer_set_fn")
    parser.add_argument("lamplicons_fn")
    parser.add_argument("--output_dir", type=str, default="results")
    parser.add_argument("--target_vcf_fn", type=str, default='')
    parser.add_argument("--ref_vcf_fn", type=str, default='')
    parser.add_argument("--max_lamplicons", type=int, default=-1)
    parser.add_argument("--threshold", type=float, default=0.80,
            help="Primer identity threshold for successful alignment")
    parser.add_argument("--debug_print", action="store_true", default=False)
    parser.add_argument("--save_target_reads", action="store_true", default=False)
    #parser.add_argument("--skip_consensus_polishing", action="store_true", default=False)
    parser.add_argument("--print_summary_stats", action="store_true", default=False)
    parser.add_argument("--save_concatemers", action="store_true", default=False)
    parser.add_argument("--high_confidence", action="store_true", default=False)
    
    return parser



if __name__ == "__main__":
    parser = argparser()
    args = parser.parse_args()

    
    #import cProfile, pstats
    #profiler = cProfile.Profile()
    #profiler.enable()

    main(args)

    #profiler.disable()
    #stats = pstats.Stats(profiler).sort_stats('tottime')
    #stats.print_stats()
    #stats.dump_stats('profile.dat')
