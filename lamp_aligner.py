#!/usr/bin/python

import sys, os, subprocess, argparse

import swalign
import mappy
import pysam

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

fwd_primer_order = ['F3', 'F2', 'F1', 'B1c', 'T', 'B2c', 'B3c']
rev_primer_order = ['B3', 'B2', 'Tc', 'B1', 'F1c', 'F2c', 'F3c']


#########################
class PrimerSet:
    '''
    Stores bases in a primer set.
    '''
    
    def __init__(self, primer_set_filename):

        # initialize primer set given file
        with open(primer_set_filename) as fp:
            for line in fp:

                # read lines in format <primer_name> <primer_seq>
                line = line.strip()
                if not line: continue
                fields = line.split(' ')
                primer_name = fields[0]
                primer_string = fields[1]

                # set each primer
                if primer_name == "F3":
                    self.F3 = primer_string
                elif primer_name == "B3":
                    self.B3 = primer_string
                elif primer_name == "FIP":
                    self.FIP = primer_string
                elif primer_name == "BIP":
                    self.BIP = primer_string
                elif primer_name == "FLP":
                    self.FLP = primer_string
                elif primer_name == "BLP":
                    self.BLP = primer_string
                elif primer_name == "F1":
                    self.F1 = primer_string
                elif primer_name == "F1c":
                    self.F1 = rc(primer_string)
                elif primer_name == "F2":
                    self.F2 = primer_string
                elif primer_name == "B1":
                    self.B1 = primer_string
                elif primer_name == "B1c":
                    self.B1 = rc(primer_string)
                elif primer_name == "B2":
                    self.B2 = primer_string
                elif primer_name == "T":
                    self.T = primer_string
                elif primer_name == "Tc":
                    self.T = rc(primer_string)
        
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

    alignment_tmp = aligner.align(seq, primer)
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
    if alignment.identity > identity_threshold and \
			alignment_len > len(primer) * identity_threshold:
        
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
def findAllPrimerAlignments(aligner, seq, primers, identity_threshold):

    alignment_list = list()
    
    alignment_list.extend(findPrimerAlignments(aligner, seq, primers.F1, "F1", identity_threshold))
    alignment_list.extend(findPrimerAlignments(aligner, seq, rc(primers.F1), "F1c", identity_threshold))
    alignment_list.extend(findPrimerAlignments(aligner, seq, primers.F2, "F2", identity_threshold))
    alignment_list.extend(findPrimerAlignments(aligner, seq, rc(primers.F2), "F2c", identity_threshold))
    alignment_list.extend(findPrimerAlignments(aligner, seq, primers.F3, "F3", identity_threshold))
    alignment_list.extend(findPrimerAlignments(aligner, seq, rc(primers.F3), "F3c", identity_threshold))
    alignment_list.extend(findPrimerAlignments(aligner, seq, primers.B1, "B1", identity_threshold))
    alignment_list.extend(findPrimerAlignments(aligner, seq, rc(primers.B1), "B1c", identity_threshold))
    alignment_list.extend(findPrimerAlignments(aligner, seq, primers.B2, "B2", identity_threshold))
    alignment_list.extend(findPrimerAlignments(aligner, seq, rc(primers.B2), "B2c", identity_threshold))
    alignment_list.extend(findPrimerAlignments(aligner, seq, primers.B3, "B3", identity_threshold))
    alignment_list.extend(findPrimerAlignments(aligner, seq, rc(primers.B3), "B3c", identity_threshold))
    #alignment_list.extend(findPrimerAlignments(aligner, seq, primers.BLP, "BLP", identity_threshold))
    #alignment_list.extend(findPrimerAlignments(aligner, seq, rc(primers.BLP), "BLPc", identity_threshold))
    #alignment_list.extend(findPrimerAlignments(aligner, seq, primers.FLP, "FLP", identity_threshold))
    #alignment_list.extend(findPrimerAlignments(aligner, seq, rc(primers.FLP), "FLPc", identity_threshold))
    alignment_list.extend(findPrimerAlignments(aligner, seq, primers.T, "T", identity_threshold))
    alignment_list.extend(findPrimerAlignments(aligner, seq, rc(primers.T), "Tc", identity_threshold))

    return(sorted(alignment_list))
    
###################
def printPrimerAlignments(seq, alignments):
    '''
    Fancy printing of all primers aligned to sequence.
    '''

    print("")
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

                # emit special char
                if found_start and found_end:
                    sys.stdout.write('X')
                    found_primer = True
                elif found_start:
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
            sys.stdout.write('\n\n')

            # increment wrap counter
            wrap_counter = wrap_counter + 1
            
    sys.stdout.write('\n')

###################
def extractAmpliconAroundTarget(alignments, target):
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
    primer_counter = fwd_primer_order.index('T') if fwd_strand else rev_primer_order.index('Tc')
    primer_counter = primer_counter + 1

    all_matched = True
    for i in range(target_index + 1, len(alignments)):

        # primer name
        primer_name = alignments[i].primer_name

        # expected primer name
        if fwd_strand:
            expected_primer_name = fwd_primer_order[primer_counter]
        else:
            expected_primer_name = rev_primer_order[primer_counter]

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
            if primer_counter == len(fwd_primer_order):
                break

    # if we matched all in the primer sequence, extend to the end of the entire read
    if all_matched:
        amplicon_end = -1

    #####
    # look to the left
    #####
    allowed_mismatches = 0
    primer_counter = fwd_primer_order.index('T') if fwd_strand else rev_primer_order.index('Tc')
    primer_counter = primer_counter - 1

    all_matched = True
    for i in range(target_index - 1, 0, -1):

        # primer name
        primer_name = alignments[i].primer_name

        # expected primer name
        if fwd_strand:
            expected_primer_name = fwd_primer_order[primer_counter]
        else:
            expected_primer_name = rev_primer_order[primer_counter]

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


#####
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


def main(args):

    # get file directory for relative script paths
    dirname = os.path.dirname(__file__)
    generate_consensus = os.path.join(dirname, 'scripts/generate_consensus.sh')
    sam_to_bam = os.path.join(dirname, 'scripts/sam_to_bam.sh')
    os.makedirs(args.output_dir, exist_ok=True)

    # parse all lamplicons from FASTQ file
    print("> parsing lamplicons")
    lamplicons = [l.seq for l in SeqIO.parse(args.lamplicons_fn, "fastq")]

    # parse primers from config file
    print("> parsing primers")
    primers = PrimerSet(args.primer_set_fn)

    # initialize aligner
    print("> initializing aligner")
    match = 2
    mismatch = -1
    scoring = swalign.NucleotideScoringMatrix(match, mismatch)
    sw = swalign.LocalAlignment(scoring)  # you can also choose gap penalties, etc...

    # iterate over all lamplicon sequences
    print("> splitting and mapping lamplicons")
    for lamp_idx, lamplicon in enumerate(lamplicons):

        # print status, stop early if we've found enough
        if lamp_idx >= args.max_lamplicons: break
        print("processing lamplicon {} of {}\r".format(lamp_idx+1, 
                min(len(lamplicons), args.max_lamplicons)), end="")
        
        # optionally print alignments of all primers
        if args.debug_print:
            alignments = findAllPrimerAlignments(sw, lamplicon, primers, args.threshold)
            printPrimerAlignments(lamplicon, alignments)

        # get just alignments of templates
        alignments = []
        alignments.extend(findPrimerAlignments(
                sw, lamplicon, primers.T, "T", args.threshold))
        alignments.extend(findPrimerAlignments(
                sw, lamplicon, rc(primers.T), "Tc", args.threshold))

        # gather targets into separate SeqIO Records
        target_idx = 0
        targets = []
        for alignment in alignments:
            if alignment.primer_name in ['T', 'Tc']:
                
                # extend sequence around target for better mapping
                seq_start, seq_end = extractAmpliconAroundTarget(alignments, alignment)

                # if we get a -1 for seq end, we never matched to the right, 
                # so just grab the whole read
                if seq_end == -1:
                    seq_end = len(lamplicon) - 1
                
                # generate new read from lamplicon substring
                record = SeqRecord(
                    lamplicon[seq_start:seq_end],
                    "{}_{}".format(lamp_idx, target_idx),
                    "{}_{}".format(lamp_idx, target_idx),
                    "lamplicon {} target {}".format(lamp_idx, target_idx)
                )
                targets.append(record)
                target_idx += 1

        # skip lamplicon if it didn't contain any targets
        if not target_idx: continue
        
        # write targets to new FASTA file
        fasta_fn = "{}/{}.fasta".format(args.output_dir, lamp_idx)
        with open(fasta_fn, "w") as output_handle:
            SeqIO.write(targets, output_handle, "fasta")
        
        # map all reads, including secondary mappings 
        out_sam_all = "{}/{}_all.sam".format(args.output_dir, lamp_idx)
        subprocess.run(["minimap2","-a","-xmap-ont","--eqx","-t 1","-w 1",
            "-n 1", "-N 5", "--secondary=no", "-m 0","-p 0.6","-s 30", "-o{}".
            format(out_sam_all),args.target_ref_fn,fasta_fn], 
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        # filter secondary/supplementary alignments from sam file
        out_sam = "{}/{}.sam".format(args.output_dir, lamp_idx)
        subprocess.run(["samtools", "view", "-h", "-F 0x900", out_sam_all, "-o{}".format(out_sam)])
        
        # convert sam to bam
        subprocess.run([sam_to_bam, out_sam], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        # no need to create consensus if we have one or fewer mapped targets
        out_bam = "{}/{}.bam".format(args.output_dir, lamp_idx)
        mapped_reads = pysam.AlignmentFile(out_bam, 'rb').mapped
        if mapped_reads <= 1: continue

        # generate consensus sequence if multiple aligned targets
        subprocess.run([generate_consensus, args.target_ref_fn, out_bam], 
                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        fastq_cns_fn = "{}/{}_cns.fastq".format(args.output_dir, lamp_idx)
        out_cns_sam = "{}/{}_cns.sam".format(args.output_dir, lamp_idx)

        # remove leading/trailing ambiguous characters from consensus sequence
        for consensus_seq in SeqIO.parse(fastq_cns_fn, "fastq"):

            # find location of last n from the start
            n_end = 0
            while consensus_seq.seq[n_end].upper() == 'N':
                n_end += 1
            n_start = len(consensus_seq.seq) - 1
            while consensus_seq.seq[n_start].upper() == 'N':
                n_start -= 1

            consensus_seq = consensus_seq[n_end:n_start+1]

            # re-write consensus fastq
            with open(fastq_cns_fn, "w") as output_handle:
                SeqIO.write(consensus_seq, output_handle, "fastq")

            break # should be a single consensus sequence
        
        # align consensus sequence
        subprocess.run(["minimap2","-a","-xmap-ont","--eqx","-t 1","-w 1", "-o{}"\
                .format(out_cns_sam), args.target_ref_fn, fastq_cns_fn], 
                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        # convert sam to bam
        subprocess.run([sam_to_bam, out_cns_sam], 
                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        # computing error from each alignment (ignoring leading/trailing clips and accounting for mutations
        orig_err = calcErrorRate(args.target_ref_fn, out_sam, 157, 'T')
        cns_err = calcErrorRate(args.target_ref_fn, out_cns_sam, 157, 'T')

    print("")



def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("target_ref_fn")
    parser.add_argument("primer_set_fn")
    parser.add_argument("lamplicons_fn")
    parser.add_argument("--output_dir", type=str, default="results")
    parser.add_argument("--max_lamplicons", type=int, default=50)
    parser.add_argument("--threshold", type=float, default=0.75,
            help="Primer identity threshold for successful alignment")
    parser.add_argument("--debug_print", action="store_true", default=False)

    return parser



if __name__ == "__main__":
    parser = argparser()
    args = parser.parse_args()
    main(args)

