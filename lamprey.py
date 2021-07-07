#
import sys, os, subprocess, argparse, collections
import cProfile
import random
import heapq
import time
import multiprocessing
from progressbar import progressbar

import random as rng

# use local swalign?
sys.path.insert(0, './lib')
import swalign 
from skbio.alignment import StripedSmithWaterman


import mappy
import bwapy
import pysam
import vcfpy
import numpy as np
import statsmodels.api as sm
import datetime

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Sequencing.Applications import SamtoolsViewCommandline

#########################
cigarOpMap = {'M' : 0,
              'I' : 1,
              'D' : 2,
              'N' : 3,
              'S' : 4,
              'H' : 5,
              'P' : 6,
              '=' : 7,
              'X' : 8,
              'B' : 9}

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
    sys.stdout.write(bcolors.BLUE + text + bcolors.ENDC)

def print_cyan(text):
    sys.stdout.write(bcolors.CYAN + text + bcolors.ENDC)

def print_red(text):
    sys.stdout.write(bcolors.FAIL + text + bcolors.ENDC)

def print_green(text):
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


#########################
class Result:
    '''
    Stores lamplicon classification, pileup info, and other post processing information.
    '''

    def __init__(self, classification='unknown', mut_count=0, wt_count=0, pileup_str='', target_depth=0, alignments=[], seq=None, idx=0, read_id = '', timestamp=None):
        self.classification = classification
        self.mut_count = mut_count
        self.wt_count = wt_count
        self.plural_base = ''
        self.plural_base_support = 0
        self.pileup_str = pileup_str
        self.target_depth = target_depth
        self.target_seq_accuracy = (0,0)
        self.bad_calls = None
        self.polished_bad_calls = None
        self.alignments = alignments
        self.seq = seq
        self.polished_seq = None
        self.idx = idx
        self.read_id = read_id
        self.timestamp = timestamp

    def __lt__(self, other):
        return self.timestamp < other.timestamp

##########################
def printHeader():
    s = ''
    s += '    __    ___    __  _______\n'
    s += '   / /   /   |  /  |/  / __ \________  __  __\n'
    s += '  / /   / /| | / /|_/ / /_/ / ___/ _ \/ / / /\n'
    s += ' / /___/ ___ |/ /  / / ____/ /  /  __/ /_/ /\n'
    s += '/_____/_/  |_/_/  /_/_/   /_/   \___/\__, /\n'
    s += '                                    /____/   \n'
    s += ' by Jack Wadden\n'
    s += ' Version 0.1\n'
    
    print(s)

#######
def cigarStringToTuples(cigar):
    '''
    Converts a cigar string to a pysam endocded list of cigar tuples.
    '''
    tuples = list()

    # walk the cigar string
    is_numeric = True
    count_str = ''
    for i in range(len(cigar)):
        if cigar[i].isnumeric():
            count_str += cigar[i]
        else:
            op = cigarOpMap[cigar[i]]
            count = int(count_str)
            tuples.append((op,count))
            count_str = ''

    return tuples

#######
def alignSequence_swalign(aligner, seq, primer, primer_name):
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
def alignSequence_skbio(seq, primer, primer_name):
    '''
    Aligns two sequences using the skbio Smith-Waterman library.
    '''

    query = StripedSmithWaterman(primer)
    alignment = query(str(seq))

    # compute matches between aligned query/target
    matches = 0
    for i in range(0, len(alignment.aligned_query_sequence)):
        if alignment.aligned_query_sequence[i] == alignment.aligned_target_sequence[i] :
            matches = matches + 1

    identity = float(matches)/float(len(primer))

    ret_alignment = Alignment(alignment.target_begin,
                              alignment.target_end_optimal + 1,
                              identity,
                              primer_name)

    return ret_alignment

#######
def findPrimerAlignments(aligner, seq, primer, primer_name, identity_threshold, args):
    ''' 
    Greedy approach which finds the best primer alignment for a long sequence, 
    then uses left/right recursion to find all other possible (non-overlapping) 
    alignments. 
    '''

    # find optimal alignment
    #print("*")
    #
    #print(alignment)
    if args.swalign:
        alignment = alignSequence_swalign(aligner, seq, primer, primer_name)
    else:
        alignment = alignSequence_skbio(seq, primer, primer_name)
    #print(alignment)
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
                aligner, left_seq, primer, primer_name, identity_threshold, args) 
            alignment_list.extend(left_alignments)
        
        # recurse right
        if len(right_seq) >= len(primer):
            right_alignments = findPrimerAlignments(
                    aligner, right_seq, primer, primer_name, identity_threshold, args)

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

    alignment_list.extend(findPrimerAlignments(aligner, seq, primers.primer_dict["T"], "T", identity_threshold, args))
    alignment_list.extend(findPrimerAlignments(aligner, seq, rc(primers.primer_dict["T"]), "Tc", identity_threshold, args))

    # bail if we didn't find a target in this read
    # this is a lazy shortcut optimization
    if len(alignment_list) == 0:
        return(sorted(alignment_list), 0.0)

    if args.high_confidence and len(alignment_list) < 3:
        return(sorted(alignment_list), 0.0)
    
    for primer_name, primer_seq in primers.primer_dict.items():
        if primer_name == 'T':
            continue
        else:
            # fwd
            alignment_list.extend(findPrimerAlignments(aligner, seq, primer_seq, primer_name, identity_threshold, args))
            # rc
            alignment_list.extend(findPrimerAlignments(aligner, seq, rc(primer_seq), primer_name + "c", identity_threshold, args))


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
def removeOverlappingPrimerAlignments(debug_print, alignments, allowed_overlap):

    if debug_print:
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

            if debug_print:
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
def countCallsFromPileup(ref_char, pileup_str):

    pileup_calls = dict()
    pileup_calls['a'] = 0
    pileup_calls['t'] = 0
    pileup_calls['g'] = 0
    pileup_calls['c'] = 0
    pileup_calls['*'] = 0 # deletion
    
    #
    total_count = 0
    skip_counter = 0
    set_skip = False

    for c in pileup_str:

        if(skip_counter > 0):
            skip_counter = skip_counter - 1
            continue

        if(set_skip):
            skip_counter = int(c)
            set_skip = False
            continue

        total_count = total_count + 1

        if c == ',' or c == '.':
            pileup_calls[ref_char.lower()] = pileup_calls[ref_char.lower()] + 1

        elif c.lower() in pileup_calls.keys():
            pileup_calls[c.lower()] = pileup_calls[c.lower()] + 1

        elif c == '-' or c == '+':
            set_skip = True

        elif c == '^':
            skip_counter = 1
            
    return pileup_calls

################################################
def getPluralityVoter(calls):

    total_calls = sum(calls.values())

    # is there a plurality voter?
    plural_bases = set()
    majority_base = ''
    majority_vote = False
    max_count = 0
    for base, count in calls.items():

        if count > max_count:
            max_count = count
            majority_base = base
            majority_vote = True
                                                                
        if count == max_count:
            plurality_vote = False

    for base, count in calls.items():
        if count == max_count:
            plural_bases.add(base)
                                
    # if there was a plurality vote, pick that
    # else choose, random call among the tied plurality
    if majority_vote == False and plurality_vote == False:
        majority_base = random.choice(plural_bases)

    return majority_base, max_count

################################################
def extractPileupInfo(pileup_fn, contig, target_locus, limit):
    # parse lines in the pileup, extract the depth at the locus

    ref = ''
    depth = 0
    code_string = ''
    covers_roi = False
    within_roi = False

    vote_count = 0
    ref_match_count = 0

    bad_calls = dict({'a' : 0, 't' : 0, 'c' : 0, 'g' : 0, '*' : 0, 'ref' : 0})
    polished_bad_calls = dict({'a' : 0, 't' : 0, 'c' : 0, 'g' : 0, '*' : 0, 'ref' : 0})
    
    with open(pileup_fn) as fp:
        for line in fp:
            fields = line.split('\t')
            if fields[0] == contig:

                pileup_locus = int(fields[1])

                # retrieve target specific information
                if pileup_locus == target_locus:
                    ref = fields[2]
                    depth = int(fields[3])
                    code_string = fields[4]
                    covers_roi = True
                elif pileup_locus < (target_locus + limit) and pileup_locus > (target_locus - limit):
                    within_roi = True

                # consider entire read accuracy if the number of votes is 2+
                locus_depth = int(fields[3])

                if locus_depth > 0 and pileup_locus != target_locus:

                    locus_code_string = fields[4]
                    locus_ref = fields[2].lower()
                    
                    calls = countCallsFromPileup(locus_ref, locus_code_string)
                    majority_base, count = getPluralityVoter(calls)
                    

                    # is majority/random base right or wrong?
                    if majority_base == locus_ref.lower():
                        ref_match_count = ref_match_count + 1
                        polished_bad_calls['ref'] = polished_bad_calls['ref'] + 1
                    else:
                        polished_bad_calls[majority_base] = polished_bad_calls[majority_base] + 1
                        
                    vote_count = vote_count + 1

                    # collect all incorrect calls
                    
                    for key in bad_calls:
                        if key in calls:
                            if key == locus_ref.lower():
                                bad_calls['ref'] = bad_calls['ref'] + calls[key]
                            else:
                                bad_calls[key] = bad_calls[key] + calls[key]
                    
                        
    return ref, depth, code_string, covers_roi, within_roi, ref_match_count, vote_count, bad_calls, polished_bad_calls

################################################
def proportionConfidenceInterval(mut, wt, conf):


    conf_int = sm.stats.proportion_confint(count=mut,
                                           nobs=mut+wt,
                                           alpha=(1 - conf))

    if (mut + wt) > 0:
        samp_mean = float(mut)/float(mut + wt) * 100.0
    else:
        samp_mean = 0.0
        
    if conf_int[0] == 0.0:
        low = 0.0
    else:
        low = samp_mean - conf_int[0] * 100.0
    high = conf_int[1] * 100.0 - samp_mean

    #print("conf:", conf, "%.2f" % low, " -- ", "%.2f" % samp_mean, " -- ", "%.2f" % high, ")")
    #print("conf:", conf, "%.4f" % (conf_int[0]*100.0), " -- ", "%.2f" % samp_mean, " -- ", "%.4f" % (conf_int[1]*100.0), ")")

    return conf_int[0], samp_mean, conf_int[1]

def isStatisticallySignificant(mut, wt, conf, err):

    total = mut + wt
    err_low, err_mean, err_high = proportionConfidenceInterval(err * float(total), (1.0 - err) * float(total), conf)
    obs_low, obs_mean, obs_high = proportionConfidenceInterval(mut, wt, conf)

    if (obs_low > err_high):
        return True
    else:
        return False
    
    
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

######################################
def parseTimestamp(timestamp):
    # format: 2019-07-10T20:53:45Z
    parsed_time = datetime.datetime.strptime(timestamp, "%Y-%m-%dT%H:%M:%SZ")
    return parsed_time

def timestampDelta(start, end):
    delta = end - start
    return delta.total_seconds()

##############################################
def polishLamplicon(lamplicon, subread_alignment_file, alignment_interval_dict):

    # generate new read from lamplicon.seq substring
    polished_lamplicon = SeqRecord(lamplicon.seq,
                                   id = lamplicon.id,
                                   name = lamplicon.name,
                                   description = lamplicon.description
    )

    
    # parse the samfile into pysam
    pileup_alignments = pysam.AlignmentFile(subread_alignment_file, "r")

    # generate pileup over target region
    pileup_iter = pileup_alignments.pileup(truncate=True,
                                           stepper="samtools",
                                           min_base_quality=0,
                                           ignore_overlaps=False,
                                           ignore_orphans=False)

    deletion_removal_list = list()
    
    # iterate over the pileup
    deletion_streak_counters = dict()

    # for each pileup position
    for column in pileup_iter:

        calls = dict({'a' : 0, 't' : 0, 'c' : 0, 'g' : 0, '*' : 0, 'ref' : 0})
        call_count = 0

        # get the position in the sub-read for this column
        # NOTE: deletions return the location "before" the deletion position in the alignment
        query_positions = column.get_query_positions()
        
        # gather all votes and positions in mother read
        for read in column.pileups:
            
            # deletion
            if read.is_del:
                calls['*'] = calls['*'] + 1
                call_count = call_count + 1
            # insertion (skip)
            elif read.is_refskip:
                continue
            # mismatch
            else:
                #print(read.query_position)
                #print(read.alignment.query_sequence)
                base = read.alignment.query_sequence[read.query_position].lower() 
                #print("BASE: {}".format(base))
                #print(read.alignment.query_sequence)
                calls[base] = calls[base] + 1
                call_count = call_count + 1

        # compute plurality vote
        vote_base, plurality_count = getPluralityVoter(calls)
        #print(vote_base, plurality_count)


        # if we end up with zero votes, skip this location
        if plurality_count < 1:
            continue
                
        # polish bases in mother read
        read_idx = 0
        for read in column.pileups:

            # get sub-read direction info
            read_start, read_end, direction = alignment_interval_dict[read.alignment.query_name]

            # where does this base lie?
            offset = query_positions[read_idx]            

            # if the sub-read aligns to the mother in the forard direction
            if direction == 'fwd':
                base = vote_base.upper()
                # if the sub-read aligns to the reference in the forward
                if not read.alignment.is_reverse:
                    mother_pos = read_start + offset
                else:
                    mother_pos = read_end - 1 - offset

            # if the sub-read aligns to the mother in the reverse direction
            else:
                base = rc(vote_base.upper())
                # if the sub-read aligns to the reference in the forward
                if not read.alignment.is_reverse:
                    mother_pos = read_start + offset
                # if the sub-read aligns to the reference in the reverse direction
                else:
                    mother_pos = read_end - 1 - offset

            # if we need to polish a deletion, save it and we'll handle it later
            # if we are a deletion, and the polished base is a deletion, skip/ignore
            if read.is_del and base != '*':
                if len(deletion_removal_list) > 0 and deletion_removal_list[-1][0] == mother_pos:
                    deletion_removal_list[-1][1].append(base)
                else:
                    deletion_removal_list.append((mother_pos, [base], direction))
                    
            # polish the base in the mother read
            else:
                    
                # adjust mother read
                polished_lamplicon.seq = polished_lamplicon.seq[:mother_pos] + base + polished_lamplicon.seq[mother_pos + 1:]
            
            # cleanup
            read_idx = read_idx + 1

    # for each deletion, make sure to insert the voted base in the correct
    # order at the right index. Offsets/coordinates are always in the read direction

    # sort so we add left->right in the mother read
    deletion_removal_list = sorted(deletion_removal_list)
    #print(deletion_removal_list)

    # offset counter accounts for inserted bases into the mother read that weren't there before
    # anytime we insert, we adjust the offset counter
    offset_counter = 0
    for position, base_list, direction in deletion_removal_list:

        # if we're a reverse aligned read, deletions are in the opposite order
        # add in 
        #  starting at the beginning (reverse end) if there's a run of deletions
        if direction == 'rev':

            base_list.reverse()
            for base in base_list:

                # add 1 to offset because we flipped directions
                final_position = position + offset_counter + 1            
                polished_lamplicon.seq = polished_lamplicon.seq[:final_position] + base + polished_lamplicon.seq[final_position:]
                
                # update offset counter because we added a base
                offset_counter = offset_counter + 1

        # if we're a forward aligned read, add in forward order as encountered
        else:

            for base in base_list:
            
                final_position = position + offset_counter
                polished_lamplicon.seq = polished_lamplicon.seq[:final_position] + base + polished_lamplicon.seq[final_position:]
                
                # update offset counter because we added a base
                offset_counter = offset_counter + 1
                
            
    # actually delete any '*' deletions from the polished read
    polished_lamplicon.seq = Seq(str(polished_lamplicon.seq).replace("*",""))

    return polished_lamplicon

################################################
def processLamplicon(sw, process_candidates_path, generate_consensus_path, minimap2, bwamem, lamp_idx, lamplicon, primers, ref_vcf_record, target_vcf_record, args):

    if args.debug_print:
        print_green("Processing Lamplicon {}".format(lamp_idx))
    
    result = Result(seq=lamplicon)

    # fragment limit
    # the distance between the locus and th start or end of the alignment to be considered a fragment
    fragment_limit = 1000

    alignments, alignment_coverage = findAllPrimerAlignments(sw, lamplicon.seq, primers, args.threshold, args)
    if len(alignments) > 0:
        # remove alignments that overlap by more than 4bp
        alignments = removeOverlappingPrimerAlignments(args.debug_print, alignments, 4)

    # save alignments
    result.alignments = alignments
        
    # optionally print alignments of all primers
    if args.debug_print:
        printPrimerAlignments(lamplicon.seq, alignments)

    # does this lamplicon possibly cover the target?
    num_ont_seqs = 0

    for alignment in alignments:
        if alignment.primer_name in ['T', 'Tc']:
            result.target_depth = result.target_depth + 1

        # track high quality ont sequences
        if alignment.primer_name.startswith('ont') and alignment.identity >= args.threshold:
            num_ont_seqs = num_ont_seqs + 1

    # if high confidence mode, only proceed if there are 2+ targets
    if args.high_confidence and result.target_depth < 3:
        return Result(target_depth=0, classification='unknown', mut_count=0, wt_count=0, pileup_str='')

                
    ####################################
    # if we found a suspected target seed
    # attempt to extend seed, and align to target region
    if result.target_depth > 0:

        ###############################################################################
        #  STAGE 1: extend target seed to left/right to identify candidate sub-reads  #
        ###############################################################################
        # extend each target into a single amplicon and emit as fasta
        target_idx = 0
        targets = []
        oriented_targets = []
        alignment_intervals = list()
        alignment_interval_dict = dict()
        for alignment in alignments:

            if alignment.primer_name == 'T' or alignment.primer_name == 'Tc':
        
                # extend sequence around target for better mapping
                seq_start, seq_end = extractAmpliconAroundTarget(primers, alignments, alignment)
        
                # if we get a -1 for seq end, we never matched to the right, 
                # so just grab the whole read
                if seq_end == -1:
                    seq_end = len(lamplicon.seq) - 1

                subread_id = "{}_{}".format(lamp_idx, target_idx)
                if alignment.primer_name == 'T':
                    direction = 'fwd'
                else:
                    direction = 'rev'
                alignment_intervals.append((seq_start, seq_end, subread_id))
                alignment_interval_dict[subread_id] = (seq_start, seq_end, direction)
            
                # generate new read from lamplicon.seq substring
                record = SeqRecord(
                    lamplicon.seq[seq_start:seq_end],
                    subread_id,
                    subread_id,
                    "lamplicon {} target {}".format(lamp_idx, target_idx)
                )

                targets.append(record)
                oriented_targets.append((record, direction))
                target_idx += 1

        # Print intervals
        if args.debug_print:
            print("* Found {} candidate sub-reads at intervals:".format(len(alignment_intervals)))
            for interval in alignment_intervals:
                print(interval)
            print("", flush=True)

        # write separated targets to new FASTA file
        fasta_fn = "{}/{}.fasta".format(args.output_dir, lamp_idx)
        with open(fasta_fn, "w") as output_handle:
            SeqIO.write(targets, output_handle, "fasta")


        ######################################################################
        #  STAGE 2: Map each sub-read candidate against the target reference #
        ######################################################################
        # map each sub read to the target sequence; ignore secondary mappings 
        if args.debug_print:
            print("* Aligning candidate sub-reads...")

        # read in target ref string (TODO THIS IS SLOW)
        target_ref_str = ''
        target_ref_id = ''
        with open (args.target_ref_fn, 'r') as target_ref_fasta:
            data = target_ref_fasta.readlines()
            target_ref_id = data[0][1:].strip()
            target_ref_str = data[1]
            
        # align the reads using skbio alignment
        # BUG this doesn't align revcomp
        # really should just rip this out for bwamem
        target_alignments = list()
        for target, direction in oriented_targets:

            # get string and align appropriate fwd version
            if direction == 'rev':
                query = StripedSmithWaterman(rc(str(target.seq)))
            else:
                query = StripedSmithWaterman(str(target.seq))
            target_alignment = query(target_ref_str)

            cigar_tuples = cigarStringToTuples(target_alignment.cigar)
            count = 0
            for letter, value in cigar_tuples:
                if letter == 2:
                    count -= value
                else:
                    count += value
           
            aligned_query_string = target_alignment.aligned_query_sequence.replace('-','')

            # generate samfile from alignments
            sam_alignment = pysam.AlignedSegment()
            sam_alignment.query_name = target.id
            sam_alignment.query_sequence = aligned_query_string
            sam_alignment.flag = 0
            sam_alignment.reference_id = 0
            sam_alignment.reference_start = target_alignment.target_begin
            sam_alignment.mapping_quality = 60
            sam_alignment.cigar = cigarStringToTuples(target_alignment.cigar)
            sam_alignment.next_reference_id = 0
            sam_alignment.next_reference_start = 0
            sam_alignment.template_length = len(target_ref_str)
            qual_str = "<" * len(aligned_query_string)
            sam_alignment.query_qualities = pysam.qualitystring_to_array(qual_str)
            sam_alignment.tags = (("NM", 1),
                                  ("RG", "L1"))

            target_alignments.append(sam_alignment)

        # write alignments out to sam file
        ref_len = len(target_ref_str)
        header = { 'HD': {'VN': '1.0'},
                   'SQ': [{'LN': ref_len, 'SN': target_ref_id}]
        }

        sam_fn = "{}/{}_all.sam".format(args.output_dir, lamp_idx)
        unsorted_bam_fn = "{}/{}_all.unsorted.bam".format(args.output_dir, lamp_idx)
        bam_fn = "{}/{}_all.bam".format(args.output_dir, lamp_idx)
        # output pileup generated by process_candidates.sh
        pileup_fn = "{}/{}.pileup".format(args.output_dir, lamp_idx)
        
        with pysam.AlignmentFile(sam_fn, "w", header=header) as outf:
            for alignment in target_alignments:
                outf.write(alignment)


        # filter
        pysam.view("-F", "0x900", "-S", "-b", "-o", unsorted_bam_fn, sam_fn, catch_stdout=False)
        # sort
        pysam.sort("-o", bam_fn, unsorted_bam_fn, catch_stdout=False)
        # index
        pysam.index(bam_fn, catch_stdout=False)
        # generate a pileup file
        pysam.mpileup("-Q", "0", "-f", args.target_ref_fn, bam_fn, "-o", pileup_fn, catch_stdout=False)
            
        ###############################
        #  Polish orig read and emit  #
        ###############################
        result.polished_seq = polishLamplicon(lamplicon, bam_fn, alignment_interval_dict)

        ###############################
        # extract pileup info
        limit = 1000
        ref, depth, code_string, covers_roi, within_roi, correct_voted_bases, total_voted_bases, bad_calls, polished_bad_calls = extractPileupInfo(pileup_fn, target_vcf_record.CHROM, target_vcf_record.POS, limit)

        result.target_seq_accuracy = (correct_voted_bases, total_voted_bases)
        result.bad_calls = bad_calls
        result.polished_bad_calls = polished_bad_calls

        # update target depth to be aligned target depth
        result.target_depth = depth
        
        for alt in target_vcf_record.ALT:
            variant = alt.value.upper()

        if args.debug_print:
            print("  - Found {} reads that cover the target".format(result.target_depth))
            print("  - Pileup result: ",ref, result.target_depth, code_string, covers_roi)

            
        # in high confidence mode, ignore pileups that have less than 2 entries
        if args.high_confidence and result.target_depth < 2:
            #0, 'unknown', 0, 0
            # TODO fix result assignment to zero here
            return result 

        # a sub-read that successfully aligns to and covers the target locus is confirmed as a target sequence
        if result.target_depth > 0:
            result.classification = 'target'
        
            # parse code string and set mut/wt counts
            wt, mut = parsePileupCodeString(code_string, variant)
            calls = countCallsFromPileup(ref, code_string)
            plural_base, plural_base_support = getPluralityVoter(calls)
            #

            result.plural_base = plural_base
            result.plural_base_support = plural_base_support
            result.mut_count = mut
            result.wt_count = wt
            result.pileup_str = code_string
            result.target_depth = depth

            if args.debug_print:
                print("  - Found {} mut and {} wt calls".format(result.mut_count, result.wt_count))
                print("  - Plural basecall: {}".format(result.plural_base))


    #######################################
    # if we never found an alignable target sequence
    # attempt to diagnose the issue and classify the read as something other than target
    if result.target_depth == 0:        

        if args.debug_print:
            print("* Didn't find any alignable targets.... diagnosing...")
            print(" - Attempting to align to reference genome...")
        # this is maybe a background genomic read or a lamplicon "fragment"
        # fragments occur due to fragmentation-based library prep or shearing
        # try to align the read to see if it's a background genomic read
        # if we have any hits, classify as background
        
        for hit in minimap2.map(lamplicon.seq):

            if args.debug_print:
                print(hit)
            
            # skip non-primary alignments
            if not hit.is_primary:
                continue
            
            # if the hit is within fragment_limit of the target, mark it as a fragment
            if is_hit_near_target(hit, ref_vcf_record, fragment_limit):
                result.classification = 'fragment'
                if args.debug_print:
                    print(" - Hit is near the target region but does not contain target. Diagnosing as fragment.")
                return result
            
            # else, just mark it as background
            result.classification = 'background'
            if args.debug_print:
                print(" - Hit is not near target. Diagnosing as background.")
            return result

        if args.debug_print:
            print(" - Could not find a genomic hit...")
            
        # classify read as "ont" if it contains > 1 ONT sequence and fewer than 2 other suspected primer seqs and greater than 50% coverage
        if num_ont_seqs > 0 and (len(alignments) - 2) <= num_ont_seqs :

            # get coverage
            if alignment_coverage >= 0.15:
                if args.debug_print:
                    print(" - ONT seqs dominate alignments. Diagnosing as ONT")
                
                result.classification = 'ont'
            else:
                # assume unknown
                if args.debug_print:
                    print(" - Diagnosing as unknown.")

                result.classification = 'unknown'

                # try to align the read to see if it's a background genomic read
                for hit in minimap2.map(lamplicon.seq):
                    if args.debug_print:
                        print(" - Read actually maps. Adjusting diagnosis to background.")

                    result.classification = 'background'

            return result
                    
            
        # classify the read as spurious (bad/erroneous primer amplification) if it has at least 3 primer sequences, didn't align at the locus, and no alignable target
        if len(alignments) >= 3:
            if args.debug_print:
                print(" - Diagnosing as spurious.")

            result.classification = 'spurious'
            return result

        # suspiciously short sequences might just be a fragment that can't align with enough identity
        if len(lamplicon.seq) < 60:
            if args.debug_print:
                print(" - Diagnosing as short.")

            result.classification = 'short'
        else:
            if args.debug_print:
                print(" - Diagnosing as unknown.")

            # final default to unknown
            result.classification = 'unknown'

        # Return the final classification
        return result


    return result


#######
def initAligners(args):
    match = 2
    mismatch = -1
    
    scoring = swalign.NucleotideScoringMatrix(match, mismatch)

    sw = swalign.LocalAlignment(scoring)  # you can also choose gap penalties, etc...

    # set up a mappy instance to align reads
    minimap2 = mappy.Aligner(args.human_ref_fn)  # load or build index
    if not minimap2: raise Exception("ERROR: failed to load/build index")

    # set up bwamem
    #index = "/data/human_reference/H3F3A.fa"
    #options = "-k {}".format(15)
    #bwamem = bwapy.BwaAligner(index, options)

    
    return sw, minimap2, None


#####
def processLamplicons(args, primers, ref_vcf_record, target_vcf_record, lamplicon_batch, q):

    # get file directory for relative script paths
    dirname = os.path.dirname(__file__)
    process_candidates_path = os.path.join(dirname, 'scripts/process_candidates.sh')
    generate_consensus_path = os.path.join(dirname, 'scripts/generate_consensus.sh')
    os.makedirs(args.output_dir, exist_ok=True)

    # initialize aligner
    if args.debug_print:
        print("> initializing aligners")
    sw, minimap2, bwamem = initAligners(args)

    results = []
    
    for lamp_idx,lamplicon in lamplicon_batch:
        
        #print_green("Processing Lamplicon {} \r".format(lamp_idx + 1))

        if args.time_sort:
            timestamp = parseTimestamp(lamplicon.description.split()[5].split('=')[1])
        else:
            timestamp = 0
            
        result = processLamplicon(sw,
                                  process_candidates_path,
                                  generate_consensus_path,
                                  minimap2,
                                  bwamem,
                                  lamp_idx,
                                  lamplicon,
                                  primers,
                                  ref_vcf_record,
                                  target_vcf_record,
                                  args);

        # add to list of results
        result.read_id = lamplicon.name
        result.idx = lamp_idx
        result.timestamp = timestamp
        results.append(result)

        q.put(result)
        
            
#####
def runProcess(thread_id, q, args, primers, rev_vcf_record, target_vcf_record, lamplicon_batch):

    processLamplicons(args, primers, rev_vcf_record, target_vcf_record, lamplicon_batch, q)


def main(args):

    # print header
    printHeader()
    
    # load and parse files
    # parse all lamplicons from FASTQ file
    print("> Parsing lamplicon file: {}".format(args.lamplicons_fn))
    #lamplicon_seqs = [l.seq for l in SeqIO.parse(args.lamplicons_fn, "fastq")]

    if args.lamplicons_fn.endswith(".fastq"):
        raw_lamplicons = SeqIO.parse(args.lamplicons_fn, "fastq")
    elif args.lamplicons_fn.endswith(".fasta"):
        raw_lamplicons = SeqIO.parse(args.lamplicons_fn, "fasta")

    # give each lamplicon an id
    lamplicons = []
    id_counter = 0
    for lamplicon in raw_lamplicons:
        lamplicons.append((id_counter, lamplicon))
        id_counter = id_counter + 1

    print("  * found {} reads.".format(len(lamplicons)))
        
    # parse primers from config file
    print("> Parsing LAMP assay schema file: {}".format(args.primer_set_fn))
    primers = PrimerSet(args.primer_set_fn)

    print("> Parsing target mutation VCF files:")
    # parse ref VCF file if there is one
    if args.ref_vcf_fn != '':
        print("  * {}".format(args.ref_vcf_fn))
        ref_vcf_record = parseVCF(args.ref_vcf_fn)

    # parse VCF file if there is one
    if args.target_vcf_fn != '':
        print("  * {}".format(args.target_vcf_fn))
        target_vcf_record = parseVCF(args.target_vcf_fn)
    

    # 
    q = multiprocessing.Queue()
    processes = []
    results = []
    num_processes = args.num_threads

    # divvy up work
    lamplicon_batch = []
    for i in range(num_processes):
        lamplicon_batch.append(list())
        
    for i in range(len(lamplicons)):
        lamplicon_batch[i % num_processes].append(lamplicons[i])

    # allocate processes
    print("> Launching {} processes, initializing aligners, and processing reads...".format(num_processes))
    for i in range(num_processes):
        processes.append(multiprocessing.Process(target = runProcess, args = (i, q, args, primers, ref_vcf_record, target_vcf_record, lamplicon_batch[i])))

    # launch processes
    for i in range(num_processes):
        processes[i].start()

    # collect results
    # priority queue will sort based on timestamp
    # reads will always have a timestamp set to 0 if time_sort flag is not set
    if not args.debug_print:
        for i in progressbar(range(len(lamplicons))):
            heapq.heappush(results, q.get(True))
            i = i + 1
    else:
        for i in range(len(lamplicons)):
            heapq.heappush(results, q.get(True))
        
    # wait for jobs to complete
    for i in range(num_processes):
        processes[i].join()

    print("> Finished processing all reads.")
        
    # process results
    if args.time_sort:
        results.sort()
    
    classification_counters = dict({'target' : 0, 'fragment' : 0, 'background' : 0, 'ont' : 0, 'spurious' : 0, 'short' : 0, 'unknown' : 0})
    target_counter = 0
    mut_count = 0
    wt_count = 0

    statistical_significance_map = dict()
    num_targets_map = dict()
    all_bad_calls = dict({'a' : 0, 't' : 0, 'c' : 0, 'g' : 0, '*' : 0, 'ref' : 0})
    all_polished_bad_calls = dict({'a' : 0, 't' : 0, 'c' : 0, 'g' : 0, '*' : 0, 'ref' : 0})
    plural_target_calls = dict({'a' : 0, 't' : 0, 'c' : 0, 'g' : 0, '*' : 0, 'ref' : 0})
    plural_mut_count = 0
    plural_wt_calls = 0
    
    polished_lamplicons = list()
    
    total_correct_calls = 0
    total_calls = 0
    
    time_ordered_index = 0
    for result in results:

        classification_counters[result.classification] = classification_counters[result.classification] + 1
        
        if result.classification == 'target':
            target_counter = target_counter + 1

        # build histogram of target count
        if result.target_depth not in num_targets_map.keys():
            num_targets_map[result.target_depth] = 1
        else:
            num_targets_map[result.target_depth] = num_targets_map[result.target_depth] + 1

            
        # mut/wt counting. majority vote among sections.
        if result.target_depth > 0:
            if float(result.mut_count)/float(result.target_depth) > 0.5:
                mut_count = mut_count + 1

            if float(result.wt_count)/float(result.target_depth) > 0.5:
                wt_count = wt_count + 1

        if result.plural_base.lower() in plural_target_calls.keys():
            if args.high_confidence:
                if result.plural_base_support > 2:
                    plural_target_calls[result.plural_base.lower()] = plural_target_calls[result.plural_base.lower()] + 1
            else:
                plural_target_calls[result.plural_base.lower()] = plural_target_calls[result.plural_base.lower()] + 1
                        
        total_correct_calls = total_correct_calls + result.target_seq_accuracy[0]
        total_calls = total_calls + result.target_seq_accuracy[1]

        # collect all bad calls to identify basecall error trends
        if result.bad_calls is not None:
            for key in all_bad_calls:
                if key in result.bad_calls:
                    all_bad_calls[key] = all_bad_calls[key] + result.bad_calls[key]

        # collect all bad calls to identify basecall error trends
        if result.polished_bad_calls is not None:
            for key in all_polished_bad_calls:
                if key in result.polished_bad_calls:
                    all_polished_bad_calls[key] = all_polished_bad_calls[key] + result.polished_bad_calls[key]


        # collect polished lamplicons to emit
        if result.polished_seq is not None:
            polished_lamplicons.append(result.polished_seq)
        else:
            polished_lamplicons.append(result.seq)
        #
        err = 0.015
        if mut_count + wt_count > 0:
            
            # confidence intervals
    
            conf = 0.95
        
            if isStatisticallySignificant(mut_count, wt_count, conf, err) and (conf not in statistical_significance_map.keys()):
                statistical_significance_map[conf] = time_ordered_index + 1

            conf = 0.99
            if isStatisticallySignificant(mut_count, wt_count, conf, err) and (conf not in statistical_significance_map.keys()):
                statistical_significance_map[conf] = time_ordered_index + 1

            conf = 0.999
            if isStatisticallySignificant(mut_count, wt_count, conf, err) and (conf not in statistical_significance_map.keys()):
                statistical_significance_map[conf] = time_ordered_index + 1

            conf = 0.9999
            if isStatisticallySignificant(mut_count, wt_count, conf, err) and (conf not in statistical_significance_map.keys()):
                statistical_significance_map[conf] = time_ordered_index + 1

        time_ordered_index = time_ordered_index + 1

        
    # Summary results
    low, samp_mean, high = proportionConfidenceInterval(mut_count, wt_count, 0.95)
        
    # print results
    s = ""
    s += "|---------------------|\n"
    s += "|   Summary Results   |\n"
    s += "|---------------------|\n"
    s += " Target Fraction: {:.2f}%\n".format(float(target_counter)/float(len(lamplicons)) * 100.0)
    s += " Classifications: "
    s += str(classification_counters) + "\n"

    # VAF calculations
    mut_base = ''
    for alt in target_vcf_record.ALT:
        mut_base = alt.value.lower()
    wt_base = target_vcf_record.REF.lower()

    
    if mut_count + wt_count > 0:
        VAF = float(mut_count)/float(mut_count + wt_count) * 100.0
    else:
        VAF = 0.0

    plural_mut_count = plural_target_calls[mut_base]
    plural_wt_count = plural_target_calls[wt_base]

    if plural_mut_count + plural_wt_count > 0:
        plural_VAF = float(plural_mut_count)/float(plural_mut_count + plural_wt_count) * 100.0
    else:
        plural_VAF = 0.0

    #plural_mut_count = plural_target_calls[
        
    s += " Majority VAF: {:.2f}% ({}/{})\n".format(VAF, mut_count, (mut_count + wt_count))
    s += " Plurality VAF: {:.2f}% ({}/{})\n".format(plural_VAF, plural_mut_count, (plural_mut_count + plural_wt_count))
    s += " Plural target calls: {}\n".format(plural_target_calls)
    low, samp_mean, high = proportionConfidenceInterval(mut_count, wt_count, 0.95)
    s += "   95%    CI: {:.2f}% -- {:.2f}%\n".format(low * 100.0, high * 100.0)
    low, samp_mean, high = proportionConfidenceInterval(mut_count, wt_count, 0.99)
    s += "   99%    CI: {:.2f}% -- {:.2f}%\n".format(low * 100.0, high * 100.0)
    low, samp_mean, high = proportionConfidenceInterval(mut_count, wt_count, 0.999)
    s += "   99.9%  CI: {:.2f}% -- {:.2f}%\n".format(low * 100.0, high * 100.0)
    low, samp_mean, high = proportionConfidenceInterval(mut_count, wt_count, 0.9999)
    s += "   99.99% CI: {:.2f}% -- {:.2f}%\n".format(low * 100.0, high * 100.0)

    
    s += "   Read # when variant call statistically significant vs. {}% error: {}\n".format(err * 100.0, statistical_significance_map)


    s += " LAMP Concatemer Histogram:\n"
    longest_concatemer = max(num_targets_map.keys())
    for i in range(longest_concatemer + 1):
        if i in num_targets_map.keys():
            s += "   {} : {}\n".format(i, num_targets_map[i])
        else:
            s += "   {} : {}\n".format(i, 0)

    s += " Average Aligned Read Accuracy:\n"
    if total_calls > 0:
        s += "    {}/{} {}%\n".format(total_correct_calls, total_calls, float(total_correct_calls)/float(total_calls)*100.0)
    else:
        s += "    {}/{} {}%\n".format(total_correct_calls, total_calls, 'NA')

    s += " Error call distribution:\n"
    s += "    {}\n".format(all_bad_calls)
    s += " Polished error call distribution:\n"
    s += "    {}\n".format(all_polished_bad_calls)

    if args.print_summary_stats:
        print("\n" + s)

    # save summary stats
    with open("{}/summary_stats.txt".format(args.output_dir), "w") as sum_stat_fh:
        sum_stat_fh.write(s)


    # save read classifications, basecalls, and time stamps
    if args.save_read_classifications:
        read_classification_list = list()
        first = True
        first_timestamp = results[0].timestamp
        for result in results:

            seconds_elapsed = timestampDelta(first_timestamp, result.timestamp)
            read_classification_list.append((result.read_id, result.classification, result.plural_base, seconds_elapsed))


        #print(read_classification_list)
        class_fn = "{}/lamprey_read_classifications.csv".format(args.output_dir)

        with open(class_fn, 'w') as filehandle:
            filehandle.writelines("{},{},{},{}\n".format(read_id,classification,plural_base,seconds_elapsed) for read_id,classification,plural_base,seconds_elapsed in read_classification_list)
        
    # Write output files
    # polished lamplicons
    fasta_fn = "{}/polished.fasta".format(args.output_dir)
    #TODO: Nonetypes gett added to the polished list sometimes....
    if len(polished_lamplicons) > 0:
        with open(fasta_fn, "w") as output_handle:
            SeqIO.write(polished_lamplicons, output_handle, "fasta")


    
def argparser():
    parser = argparse.ArgumentParser()

    # input/output files
    parser.add_argument("target_ref_fn")
    parser.add_argument("human_ref_fn")
    parser.add_argument("primer_set_fn")
    parser.add_argument("lamplicons_fn")
    parser.add_argument("--output_dir", type=str, default="results")
    parser.add_argument("--target_vcf_fn", type=str, default='')
    parser.add_argument("--ref_vcf_fn", type=str, default='')

    # algorithm options
    parser.add_argument("--threshold", type=float, default=0.75,
            help="Primer identity threshold for successful alignment")
    parser.add_argument("--time_sort", action="store_true", default=False)
    parser.add_argument("--high_confidence", action="store_true", default=False)
    parser.add_argument("--swalign", action="store_true", default=False)
    
    # run options
    parser.add_argument("--num_threads", type=int, default=1)
    
    # print options
    parser.add_argument("--debug_print", action="store_true", default=False)
    parser.add_argument("--print_summary_stats", action="store_true", default=False)
    
    # output options
    parser.add_argument("--save_target_reads", action="store_true", default=False)
    parser.add_argument("--save_concatemers", action="store_true", default=False)
    parser.add_argument("--save_read_classifications", action="store_true", default=False)


    
    
    return parser



if __name__ == "__main__":
    parser = argparser()
    args = parser.parse_args()
    main(args)

