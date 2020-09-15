#!/usr/bin/pyAthon

import sys
import os
import subprocess

import swalign
import mappy
import pysam

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

#########################
#fwd_primer_order = ['F3', 'F2', 'FLPc', 'F1', 'B1c', 'BLP', 'T', 'B2c', 'B3c']
#rev_primer_order = ['B3', 'B2', 'Tc', 'BLPc', 'B1', 'F1c', 'FLP', 'F2c', 'F3c']
fwd_primer_order = ['F3', 'F2', 'F1', 'B1c', 'T', 'B2c', 'B3c']
rev_primer_order = ['B3', 'B2', 'Tc', 'B1', 'F1c', 'F2c', 'F3c']


#########################
class PrimerSet:
    
    def __init__(self, F3='', B3='', FIP='', BIP='', FLP='', BLP='', F1='', F2='', B1='', B2='',T=''):
        self.F3 = F3
        self.B3 = B3
        self.FIP = FIP
        self.BIP = BIP
        self.FLP = FLP
        self.BLP = BLP
        self.F1 = F1
        self.F2 = F2
        self.B1 = B1
        self.B2 = B2
        self.T = T
        
#########################
class Alignment:

    def __init__(self, start=0, end=0, identity=0.0, name=''):
        self.start = start
        self.end = end
        self.identity = identity
        self.name = name

    def __repr__(self):
        return "{} at ({}, {}) identity: {}".format(self.name, self.start, self.end, self.identity)
        
    def __str__(self):
        return "{} at ({}, {}) identity: {}".format(self.name, self.start, self.end, self.identity)

    def __lt__(self, other):
        # p1 < p2 calls p1.__lt__(p2)
        return self.start < other.start
    
    def __eq__(self, other):
        self.start == other.start
    
#########################
def getPrimersFromFile(fn):

    primerSet = PrimerSet()
    
    with open(fn) as fp:
        for line in fp:

            line = line.strip()
            if not line:
                continue
            fields = line.split(' ')
            primer_name = fields[0]
            primer_string = fields[1]
            
            if primer_name == "F3":
                primerSet.F3 = primer_string
            elif primer_name == "B3":
                primerSet.B3 = primer_string
            elif primer_name == "FIP":
                primerSet.FIP = primer_string
            elif primer_name == "BIP":
                primerSet.BIP = primer_string
            elif primer_name == "FLP":
                primerSet.FLP = primer_string
            elif primer_name == "BLP":
                primerSet.BLP = primer_string
            elif primer_name == "F1":
                primerSet.F1 = primer_string
            elif primer_name == "F1c":
                primerSet.F1 = rc(primer_string)
            elif primer_name == "F2":
                primerSet.F2 = primer_string
            elif primer_name == "B1":
                primerSet.B1 = primer_string
            elif primer_name == "B1c":
                primerSet.B1 = rc(primer_string)
            elif primer_name == "B2":
                primerSet.B2 = primer_string
            elif primer_name == "T":
                primerSet.T = primer_string
            elif primer_name == "Tc":
                primerSet.T = rc(primer_string)
                
    return primerSet

#######
def alignPrimer(aligner, seq, primer, primer_name):

    # 
    alignment_tmp = aligner.align(seq, primer)
    alignment = Alignment(alignment_tmp.r_pos,
                          alignment_tmp.r_end,
                          alignment_tmp.identity,
                          primer_name)

    return alignment

#######
def findPrimerAlignments(aligner, seq, primer, primer_name, identity_threshold):

    # recursively find the best alignment and then reserve that sequence

    # find alignment
    alignment = alignPrimer(aligner, seq, primer, primer_name)
    alignment_len = alignment.end - alignment.start
    
    #
    alignment_list = list()
    
    if alignment.identity > identity_threshold and alignment_len > len(primer) * identity_threshold:
        
        # add to interval list
        alignment_list.append(alignment)
        
        # split along the alignment and recurse
        #r_pos (ref (lamplicon) start pos) r_end (lamplicon end)
                
        # recurse left
        left_seq = seq[:alignment.start]
        if len(left_seq) >= len(primer):
            left_alignments = findPrimerAlignments(aligner,
                                                   left_seq,
                                                   primer,
                                                   primer_name,
                                                   identity_threshold) 

            alignment_list.extend(left_alignments)
        
        # recurse right
        #print("Recurse right: {}".format(seq[alignment.r_end:]))
        right_seq = seq[alignment.end:]
        if len(right_seq) >= len(primer):
            right_alignments = findPrimerAlignments(aligner,
                                                    seq[alignment.end:],
                                                    primer,
                                                    primer_name,
                                                    identity_threshold)
            # adjust right alignments
            for right_alignment in right_alignments:
                right_alignment.start = right_alignment.start + alignment.end
                right_alignment.end = right_alignment.end + alignment.end

            # add all right alignments
            alignment_list.extend(right_alignments)

    return alignment_list

#########################
# from stack overflow: https://stackoverflow.com/questions/25188968/reverse-complement-of-dna-strand-using-python
def rc(seq):

    alt_map = {'ins':'0'}
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

    
    for k,v in alt_map.items():
        seq = seq.replace(k,v)

    bases = list(seq)
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)

    for k,v in alt_map.items():
        bases = bases.replace(v,k)

    return bases

############################
def pruneRapidAdapters(seq):

    adapter_ytop='GGCGTCTGCTTGGGTGTTTAACCTTTTTTTTTTAATGTACTTCGTTCAGTTACGTATTGCT'
    adapter_ybot='GCAATACGTAACTGAACGAAGT'

    # prep aligner
    match = 2
    mismatch = -1
    scoring = swalign.NucleotideScoringMatrix(match, mismatch)

    sw = swalign.LocalAlignment(scoring)  # you can also choose gap penalties, etc...

    #
    orientation = 0
    
    # find adapter_ytop forward
    alignment = sw.align(adapter_ytop, read_string)
    print(alignment.identity)
    print(alignment.q_pos)
    print(alignment.q_end)

    # check ytop
    if alignment.q_end < 100 and alignment.identity > 0.55:
        print("Found possible adapter top of fwd:")
        alignment.dump()


    # check rc ytop
    alignment = sw.align(adapter_ytop, rc(read_string))
    if alignment.q_end < 100 and alignment.identity > 0.55:
        print("Found possible adapter top on rc:")
        alignment.dump()

    # check ybot
    alignment = sw.align(adapter_ybot, read_string)
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
    
    alignment_list.extend(findPrimerAlignments(sw, lamplicon, primers.F1, "F1", identity_threshold))
    alignment_list.extend(findPrimerAlignments(sw, lamplicon, rc(primers.F1), "F1c", identity_threshold))
    alignment_list.extend(findPrimerAlignments(sw, lamplicon, primers.F2, "F2", identity_threshold))
    alignment_list.extend(findPrimerAlignments(sw, lamplicon, rc(primers.F2), "F2c", identity_threshold))
    alignment_list.extend(findPrimerAlignments(sw, lamplicon, primers.F3, "F3", identity_threshold))
    alignment_list.extend(findPrimerAlignments(sw, lamplicon, rc(primers.F3), "F3c", identity_threshold))
    alignment_list.extend(findPrimerAlignments(sw, lamplicon, primers.B1, "B1", identity_threshold))
    alignment_list.extend(findPrimerAlignments(sw, lamplicon, rc(primers.B1), "B1c", identity_threshold))
    alignment_list.extend(findPrimerAlignments(sw, lamplicon, primers.B2, "B2", identity_threshold))
    alignment_list.extend(findPrimerAlignments(sw, lamplicon, rc(primers.B2), "B2c", identity_threshold))
    alignment_list.extend(findPrimerAlignments(sw, lamplicon, primers.B3, "B3", identity_threshold))
    alignment_list.extend(findPrimerAlignments(sw, lamplicon, rc(primers.B3), "B3c", identity_threshold))
    #alignment_list.extend(findPrimerAlignments(sw, lamplicon, primers.BLP, "BLP", identity_threshold))
    #alignment_list.extend(findPrimerAlignments(sw, lamplicon, rc(primers.BLP), "BLPc", identity_threshold))
    #alignment_list.extend(findPrimerAlignments(sw, lamplicon, primers.FLP, "FLP", identity_threshold))
    #alignment_list.extend(findPrimerAlignments(sw, lamplicon, rc(primers.FLP), "FLPc", identity_threshold))
    alignment_list.extend(findPrimerAlignments(sw, lamplicon, primers.T, "T", identity_threshold))
    alignment_list.extend(findPrimerAlignments(sw, lamplicon, rc(primers.T), "Tc", identity_threshold))

    return(sorted(alignment_list))
    
###################
def printPrimerAlignments(seq, alignments):

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


        # print the alignment if we're at the wrap factor OR if we've reached the end of the string
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
                        primer_string = alignment.name
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

                    # if no special char, print either a primer string, -, or space if not currently in a primer
                    if found_primer:
                        if primer_string_counter < len(primer_string):
                            sys.stdout.write(primer_string[primer_string_counter])
                            primer_string_counter = primer_string_counter + 1
                        else:
                            sys.stdout.write('-')
                    else:
                        sys.stdout.write(' ')

            # end alignment line
            sys.stdout.write('\n')

            # increment wrap counter
            wrap_counter = wrap_counter + 1

            
            
    sys.stdout.write('\n')

###################
def extractAmpliconAroundTarget(alignments, target):

    amplicon_start = 0
    amplicon_end = 0
    
    # is target forward or reverse?
    fwd_strand = True if target.name == "T" else False

    #print("Forward strand?: ",fwd_strand)

    # find target alignment index
    target_index = alignments.index(target)
    #print("Target index: ", target_index)

    #####
    # look to the right
    #####
    allowed_mismatches = 0
    #print("LOOKING RIGHT")
    primer_counter = fwd_primer_order.index('T') if fwd_strand else rev_primer_order.index('Tc')
    primer_counter = primer_counter + 1

    #print(primer_counter)
    
    all_matched = True
    for i in range(target_index + 1, len(alignments)):

        # primer name
        primer_name = alignments[i].name

        # expected primer name
        if fwd_strand:
            expected_primer_name = fwd_primer_order[primer_counter]
        else:
            expected_primer_name = rev_primer_order[primer_counter]

        #print("primer name: ", primer_name)
        #print("expected primer name: ", expected_primer_name)

        # on a mismatch, end the search
        if primer_name != expected_primer_name:

            if allowed_mismatches > 0:
                allowed_mismatches = allowed_mismatches - 1
            else:
                # we've found our end, point, so make our end the last match end
                amplicon_end = alignments[i-1].end
                #print("MISMATCH!")
                all_matched = False
                break
        else:
            #print("MATCH")
            primer_counter = primer_counter + 1
            # if we matched all primers, we have to quit no matter what
            if primer_counter == len(fwd_primer_order):
                break


    # if we matched all in the primer sequence, extend to the end of the entire read
    if all_matched:
        amplicon_end = -1

    #print("Trim end: ", amplicon_end)
        
    #####
    # look to the left
    #####
    allowed_mismatches = 0
    #print("LOOKING LEFT")
    primer_counter = fwd_primer_order.index('T') if fwd_strand else rev_primer_order.index('Tc')
    primer_counter = primer_counter - 1

    all_matched = True
    for i in range(target_index - 1, 0, -1):

        #print(i)
        
        # primer name
        primer_name = alignments[i].name

        # expected primer name
        if fwd_strand:
            expected_primer_name = fwd_primer_order[primer_counter]
        else:
            expected_primer_name = rev_primer_order[primer_counter]

        #print("primer name: ", primer_name)
        #print("expected primer name: ", expected_primer_name)

        # on a mismatch, end the search
        if primer_name != expected_primer_name:
            if allowed_mismatches > 0:
                allowed_mismatches = allowed_mismatches - 1
                
            else:
                # we've found our end, point, so make our end the last match end
                amplicon_start = alignments[i+1].start
                all_matched = False
                #print("MISMATCH!")
                break
        else:
            #print("MATCH")
            primer_counter = primer_counter - 1

    if all_matched:
        amplicon_start = 0

    #print("Trim start: ", amplicon_start)
        
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

def getAccuracy(read, alt_pos, alt_base):

    Eq = getEq(read)
    X = getX(read)
    I = getI(read)
    D = getD(read)

    # adjust Eq and X if this is a suspected alt read
    # find the alignment pair for the ref position
    # query the query pos that corresponds, and check if the base is the alt base
    for pair in read.get_aligned_pairs():
        if pair[0] == None or pair[1] == None:
            continue
        if pair[1] == (alt_pos - 1):
            print(pair)
            print(len(read.query_alignment_sequence))
            if read.query_alignment_sequence[pair[0]-read.query_alignment_start] == alt_base:
                print("TUMOR")
                Eq = Eq + 1
                X = X - 1
            else:
                print("NORMAL")
                
    return (float(Eq) / float(Eq + I + X + D))

#### parses CIGAR string and calc total errors.
# ignores leading/trailing errors
# can accept one alt base possibility/position pair to account for tumor reads
# if an alt base is provided, a consensus vote for the alt at that position will
# not be counted as an error
def calcErrorRate(ref_fn, sam_fn, alt_pos = -1, alt_base ='N'):

    samfile = pysam.AlignmentFile(sam_fn, 'r')

    cumulative_err_rate = 0.0
    cumulative_I_rate = 0.0
    cumulative_D_rate = 0.0
    cumulative_X_rate = 0.0
    mapped_read_count = 0.0
        
    for read in samfile:
        if not read.is_unmapped:
            mapped_read_count += 1.0
            #cigar = read.cigarstring
            #print(cigar)
            #print(stats[0])
            #print(getAccuracy(read))
            cumulative_err_rate += getAccuracy(read, alt_pos, alt_base)
            #cumulative_I_rate += getI(read) / (getI(read) + getEq(read) + getD(read) + getX(read))
            #cumulative_D_rate += getD(read) / (getI(read) + getEq(read) + getD(read) + getX(read))
            #cumulative_X_rate += getX(read) / (getI(read) + getEq(read) + getD(read) + getX(read))
            
            
    avg_err_rate = cumulative_err_rate / mapped_read_count
            
    #avg_I_rate = cumulative_I_rate / mapped_read_count
    #avg_D_rate = cumulative_D_rate / mapped_read_count
    #avg_X_rate = cumulative_X_rate / mapped_read_count

    #print(avg_err_rate)
    #print("Error {}% = {}%I {}%D {}%X".format(1.0-avg_err_rate, avg_I_rate, avg_D_rate, avg_X_rate))

    return avg_err_rate
    

################################################

# get file directory for relative script paths
dirname = os.path.dirname(__file__)
generate_consensus = os.path.join(dirname, 'scripts/generate_consensus.sh')
sam_to_bam = os.path.join(dirname, 'scripts/sam_to_bam.sh')


# parse args
args = sys.argv

# print usage
if len(sys.argv) != 5:
    print("usage: python3 lamp_aligner.py <target_ref.fa> <primer_set> <lamplicons.fq> <lamplicon_target>")
    sys.exit(1)


target_ref_fn = args[1]
primer_set_fn = args[2]
lamplicon_fn = args[3]
lamplicon_target = int(args[4])


# parse reads
#print("Parsing fastq file...")
lamplicons = SeqIO.parse(lamplicon_fn, "fastq")
#print("    Done.")

# parse lamp primers
#print("Parsing lamp primers...")
primers = getPrimersFromFile(primer_set_fn)
#print("    Done.")

# prep aligner
match = 2
mismatch = -1
scoring = swalign.NucleotideScoringMatrix(match, mismatch)
sw = swalign.LocalAlignment(scoring)  # you can also choose gap penalties, etc...


lamplicon = ''
lamplicon_counter = 0

alignment_identity_threshold = 0.75

single_lamplicon = True

for lampli in lamplicons:

    #if single_lamplicon and lamplicon_counter < lamplicon_target:
    #    lamplicon_counter = lamplicon_counter + 1
    #    continue
    
    #if lamplicon_counter > lamplicon_target:
    #    break
    
    lamplicon = lampli.seq

    alignments = findAllPrimerAlignments(sw, lamplicon, primers, alignment_identity_threshold)

    alignments = sorted(alignments)

    printPrimerAlignments(lamplicon, alignments)

    target_counter = 0
    targets = []
    buf = 60

    # gather targets into separate SeqIO Records
    for alignment in alignments:
        if alignment.name == 'T' or alignment.name == 'Tc':
            target_counter = target_counter + 1

            
            # reach left and right to find concatemer cut points

            # Naive implementation where we just consider target + buf and target - buf region
            # didn't work nearly as well in practice as extending the chain to encompass as much of a full
            # amplicon as possible
            #seq_start = alignment.start - buf if (alignment.start - buf) > 0 else 0
            #seq_end = alignment.end + buf if (alignment.end + buf) < len(lamplicon) else len(lamplicon) - 1

            # primer token chain extension
            seq_start, seq_end = extractAmpliconAroundTarget(alignments, alignment)

            # if we get a -1 for seq end, it means we never matched to the right, so just grab the whole read
            if seq_end == -1:
                seq_end = len(lamplicon) - 1
            
            # generate new read from lamplicon substring
            record = SeqRecord(
                lamplicon[seq_start:seq_end],
                "{}_{}".format(lamplicon_counter, target_counter),
                "{}_{}".format(lamplicon_counter, target_counter),
                "Lamplicon {} target {}".format(lamplicon_counter, target_counter)
            )

            #
            print("  {}:{}".format(seq_start, seq_end))
            
            #print(record)
            targets.append(record)
            
    print("Found Targets: {}".format(target_counter))
    
    # write targets to new fasta file
    print("Writing {} target sequences to {}_{}.fasta".format(target_counter, lamplicon_counter, target_counter))
    fasta_fn = "{}_{}.fasta".format(lamplicon_counter, target_counter)
    with open(fasta_fn, "w") as output_handle:
        SeqIO.write(targets, output_handle, "fasta")


    print("Generating alignments and consensus sequence...")
    # align targets to ref
    ref_fn = target_ref_fn

    out_sam_temp = "{}_{}_temp.sam".format(lamplicon_counter,target_counter)
    out_sam = "{}_{}.sam".format(lamplicon_counter,target_counter)
    in_bam = "{}_{}.bam".format(lamplicon_counter,target_counter)
    
    subprocess.run(["minimap2","-a","-xmap-ont","--eqx","-t 1","-w 1","-n 1", "-N 0", "-m 0","-p 0.6","-s 30", "-o{}".format(out_sam_temp),ref_fn,fasta_fn], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    # filter secondary/supplementary alignments from sam file
    # samtools view -h -F 0x900 filename.bam
    subprocess.run(["samtools", "view", "-h", "-F 0x900", out_sam_temp, "-o{}".format(out_sam)])
    
    # convert sam to bam
    subprocess.run([sam_to_bam, out_sam], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    # generate consensus sequence if we had at least one valid alignment
    mapped_reads = pysam.AlignmentFile(in_bam, 'rb').mapped
    print("Mapped {} of {} targets".format(mapped_reads, target_counter))
    if mapped_reads <= 1:
        # consider next lamplicon candidate
        lamplicon_counter = lamplicon_counter + 1
        continue

    subprocess.run([generate_consensus, ref_fn, in_bam], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    fastq_cns_fn = "{}_{}_cns.fastq".format(lamplicon_counter, target_counter)
    out_cns_sam = "{}_{}_cns.sam".format(lamplicon_counter, target_counter)

    # remove leading/trailing ambiguous characters from consensus sequence
    for consensus_seq in SeqIO.parse(fastq_cns_fn, "fastq"):

        # find location of last n from the start
        n_end = 0
        while consensus_seq.seq[n_end] == 'n' or consensus_seq.seq[n_end] == 'N':
            n_end = n_end + 1

        n_start = len(consensus_seq.seq) - 1
        while consensus_seq.seq[n_start] == 'n' or consensus_seq.seq[n_start] == 'N':
            n_start = n_start - 1

        consensus_seq = consensus_seq[n_end:n_start+1]

    # re-write consensus fastq
    with open(fastq_cns_fn, "w") as output_handle:
        SeqIO.write(consensus_seq, output_handle, "fastq")
    
    # align consensus sequence
    subprocess.run(["minimap2","-a","-xmap-ont","--eqx","-t 1","-w 1", "-o{}".format(out_cns_sam),ref_fn,fastq_cns_fn], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    # convert sam to bam
    subprocess.run([sam_to_bam, out_cns_sam], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    
    print("  DONE")

    # computing error from each alignment (ignoring leading/trailing clips and accounting for mutations
    orig_err = calcErrorRate(ref_fn, out_sam, 157, 'T')
    cns_err = calcErrorRate(ref_fn, out_cns_sam, 157, 'T')

    print("Orig Acc: {}%".format(orig_err*100))
    print("Cons Acc: {}%".format(cns_err*100))
    
    # consider next lamplicon candidate
    lamplicon_counter = lamplicon_counter + 1

