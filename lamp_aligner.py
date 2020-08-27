#!/usr/bin/pyAthon

import sys


import swalign
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

#########################
fwd_primer_order = ['F3', 'F2', 'FLPc', 'F1', 'B1c', 'BLP', 'T', 'B2c', 'B3c']
rev_primer_order = ['B3', 'B2', 'Tc', 'BLPc', 'B1', 'F1c', 'FLP', 'F2c', 'F3c']


#########################
class PrimerSet:
    
    def __init__(self, F3='', B3='', FIP='', BIP='', FLP='', BLP='', F1='', F2='', B1='', B2=''):
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
    alignment_list.extend(findPrimerAlignments(sw, lamplicon, primers.BLP, "BLP", identity_threshold))
    alignment_list.extend(findPrimerAlignments(sw, lamplicon, rc(primers.BLP), "BLPc", identity_threshold))
    alignment_list.extend(findPrimerAlignments(sw, lamplicon, primers.FLP, "FLP", identity_threshold))
    alignment_list.extend(findPrimerAlignments(sw, lamplicon, rc(primers.FLP), "FLPc", identity_threshold))
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

    print("Forward strand?: ",fwd_strand)

    # find target alignment index
    target_index = alignments.index(target)
    print("Target index: ", target_index)

    #####
    # look to the right
    #####
    allowed_mismatches = 1
    #print("LOOKING RIGHT")
    primer_counter = fwd_primer_order.index('T') if fwd_strand else rev_primer_order.index('Tc')
    primer_counter = primer_counter + 1

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


    if all_matched:
        amplicon_end = alignments[-1].end

    print("Trim end: ", amplicon_end)
        
    #####
    # look to the left
    #####
    allowed_mismatches = 1
    #print("LOOKING LEFT")
    primer_counter = fwd_primer_order.index('T') if fwd_strand else rev_primer_order.index('Tc')
    primer_counter = primer_counter - 1

    all_matched = True
    for i in range(target_index - 1, 0, -1):

        print(i)
        
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

    print("Trim start: ", amplicon_start)
        
    return amplicon_start, amplicon_end
    
################################################

# parse args
args = sys.argv
lamplicon_target = int(args[1])

# parse reads
#print("Parsing fastq file...")
lamplicons = SeqIO.parse("data/lamp_H33H31mux/UMPED18B_H33.fastq", "fastq")
#print("    Done.")

# parse lamp primers
#print("Parsing lamp primers...")
primers = getPrimersFromFile("data/primers/H3.3_set1.primers")
#print("    Done.")

# prep aligner
match = 2
mismatch = -1
scoring = swalign.NucleotideScoringMatrix(match, mismatch)
sw = swalign.LocalAlignment(scoring)  # you can also choose gap penalties, etc...


lamplicon = ''
lamplicon_counter = 0

alignment_identity_threshold = 0.85

single_lamplicon = True

for lampli in lamplicons:

    if single_lamplicon and lamplicon_counter < lamplicon_target:
        lamplicon_counter = lamplicon_counter + 1
        continue
    
    if lamplicon_counter > lamplicon_target:
        break
    
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
            seq_start = alignment.start - buf if (alignment.start - buf) > 0 else 0
            seq_end = alignment.end + buf if (alignment.end + buf) < len(lamplicon) else len(lamplicon) - 1

            seq_start, seq_end = extractAmpliconAroundTarget(alignments, alignment)
            
            record = SeqRecord(
                lamplicon[seq_start:seq_end],
                "{}_{}".format(lamplicon_counter, target_counter),
                "{}_{}".format(lamplicon_counter, target_counter),
                "Lamplicon {} target {}".format(lamplicon_counter, target_counter)
            )
            #print(record)
            targets.append(record)
            
    print("NUM TARGETS: {}".format(target_counter))

    # write targets to new fasta file
    fn = "alignments/{}_{}.fasta".format(lamplicon_counter, target_counter)
    with open(fn, "w") as output_handle:
        SeqIO.write(targets, output_handle, "fasta")
    
    lamplicon_counter = lamplicon_counter + 1

