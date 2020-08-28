#!/usr/bin/python

import sys

import pysam

#accuracy = counts['='] / (counts['='] + counts['I'] + counts['X'] + counts['D'])

def getM(cigar_stats):
    return cigar_stats[0][0]

def getI(cigar_stats):
    return cigar_stats[0][1]

def getD(cigar_stats):
    return cigar_stats[0][2]

def getEq(cigar_stats):
    return cigar_stats[0][7]

def getX(cigar_stats):
    return cigar_stats[0][8]

def getAccuracy(cigar_stats):

    Eq = getEq(cigar_stats)
    I = getI(cigar_stats)
    D = getD(cigar_stats)
    X = getX(cigar_stats)
    
    return (float(Eq) / float(Eq + I + X + D))

    
    
##########################
args = sys.argv

samfile = args[1]

# M BAM_CMATCH 0
# I BAM_CINS 1
# D BAM_CDEL 2
# N BAM_CREF_SKIP 3
# S BAM_CSOFT_CLIP 4
# H BAM_CHARD_CLIP 5
# P BAM_CPAD 6
# = BAM_CEQUAL 7
# X BAM_CDIFF 8
# B BAM_CBACK 9

samfile = pysam.AlignmentFile(samfile, 'r')

cumulative_err_rate = 0.0
cumulative_I_rate = 0.0
cumulative_D_rate = 0.0
cumulative_X_rate = 0.0
mapped_read_count = 0.0


for read in samfile:
    if not read.is_unmapped:
        mapped_read_count += 1.0
        cigar = read.cigarstring
        #print(cigar)
        stats = read.get_cigar_stats()
        #print(stats[0])
        #print(getAccuracy(stats))
        cumulative_err_rate += getAccuracy(stats)
        cumulative_I_rate += getI(stats) / (getI(stats) + getEq(stats) + getD(stats) + getX(stats))
        cumulative_D_rate += getD(stats) / (getI(stats) + getEq(stats) + getD(stats) + getX(stats))
        cumulative_X_rate += getX(stats) / (getI(stats) + getEq(stats) + getD(stats) + getX(stats))


avg_err_rate = cumulative_err_rate / mapped_read_count

avg_I_rate = cumulative_I_rate / mapped_read_count
avg_D_rate = cumulative_D_rate / mapped_read_count
avg_X_rate = cumulative_X_rate / mapped_read_count



print(avg_err_rate)
print("Error {}% = {}%I {}%D {}%X".format(1.0-avg_err_rate, avg_I_rate, avg_D_rate, avg_X_rate))
