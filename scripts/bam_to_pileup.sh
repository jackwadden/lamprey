#!/bin/bash

REF=${1}
BAM=${2}

key=${BAM%.bam}

echo "Generating pileup..."
samtools mpileup \
         -Q 0 \
         -f ${REF} \
         ${BAM} \
         > ${key}.pileup
