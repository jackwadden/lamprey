#!/bin/bash

if [ $# -ne 2 ]
  then
      echo "./generate_consensus.sh <reference.fa> <bamfile.bam>"
      exit
fi

REF=${1}
BAM=${2}

KEY=${BAM%.bam}

echo $KEY

bcftools mpileup -f ${REF} ${BAM} | \
    bcftools call -c --ploidy 1 | \
    vcfutils.pl vcf2fq > \
    ${KEY}_cns.fastq \
