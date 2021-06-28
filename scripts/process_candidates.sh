#!/bin/bash

FASTA=${1}
SAM_FN=${2}
REF=${3}

# align each extracted concatemer piece
#minimap2 \
#    -a \
#    -xmap-ont \
#    --eqx \
#    -t 1 \
#    -w 1 \
#    -n 1 \
#    -N 5 \
#    --secondary=no \
#    -m 0 \
#3    -p 0.6 \
#    -s 30 \
#    -o ${SAM_FN} \
#    ${REF} \
#    ${FASTA}

KEY=${SAM_FN%_all.sam}

# filter, sort, and generate a pileup of the subs; piping for efficiency
full_path=$(realpath $0)
dir_path=$(dirname $full_path)
${dir_path}/sam_to_bam.sh ${SAM_FN}
samtools mpileup \
         -Q 0 \
         -f ${REF} \
	 ${KEY}_all.bam \
         > ${KEY}.pileup

#samtools view -Sb -h -F 0x900 ${SAM_FN} | \
#samtools sort | \
#samtools mpileup \
#         -Q 0 \
#         -f ${REF} \
#	 - \
#         > ${KEY}.pileup
