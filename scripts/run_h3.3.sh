#!/bin/bash

if [ $# -eq 0 ]
then
    echo "USAGE: run_h3.3.sh <FASTQ>"
    exit
fi

# 	--high_confidence \
#	--debug_print \

REF=${1}
FASTQ=${2}


python3 ../lamprey.py \
     --time_sort \
     --save_target_reads \
     --save_concatemers \
     --save_read_classifications \
     --print_summary_stats \
     --num_threads 2 \
     --threshold 0.75 \
     --ref_vcf_fn ../data/H33/H3.3k27m_ref.vcf \
     --target_vcf_fn ../data/H33/H3.3k27m_target.vcf \
     ../data/references/H3F3A.fa \
     ${REF} \
     ~/lamprey/data/H33/H3.3_set1.primers \
     ${FASTQ}


