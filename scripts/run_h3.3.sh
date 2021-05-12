#!/bin/bash

if [ $# -eq 0 ]
then
    echo "USAGE: run_h3.3.sh <FASTQ>"
    exit
fi

REF=${1}
FASTQ=${2}

full_path=$(realpath $0)
dir_path=$(dirname $full_path)

python3 ${dir_path}/../lamprey.py \
	--debug_print \
	--time_sort \
	--save_target_reads \
	--save_concatemers \
	--save_read_classifications \
	--print_summary_stats \
	--num_threads 1 \
	--threshold 0.75 \
	--ref_vcf_fn ${dir_path}/../data/H33/H3.3k27m_ref.vcf \
	--target_vcf_fn ${dir_path}/../data/H33/H3.3k27m_target.vcf \
	${dir_path}/../data/references/H3F3A.fa \
	${REF} \
	${dir_path}/../data/H33/H3.3_set1.primers \
	${FASTQ}


