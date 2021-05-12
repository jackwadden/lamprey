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
	--ref_vcf_fn ${dir_path}/../data/H31/H3.1k27m_ref.vcf \
	--target_vcf_fn ${dir_path}/../data/H31/H3.1k27m_target.vcf \
	${dir_path}/../data/references/HIST1H3B.fa \
	${REF} \
	${dir_path}/../data/H31/H3.1_set1.primers \
	${FASTQ}


