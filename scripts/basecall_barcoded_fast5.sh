#!/bin/bash

#    --barcode_kits "EXP-NBD104 EXP-NBD114" \

#DATA_ROOT=.
FAST5_SRC=${1}
FASTQ_DST=${2}
INPUT_FILE_LIST=${3}
OUTPUT_SAVE_FILE=${4}


guppy_basecaller \
    -i ${FAST5_SRC} \
    --chunks_per_runner 40 \
    -s ${FASTQ_DST} \
    --input_file_list ${INPUT_FILE_LIST} \
    --disable_qscore_filtering \
    --quiet \
    --trim_barcodes \
    --kit SQK-RBK004 \
    --barcode_kits "SQK-RBK004" \
    --flowcell FLO-MIN106 \
    --num_callers 4 \
    -x "cuda:0"

cp -r ${FASTQ_DST}/ ./fastq_saved
