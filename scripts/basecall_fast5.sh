#!/bin/bash

#    --barcode_kits "EXP-NBD104 EXP-NBD114" \

#DATA_ROOT=.
FAST5_SRC=${1}
FASTQ_DST=${2}
INPUT_FILE_LIST=${3}


guppy_basecaller \
    -i ${FAST5_SRC} \
    --chunks_per_runner 70 \
    -s ${FASTQ_DST} \
    --input_file_list ${INPUT_FILE_LIST} \
    --disable_qscore_filtering \
    --kit SQK-RAD004 \
    --flowcell FLO-MIN106 \
    --num_callers 4 \
    -x "cuda:0"
