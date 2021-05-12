#!/bin/bash

fn=${1}

# validate input is .sam format
if [[ ${fn} == *.sam ]]
then

    key=${fn::-4}
    
    
    echo "Converting Sam to Bam...."
    # filters secondary 0x100 and supplementary 0x800
    samtools view -F 0x900 -Sb ${fn} > ${key}.unsorted.bam
    
    echo "Sorting Bam file...."
    samtools sort ${key}.unsorted.bam -o ${key}.bam
    
    echo "Indexing sorted Bam...."
    samtools index ${key}.bam

    echo "Cleaning up..."
    rm ${key}.unsorted.bam

    echo "Done!"
    
    #echo "Generating pileup..."
    #samtools mpileup \
#	     -Q 0 \
#	     -f ~/wadden/datasets/human_reference/human_g1k_v37.fasta \
#	     ./${key}.sorted.bam \
#	     > pilup.txt
    
fi
