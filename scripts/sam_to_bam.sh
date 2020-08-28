#!/bin/bash

fn=${1}

# validate input is .sam format
if [[ ${fn} == *.sam ]]
then

    key=${fn%.sam}
    echo ${key}    
    
    echo "Converting Sam to Bam...."
    samtools view -Sb ${fn} > ${key}.unsorted.bam
    
    echo "Sorting Bam file...."
    samtools sort ${key}.unsorted.bam -o ${key}.bam
    
    echo "Indexing sorted Bam...."
    samtools index ${key}.bam

    echo "Cleaning up..."
    rm ${key}.unsorted.bam

    echo "Done!"
    
fi
