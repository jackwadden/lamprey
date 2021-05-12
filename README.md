# lamprey

A LAMP concatemer (lamplicon) analysis, variant calling, and polishing tool. LAMPrey follows a seed, chain, align paradigm to improve recovery of relevent information from lamplicons and aid in assay debugging.

## Installation
First, install `samtools-1.10`, `bcftools`, and `minimap2`.

Next, install `lamprey`:
```bash
$ git clone https://github.com/jackwadden/lamprey
$ cd lamprey
$ python3.7 -m venv venv3 --prompt lamprey
$ source venv3/bin/activate
(lamprey) $ pip install --upgrade pip
(lamprey) $ pip install -r requirements.txt
```

## Quick Start

An example run script is provided to help users get started and play with various options.
```
$ data/H33/run_h3.3.sh <reference_minimap_idx> <fastq>
```

## Example LAMP Schema File
LAMP schema files outline the ordering of expected LAMP sequences and also the sequences themselves. Sequences can be listed as forward or reverse complemented. Reverse complemented sequence names should end in a lower-case 'c'. The target sequence (```T```) is a >15bp region that contains the target information of interest to the assay. In our example, this target sequence covers the H3F3A p.K27M mutation.

The expected ordering of LAMP sequences should be put into both ```fwd_primer_order``` and ```rev_primer_order``` entries. Each entry lists the LAMP sequences as expected from a properly formed concatemer and the location of the target sequence (which usually exists between the F2/B2 sequences).

LAMPrey can identify arbitrary sequences of interest and is specifically designed to identify Oxford Nanopore Technologies fragmentation chemistry sequences, adapters, and barcodes. These sequences should start with lower case ```ont```.

An example schema file containing a real primer set is shown below. A full example is located at ```data/H33/H3.3_set1.primers```.

```
fwd_primer_order F3 F2 F1 B1c T B2c B3c
rev_primer_order B3 B2 Tc B1 F1c F2c F3c
F3 GTTTGGTAGTTGCATATGGTG
B3 ATACCTGTAACGATGAGGTTTC
F1c GCGGGCAGTCTGCTTTGTA
F2 ATGCTGGTAGGTAAGTAAGGA
B1c CGACCGGTGGTAAAGCACC
B2 CACCCCTCCAGTAGAG
# FLP CGAGCCATGGTACAGAGAC
# BLP CAGGAAGCAACTGGCTACA
T GCTCGCAAGAGTGCG
ontFRA3p GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA
ontFRA5p GCTTGGGTGTTTAACC
```

