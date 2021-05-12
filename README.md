# LAMPrey

LAMPrey is a Loop-Mediated Isothermap Amplification (LAMP) concatemer analysis, variant calling, and polishing tool. LAMPrey follows a seed, chain, align paradigm to improve recovery of relevent information from lamplicons and aid in assay debugging.

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

LAMPrey currently relies heavily on the ```swalign``` library. This library is a pure python implementation of the Smith-Waterman local alignment algorithm and so is.... not fast. LAMPrey hugely benefits from a few loop optimizations and cython compilation of this package. An upgraded version of ```swalign``` is located at [jackwadden/swalign](https://www.github.com/jackwadden/swalign). To compile and install the accelerated version of swalign, clone the repo, run the build script, and move the resulting .so file to LAMPrey's lib/ directory.
```
(lamprey) $ git clone https://github.com/jackwadden/swalign.git
(lamprey) $ cd swalign
(lamprey) $ python3 setup.py build_ext --inplace
(lamprey) $ mkdir path-to-lamprey/lib
(lamprey) $ cp swalign.cpython#####.so path-to-lamprey/lib
```

## Quick Start

An example run script is provided to help users get started and play with various options.
```
$ scripts/run_h3.3.sh <reference_minimap_idx> <fastq>
```
## General LAMPrey Usage

### Reference/VCF Files
LAMPrey uses alignment to a local gene target as well as alignment to the human reference to diagnose LAMP assay performance. Therefore, you will need to supply four files (on top of a FASTQ input and LAMP schema file) for LAMPrey to run:
1. A target reference (target_ref_fn). this is a small FASTA sequence representing the gene or genomic region the LAMP assay targets.
2. A VCF file outlining the mutation according to the target reference (--target_vcf_fn).
3. A genomic reference (human_ref_fn). This is a genomic reference (or ideally minimap2 index) of the host. The reference used for the examples should be HG19.
4. A VCF file outlining the mutation according to the genomic reference (--ref_vcf_fn).

1/2 are technically not required and will be removed in the future in favor of an automatic extraction of the genomic region based on the ref_vcf mutation.

### LAMP Schema File
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

