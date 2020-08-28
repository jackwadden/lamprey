# lamprey

A LAMP concatemer (lamplicon) polisher to improve single read accuracy.

## Installation

## Usage

### Breaking up lamplicons into separate reads
The following will convert all lamplicon sequences to files of separated concatemer sequences.
```
$ python3 lamprey.py primer_set.primers lamplicons.fastq
```

### Aligning separate reads to the target reference
```
$ minimap2 -a \
     -x map-ont \
     --eqx \
     -t 1 \
     -w 1 \
     -o out.sam \
     ref.fa \
     mers.fa
```

### Converting SAM to sorted BAM
```
$ scripts/sam_to_bam.sh out.sam
```

### Generating a consensus from aligned mers
```
$ scripts/generate_consensus.sh ref.fa out.bam
```

### Aligning consensus sequence
```
$ minimap2 -a \
     -x map-ont \
     --eqx \
     -t 1 \
     -w 1 \
     -o out.sam \
     ref.fa \
     mers_cns.fa
```

## Example Primer Set File
```
F3 ATGTCGTGCCGTTG
F2 ATGCGTATAGATGTATAGA
F1c AAAATTTTGGGGCCCC
FLP ATGTTTGGTTGAGATG
T ATATATATATATATAAT
```
