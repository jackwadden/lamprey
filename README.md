# lamprey

A LAMP concatemer (lamplicon) polisher to improve single read accuracy.

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

## Usage

### Breaking up lamplicons into separate reads
The following will convert all lamplicon sequences to files of separated concatemer sequences.
```
$ python3 lamprey.py primer_set.primers lamplicons.fastq
```

### Aligning separate reads to the target reference
The following will align each separate mer into out.sam.

```
$ minimap2 -a \
     -x map-ont \
     --eqx \
     -t 1 \
     -w 1 \
     -o mers.sam \
     ref.fa \
     mers.fa
```

### Converting SAM to sorted BAM
The following will convert the aligned mers to a sorted BAM file.
```
$ scripts/sam_to_bam.sh mers.sam
```

### Generating a consensus from aligned mers
The following will generate a consensus sequence (out_cns.fa) from the aligned mers.
```
$ scripts/generate_consensus.sh ref.fa out.bam
```

### Aligning consensus sequence
The following will align the consensus sequence to a reference.
```
$ minimap2 -a \
     -x map-ont \
     --eqx \
     -t 1 \
     -w 1 \
     -o mers_cns.sam \
     ref.fa \
     mers_cns.fa
```

### Measuring Error
The following will compute the average error of the original aligned mers and the polished concensus sequence.
```
$ python3 scripts/error_rate.py mers.sam
$ python3 scripts/error_rate.py mers_cns.sam
```

## Example Primer Set File
```
F3 ATGTCGTGCCGTTG
F2 ATGCGTATAGATGTATAGA
F1c AAAATTTTGGGGCCCC
FLP ATGTTTGGTTGAGATG
T ATATATATATATATAAT
```
