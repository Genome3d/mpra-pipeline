Massively parallel reporter assay (MPRA) data analysis pipeline
===============================================================

This repository contain scripts used to analyse MPRA sequencing 
data
     
Requirements
============

* python >3.0
* Biopython
* samtools
* Burrows-Wheeler Aligner (BWA)
* cutadapt
* seqtk

Usage
=====

```
python mpra_pipeline.py -h

usage: mpra_pipeline.py [-h] [-r1 READ1] [-r2 READ2] [-ds1 DNA_SEQ1 [DNA_SEQ1 ...]] [-ds2 DNA_SEQ2 [DNA_SEQ2 ...]] -o OUTDIR [-i INDEX_DIR] -m {enhancer,cdna,gdna}

MPRA sequencing data analysis pipeline.

optional arguments:
  -h, --help            show this help message and exit
  -r1 READ1, --read1 READ1
                        Path to a gzipped fastq file. Sequences in this file will be mapped against the reference MPRA library to find refSNP containing reads. Required if run --mode='enhancer'.
  -r2 READ2, --read2 READ2
                        Path to a gzipped fastq file. Sequences in this file will be used to search for barcodes of refSNPs. Required if run --mode='enhancer'.
  -ds1 DNA_SEQ1 [DNA_SEQ1 ...], --dna-seq1 DNA_SEQ1 [DNA_SEQ1 ...]
                        Provide a fastq file or a list of space separated fastq files (forward reads) containing 1) cDNA sequences (required if run --mode='cdna') 2) gDNA sequences (required if run
                        --mode='gdna').
  -ds2 DNA_SEQ2 [DNA_SEQ2 ...], --dna-seq2 DNA_SEQ2 [DNA_SEQ2 ...]
                        Provide a fastq file or a list of space separated fastq files (reverse reads) containing 1) cDNA sequences (required if run --mode='cdna') or 2) gDNA sequences (required if run
                        --mode='gdna').
  -o OUTDIR, --outdir OUTDIR
                        Provide a path to save output files.
  -i INDEX_DIR, --index-dir INDEX_DIR
                        Provide path to the directory where reference index (generated using bwa) files are located and also include prefix of the file
  -m {enhancer,cdna,gdna}, --mode {enhancer,cdna,gdna}
                        Choose a run mode to find the barcodes from enhancer, cDNA or gDNA sequences

```

