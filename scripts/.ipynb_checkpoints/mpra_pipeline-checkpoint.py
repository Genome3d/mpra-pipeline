import pandas as pd
import os
import sys
import argparse
import gzip
import multiprocessing as mp
import multiprocessing.pool as mpp
import pyfastx
from Bio import SeqIO

def parse_fastq(read1):
    rd1_ls = ()
    with gzip.open(read1, 'rt') as rd1:
        for record in rd1:
            name = list(record.rstrip().split('\n'))[0]
            #print(name)
            #seq = list(record.split('\n'))
            print(name)
            print(name[1])
      
        
    
def parse_args():
    parser = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description= '''MPRA sequencing data analysis pipeline.''')
    parser.add_argument(
            '-r1', '--read1', required=True, 
            help='''Read 1 fastq file. Takes single fastq file.''')
    parser.add_argument(
            '-r2', '--read2', required=True,
            help='''Read 2 fastq file. Takes single fastq file.''')
    parser.add_argument(
            '-l', '--library', required=True,
            help='''The 200bp sequence library. Should also contain Unique IDs''')
    parser.add_argument(
            '-o', '--outdir', required=True,
            help='''Provide a path to the output directory''')
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()

    # validate arguments

    if not (args.read1 and args.read2):
        sys.exit('''Missing Read1 and Read2 files.
                See 'mpra_pipeline.py -h' for more details''')
    if not args.library:
        sys.exit('''Missing reference library file.
                See 'mpra_pipeline.py -h' for more details''')
    if not args.outdir:
        sys.exit('''Missing argument: output directory.
                See 'mpra_pipeline.py -h' for more details''')
    
    # Currently only accept gzipped fastq files 
    if not (args.read1.endswith('.gz') and args.read2.endswith('.gz')):
        sys.exit('''Input file format not recognised.
                Make sure to provide gzipped fastq reads.''')  
    else:
        read1 = args.read1
        read2 = args.read2
    
    parse_fastq(read1)
        