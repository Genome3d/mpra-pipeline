import pandas as pd
import os
import sys
import argparse


def write_fasta(mpra_lib, outdir):
    out_file = []
    df = pd.read_csv(mpra_lib, sep = '\t', header=None)
    for _,line in df.iterrows():
        ids = '> '+ line[0]
        out_file.append(ids)
        seq = line[1]
        out_file.append(seq)
    out_df = pd.DataFrame(out_file)
    out_df.to_csv(os.path.join(outdir, 'mpra_library_fasta.fa'), 
            sep = '\n', 
            header=False, 
            index=False)

def parse_args():
    parser = argparse.ArgumentParser(
        description= "Make your own fasta file.")
    parser.add_argument(
        '-l', '--mpra-lib', required=True, 
        help='''Text file containing two columns without header 
        (column1: sequence description/identifier,
        column2: sequence''')
    parser.add_argument(
        '-o', '--outdir', required=True,
        help='''Path to the output directory''')
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir, exist_ok=True)

    lib_fasta = write_fasta(args.mpra_lib, args.outdir)




