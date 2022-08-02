import pandas as pd
import numpy as np
import os
import sys
import argparse
import gzip
import itertools
import subprocess as sp
import time
import datetime
import pysam
import tqdm
from Bio.Seq import Seq

def parse_fastq(fastq):
    processed = []
    col = ['name', 'seq', 'optional', 'qual']
    with gzip.open(fastq, 'rt') as fh:
        fh_iterator = (l.rstrip() for l in fh)
        for read in itertools.zip_longest(*[fh_iterator] * 4):
            read_dict = {key: val for key, val in zip(col, read)}
            keys_to_keep = ('name', 'seq')
            sub_reads_dict = {x: read_dict[x] for x in keys_to_keep
                    if x in read_dict}
            processed.append(sub_reads_dict)
    fastq_df = pd.DataFrame.from_dict(processed)
    fastq_df['name'] = fastq_df['name'].str.split(' ').str[0].str.strip('@').str.strip()
    return fastq_df

def find_barcodes(df0, outdir):
    logger.write("\tFinding barcodes...")
    extract = 'CGCCGAGGCCCGACGCTCTTCCGATCT(.*?)TCTAGAGGTACCGCAGGAGCCGCAGTG'
    df0['barcode'] = df0.seq.astype(str).str.extract(r'{}'.format(extract))
    df0 = df0.dropna()
    df = df0.copy()
    df['barcode_len'] = df['barcode'].str.len() 
    df.to_csv(os.path.join(outdir, "barcodes.txt"), 
            sep="\t", 
            header=True, 
            index=False)
    df_20 = df[df['barcode_len'] == 20]
    df_20.to_csv(os.path.join(outdir, "barcodes_len20.txt"), 
            sep="\t", 
            header=True, 
            index=False)
    return

def run_shell_commands(cmd):
    try:
        shell_cmd = sp.check_output(cmd, 
                            shell=True, 
                            stderr=sp.STDOUT)
    except sp.CalledProcessError as exception:
        logger.write ("\tShell command exited with an error.")
        logger.write(exception.output)
        return None

def run_shell_pipe(cmd):
    try:
        sp.Popen(cmd, shell=True,
                stdin=sp.PIPE,
                stdout=sp.PIPE,
                stderr=sp.PIPE).wait()
    except sp.CalledProcessError as exception:
        logger.write ("\tShell command exited with an error.")
        logger.write(exception.output)
        return None

def find_proper_alignments(bamfile, outdir):
    save = pysam.set_verbosity(0) # to suppress index not found warning
    bam_fh = pysam.AlignmentFile(bamfile, "rb")
    pysam.set_verbosity(save)
    cigar_matched = pysam.AlignmentFile(os.path.join(outdir, 
                                        "matched_reads.bam"), 
                                        "wb",
                                        template = bam_fh)
    alt_reads_file = pysam.AlignmentFile(os.path.join(outdir, 
                                        "indel_reads.bam"), 
                                        "wb",
                                        template = bam_fh)
    cigar_not_complete_match = pysam.AlignmentFile(os.path.join(outdir,
                                                   "not_complete_match_reads.bam"),
                                                   "wb",
                                                   template = bam_fh)

    read_ref_ls = []
    reads = []
    logger.write("\tFinding sequences that map to ref library...")
    for read in bam_fh:
        if read.cigarstring != None:
            cigarline = read.cigar
            tags = read.tags
            cigarstr = read.cigarstring
            name=read.qname
            len_on_ref = read.alen
            start = read.reference_start
            end = read.reference_end
            flag = read.flag
            rname = bam_fh.get_reference_name(read.tid)
            uniq_id = '_'.join([rname,str(start),str(end),str(len_on_ref)])
            # get cigar codes '0'= match/mismatch
            c_type = [1,2,3,4,5,6,8]
            codes = [cigar[0] for _, cigar in enumerate(cigarline) 
                    if cigar[0] in c_type]
            #reads perfectly mapped onto ref (cogar code = '0' and edit distance = 0)
            if not codes:
                for _,tagtype in enumerate(tags):
                    if tagtype[0] == "NM" and tagtype[1] == 0:
                        cigar_matched.write(read)
                        read_ref_ls.append([name, rname])
                        reads.append(name)
                    elif tagtype[0] == "NM" and tagtype[1] > 0:
                        cigar_not_complete_match.write(read)
                        continue
            else:
                alt_reads_file.write(read)
    read_ref_out = pd.DataFrame(read_ref_ls)
    read_ref_out.to_csv(os.path.join(outdir, "reads_refs_ids.txt"), 
            sep = "\t",
            header = False, 
            index = False)
    reads_out = pd.DataFrame(reads)
    reads_out.to_csv(os.path.join(outdir, "reads_ids.txt"), 
            header = False,
            index = False)
    logger.write("\tDone.")

def find_dna_barcodes(seq1, seq2, outdir, mode):
    n = 500000
    logger.write("  * Parsing file {}...".format(os.path.basename(seq1)))
    seq1_fname = os.path.basename(('.').join(seq1.split('.')[:-2]))
    seq1_df = parse_fastq(seq1)
    seq1_df_ls = [seq1_df[i:i+n] for i in range(0, seq1_df.shape[0], n)]
    seq1_parsed_full = []
    seq1_parsed_20bbc = []
    bar_format = '        |{bar}| {percentage:3.0f}% {n_fmt}/{total_fmt} {unit}'
    gdna_seq1_ext2 = 'TCGCCGTGTAATAATTCTAGA(.*?)AGATCGGAAGAGCG' # convert
    gdna_seq1_ext1 = 'CCGACGCTCTTCCGATCT(.*?)TCTAGAATTATTACACGG' #no convert
    gdna_seq1_ext = 'CCGACGCTCTTCCGATCT(.*?)TCTAGAATTATTACACGG|TCGCCGTGTAATAATTCTAGA(.*?)AGATCGGAAGAGCG'
    logger.write("\tSearching for barcodes...")
    for df in tqdm.tqdm(seq1_df_ls,
            total = len(seq1_df_ls),
            unit = ' batches',
            ncols = 75,
            bar_format = bar_format):
        df1 = df.copy()
        if mode == 'cdna':# no rev_com convertion
            df1['barcode'] = df1['seq'].apply(lambda x: x.split('TCTAGAATTATTACACGG', 1)[0])
            df1 = df1.dropna()
            df1['barcode_len'] = df1['barcode'].str.len()
            df1_20bbc = df1[df1['barcode_len'] == 20]
            seq1_parsed_full.append(df1)
            seq1_parsed_20bbc.append(df1_20bbc)
        else:
            df1_ext_res = df1.seq.astype(str).str.extract(r'{}'.format(gdna_seq1_ext))
            df1 = pd.concat([df1,df1_ext_res],axis=1)
            df1.columns = ['name', 'seq', 'barcodes_keep', 'barcodes_convert']
            df1['barcodes_convert'] = np.where(df1['barcodes_convert'].isna(), '0',
                    df1['barcodes_convert'])
            df1['barcodes_convert'] = df1['barcodes_convert'].apply(lambda x:
                    ''.join(Seq(x).reverse_complement()))
            df1['barcode'] = np.where(df1['barcodes_convert'] == '0', df1['barcodes_keep'],
                    df1['barcodes_convert'])
            df1['barcode_len'] = df1['barcode'].str.len()
            df1_20bbc = df1[df1['barcode_len'] == 20]
            #df1_20bbc = df1_20bbc[df1_20bbc['barcode'].notnull()].copy()
            df1_bc = df1[['name', 'seq','barcode','barcode_len']]
            df1_bc = df1_bc.dropna()
            seq1_parsed_full.append(df1_bc)
            seq1_parsed_20bbc.append(df1_20bbc)
    seq1_full_out = pd.concat(seq1_parsed_full)
    logger.write("\tBarcodes found in {}/{} reads.".format(len(seq1_full_out), 
        len(seq1_df)))
    seq1_20bbc_out = pd.concat(seq1_parsed_20bbc)
    logger.write("\t20b barcodes found in {}/{} reads.".format(len(seq1_20bbc_out),
        len(seq1_df)))
    logger.write("\tWriting output files...")
    seq1_full_ofile = seq1_fname+'_bc'+'.txt.gz'
    seq1_full_out.to_csv(os.path.join(outdir,seq1_full_ofile),
            sep="\t",
            header = True,
            index = False,
            compression = 'gzip')
    seq1_20bbc_ofile = seq1_fname+'_20b_bc'+'.txt.gz'
    seq1_20bbc_out.to_csv(os.path.join(outdir,seq1_20bbc_ofile),
            sep="\t",
            header = True,
            index = False,
            compression = 'gzip')
    
    logger.write("  * Parsing file {}...".format(os.path.basename(seq2)))
    seq2_fname = os.path.basename(('.').join(seq2.split('.')[:-2]))
    seq2_df = parse_fastq(seq2)
    seq2_df_ls = [seq2_df[i:i+n] for i in range(0, seq2_df.shape[0], n)]
    seq2_parsed_full = []
    seq2_parsed_20bbc = []
    cdna_seq2_ext = 'GTAATAATTCTAGA(.*?)AGATCGGAAGAGCGTC' # convert
    gdna_seq2_ext1 = 'TCGCCGTGTAATAATTCTAGA(.*?)AGATCGGAAGAGCG' # convert 
    gdna_seq2_ext2 = 'CCGACGCTCTTCCGATCT(.*?)TCTAGAATTATTACACGG' #no convert
    gdna_seq2_ext = 'CCGACGCTCTTCCGATCT(.*?)TCTAGAATTATTACACGG|TCGCCGTGTAATAATTCTAGA(.*?)AGATCGGAAGAGCG'
    for df in tqdm.tqdm(seq2_df_ls,
            total = len(seq2_df_ls),
            unit = ' batches',
            ncols = 75,
            bar_format = bar_format):
        df1 = df.copy()
        if mode == 'cdna': # with rev_com conversion
            df1['barcode'] = df1.seq.astype(str).str.extract(r'{}'.format(cdna_seq2_ext))
            df1 = df1.dropna()
            df1['barcode_len'] = df1['barcode'].str.len()
            df1_20bbc = df1[df1['barcode_len'] == 20]
            df1_20bbc = df1_20bbc[df1_20bbc['barcode'].notnull()].copy()
            df1_20bbc['barcode'] = df1_20bbc['barcode'].apply(lambda x:
                    ''.join(Seq(x).reverse_complement()))
            seq2_parsed_full.append(df1)
            seq2_parsed_20bbc.append(df1_20bbc)
        else: # mode == 'gdna' 
            df1_ext_res = df1.seq.astype(str).str.extract(r'{}'.format(gdna_seq2_ext))
            df1 = pd.concat([df1,df1_ext_res],axis=1)
            df1.columns = ['name', 'seq', 'barcodes_keep', 'barcodes_convert']
            df1['barcodes_convert'] = np.where(df1['barcodes_convert'].isna(), '0',
                    df1['barcodes_convert'])
            df1['barcodes_convert'] = df1['barcodes_convert'].apply(lambda x:
                    ''.join(Seq(x).reverse_complement()))
            df1['barcode'] = np.where(df1['barcodes_convert'] == '0', df1['barcodes_keep'],
                    df1['barcodes_convert'])
            df1['barcode_len'] = df1['barcode'].str.len()
            df1_20bbc = df1[df1['barcode_len'] == 20]
            df1_bc = df1[['name', 'seq','barcode','barcode_len']]
            df1_bc = df1_bc.dropna()
            seq2_parsed_full.append(df1_bc)
            seq2_parsed_20bbc.append(df1_20bbc)
    seq2_full_out = pd.concat(seq2_parsed_full)
    logger.write("\tBarcodes found in {}/{} reads.".format(len(seq2_full_out),
        len(seq2_df)))
    seq2_20bbc_out = pd.concat(seq2_parsed_20bbc)
    logger.write("\t20b barcodes found in {}/{} reads.".format(len(seq2_20bbc_out),
        len(seq2_df)))
    logger.write("\tWriting output files...")
    seq2_full_ofile = seq2_fname+'_bc'+'.txt.gz'
    seq2_full_out.to_csv(os.path.join(outdir,seq2_full_ofile),
            sep="\t",
            header = True,
            index = False,
            compression = 'gzip')
    seq2_20bbc_ofile = seq2_fname+'_20b_bc'+'.txt.gz'
    seq2_20bbc_out.to_csv(os.path.join(outdir,seq2_20bbc_ofile),
            sep="\t",
            header = True,
            index = False,
            compression = 'gzip')

class Logger(object):
    def __init__(self, logfile=None, verbose=True):
        self.console = sys.stdout
        self.verbose = verbose
        if logfile is not None:
            self.log = open(logfile, 'w')
        else:
            self.log = None
    def write(self, message):
        if self.verbose:
            self.console.write(message+'\n')
        if self.log is not None:
            self.log.write(message+'\n')
            self.log.flush()
    def verbose(self, verbose=True):
        self.verbose = verbose

def log_settings(args, logger):
    now = datetime.datetime.now()
    logger.write('\n')
    logger.write(f'{now.strftime("%d/%m/%Y %H:%M:%S")}')
    if args.mode == 'enhancer':
        logger.write("#~*~*~*~ RUN MODE: Enhancer ~*~*~*~#")
    elif args.mode == 'cdna':
        logger.write("#~*~*~*~ RUN MODE: cdna ~*~*~*~#")
    elif args.mode == 'gdna':
        logger.write("#~*~*~*~ RUN MODE: gdna ~*~*~*~#")

def parse_args():
    parser = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description= '''MPRA sequencing data analysis pipeline.''')
    parser.add_argument(
            '-r1', '--read1', 
            help='''Path to a gzipped fastq file. Sequences in this file will 
            be mapped against the reference MPRA library to find refSNP containing
            reads. Required if run --mode='enhancer'.''')
    parser.add_argument(
            '-r2', '--read2',
            help='''Path to a gzipped fastq file. Sequences in this file will
            be used to search for barcodes of refSNPs. Required if run --mode='enhancer'.''')
    parser.add_argument(
            '-ds1', '--dna-seq1', nargs='+',
            help='''Provide a fastq file or a list of space separated fastq files 
            (forward reads) containing 
            1) cDNA sequences (required if run --mode='cdna')
            2) gDNA sequences (required if run --mode='gdna').''')
    parser.add_argument(
            '-ds2', '--dna-seq2', nargs='+',
            help='''Provide a fastq file or a list of space separated fastq files
            (reverse reads) containing
            1) cDNA sequences (required if run --mode='cdna') or
            2) gDNA sequences (required if run --mode='gdna').''')
    parser.add_argument(
            '-o', '--outdir', required = True,
            help='''Provide a path to save output files.''')
    parser.add_argument(
            '-i', '--index-dir', 
            help='''Provide path to the directory where reference index (generated 
            using bwa) files are located and also include prefix of the file''')
    parser.add_argument(
            '-m', '--mode', required=True, type=str, choices=['enhancer', 'cdna', 'gdna'],
            help='''Choose a run mode to find the barcodes from enhancer, cDNA or 
            gDNA sequences''')
    return parser.parse_args()

if __name__ == '__main__':

    args = parse_args()
    
    if args.mode == 'enhancer':
        if not (args.read1 and args.read2 and args.index_dir):
            sys.exit('''\n\tMissing arguments.
            See 'mpra_pipeline.py -h' for usage info.\n''')
        if not (args.read1.endswith('.gz') and args.read2.endswith('.gz')):
            sys.exit('''\n\tInput file format not recognised.
            Make sure to provide gzipped fastq file.\n''')

    if (args.mode == 'cdna' or args.mode == 'gdna') and not args.mode == 'enhancer':
        if not (args.dna_seq1 and args.dna_seq2):
            sys.exit('''\nMissing input fastq read files. '-ds1' and/or '-ds2'.
            See 'mpra_pipeline.py -h' for usage info.\n''')
        
        for (seq1file,seq2file)  in zip(args.dna_seq1,args.dna_seq2):
            if not (seq1file.endswith('.gz') and seq2file.endswith('.gz')):
                sys.exit('''\n\tInput file format not recognised.
                Make sure to provide gzipped fastq file.\n''')
         
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir, exist_ok=True)

    logger = Logger(logfile = os.path.join(args.outdir, 'mpra_pipeline.log'))
    start_time = time.time()
    log_settings(args, logger)
    
    if args.mode == 'enhancer':
        fq_read1 = args.read1
        fq_read2 = args.read2

        # Trim constant sequences 
        logger.write(''' *** Trimming CS...''')
        cutadapt_out = os.path.join(args.outdir, 'valid_reads_post_trimming_CS.fq.gz')
        cutadapt_log = os.path.join(args.outdir, 'cutadapt_CS_trim_log.txt')
        cutadapt_discard = os.path.join(args.outdir, 'discarded_reads_post_trim.fq.gz')
        cut_cmd ="cutadapt -j 16 --no-indels --discard-untrimmed -e 0.01 -O 21 -m 125\
                --too-short-output "+cutadapt_discard+"\
                -g GGCCTAACTGGCCGCTTGACG -o "+cutadapt_out+" "+fq_read1+" > "+cutadapt_log+"" 
        run_shell_commands(cut_cmd)
        logger.write("\tTrimming done.")

        # Align trimmed reads to custom MPRA library using bwa
        logger.write(" *** Aligning trimmed reads to custom MPRA library...")
        index = args.index_dir
        R1_trimmed = os.path.join(args.outdir, 'valid_reads_post_trimming_CS.fq.gz')
        aligned_out = os.path.join(args.outdir, 'aligned_reads.sam')
        align_cmd = "bwa mem -t 16 -K 100000000 -M -B 40 -O 60 -E 10 -L 150 -v 0\
                "+index+" "+R1_trimmed+" > "+aligned_out+""
        run_shell_commands(align_cmd)
        logger.write("\tRead alignment done.")

        logger.write(" *** Converting SAM into BAM file...")
        bam_out = os.path.join(args.outdir, "aligned_reads.bam")
        s2b_cmd = "samtools view -Sb "+aligned_out+" > "+bam_out+""
        run_shell_commands(s2b_cmd)
    
        logger.write("\tFiltering out unmapped/supplementary/secondary alignment reads...")
        sam_mpri = os.path.join(args.outdir, "mapped_primary_alignment.sam")
        filter_cmd = "samtools view -h -F 4 "+aligned_out+"\
                | grep -v 'SA:Z:' | grep -v 'XA:Z:' > "+sam_mpri+""
        run_shell_pipe(filter_cmd)

        bam_mpri = os.path.join(args.outdir, "mapped_primary_alignment.bam")
        bam_cmd = "samtools view -Sb "+sam_mpri+" > "+bam_mpri+""
        sp.Popen(bam_cmd, shell=True).wait()
    
        #find reads that map perfectly to reference library 
        find_proper_alignments(bam_mpri, args.outdir)
    
        # subset reverse reads based on read names
        logger.write(" *** Subsetting reads containing barcodes...")
        q_names =  os.path.join(args.outdir, "reads_ids.txt")
        R2_barcode = os.path.join(args.outdir, "subset.barcoded.reads.fastq.gz")
        seqtk_cmd = "seqtk subseq "+fq_read2+" "+q_names+" | gzip -c > "+R2_barcode+""
        run_shell_pipe(seqtk_cmd)
        barcode_df = parse_fastq(R2_barcode)
        ref_ids = pd.read_csv(os.path.join(args.outdir, "reads_refs_ids.txt"),
                sep="\t", names=['name', 'ref'])
        rname_ref_barc = pd.merge(ref_ids, barcode_df, how='inner',
                on='name', sort=False)
        find_barcodes(rname_ref_barc, args.outdir)
        logger.write("\tWriting outputs.")

    if args.mode == 'cdna' or args.mode == 'gdna':
        for seq1,seq2 in zip(args.dna_seq1, args.dna_seq2):
            find_dna_barcodes(seq1,seq2,args.outdir,args.mode)

    msg = ' *** Program completed successfully ***\nTotal time elasped: {:.2f} mins.'.format(
            (time.time() - start_time)/60)
    logger.write('{}\n'.format(msg))
