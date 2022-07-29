import pandas as pd
import os
import sys
import argparse
import gzip
import itertools
import subprocess as sp
import time
import datetime
import pysam
import re

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
    print("\tFinding barcodes...")
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
        print ("\tShell command exited with an error.")
        print(exception.output)
        return None

def run_shell_pipe(cmd):
    try:
        sp.Popen(cmd, shell=True,
                stdin=sp.PIPE,
                stdout=sp.PIPE,
                stderr=sp.PIPE).wait()
    except sp.CalledProcessError as exception:
        print ("\tShell command exited with an error.")
        print(exception.output)
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
    read_ref_ls = []
    reads = []
    print("\tFinding sequences that map to ref library...")
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
            # get seq alignment record
            c_type = [1,2,3,4,5,6,8]
            codes = [cigar[0] for _, cigar in enumerate(cigarline) 
                    if cigar[0] in c_type]
            #alignment match (can be a sequence match or mismatch)
            if not codes:
                for _,tagtype in enumerate(tags):
                    if tagtype[0] == "NM" and tagtype[1] == 0:
                        cigar_matched.write(read)
                        read_ref_ls.append([name, rname])
                        reads.append(name)
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
    print("\tDone.")

def parse_args():
    parser = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description= '''MPRA sequencing data analysis pipeline.''')
    parser.add_argument(
            '-r1', '--read1', required=True, 
            help='''Path to a gzipped fastq file. Sequences in this file will 
            be mapped against the reference MPRA library to find SNPs.''')
    parser.add_argument(
            '-r2', '--read2', required=True,
            help='''Path to a gzipped fastq file. Sequences in this file will
            be used to search for barcodes.''')
    parser.add_argument(
            '-o', '--outdir', required=True,
            help='''Provide a path to save output files.''')
    parser.add_argument(
            '-i', '--index-dir', required=True,
            help='''Provide a path to the directory where reference index (generated 
            using bwa) files are located along with the file prefix''')
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()

    if not (args.read1 and args.read2):
        sys.exit('''Missing input fastq files.
                See 'mpra_pipeline.py -h' for more details''')
    if not args.index_dir:
        sys.exit('''Missing reference index files.
                See 'mpra_pipeline.py -h' for more details''')
    if not args.outdir:
        sys.exit('''Missing argument: output directory.
                See 'mpra_pipeline.py -h' for more details''')
    
    if not (args.read1.endswith('.gz') and args.read2.endswith('.gz')):
        sys.exit('''Input file format not recognised.
                Make sure to provide gzipped fastq file.''')  
    else:
        fq_read1 = args.read1
        fq_read2 = args.read2

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir, exist_ok=True)

    start_time = time.time()
    
    # Trim constant sequences 
    print(''' *** Trimming CS...''')
    cutadapt_out = os.path.join(args.outdir, 'valid_reads_post_trimming_CS.fq.gz')
    cutadapt_log = os.path.join(args.outdir, 'cutadapt_CS_trim_log.txt')
    cutadapt_discard = os.path.join(args.outdir, 'discarded_reads_post_trim.fq.gz')
    cut_cmd ="cutadapt -j 16 --no-indels --discard-untrimmed -e 0.01 -O 21 -m 125\
            --too-short-output "+cutadapt_discard+"\
            -g GGCCTAACTGGCCGCTTGACG -o "+cutadapt_out+" "+fq_read1+" > "+cutadapt_log+"" 
    run_shell_commands(cut_cmd)
    print("\tTrimming done.")

    # Align trimmed reads to custom MPRA library using bwa
    print(" *** Aligning trimmed reads to custom MPRA library...")
    index = args.index_dir
    R1_trimmed = os.path.join(args.outdir, 'valid_reads_post_trimming_CS.fq.gz')
    aligned_out = os.path.join(args.outdir, 'aligned_reads.sam')
    align_cmd = "bwa mem -t 16 -K 100000000 -M -B 40 -O 60 -E 10 -L 150 -v 0\
            "+index+" "+R1_trimmed+" > "+aligned_out+""
    run_shell_commands(align_cmd)
    print("\tRead alignment done.")

    print(" *** Converting SAM into BAM file...")
    bam_out = os.path.join(args.outdir, "aligned_reads.bam")
    s2b_cmd = "samtools view -Sb "+aligned_out+" > "+bam_out+""
    run_shell_commands(s2b_cmd)
    
    unM_count_cmd = "samtools view -f4 -c "+bam_out+"" # get the count of unmapped reads
    unMapped = sp.check_output(unM_count_cmd, shell=True,
            stderr=sp.STDOUT)
    unMapped = int(unMapped.decode("utf-8"))
    print("\tNumber of reads unmapped: {}".format(unMapped))

    ## Remove unmapped/supplementary/secondary alignment reads from BAM file
    # TODO: Keep discarded reads?
    print("\tFiltering out unmapped/supplementary/secondary alignment reads...")
    sam_mpri = os.path.join(args.outdir, "mapped_primary_alignment.sam")
    filter_cmd = "samtools view -h -F 4 "+aligned_out+"\
            | grep -v 'SA:Z:' | grep -v 'XA:Z:' > "+sam_mpri+""
    run_shell_pipe(filter_cmd)

    bam_mpri = os.path.join(args.outdir, "mapped_primary_alignment.bam")
    bam_cmd = "samtools view -Sb "+sam_mpri+" > "+bam_mpri+""
    sp.Popen(bam_cmd, shell=True).wait()
    
    #find reads that map perfectly(?) to reference library 
    find_proper_alignments(bam_mpri, args.outdir)
    
    # subset R2 reads based on read names
    print(" *** Subsetting reads containing barcodes...")
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
    print("\tDone.")
    msg = ' *** Program completed successfully ***\nTotal time elasped: {:.2f} mins.'.format(
            (time.time() - start_time)/60)
    print(msg)
