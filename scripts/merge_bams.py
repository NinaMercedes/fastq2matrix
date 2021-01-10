import fastq2matrix as fm
import argparse
import sys
from uuid import uuid4
import os


def main(args):
    if args.prefix and not args.bams:
        individual_bams = ["%s/%s%s" % (args.dir,run,args.suffix) for run in args.prefix.split("_")]
        new_id = args.new_id if args.new_id else args.prefix
    elif args.bams:
        individual_bams = args.bams.split(",")
        if args.prefix:
            new_id = args.prefix
        else:
            new_id = args.new_id if args.new_id else "_".join([bam.split("/")[-1].replace(args.suffix,"") for bam in individual_bams])
    if not args.prefix and not args.bams:
        quit("\nExpected at least --prefix or --bams... Exiting!\n")

    if len(individual_bams)==1:
        sys.stderr.write("Need more than one bam... Exiting!\n")
        quit()
    for bam in individual_bams:
        fm.filecheck(bam)
    new_bamfile = "%s/%s%s" % (args.dir,new_id,args.suffix)
    tmp_file = str(uuid4())
    with open(tmp_file,"w") as O:
        for l in fm.cmd_out("samtools view -H %s" % individual_bams[0]):
            row = l.strip().split("\t")
            if args.format=="cram" and row[0]=="@SQ":
                row[4] = "UR:%s" % os.path.abspath(os.path.expanduser(args.ref))
            if row[0]=="@RG":
                # continue
                row[1] = "ID:%s" % new_id
                row[2] = "SM:%s" % new_id
            O.write("%s\n" % "\t".join(row))
    if args.format=="cram":
        if not args.ref:
            quit("\nOutput in cram format requires referene with --ref\n")
        cram_arg = "| samtools view -CT %s " % args.ref
    else:
        cram_arg = ""
    # fm.run_cmd("samtools merge %s -@ %s - %s | samtools reheader -i %s - | samtools addreplacerg -@ %s - -r 'ID:%s\\tSM:%s\\tPL:Illumina' -o %s" % (
    tmp_bam = str(uuid4())
    fm.run_cmd("samtools merge -@ %s - %s | samtools addreplacerg -@ %s - -r 'ID:%s\\tSM:%s\\tPL:Illumina' -O BAM | samtools reheader %s - %s > %s" % (
        args.threads," ".join(individual_bams), args.threads,new_id,new_id,tmp_file,cram_arg,new_bamfile)
    )
    # fm.run_cmd("samtools reheader %s %s %s  >  %s" % (tmp_file, tmp_bam, cram_arg,new_bamfile))
    fm.run_cmd("samtools index %s" % new_bamfile)
    # fm.rm_files([tmp_file])

parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--bams',help='Comma seperated BAM file list')
parser.add_argument('--new-id',help='New id for BAM')
parser.add_argument('--prefix',help='Prefix')
parser.add_argument('--suffix',default=".bqsr.bam",help='BAM file suffix')
parser.add_argument('--dir',default=".",help='Directory containing bams')
parser.add_argument('--ref',help='Reference for cram output')
parser.add_argument('--format',default="bam",choices=["bam","cram"],help='Output format')
parser.add_argument('--threads',default=4, type=int, help='Number of threads for parallel operations')
parser.set_defaults(func=main)
args = parser.parse_args()
args.func(args)
#
