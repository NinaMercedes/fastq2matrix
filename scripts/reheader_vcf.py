import sys
import argparse
from uuid import uuid4
import fastq2matrix as fm

def main(args):
    args.uuid = str(uuid4())
    with open(args.uuid,"w") as O:
        for l in fm.cmd_out("bcftools view -h %(vcf)s" % vars(args)):
            row = l.strip().split("\t")
            if row[0]=="#CHROM":
                # continue
                row[9] = args.new_id
            O.write("%s\n" % "\t".join(row))
    fm.run_cmd("bcftools reheader -h %(uuid)s %(vcf)s -o %(new_id)s%(extension)s" % vars(args))
    fm.rm_files([args.uuid])

parser = argparse.ArgumentParser(description='XXX pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--vcf',help='Bam/Cram file',required=True)
parser.add_argument('--new-id',help='New ID file',required=True)
parser.add_argument('--extension',help='Extension file',required=True)
parser.set_defaults(func=main)
args = parser.parse_args()
args.func(args)
