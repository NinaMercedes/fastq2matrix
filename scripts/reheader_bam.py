import sys
import argparse
from uuid import uuid4
import fastq2matrix as fm

def main(args):
    args.uuid = str(uuid4())
    rg_done = False
    with open(args.uuid,"w") as O:
        for l in fm.cmd_out("samtools view -H %(bam)s" % vars(args)):
            row = l.strip().split("\t")
            if row[0]=="@RG":
                # continue
                if rg_done: continue
                row[1] = "ID:%s" % args.new_id
                row[2] = "SM:%s" % args.new_id
                rg_done = True
            O.write("%s\n" % "\t".join(row))
    fm.run_cmd("samtools reheader %(uuid)s %(bam)s > %(new_id)s%(extension)s" % vars(args))
    fm.rm_files([args.uuid])

parser = argparse.ArgumentParser(description='XXX pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--bam',help='Bam/Cram file',required=True)
parser.add_argument('--new-id',help='New ID file',required=True)
parser.add_argument('--extension',help='Extension file',required=True)
parser.set_defaults(func=main)
args = parser.parse_args()
args.func(args)
