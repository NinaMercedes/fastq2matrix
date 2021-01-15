import sys
import argparse
import os
from tqdm import tqdm


#15650507 + 0 in total (QC-passed reads + QC-failed reads)
#0 + 0 secondary
#34633 + 0 supplementary
#166695 + 0 duplicates
#15546150 + 0 mapped (99.33% : N/A)
#15615874 + 0 paired in sequencing
#7807937 + 0 read1
#7807937 + 0 read2
#14991678 + 0 properly paired (96.00% : N/A)
#15499656 + 0 with itself and mate mapped
#11861 + 0 singletons (0.08% : N/A)
#0 + 0 with mate mapped to a different chr
#0 + 0 with mate mapped to a different chr (mapQ>=5)

def main(args):
    args.flagstat_extension = args.bam_extension+".flagstat.txt"
    args.coverage_extension = args.bam_extension+".genomecov.txt"
    if args.sample_file:
        samples = [x.rstrip() for x in open(args.sample_file).readlines()]
    else:
        samples = [x.replace(args.bam_extension,"") for x in os.listdir("%s/" % args.dir) if x[- len(args.bam_extension):]==args.bam_extension]
    with open(args.out,"w") as O:
        for s in tqdm(samples):
            res = []
            if args.allow_missing and (not os.path.isfile("%s/%s%s" % (args.dir,s,args.coverage_extension)) or not os.path.isfile("%s/%s%s" % (args.dir,s,args.flagstat_extension))):
                O.write("%s\tNA\tNA\tNA\tNA\n" % (s))
                continue

            for i,l in enumerate(open("%s/%s%s" % (args.dir,s,args.flagstat_extension))):
                row = l.rstrip().split()
                if i==4:
                    res.append(str(row[0]))
                    res.append(str(row[4][1:-1]))


            cutoff=10
            cumfrac = 0.0
            tmp_positions = 0
            for l in open("%s/%s%s" % (args.dir,s,args.coverage_extension)):
                row = l.strip().split()
                chr = row[0]
                if row[0]!="genome": continue
                dp = int(row[1])
                freq = int(row[2])
                seqsize = int(row[3])
                fraction = float(row[4])
                if dp<args.depth_cutoff:
                    cumfrac+=fraction
                if dp==0:
                    genome_size = seqsize
                num_positions = tmp_positions+freq
                if num_positions>(genome_size)/2 and tmp_positions<=(genome_size)/2:
                    if genome_size%2==0 and genome_size//2==tmp_positions:
                        median_dp = dp-0.5
                    else:
                        median_dp = dp
                tmp_positions = num_positions
            res.append(str(median_dp))
            res.append(str(cumfrac))
            O.write("%s\t%s\n" % (s,"\t".join(res)))

parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--out',type=str,help='Sample file',required=True)
parser.add_argument('--sample-file',type=str,help='Sample file')
parser.add_argument('--dir',default=".",type=str,help='Directory')
parser.add_argument('--depth-cutoff',default=10,type=int,help='Add depth info')
parser.add_argument('--bam-extension',default=".bqsr.cram",type=str,help='Extension of bam files')
parser.add_argument('--allow-missing',action="store_true",help='Extension of bam files')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
