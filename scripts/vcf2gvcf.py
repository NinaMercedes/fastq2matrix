import argparse
import sys
import pysam
import numpy as np
import subprocess
from uuid import uuid4
import os
from collections import defaultdict

#Chromosome      1       .       T       <NON_REF>       .       .       END=1976        GT:DP:GQ:MIN_DP:PL      0/0:80:99:51:0,108,1800
#Chromosome      1977    .       A       G,<NON_REF>     2223.06 .       DP=68;ExcessHet=0.0000;MLEAC=2,0;MLEAF=1.00,0.00;RAW_MQandDP=244800,68  GT:AD:DP:GQ:PL:SB       1/1:0,64,0:64:99:2237,192,0,2237,192,2237:0,0,27,37

# Required info lines:
##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">
##INFO=<ID=ExcessHet,Number=1,Type=Float,Description="Phred-scaled p-value for exact test of excess heterozygosity">
##INFO=<ID=InbreedingCoeff,Number=1,Type=Float,Description="Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation">
##INFO=<ID=MLEAC,Number=A,Type=Integer,Description="Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed">
##INFO=<ID=MLEAF,Number=A,Type=Float,Description="Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed">
##INFO=<ID=MQRankSum,Number=1,Type=Float,Description="Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities">
##INFO=<ID=RAW_MQandDP,Number=2,Type=Integer,Description="Raw data (sum of squared MQ and total depth) for improved RMS Mapping Quality calculation. Incompatible with deprecated RAW_MQ formulation.">
##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias">

# Required format lines:
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum DP observed within the GVCF block">
##FORMAT=<ID=PGT,Number=1,Type=String,Description="Physical phasing haplotype information, describing how the alternate alleles are phased in relation to one another; will always be heterozygous and is not intended to describe called alleles">
##FORMAT=<ID=PID,Number=1,Type=String,Description="Physical phasing ID information, where each unique ID within a given sample (but not across samples) connects records within a phasing group">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phasing set (typically the position of the first variant in the set)">
##FORMAT=<ID=SB,Number=4,Type=Integer,Description="Per-sample component statistics which comprise the Fisher's Exact Test to detect strand bias.">


info_tags = [
    {'ID': 'BaseQRankSum', 'Number': '1', 'Type':'Float', 'Description': 'Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities'},
    {'ID': 'DP', 'Number': '1', 'Type':'Integer', 'Description': 'Approximate read depth; some reads may have been filtered'},
    {'ID': 'END', 'Number': '1', 'Type':'Integer', 'Description': 'Stop position of the interval'},
    {'ID': 'ExcessHet', 'Number': '1', 'Type':'Float', 'Description': 'Phred-scaled p-value for exact test of excess heterozygosity'},
    {'ID': 'InbreedingCoeff', 'Number': '1', 'Type':'Float', 'Description': 'Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation'},
    {'ID': 'MLEAC', 'Number': 'A', 'Type':'Integer', 'Description': 'Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed'},
    {'ID': 'MLEAF', 'Number': 'A', 'Type':'Float', 'Description': 'Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed'},
    {'ID': 'MQRankSum', 'Number': '1', 'Type':'Float', 'Description': 'Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities'},
    {'ID': 'RAW_MQandDP', 'Number': '2', 'Type':'Integer', 'Description': 'Raw data (sum of squared MQ and total depth) for improved RMS Mapping Quality calculation. Incompatible with deprecated RAW_MQ formulation.'},
    {'ID': 'ReadPosRankSum', 'Number': '1', 'Type':'Float', 'Description': 'Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias'}
]

format_tags = [
    {'ID': 'AD', 'Number': 'R', 'Type':'Integer', 'Description': 'Allelic depths for the ref and alt alleles in the order listed'},
    {'ID': 'DP', 'Number': '1', 'Type':'Integer', 'Description': 'Approximate read depth (reads with MQ=255 or with bad mates are filtered)'},
    {'ID': 'GQ', 'Number': '1', 'Type':'Integer', 'Description': 'Genotype Quality'},
    {'ID': 'GT', 'Number': '1', 'Type':'String', 'Description': 'Genotype'},
    {'ID': 'MIN_DP', 'Number': '1', 'Type':'Integer', 'Description': 'Minimum DP observed within the GVCF block'},
    {'ID': 'PGT', 'Number': '1', 'Type':'String', 'Description': 'Physical phasing haplotype information, describing how the alternate alleles are phased in relation to one another; will always be heterozygous and is not intended to describe called alleles'},
    {'ID': 'PID', 'Number': '1', 'Type':'String', 'Description': 'Physical phasing ID information, where each unique ID within a given sample (but not across samples) connects records within a phasing group'},
    {'ID': 'PL', 'Number': 'G', 'Type':'Integer', 'Description': 'Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification'},
    {'ID': 'PS', 'Number': '1', 'Type':'Integer', 'Description': 'Phasing set (typically the position of the first variant in the set)'},
    {'ID': 'SB', 'Number': '4', 'Type':'Integer', 'Description': 'Per-sample component statistics which comprise the Fisher\'s Exact Test to detect strand bias.'}
]


def consecutive(data, stepsize=1):
    return np.split(data, np.where(np.diff(data) != stepsize)[0]+1)

def get_blocks(dp,variant_positions,min_depth,chrom):
    dp = np.array(dp)
    good_blocks = consecutive(np.setdiff1d(np.where(dp>=min_depth)[0],variant_positions))
    good_blocks = [(chrom,d[0],d[-1]+1,True) for d in good_blocks]

    bad_blocks = consecutive(np.where(dp<min_depth)[0])
    bad_blocks = [(chrom,d[0],d[-1]+1,False) for d in bad_blocks]

    return sorted(good_blocks + bad_blocks,key=lambda x: x[1])

def load_fasta(fasta):
    refseq = {}
    for line in open(fasta):
        if line.startswith(">"):
            chrom = line.strip().split()[0][1:]
            refseq[chrom] = []
        else:
            refseq[chrom].append(line.strip())

    for chrom in refseq:
        refseq[chrom] = "".join(refseq[chrom])

    return refseq

def get_nonref_blocks_from_depth(bam,variants,vcf,refseq,min_depth=10):
    tmpfile = str(uuid4())
    # tmpfile = "dp.txt"
    cmd = "samtools depth -a %s > %s" % (bam,tmpfile)
    sys.stderr.write(cmd+"\n")
    subprocess.check_call(cmd,shell=True)
    blocks = []
    dp = []
    current_chrom = None
    for line in open(tmpfile):
        chrom,pos,depth = line.strip().split()
        pos = int(pos)
        depth = int(depth)
        if chrom != current_chrom:
            if current_chrom is not None:
                variant_positions = [v.pos-1 for v in variants if v.chrom == current_chrom]
                blocks = blocks + get_blocks(dp,variant_positions,min_depth,chrom)
            dp = [depth]
            current_chrom = chrom
        else:
            dp.append(depth)
    variant_positions = [v.pos-1 for v in variants if v.chrom == current_chrom]
    blocks = blocks + get_blocks(dp,variant_positions,min_depth,chrom)

    # print(blocks)
    vcf_lines = []
    for block in blocks:
        # print(block)
        v = vcf.new_record(
            contig=block[0], 
            start=block[1],
            stop=block[2],
            alleles=(refseq[block[0]][block[1]],'<NON_REF>')
        )
        if block[3]==True:
            v.samples[0]['GT'] = (0,0)
            v.samples[0]['DP'] = 10
            v.samples[0]['GQ'] = 10
            v.samples[0]['MIN_DP'] = 10
            v.samples[0]['PL'] = (0,0,0)
        else:
            v.samples[0]['GT'] = (0,0)
            v.samples[0]['DP'] = 0
            v.samples[0]['GQ'] = 10
            v.samples[0]['MIN_DP'] = 0
            v.samples[0]['PL'] = (0,0,0)

        vcf_lines.append(v)
        
    os.remove(tmpfile)
    return vcf_lines
    

def main(args):

    refseq = load_fasta(args.fasta)



    vcf = pysam.VariantFile(args.vcf)
    
    vcfh = pysam.VariantHeader()
    for tag in info_tags:
        vcfh.add_meta('INFO',items=tag.items())
    for tag in format_tags:
        vcfh.add_meta('FORMAT',items=tag.items())
    # add contigs lines
    for chrom in refseq:
        vcfh.contigs.add(chrom,length=len(refseq[chrom]))
    # add sample line
    vcfh.add_sample('sample')
    new_vcf = pysam.VariantFile(args.out,'w',header=vcfh)

    variants = []
    for variant in vcf:
        v = new_vcf.new_record(
            contig=variant.chrom,
            start=variant.pos-1,
            alleles=variant.alleles
        )
        v.alts += ('<NON_REF>',)
        # print(str(v))
        # print(str(variant))

        v.info['DP'] = variant.info['DP']
        v.info['ExcessHet'] = (0,)
        v.info['MLEAC'] = (0,0)
        v.info['MLEAF'] = (0,0)
        v.info['RAW_MQandDP'] = (0,0)
        
        v.samples[0]['GT'] = variant.samples[0]['GT']
        v.samples[0]['DP'] = variant.samples[0]['DP']
        v.samples[0]['AD'] = (0,1,0)
        v.samples[0]['GQ'] = (99,)
        v.samples[0]['PL'] = (2237,  192,    0, 2237,  192, 2237)
        v.samples[0]['SB'] = (variant.info['SRF'], variant.info['SRR'], variant.info['SAF'][0], variant.info['SAR'][0])
        
        variants.append(v)

    nonref_blocks = get_nonref_blocks_from_depth(args.bam,variants,new_vcf,refseq)

    for v in sorted(variants + nonref_blocks,key=lambda x: (x.chrom,x.pos)):

        new_vcf.write(v)
        


parser = argparse.ArgumentParser(description='fastq2matrix pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--vcf','-v',help='VCF file',required=True)
parser.add_argument('--bam','-b',help='BAM file',required=True)
parser.add_argument('--fasta','-f',help='FASTA file',required=True)
parser.add_argument('--out','-o',help='BAM file',required=True)

parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)