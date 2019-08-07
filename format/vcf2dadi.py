#coding:utf-8

"""
Author : Sitao Zhu
Date: 2019-4-18
Desc:
# genome file is the location of genome file;
# sample list file is like this:
# BEGIN
# sample1    population1
# sample2    population1
# sample3    population2
# END
# the sample name should be same as which in vcf file.
# Usage: python vcf2dadi.py -r <genome file> -f <vcf file> -l <sample list file>

"""


import os
import re
import getopt
from format.fasta import Fasta



def parse_sample(sample_list):
    """ put sample list file into a dict """
    samples={}
    with open(sample_list,'r') as F:
        for line in F.readline():
            if line.startswith('#'):
                continue
            else:
                fields = line.strip().split()
                samples[fields[0]] = fields[1]
    return samples

def parse_vcf(vcf,samples):
    """ put vcf file into a dict """
    pattern = re.compile(r'/^(.)\/(.):/')
    record = {}
    vcf = {}
    pop = {}
    with open(vcf,'r') as F:
        for line in F.readline():
            line = line.strip()
            if line.startswith('##'):continue
            if line.startswith('#CHROM'):
                fields = line.split()
                for i in range(9,len(fields)):
                    if fields[i] not in sample.keys():
                        continue
                    else:
                        record[i] = samples[fields[i]]
            a = line.split()
            chrom,pos,ref,alt = a[0],a[1],a[3],a[4]
            if ',' in alt:continue
            vcf[chrom] = {pos: {'ref': ref}}
            vcf[chrom] = {pos: {'alt': alt}}


            for key in record.keys():
                indv = record[key]
                # vcf[chrom][pos] = {indv: {'a1', 0}}
                # vcf[chrom][pos] = {indv: {'a2', 0}}
                pop[indv] = 1
                geno = a[key]
                match = pattern.match(geno)
                if match:
                    tempa = match.group(1)
                    tempb = match.group(2)
                    a1,a2 = 0,0
                    if tempa == '.' or tempb == '.':
                        a1 = 0
                        a2 = 0
                    elif int(tempa) + int(tempb) == 1:
                        a1 = 1
                        a2 = 1
                    elif int(tempa) + int(tempb) == 2:
                        a1 = 0
                        a2 = 2
                    elif int(tempa) + int(tempb) == 0 :
                        a1 = 2
                        a2 = 0
                    if not vcf[chrom][pos][indv]:
                        vcf[chrom][pos] = {indv:{'a1':a1}}
                        vcf[chrom][pos] = {indv:{'a2':a2}}
                    else:
                        vcf[chrom][pos][indv]['a1'] += a1
                        vcf[chrom][pos][indv]['a2'] += a2

    with open("list.data",'w') as F:
        title = "NAME\tOUT\tAllele1"
        for key in sorted(pop.keys()):
            title += "\t%s"%(key)
        title += "\tAllele2"
        for key in sorted(pop.keys()):
            title += "\t%s"%(key)
        title += "\tGene\tPostion\n"
        F.writelines(title+"\n")
        fasta = Fasta(fasta_path)
