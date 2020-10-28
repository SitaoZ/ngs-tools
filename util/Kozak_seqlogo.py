#coding:utf-8 

'''
Author: Sitao Zhu
Date: 2020-10-22
Desc: plot a Position Weight Matrix for fasta file 
      using weblogo form http://weblogo.threeplusone.com/manual.html
'''
import os
import re
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse
parser = argparse.ArgumentParser(description='Process input file')
parser.add_argument('-i','--input', type=str,
                    help='input for plot')

pattern = re.compile(r'CDS=(\d+)-(\d+)')
args = parser.parse_args()
transcript = args.input
output = []
i = 0
for record in SeqIO.parse(transcript, 'fasta'):
    i += 1
    chrom = record.id
    seq = record.seq
    seq = seq.upper()
    matched = pattern.search(record.description)
    if matched and len(seq)>9 :
        start = int(matched.group(1))
        s = start - 4
        if s < 0 :
            continue
        e = start + 5
        seq10 = seq[s:e]
        rec = SeqRecord(
                    seq10,
                    id=chrom,
                    description="5utr_length=%d" % (8))
        output.append(rec)

SeqIO.write(output, 'utr_8_AUG.fa', 'fasta')
os.system('weblogo -f utr_8_AUG.fa -A dna -S 2 --errorbars NO --fineprint weglogo --annotate -3,-2,-1,1,2,3,4,5,6 -o kozak.eps')


