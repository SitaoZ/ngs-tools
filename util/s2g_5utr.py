#coding:utf-8 

'''
Author: Sitao Zhu
Date: 2020-10-22
Desc: plot a Position Weight Matrix for fasta file 
      using weblogo form http://weblogo.threeplusone.com/manual.html
'''
import os
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse
parser = argparse.ArgumentParser(description='Process input file')
parser.add_argument('-i','--input', type=str,
                    help='input for plot')

#transcript = sys.argv[1]
args = parser.parse_args()
transcript = args.input
output = []
i = 0
for record in SeqIO.parse(transcript, 'fasta'):
    i += 1
    chrom = record.id
    seq = record.seq
    seq = seq.upper()
    if len(seq) > 10 :
        seq10 = seq[-10:]
        rec = SeqRecord(
                    seq10,
                    id=chrom,
                    description="5utr_length=%d" % (10))
        output.append(rec)

SeqIO.write(output, 'utr_10_AUG.fa', 'fasta')
os.system('weblogo -f utr_10_AUG.fa -A dna -S 0.5 --errorbars NO --fineprint weglogo -o cap.eps')


