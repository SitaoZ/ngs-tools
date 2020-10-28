# coding:utf-8


import argparse
import pandas as pd
import numpy as np
from Bio import SeqIO


parser = argparse.ArgumentParser(description='Process some files.')
parser.add_argument('-i','--input',
                    help='an fasta file for input')
parser.add_argument('-o','--output',
                    help='an fasta file for output')

args = parser.parse_args()

df = pd.DataFrame()
fasta = args.input

j = 0
for record in SeqIO.parse(fasta, 'fasta'):
    j += 1
    trans_id = record.id
    desc = record.description
    seq = record.seq
    gene_id = trans_id.split(',')
    if len(gene_id) == 1:
        gene_id = gene_id[0].split('.')[0]
    else:
        A = []
        for i in gene_id:
            i = i.split('.')[0]
            A.append(i)
        gene_id = np.unique(A)[0]
    data = {
        'geneID': gene_id,
        'transcriptID': trans_id,
        '5UTR_primer': seq
    }
    # df2 = pd.DataFrame([data],index=[str(j)])
    series = pd.Series(data, name=str(j))

    df = df.append(series)

df.to_excel(args.output, header=['geneID', 'transcriptID', '5UTR_primer'])






