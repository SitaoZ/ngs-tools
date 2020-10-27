# coding:utf-8
import re
import argparse
from collections import defaultdict
from Bio import SeqIO

parser = argparse.ArgumentParser(description="Process some files")

parser.add_argument('-g', '--gff', help='an fasta file for gff3')
parser.add_argument('-t', '--transcript', help='an fasta file for transcript')

args = parser.parse_args()
Dict = dict()

pattern = re.compile(r'ID=gene:(\S+);Name=(\S+);biotype')

# parser gff3 to stare transcript id and gene symbol
# sometimes a specific gene do not have a symbol
with open(args.gff, 'r') as f:
    for line in f.readlines():
        line = line.strip()
        if line.startswith('#'):
            continue
        fields = line.split()
        feature_type = fields[2]
        start, end, strand = fields[3], fields[4], fields[6]
        info = fields[8]

        if feature_type == 'gene':
            matched = pattern.search(info)
            if matched:
                gene_id = matched.group(1)
                symbol_id = matched.group(2).lower()  # lower
                Dict[symbol_id] = gene_id

Fasta = defaultdict(dict)
pattern_c = re.compile(r'CDS=(\d+?)-(\d+?)$')

# info save gene:transcript:[utr5,utr3]
for record in SeqIO.parse(args.transcript, 'fasta'):
    trans_id = record.id.replace('transcript:', '')
    desc = record.description
    seq = record.seq
    matched_c = pattern_c.search(desc)
    if matched_c:
        start, end = map(int, [matched_c.group(1), matched_c.group(2)])
        start = start - 1
        gene_id = trans_id.split('.')[0]
        utr5, utr3 = seq[:start], seq[end:]
        Fasta[gene_id][trans_id] = [utr5, utr3]

with open('sym_list.txt') as f:
    for gene_id in f.readlines():
        gene_id = gene_id.strip()  # lower
        for trans_id in sorted(Fasta[gene_id].keys(), key=lambda x: x.split('.')[1]):
            utr5 = Fasta[gene_id][trans_id][0]
            utr3 = Fasta[gene_id][trans_id][1]
            print(">{}\tgene:{}|type:utr5\n{}".format(trans_id, gene_id, utr5))

with open('./gene', 'r') as f:
    for line in f.readlines():
        gene_id = line.strip()
        for trans_id in sorted(Fasta[gene_id].keys(), key=lambda x: x.split('.')[1]):
            utr5 = Fasta[gene_id][trans_id][0]
            utr3 = Fasta[gene_id][trans_id][1]
            # print(gene_id, trans_id, utr5, utr3)
            # print(">{}\tgene:{}|type:utr5\n{}".format(trans_id, gene_id, utr5))
