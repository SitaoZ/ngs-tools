import re
import argparse
import pandas as pd
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(description='get gene transcript mapping file')

    parser.add_argument('--gff', '-g', type=str,
                        help='ggnome annotion file gff')
    parser.add_argument('--output', '-o', type=str, help='output file')

    return parser.parse_args()

def main(gff, output):
    tp = re.compile(r'ID=transcript:(ENST\d+);')
    gp = re.compile(r'Parent=gene:(ENSG\d+);')
    gene2tr = defaultdict(list)
    assert('gz' not in gff)
    with open(gff,'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#'):
                continue
            a = line.split('\t')
            ftype, desc = a[2], a[8]
            if ftype == 'mRNA':
                matcht = tp.search(desc)
                matchg = gp.search(desc)
                if matcht and matchg :
                    transcript = matcht.group(1)
                    gene = matchg.group(1)
                    #print(transcript, gene)
                    if gene2tr.get(gene):
                        gene2tr[gene].append(transcript)
                    else:
                        gene2tr[gene] = [transcript]
    for gene in gene2tr.keys():
        print(*[gene,','.join(gene2tr[gene])], sep='\t', file=output)


if __name__ == '__main__':
    args = parse_args()
    main(args.gff, args.output)
