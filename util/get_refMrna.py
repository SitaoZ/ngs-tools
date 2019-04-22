import re
import sys
from itertools import groupby

def read_fasta(genome_path):
        fasta = {}
        try:
                fh = open(genome_path)
        except KeyError:
                print("The path is a compressed file")
        faiter =(x[1] for x in groupby(fh,lambda line : line[0] == ">"))
        for header in faiter:
                header = header.__next__()[1:].strip() # [1:] 为了去除 > 符号
                header = header # header.split()[0] 对名称进行简化，当前的做法是保存全部名称
                seq ="".join(s.strip() for s in faiter.__next__())
                fasta[header] = seq.upper()
        return fasta
def reverse(seq):
        return seq[::-1]
def complement(seq):
        """ complement seq """
        COMPLEMENT_TRANS = str.maketrans('TAGCtagc', 'ATCGATCG')
        return seq.translate(COMPLEMENT_TRANS)
def read_gff(gff_path):
        pattern = re.compile(r'ID=(?P<name>\S+);')
        mRNA={}
        with open(gff_path,'r') as F:
                for line in F.readlines():
                        fields = line.strip().split()
                        chrom,c_type,start,end,strand,info= fields[0],fields[2],int(fields[3]),int(fields[4]),fields[6],fields[8]
                        match = pattern.search(info)
                        if c_type == 'mRNA' and match:
                                mRNA[match.group('name')] = [chrom,strand,start,end]
        return mRNA


if __name__=="__main__":
        genome,gff = sys.argv[1:]
        fasta = read_fasta(genome)
        mRNA = read_gff(gff)
        for key in mRNA.keys():
                chrom,strand,start,end = mRNA[key]
                if strand == '-':
                        seq = fasta[chrom][start-1:end]
                        print(">%s"%(key))
                        print(complement(reverse(seq)))
                else:
                        print(">%s"%(key))
                        print(fasta[chrom][start-1:end])
