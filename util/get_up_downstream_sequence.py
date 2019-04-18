#encoding:utf-8
"""
Author:Sitao Zhu
Date:2019-4-18
Desc: find promoter sequence,running on python3 platform
"""
import re
import sys
import getopt
from itertools import groupby


def usage():
	print (""" usage: python %s -r genome_path -g genelist -f gff3 """%(sys.argv[0]))
	
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
def get_position(gff,geneid):
	""" get gene position """
	pattern = re.compile(r'ID=evm.TU.(?P<name>\S+);')
	with open(gff,'r') as F:
		for line in F.readlines():
			if line.startswith('#'):
				continue
			else:
				fileds = line.strip().split()
				chrom,c_type,start,end,geneinfo = fileds[0],fileds[2],fileds[3],fileds[4],fileds[8]	
				if c_type == 'gene':
					if pattern.search(geneinfo).group('name') == geneid:
						return pattern.search(geneinfo).group('name'),chrom,int(start),int(end)
if __name__=="__main__":
	opts,args = getopt.getopt(sys.argv[1:],'-hr:g:f:',["help","reference","genelist","gff3"])
	for name, value in opts:
		if name in ("-h", "--help"):
			usage()
			sys.exit()
		if name in ("-r","--reference"):
			refer = value
		if name in ("-g","--genelist"):
			genelist = value
		if name in ("-f","--gff3"):
			gff3 = value
	fasta = read_fasta(refer)
	with open(genelist,'r') as F:
		for line in F:
			genename = line.strip()
			geneid,chrom,start,end = get_position(gff3,genename)
			point1 = start - 1001
			point2= start - 1
			upstream = fasta[chrom][point1:point2]
			print (">{geneid}\t{chrom}\t{start}\t{end}\tupstream\n{pro}".format(geneid=geneid,chrom=chrom,start=point1,end=point2,pro=upstream))
			point3 = end + 1001
			point4 = end + 1
			downstream = fasta[chrom][point4:point3]
			print (">{geneid}\t{chrom}\t{start}\t{end}\tdownstream\n{pro}".format(geneid=geneid,chrom=chrom,start=point4,end=point3,pro=downstream))	
