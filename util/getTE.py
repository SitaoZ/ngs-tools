#!/usr/bin/python 
import sys,re
from itertools import groupby
from collections import defaultdict

"""
Author : Zhu Sitao, zhusitao@genomics.cn
Version: 1.0, Date : 2018-1-31
Python version: python2.7
"""

def tree(): 
	""" create a tree dataSturcture, it's smarter than dict """
	return defaultdict(tree)

def readGff(fileIn,ref):
	""" read gff file , return a tree containing total gff information """
	pattern = re.compile("^ID=([^;]+);*")
	with open(fileIn,'r') as F:
		for line in F.readlines():
			line = line.strip()
			if not len(line) or line.startswith("#"):
				continue
			lineArray = line.split()
			tname = lineArray[0]
			match = pattern.search(lineArray[8])
			if match:
				qname = match.group(1)
				ref[tname][qname]=[lineArray[3],lineArray[4],lineArray[6]]
	return ref

def complementReverse(sequence):
	""" complement and reverse a special DNA sequence """
	sequence = sequence.upper()
	sequence = sequence.replace('A', 't')
	sequence = sequence.replace('T', 'a')
	sequence = sequence.replace('C', 'g')
	sequence = sequence.replace('G', 'c')
	sequence = sequence[::-1]
	sequence = sequence.upper()
	return sequence
	
def displaySeq(seqIn,number):
	seqOut = ''
	for i in range(1,len(seqIn)/number+1):
		start = (i-1)*number
		end = i*number
		seqOut += seqIn[start:end]+"\n"
	left = len(seqIn) % number
	remainder = seqIn[-left:]
	return seqOut + remainder
	

class FastaFile:
	def __init__(self,path):
		self.path = path
		self._map = {}
		self.__fasta_iter()
	def __str__(self):
		return self._map.__str__()
	def __fasta_iter(self):
		fh = open(self.path)
		faiter =(x[1] for x in groupby(fh,lambda line : line[0] == ">")) # lambda tells to use the first item(line[0] == ">") as the grouping key
		for header in faiter:
			header = header.next()[1:].strip()
			header = header.split()[0]
			seq ="".join(s.strip() for s in faiter.next())
			self._map[header] = seq

def usage():
	USAGE = """
Ddescription:
	getTE.py  --  get TE elements from sequences according to coordinates
	The program should be running in python2.7
Example:
	python getTE.py chr01.gff chr01.fa >  chr01.TE.fa	
		"""
	print USAGE

if __name__ == '__main__':
	if len(sys.argv) == 1:
		usage()
		sys.exit(1)
	pos_file = sys.argv[1]
	seq_file = sys.argv[2]
	Element = tree()
	Element = readGff(pos_file,Element)
	Fasta=FastaFile(seq_file)._map
		
	for chrom in sorted(Fasta.keys()):
		seq = Fasta[chrom]
		chr_pp = Element[chrom]
		for TE_id in sorted(chr_pp.keys()):
			TE_p = chr_pp[TE_id]
			TE_seq = seq[int(TE_p[0])-1:int(TE_p[1])]
			if TE_p[2] == "-":
				TE_seq = complementReverse(TE_seq)
			TE_seqOut = displaySeq(TE_seq,50)
			output = ">%s seq:%s:%s:%s:%s\n%s" %(TE_id,chrom,TE_p[0],TE_p[1],TE_p[2],TE_seqOut)
			print output
	

