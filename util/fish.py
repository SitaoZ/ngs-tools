#!/usr/bin/python3

import sys,os,re
import argparse
from itertools import groupby

from a import *

"""
Author: Zhu Sitao
Date: 2018-4-3
Dset: get out all the fishes which match to the given baits 
"""

USAGE = """
	python fish.py -a -b -c -d 8 bait.fa fish.fa > result.fa	
	"""
parser = argparse.ArgumentParser(description='The program use bait to get wanted fish!')
parser.add_argument('-bf','--bait_format',help='bait file formate [table,fasta,gff,fq]',action='store')
parser.add_argument('-ff','--fish_format',help='fish file formate [table,fasta,gff,fq]',action='store')
parser.add_argument('-bc','--bait_column',help='bait column for extract, default 1st column [-bc 1 or --bait_column 1]',action='store',type=int,default=1)
parser.add_argument('-fc','--fish_column',help='fish column for extract, default 1st column [-fc 1 or --fish_column 1]',action='store',type=int,default=1)
parser.add_argument('-E','--Except',help='if -e means contain in fish pool,or not',action='store_true',default=False)
parser.add_argument('bait',help='bait file [table,fasta,gff,fq]')
parser.add_argument('fish',help='fish file [table,fasta,gff,fq]')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')

args = parser.parse_args()
print (args)
print (args.bait)
print (args.Except)

Except = args.Except


def readGff(baitPath,Bait):
	""" read gff file"""
	pattern = re.compile(r'(ID|Parent)=([^;]+);*/)')
	gff = GFF(baitPath)
	for line in gff.readGFF():
		line = line.strip()
		if line.startswith("#"):continue
		match = pattern.search(line.split()[8])
		if match:
			Bait[match.group(2)] = 1




def outputGff(fishPath,Bait):
	""" output gff file from fish pool """
	idPattern = re.compile(r'(ID|Parent)=([^;]+);*')
	with open(fishPath,'r') as F:
		for line in F:
			line = line.strip()
			if line.startswith("#"):
				annotation = line
			else:
				geneInfo = line.split()[8]
				matched = idPattern.match(geneInfo)
				geneID = matched.group(2)
				if not Except:
					if Bait.get(geneID):
						print (line)
				else:
					if not Bait.get(geneID):
						print (line)


def readFasta(baitPath,Bait):
	""" read fasta file and put into a dict """
	fasta = Fasta(baitPath)
	for ID in fasta._fasta:
		pass
		
		
		

if __name__ == '__main__':
	Bait = dict()
	if args.bait_format == 'gff' or (args.bait_format == None and args.bait.endswith('.gff')):
		readGff(args.bait,Bait)
		outputGff(args.fish,Bait)
	elif args.bait_format == 'fasta' or (args.bait_format == None and args.bait.endswith('.fa')):
		readFasta(args.bait,Bait)
		outputFasta(args.fish,Bait)
	elif args.bait_format == 'fq' or (args.bait_format == None and args.bait.endswith('.fq')):
		readFq(args.bait,Bait)
		outputFq(args.fish,Bait)
	else:
		readTable(args.bait,Bait)
		outputTable(args.fish,Bait)
