#!/usr/bin/python3

import sys,os,re
import argparse
from itertools import groupby
sys.path.append('/Users/zhusitao/Flow/ngs-tools/format')
from gff import GFF
from fasta import Fasta
from fastq import Fastq

"""
Author: Zhu Sitao
Date: 2018-4-3
Dset: get out all the fishes which match to the given baits 
"""



def readGff(baitPath,Bait:dict):
	""" read gff file"""
	pattern = re.compile(r'(ID|Parent)=([^;]+);*')
	gff = GFF(baitPath)
	for line in gff.readGFF():
		line = line.strip()
		if line.startswith("#"):continue
		match = pattern.search(line.split('\t')[8])
		if match:
			Bait[match.group(2)] = 1

def outputGff(fishPath,Bait,Except):
	""" output gff file from fish pool """
	pattern = re.compile(r'(ID|Parent)=([^;]+);*')
	gff = GFF(fishPath)
	for line in gff.readGFF():
		line = line.strip()
		if line.startswith("#"):
			annotation = line
		else:
			geneInfo = line.split()[8]
			matched = pattern.match(geneInfo)
			geneID = matched.group(2)
			if not Except:
				if Bait.get(geneID):
					print (line)
			else:
				if not Bait.get(geneID):
					print (line)


def readFasta(baitPath,Bait):
	""" read fasta file and put into a dict """
	fasta = Fasta(baitPath).readFasta()
	for ID in fasta.keys():
		Bait[ID] = 1

def outputFasta(fishPath,Bait,Except):
	""" output fasta from fish pool """
	fasta = Fasta(fishPath).readFasta()
	for ID,record in fasta.items():
		line_out = record.fasta_parse()
		if not Except:
			if Bait.get(ID):
				print(line_out.strip())
		else:
			if not Bait.get(ID):
				print(line_out.strip())


def readFastq(baitPath,Bait):
	bait_fastq = Fastq(baitPath)
	fqdict = bait_fastq.fastq_to_dict()
	for key in fqdict:
		Bait[key] = 1

def outputFastq(fishPath,Bait,Except):
	fish_fastq = Fastq(fishPath)
	fqdict = fish_fastq.fastq_to_dict()
	for ID,record in fqdict.items():
		line_out = record.fastq_parse()
		if not Except:
			if Bait.get(ID):
				print(line_out.strip())
		else:
			if not Bait.get(ID):
				print(line_out.strip())


def readTable(baitPath,Bait,bait_column,separator="\t"):
	bait_table = {}
	with open(baitPath,'r') as F:
		for line in F.readlines():
			line = line.strip()
			tmp = line.split(sep=separator)
			ID = tmp[bait_column-1]
			bait_table[ID] = 1

def outputTable(fishPath,Bait,fish_column,Except,separator="\t"):
	fish_table = {}
	with open(fishPath,'r') as F:
		for line in F.readlines():
			line = line.strip()
			tmp = line.split(sep=separator)
			ID = tmp[fish_column-1]
			fish_table[ID] = line
	for ID in fish_table.keys():
		if not Except:
			if Bait.get(ID):
				print(fish_table[ID])
		else:
			if not Bait.get(ID):
				print(fish_table[ID])





def main():
	USAGE = """
		python fish.py -bf fasta -bc 1 -ff fasta -fc 1 bait.fa fish.fa > result.fa    # bait.fa in fish.fa
		python fish.py -bf fasta -ff fasta -E bait.fa fish.fa > result.fa # bait.fa not in fish.fa
		python fish.py -bf gff -ff gff 
		"""
	parser = argparse.ArgumentParser(description='The program use bait to get wanted fish!')
	parser.add_argument('-bf', '--bait_format', help='bait file formate [table,fasta,gff,fq]', action='store')
	parser.add_argument('-ff', '--fish_format', help='fish file formate [table,fasta,gff,fq]', action='store')
	parser.add_argument('-bc', '--bait_column',
						help='bait column for extract, default 1st column [-bc 1 or --bait_column 1]', action='store',
						type=int, default=1)
	parser.add_argument('-fc', '--fish_column',
						help='fish column for extract, default 1st column [-fc 1 or --fish_column 1]', action='store',
						type=int, default=1)
	parser.add_argument('-E', '--Except', help='if -e means contain in fish pool,or not', action='store_true',
						default=False)
	parser.add_argument('bait', help='bait file [table,fasta,gff,fq]')
	parser.add_argument('fish', help='fish file [table,fasta,gff,fq]')
	parser.add_argument('--version', action='version', version='%(prog)s 1.0')

	args = parser.parse_args()
	print('args', args)
	print('bait', args.bait)
	print('fish', args.fish)
	print('Except', args.Except)
	Except = args.Except

	Bait = {}
	if args.bait_format == 'gff' or (args.bait_format == None and args.bait.endswith('.gff')):
		a = readGff(args.bait, Bait)
		outputGff(args.fish, Bait,Except)

	elif args.bait_format == 'fasta' or (args.bait_format == None and args.bait.endswith('.fa')):
		readFasta(args.bait, Bait)
		outputFasta(args.fish, Bait,Except)
	elif args.bait_format == 'fastq' or (args.bait_format == None and args.bait.endswith('.fq')):
		readFastq(args.bait, Bait)
		outputFastq(args.fish, Bait,Except)
	else:
		readTable(args.bait, Bait,1)
		outputTable(args.fish, Bait,1,Except)

if __name__ == '__main__':
	main()