import sys

"""
Author : Zhu Sitao ,zhusitao@genomics.cn
Version: 1.0
Date   : 2018-2-1
Aim    : It's design to get a fasta total length 
"""
if len(sys.argv) == 1:
	print "Usage: python totalLength.py xx.fa > outFile"
	print "The program should be running in python2.7"
	sys.exit(1)
fileIn = sys.argv[1]
with open(fileIn,'r') as F:
	seq = ''
	for line in F.readlines():
		line = line.strip()
		if line.startswith(">"):
			ID = line
		else:
			seq += line
	totalLength = len(seq)
	print ">total_length\t%s"%(totalLength)
