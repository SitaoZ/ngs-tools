import sys
from format.fasta import Fasta

"""
Author:Zhu Sitao
Date : 2018-3-28
Dest : extract Fasta 
"""
InputFile = sys.argv[1]
IdList = sys.argv[2]

fastaFile = Fasta(InputFile)
ID = open(IdList,'r')
for line in ID:
	line = line.strip()
	if len(line) != 0:
		fastaFile[line]
ID.close()