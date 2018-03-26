import re
from itertools import groupby

"""
Author: Zhu Sitao
Date : 2018-3-21
Dest: a DNA fasta class; used in python 3.6
"""

class ReadFasta(object):
	def __init__(self,filePath):
		self._fasta = dict()
		self.path = filePath
		self.readFasta()
	def readFasta(self):
		""" Read Fasta file and load in a dict """
		fh = open(self.path)
		faiter =(x[1] for x in groupby(fh,lambda line : line[0] == ">"))
		for header in faiter:
			header = header.__next__()[1:].strip()
			header = header.split()[0]
			seq ="".join(s.strip() for s in faiter.__next__())
			self._fasta[header] = seq.upper()
	def __str__(self):
		""" Print the fasta dict """
		return self._fasta.__str__()
	def show(self, number = 60):
		""" Printf sequence """
		for ID in self._fasta:
			seqIn = self._fasta[ID]
			seqOut = ''
			for i in range(1,len(seqIn)//number+1):
				start = (i-1)*number
				end = i*number
				seqOut += seqIn[start:end]+"\n"
			left = len(seqIn) % number
			remainder = seqIn[-left:] if left != 0 else ""
			print(">"+ID+"\n"+seqOut + remainder)
	def GC(self):
		""" Show each fasta GC and GC rate"""
		for ID in self._fasta:
			seqIn = self._fasta[ID]
			GC = seqIn.count("G") + seqIn.count("C")
			GCrate = GC/len(seqIn)
			print ("%s\t%d\t%.2f" % (ID,GC,GCrate))
		
	def __len__(self):
		""" Total Length of the fasta file """
		totalLength = 0
		for ID in self._fasta:
			totalLength += len(self._fasta[ID])
		return totalLength
	def stat(self):
		""" Statistic of each id length """
		for ID in self._fasta:
			print (ID,len(self._fasta[ID]))
	def __iter__(self):
		""" Supports traversal with a for loop"""
		return iter(self._fasta)
	def __getitem__(self,index):
		""" Index one fasta ID """
		return index,self._fasta[index]
