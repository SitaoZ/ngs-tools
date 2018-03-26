import re

"""
Author: Zhu Sitao
Date: 2018-3-21
Dest: The readGFF.py is use to read GFF3,Python version is 3.6
"""

class ReadGFF(object):
	""" Read GFF(general feature format) """
	def __init__(self,filePath):
		self._gff = dict()
		self._filePath = filePath
		self.readGFF()
	def __str__(self):
		""" Print GFF dict """
		return self._gff.__str__()
	def __iter__(self):
		""" Supports traversal with a for loop"""
		return iter(self._gff)
	def readGFF(self):
		""" Read gff """
		fh = open(self._filePath)
		pattern = re.compile(r'ID=(.+);Shift')# only match mRNA 
		for line in fh:
			line = line.strip() 
			array = line.split('\t')
			scaff = array[0]
			program = array[1]
			structure = array[2]
			start,end = array[3],array[4]
			orient = array[6]
			matched = pattern.search(line)
			if matched:
				geneID = matched.group(1)
				self._gff[geneID] = [scaff,orient,start,end]
	def geneID(self):
		"""Get gene ID list ---> mRNA"""
		for geneid in self._gff:
			yield geneid
	def geneCount(self):
		"""Get gene total number ---> mRNA"""
		return len(self._gff)
