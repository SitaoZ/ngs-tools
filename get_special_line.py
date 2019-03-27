import re,sys
import bz2file

"""
Author: Zhu Sitao
Date: 2018-3-27
Dest : get given line from bz2 compress file
"""

class FileType(object):
	def __init__(self,filepath,fileformat):
		self.path = filepath
		self.format = fileformat

	def bz2(self):
		""" bz2 type """
		if self.format = "bz2":
			with bz2file.open(self.path,'r') as F:
				for line in F:
					line = line.decode('utf-8')
					print (line)
	def gz(self):
		pass
	
