import re,os

"""
Author :Zhu Sitao
Date: 2018-3-23
Dest: read bam file by samtools ;samtolls should be installed in your ~/.bashrc
      python version is 3.6
"""


class ReadBam(object):
	""" a readBam class """
	# class global variable , each funcation could use the variable
	patternAnnotation = re.compile(r'^@')

	def __init__(self,filePath):
		self.path = filePath
		self.sam = self.readBam()

	def readBam(self):
		""" Read bam from samtools view """
		samtools = os.popen('which samtools').read().strip()
		sam = os.popen('%s view -h %s' % (samtools,self.path)).read()
		iterm = sam.split("\n")
		itermFilter = [x for x in iterm if x !='']
		return itermFilter
		#for i in iterm:
		#	yield i

	def header(self):
		""" Bam header """
		for line in self.sam:
			if self.patternAnnotation.match(line) and line != '':
				print (line)
	def mapRate(self):
		""" Bam maprate """
		# initialize
		total,duplicate,mapped,unmapped = 0,0,0,0
		for line in self.sam:
			if not self.patternAnnotation.match(line):
				total += 1
				flag = line.split('\t')[1]
				position = line.split('\t')[3]
				if int(flag) & 0x0040: # duplicates
					duplicate += 1 
				elif position != 0:
					mapped += 1
				elif position == 0: # unmapped
					unmapped += 1
		dupRate = "%.4f" % (duplicate/total)
		mapRate = "%.4f" % (mapped/total)
		unMapRate = "%.4f" % (unmapped/total)
		return dupRate,mapRate,unMapRate

	def readLength(self):
		""" Read length """
		for line in self.sam:
			if not self.patternAnnotation.match(line) and line != '':
				readLength = len(line.split('\t')[9])
				return readLength

	def cleanReads(self):
		""" Clean reads """
		cleanReads = 0
		for line in self.sam:
			if line !='' and not self.patternAnnotation.match(line):
				cleanReads += 1
		return cleanReads

	def cleanBases(self):
		""" Clean bases """
		readLength = self.readLength()
		cleanBases = 0
		for line in self.sam:
			if line != '' and not self.patternAnnotation.match(line):
				cleanBases += readLength
		return cleanBases
				
	def depth(self):
		""" Bam deepth """
		pass

	def extract(self,chrName):
		""" extract given ID name bam file """
		for line in self.sam:
			if line != '' and not self.patternAnnotation.match(line):
				chrID = line.split("\t")[2]
				if chrName == chrID:
					print (line)
