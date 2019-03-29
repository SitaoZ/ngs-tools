import re
from fasta import Fasta
"""
Author: Zhu Sitao
Date: 2018-3-21
Dest: The readGFF.py is use to read GFF3,Python version is 3.6
GFF is a standard file format for storing genomic features in a text file. 
GFF stands for Generic Feature Format. GFF files are plain text, 9 column, tab-delimited files. 
GFF databases also exist. They use a schema custom built to represent GFF data.
GFF is frequently used in GMOD for data exchange and representation of genomic data.
http://gmod.org/wiki/GFF3
"""

class GFF(object):
	""" Read GFF(general feature format) """
	def __init__(self,gffPath):
		self._gff = dict()
		self._filePath = gffPath
		self._class = none
		self._mRNA = dict()
		self._CDS = dict()

		self.geneid_pattern = re.compile(r'geneID=(?P<id>\S+);?')
		self.mRNA_pattern = re.compile(r'ID=(?P<id>\S+);')
		self.cds_pattern = re.compile(r'Parent=(?P<id>\S+);?')


	def __iter__(self):
		""" Supports traversal with a for loop"""
		return iter(self._gff)
	def readGFF(self):
		""" Return gff file line by line through generator """
		fh = open(self._filePath)
		for line in fh:
			yield line
	def annotation_chr(self):
		""" Return a chromosome list """
		chr_list = []
		for line in self.readGFF():
			if line.startswith("#"):
				continue
			else:
				chr_list.append(line.strip().split()[0])
		chromosomes = list(set(chr_list))
		chromosomes.sort(key=chr_list.index)
		return chromosomes

	def annotation_source(self):
		""" Return annotation source class list """
		source_list = []
		for line in self.readGFF():
			if line.startswith("#"):
				continue
			else:
				source_list.append(line.strip().split()[1])
		source = list(set(source_list))
		source.sort(key=source_list.index)
		return source

	def annotation_type(self):
		""" Return annotation type list """
		type_list = []
		for line in self.readGFF():
			if line.startswith("#"):
				continue
			else:
				type_list.append(line.strip().split()[2])
		types = list(set(type_list))
		types.sort(key=type_list.index)
		return types


	def geneID(self):
		""" Return gene id list """
		geneid = []
		for line in self.readGFF():
			line = line.strip()
			match = self.geneid_pattern.search(line)
			if match:
				id = match.group("id")
				geneid.append(id)
		geneid_uniq = list(set(geneid))
		geneid_uniq_sort = sorted(geneid_uniq)
		return geneid_uniq_sort
	def geneCount(self):
		""" Return gene total number """
		geneid = self.geneID()
		return len(geneid)

	def mrnaID(self):
		""" Return transcript id list """
		mrnaid = []
		for line in self.readGFF():
			line = line.strip()
			match = self.mRNA_pattern.search(line)
			if match:
				id = match.group("id")
				mrnaid.append(id)
		mrnaid_uniq = list(set(mrnaid))
		mrnaid_uniq_sort = sorted(mrnaid_uniq)
		return mrnaid_uniq_sort
	def mrnaCount(self):
		mrnaid = self.mrnaID()
		return len(mrnaid)

	def exonCount(self):
		exon = 0
		for line in self.readGFF():
			if line.startswith("#"):
				continue
			else:
				if line.strip().split()[2].lower() == "exon":
					exon += 1
		return exon


	def mRNA_fasta(self,genomefasta_path,mRNAfasta_path):
		genome = Fasta(genomefasta_path)
		genome.readFasta()
		genomeDict = genome._fasta
		forward = {}
		for line in self.readGFF():
			if line.startswith("#"):
				continue
			array = line.strip().split('\t')
			scaff,program,structure,start,end,orient = array[0],array[1],array[2],int(array[3]),int(array[4]),array[6]
			if structure == "mRNA":
				key = self.mRNA_pattern.search(line).group("id")
				seq = genomeDict[scaff][start-1:end]
				if orient == "+":
					forward[key] = seq
				else:
					seq = seq.replace('A','{A}').replace('T','{T}').replace('C','{C}').replace('G','{G}')
					seq = seq.format(A='T', T='A', C='G', G='C')[::-1]
					forward[key] =seq

		mRNAfasta = open(mRNAfasta_path,'w')
		for key in forward:
			mRNAfasta.writelines(">"+key+"\n")
			mRNAfasta.writelines(forward[key]+"\n")
		mRNAfasta.close()





	def CDS_fasta(self,genomefasta_path,cds_outpath):
		genome = Fasta(genomefasta_path)
		genome.readFasta()
		genomeDict = genome._fasta
		forward = {}
		reverse = {}
		for line in self.readGFF():
			if line.startswith("#"):
				continue
			array = line.strip().split('\t')
			scaff,program,structure,start,end,orient = array[0],array[1],array[2],int(array[3]),int(array[4]),array[6]
			if structure == "mRNA":
				seq2 = ''
				seq3 = ''
			elif structure == "CDS":
				key = self.cds_pattern.search(line).group("id")
				out = genomeDict[scaff][start-1:end]
				seq2 += out
				seq3 += out
				if orient == "+":
					forward[key] = seq2
				else:
					reverse[key] = seq3
		for key in reverse.keys():
			seq = reverse[key]
			seq = seq.replace('A','{A}').replace('T','{T}').replace('C','{C}').replace('G','{G}')
			seq = seq.format(A='T', T='A', C='G', G='C')[::-1]
			forward[key] = seq
		CDS = open(cds_outpath,'w')
		for key in sorted(forward.keys()):
			CDS.writelines(">"+key+"\n")
			j = 0
			out = ''
			for j in range(0, len(forward[key]), 60):
				out = forward[key][j:j+60]
				CDS.write(out+"\n")
		CDS.close()




