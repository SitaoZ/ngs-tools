import re,os
from random import sample
from itertools import groupby
from collections import defaultdict

"""
Author: Zhu Sitao
Date : 2018-3-21
Dest: a DNA fasta class; used in python 3.6
"""

class Fasta(object):
	""" Fasta class """
	fastaName = "zhusitao"
	BASES = ['T', 'C', 'A', 'G']
	CODONS = [a + b + c for a in BASES for b in BASES for c in BASES]
	AMINO_ACIDS = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
	CODON_TABLE = dict(zip(CODONS, AMINO_ACIDS))
	
	def __init__(self,filePath):
		self._fasta = dict()
		if isinstance(filePath,str):
			self.path = filePath
		self.readFasta()
	def __str__(self):
		""" Print the fasta dict """
		return self._fasta.__str__()
	def __iter__(self):
		""" Supports traversal with a for loop ,for ID loop """
		return iter(self._fasta)
	def __getitem__(self,index):
		""" Index one fasta ID """
		outPut = self.fasta_parse(index,self._fasta[index],number = 50)
		return outPut
	def __len__(self):
		""" Total sequence length of the fasta file """
		totalLength = 0
		for ID in self._fasta:
			totalLength += len(self._fasta[ID])
		return totalLength

	def readFasta(self):
		""" Read Fasta file and load in a dict ,normal method
			使用groupby 将文本文件做成一个生成器，生成器没有把所有值存在内存中，而是在运行时生成的值，可以快速访问大文件。
			生成器你只能对其迭代一次。
		"""
		fh = open(self.path)
		faiter =(x[1] for x in groupby(fh,lambda line : line[0] == ">"))
		for header in faiter:
			header = header.__next__()[1:].strip() # [1:] 为了去除 > 符号
			header = header # header.split()[0] 对名称进行简化，当前的做法是保存全部名称
			seq ="".join(s.strip() for s in faiter.__next__()) 
			self._fasta[header] = seq.upper()

	def fasta_key(self):
		""" return a list of fasta name"""
		return self._fasta.keys()

	def fasta_sequence(self):
		""" return a list of sequence"""
		return self._fasta.values()

	def item_count(self):
		""" return numbers fasta names"""
		item_total_number = len(self.fasta_key())
		return item_total_number

	def std_out(self, number = 60):
		""" printf sequence to stdout """
		for ID in self.fasta_key():
			seqIn = self._fasta[ID]
			print(self.fasta_parse(ID,seqIn,number = 50))
	def gc_rate(self,output_txt):
		""" Show each fasta GC and GC rate"""
		Stat = open(output_txt,'w')
		Stat.writelines("ID\tGC\tGCrate(%)\tN\tNrate(%)\n")
		totalGC,totalN = 0,0
		for ID in self._fasta:
			seqIn = self._fasta[ID]
			GC = seqIn.count("G") + seqIn.count("C")
			N = seqIn.count("N")
			totalGC += GC
			totalN += N
			GCrate = GC/len(seqIn)
			Nrate = N/len(seqIn)
			Stat.writelines("%s\t%d\t%.4f\t%d\t%.4f\n" % (ID,GC,GCrate,N,Nrate))
		totalGCrate = totalGC/len(self)
		totalNrate = totalN/len(self)
		Stat.writelines("total\t%d\t%.4f\t%d\t%.4f"%(totalGC,totalGCrate,totalN,totalNrate))
		Stat.close()

	def stat_length(self):
		""" Statistic of each id length """
		#LenDict = defaultdict(int)
		LenDict = dict()
		for ID in self._fasta:
			LenDict[ID]=len(self._fasta[ID])
		return LenDict

	def fasta_parse(self,index,seqIn,number = 60):
		""" a parse for a fasta item """
		seqOut = ''
		for i in range(1,len(seqIn)//number+1):
			start = (i-1)*number
			end = i*number
			seqOut += seqIn[start:end]+"\n"
		left = len(seqIn) % number
		remainder = seqIn[-left:] if left != 0 else ""
		totalOut = ">"+index+"\n"+seqOut + remainder+ "\n"
		return totalOut

	def reverse(self,seq):
		"""reverse fasta seq"""
		return seq[::-1]

	def complement(self,seq):
		""" complement seq """
		COMPLEMENT_TRANS = str.maketrans('TAGCtagc', 'ATCGATCG')
		return seq.translate(COMPLEMENT_TRANS)

	def reverse_complement(self,seq):
		""" return reverse and complement seq"""
		return self.reverse(self.complement(seq))

	def translate(self,seq):
		seq = seq.lower().replace('\n', '').replace(' ', '')
		peptide = ''
		for i in range(0, len(seq), 3):
			codon = seq[i: i + 3]
			amino_acid = CODON_TABLE.get(codon, '!')
			if amino_acid != '!':  # end of seq
				peptide += amino_acid
		return peptide

	@classmethod
	def basename(cls):
		""" class method """
		return cls.fastaName
	
	@staticmethod
	def print_working_directory():
		""" static method ,no default parameter """
		return os.getcwd()


	def _max_min(self,type):
		sort_dict = sorted(self.stat_length().items(), key=lambda d: d[1], reverse=True)
		if type == 'max':
			max_item = sort_dict[0]
			return max_item
		elif type == 'min':
			min_item = sort_dict[self.item_count()-1]
			return min_item
		else:
			return "%s not max or min"%type

	def extract_item(self,type,outfile):
		""" extract max or min length fasta item"""
		if type.lower() not in ['max','min']:
			raise KeyError(str(type)+" not contain,should be min or max")
		with open(outfile,'w') as OUT:
			max_min_item = self._max_min(type)
			key, value = max_min_item
			if type == 'max':
				OUT.writelines(self.fasta_parse(key,self._fasta[key]))
			elif type == 'min':
				OUT.writelines(self.fasta_parse(key, self._fasta[key]))
			else:
				print ("%s not correct,should be max or min"%type)
	def random_sample(self,number,outpath):
		""" ramdom sample from the all fasta file"""
		OUT = open(outpath,'w')
		key_list = self.fasta_key()
		sample_list = sample(key_list,number)
		for key in sample_list:
			OUT.writelines(self.fasta_parse(key,self._fasta[key]))
		OUT.close()


