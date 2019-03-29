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
	
	def __init__(self):
		self._fasta = dict()

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
		""" return fasta item length """
		return len(self._fasta.keys())
	def __contains__(self, item):
		""" return true if a fasta class contains the item,or false"""
		if item in self._fasta.keys():
			return True
		else:
			return False
	def __eq__(self, other):
		""" Return True if self equals other or False otherwise """
		if self is other: return True
		if type(self) != type(other) or len(self) != len(other):
			return False
		for item in self:
			if item not in other:
				return False
		return True


	def __add__(self, other):
		""" return a new fasta contains the contents of self and other"""
		result = Fasta()
		for item in other:
			result._fasta[item] = other._fasta[item]
		for item in self:
			result._fasta[item] = self._fasta[item]
		return result




	def readFasta(self,path):
		""" Read Fasta file and load in a dict ,normal method
			使用groupby 将文本文件做成一个生成器，生成器没有把所有值存在内存中，而是在运行时生成的值，可以快速访问大文件。
			生成器你只能对其迭代一次。
		"""

		try:
			fh = open(path)
		except KeyError:
			print("The path is a compressed file")
		faiter =(x[1] for x in groupby(fh,lambda line : line[0] == ">"))
		for header in faiter:
			header = header.__next__()[1:].strip() # [1:] 为了去除 > 符号
			header = header # header.split()[0] 对名称进行简化，当前的做法是保存全部名称
			seq ="".join(s.strip() for s in faiter.__next__()) 
			self._fasta[header] = seq.upper()

	def fasta_key_list(self):
		""" return a list of fasta name"""
		return self._fasta.keys()

	def fasta_sequence_list(self):
		""" return a list of sequence"""
		return self._fasta.values()


	def fasta_key_count(self):
		""" return numbers fasta names"""
		item_total_number = len(self.fasta_key_list())
		return item_total_number

	def fasta_sequence_length(self):
		""" return all fasta length """
		totalLength = 0
		for ID in self._fasta:
			totalLength += len(self._fasta[ID])
		return totalLength

	def std_out(self, line_number = 60):
		""" printf sequence to stdout """
		for ID in self.fasta_key_list():
			seqIn = self._fasta[ID]
			print(self.fasta_parse(ID,seqIn,number = line_number))
	def gc_rate(self,output_txt):
		""" Write each fasta GC and GC rate to a file """
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

	def stat_item_length(self):
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

	def reverse(self):
		"""reverse fasta seq"""
		for line in self.readFasta():
			pass

	@classmethod
	def basename(cls):
		""" class method """
		return cls.fastaName
	
	@staticmethod
	def print_working_directory():
		""" static method ,no default parameter """
		return os.getcwd()


	def _max_min(self,type):
		sort_dict = sorted(self.stat_item_length().items(), key=lambda d: d[1], reverse=True)
		if type == 'max':
			max_item = sort_dict[0]
			return max_item
		elif type == 'min':
			min_item = sort_dict[self.fasta_key_count()-1]
			return min_item
		else:
			return "%s not max or min"%type

	def max_min_length(self,type,outfile):
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

	def extract_item(self,key):
		""" extract a given fasta item """
		if key in self.fasta_key_list():
			result = self.fasta_parse(key,self._fasta[key])
		else:
			raise KeyError(str(key)+" not contain,should be min or max")

	def random_sample(self,number,outpath):
		""" ramdom sample from the all fasta file"""
		OUT = open(outpath,'w')
		sample_list = sample(self.fasta_key_list(),number)
		for key in sample_list:
			OUT.writelines(self.fasta_parse(key,self._fasta[key]))
		OUT.close()

    def gc(self,binsize = 500):
        """
        %prog gc fastafile
        Plot G+C content distribution.
        """
        allbins = []
        for name, seq in self._fasta.items():
            for i in range(len(seq) // binsize):
                atcnt = gccnt = 0
                for c in seq[i * binsize: (i + 1) * binsize].upper():
                    if c in "AT":
                        atcnt += 1
                    elif c in "GC":
                        gccnt += 1
                totalcnt = atcnt + gccnt
                if totalcnt == 0:
                    continue
                gcpct = gccnt * 100 // totalcnt
                allbins.append(gcpct)

        from jcvi.graphics.base import asciiplot
        from collections import Counter

        title = "Total number of bins={}".format(len(allbins))
        c = Counter(allbins)
        x, y = zip(*sorted(c.items()))
        asciiplot(x, y, title=title)
