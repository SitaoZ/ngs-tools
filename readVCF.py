import re,os
import sys

class ReadVCF(object):
	""" VCF class """
	#patterm Annoation line
	
	def __init__(self,vcfPath):
		self.VCF = dict()
		self.path = vcfPath
		self.info_pattern = re.compile(r'''\#\#INFO=<
            		ID=(?P<id>[^,]+),\s*
            		Number=(?P<number>-?\d+|\.|[AGR])?,\s*
            		Type=(?P<type>Integer|Float|Flag|Character|String),\s*
            		Description="(?P<desc>[^"]*)"
            		(?:,\s*Source="(?P<source>[^"]*)")?
            		(?:,\s*Version="?(?P<version>[^"]*)"?)?
            		>''', re.VERBOSE)
		self.filter_pattern = re.compile(r'''\#\#FILTER=<
            		ID=(?P<id>[^,]+),\s*
            		Description="(?P<desc>[^"]*)"
            		>''', re.VERBOSE)
		self.alt_pattern = re.compile(r'''\#\#ALT=<
            		ID=(?P<id>[^,]+),\s*
            		Description="(?P<desc>[^"]*)"
            		>''', re.VERBOSE)
		self.format_pattern = re.compile(r'''\#\#FORMAT=<
  			ID=(?P<id>.+),\s*
            		Number=(?P<number>-?\d+|\.|[AGR]),\s*
            		Type=(?P<type>.+),\s*
            		Description="(?P<desc>.*)"
            		>''', re.VERBOSE)
		self.contig_pattern = re.compile(r'''\#\#contig=<
            		ID=(?P<id>[^>,]+)
            		(,.*length=(?P<length>-?\d+))?
            		.*
            		>''', re.VERBOSE)
		self.head_pattern = re.compile(r'\#CHROM\s*POS',re.VERBOSE)
		self.meta_pattern = re.compile(r'''##(?P<key>.+?)=(?P<val>.+)''')

	def readVCF(self):
		if self.path.startswith("gz"):
			sys.exit(1)
		with open(self.path,'r') as F:
			for line in F:
				yield line

	def __iter__(self):
		return iter(self.VCF)

	def __str__(self):
		return self.VCF.__str__()

	def __len__(self):
		totalLength = 0
		for chrom in self.VCF:
			totalLength += 1
		return totalLength
	
	def getChromLength(self):
		chromLenDict = dict()
		for line in self.readVCF():
			if len(line) !=0 and line.startswith("#"):
				match = self.contig_pattern.match(line)
				if match:
					chromID = match.group('id')
					chromLen = match.group('length')
					chromLenDict[chromID] = chromLen
		return chromLenDict


	def getReferPath(self):
		referPatter = re.compile(r'##reference=file:\/\/(\S+)')
		for line in self.readVCF():
			if len(line) != 0 and line.startswith("#"):
				if referPatter.search(line):
					print(referPatter.search(line).group(1))
					break
	def getReferLength(self):
		total_length = 0
		chromLenDict = self.getChromLength()
		for key in chromLenDict.keys():
			total_length += int(chromLenDict[key])
		return total_length
	
	def getSampleID(self):
		"#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  PE100-2"
		sampleList = []
		for line in self.readVCF():
			match = self.head_pattern.match(line)
			if match:
				temp = line.strip().split()
				for sample in range(9,len(temp)):
					sampleList.append(temp[sample])
		return sampleList
					
	
	def Homo_Hete_Sites_Stat(self):
		""" state variation types and print the result"""
		HomozygousSitesStat = 0
		heterozygousSitesStat = 0
		missingSitesStat = 0
		nonSitesStat = 0

		for line in self.readVCF():
			if len(line) != 0 and not line.startswith("#"):
				ale,genotype = line.strip().split("\t")[3],line.split("\t")[9]
				GT = genotype.split(":")[0].split("/")
				if GT[0] == "." and GT[1] == ".":
					nonSitesStat += 1
				elif "*" in GT[0] or "*" in GT[1]:
					missingSitesStat += 1
				## homo
				else:
					if GT[0] == GT[1]:
						HomozygousSitesStat += 1
					else:
						heterozygousSitesStat +=1
		print ("Homozygous   Sites:",HomozygousSitesStat)
		print ("heterozygous Sites:",heterozygousSitesStat)
		print ("missingSites(contain *):",missingSitesStat)
		print ("nonSites(contain .):",nonSitesStat)

	def getHomoSites(self,HomoSites_outpath):
		"""get homo sites , including header line and variation line"""
		HomoSites =open(HomoSites_outpath,'w')

		for line in self.readVCF():
			if line.startswith("#"):
				HomoSites.write(line)
			elif len(line) != 0 and not line.startswith("#"):
				temp = line.strip().split("\t")
				ale, genotype = temp[3], temp[9]
				GT = genotype.split(":")[0].split("/")
				## homo
				if GT[0] == "." and GT[1] == ".":
					continue
				elif "*" in GT[0] or "*" in GT[1]:
					continue
				else:
					if GT[0] == GT[1]:
						HomoSites.write(line)
		HomoSites.close()


	def getHeteSites(self,HeteSites_outpath):
		"""get homo sites , including header and variation line """
		HeteSites = open(HeteSites_outpath, 'w')

		for line in self.readVCF():
			if line.startswith("#"):
				HeteSites.write(line)
			elif len(line) != 0 and not line.startswith("#"):
				temp = line.strip().split("\t")
				ale, genotype = temp[3], temp[9]
				GT = genotype.split(":")[0].split("/")
				## homo
				if GT[0] == "." and GT[1] == ".":
					continue
				elif "*" in GT[0] or "*" in GT[1]:
					continue
				else:
					if GT[0] != GT[1]:
						HeteSites.write(line)
		HeteSites.close()

	def reduceNonSites(self):
		pass


	def vcf2genotype(self,genotype_outpath):
		IUPAC = {'-' : '-', 'AA' : 'A', 'TT' : 'T', 'CC' : 'C', 'GG' : 'G', 'AG' : 'R',
				'GA' : 'R', 'CT' : 'Y', 'TC' : 'Y', 'AC' : 'M', 'CA' : 'M', 'GT' : 'K',
				'TG' : 'K', 'GC' : 'S', 'CG' : 'S', 'AT' : 'W', 'TA' : 'W',
				'ATC' : 'H', 'ACT' : 'H', 'TAC' : 'H', 'TCA' : 'H', 'CAT' : 'H', 'CTA' : 'H',
				'GTC' : 'B', 'GCT' : 'B', 'CGT' : 'B', 'CTG' : 'B', 'TCG' : 'B', 'TGC' : 'B',
				'GAC' : 'V', 'GCA' : 'V', 'CGA' : 'V', 'CAG' : 'V', 'AGC' : 'V', 'ACG' : 'V',
				'GAT' : 'D', 'GTA' : 'D', 'AGT' : 'D', 'ATG' : 'D', 'TGA' : 'D', 'TAG' : 'D', 
				'ATCG' : 'N', 'ATGC' : 'N', 'ACTG' : 'N', 'ACGT' : 'N', 'AGTC' : 'N', 'AGCT' : 'N', 'TACG' : 'N',
				'TAGC' : 'N', 'TCAG' : 'N', 'TCGA' : 'N', 'TGAC' : 'N', 'TGCA' : 'N', 'CATG' : 'N', 'CAGT' : 'N',
				'CTAG' : 'N', 'CTGA' : 'N', 'CGAT' : 'N', 'CGTA' : 'N', 'GATC' : 'N', 'GACT' : 'N', 'GTAC' : 'N',
				'GTCA' : 'N', 'GCAT' : 'N', 'GCTA' : 'N'}
				#this hash was use for converting the genotype to a single letter according to th IUPAC.
		OUT = open(genotype_outpath,"w")
		header = "Chr\tPosi\tRef"
		for sample in self.getSampleID():
			sampleTab = "\t"+sample
			header += sampleTab
		OUT.write(header+"\n")

		for line in self.readVCF():
			temp,alts = [],[]
			if len(line) != 0 and not line.startswith("#"):
				temp = line.strip().split("\t")
				OUT.write(temp[0]+"\t"+temp[1]+"\t"+temp[3]+"\t") #output the chromosome number, positions and the allele of reference.
				if "," in temp[4]:
					alts = temp[4].split(',')
				else:
					alts.append(temp[4])
				alts.insert(0,temp[3])
				for num in range(9,len(temp)):
					geno = ''
					if re.compile(r'\.\/\.').match(temp[num]):
						geno = "-"
					else:
						aa = re.compile(r'(\d)\/(\d).*').match(temp[num])
						allele1 = alts[int(aa.group(1))]
						allele2 = alts[int(aa.group(2))]
						if "*" in allele1 or "*" in allele2:
							geno = "-"
						else:
							geno = allele1+allele2
					OUT.write(IUPAC[geno]+"\t")
				OUT.write("\n")
		OUT.close()


					
