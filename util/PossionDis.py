#!/usr/bin/python 
import os,re,sys
import getopt
from collections import defaultdict
from math import log,fabs

"""
Author : Zhu Sitao , zhusitao@genomics.cn
Version: v1.0 ,run in python 3.xx
Date   : 2018-2-11
Aim    : Decription:find Different Expression Genes between two samples
"""
USAGE = """
	Author : Zhu Sitao , zhusitao@genomics.cn
	Version: v1.0 ,run in python 3.xx
	Date   : 2018-2-11
	Dest   : Decription:find Different Expression Genes between two samples

	Usage  :
		python PossionDis.py [options]
		
			-l , --list <compare list of two samples,fromat:sampleName<tab>geneExpression.xls,the geneExpression file must be bed file>
			-c , --GeneIDColumn <the column of gene Id in geneExpression.xls file>
			-e , --GeneExpColumn <the column of reads counts of gene in geneExpression.xls file>
			-L , --GeneLenColumn <the column of gene length in geneExpression.xls file>
			-s , --standardColumn <the column of standard expression in geneExpression.xls file>
			-m , --standardMethod <method which is used to standardize gene expression>
			-g , --log2 <log2 cutoff of fold-change, default is 1>
			-f , --fdr  <fdr cutoff, default is 0.001>
			-o , --outdir <output directory, default is current directory ".">
			-h , --help
	"""

try:
        opts,args = getopt.getopt(sys.argv[1:],"hl:c:e:L:s:m:g:f:o:",["list=","GeneIDColumn=","GeneExpColumn=","GeneLenColumn=","standardColumn=","standardMethod=","log2=","fdr=","outdir="])
except getopt.GetoptError:
        print(USAGE)
        sys.exit(2)
for (opt,arg) in opts:
	if opt in ('-h', '--help'):
		print (USAGE)
		sys.exit()
	if opt in ("-l","--list"):
                List = arg
	if opt in ("-c","--GeneIDColumn"):
		GeneIDColumn = int(arg) 
	if opt in ("-e","--GeneExpColumn"):
		GeneExpColumn = int(arg)
	if opt in ("-L","--GeneLenColumn"):
		GeneLenColumn = int(arg)
	if opt in ("-s","--standardColumn"):
		standardColumn = int(arg)
	if opt in ("-m","--standardMethod"):
		standardMethod = int(arg)
	if opt in ("-g","--log2"):
		log2 = float(arg)
	if opt in ("-f","--fdr"):
		fdr = float(arg)
	if opt in ("-o","--outdir"):
		outdir = arg
## default set 
if "GeneIDColumn" not in locals().keys():
	GeneIDColumn = 1
if "GeneExpColumn" not in locals().keys():
	GeneExpColumn = 4
if "GeneLenColumn" not in locals().keys():
	GeneLenColumn = 3
if "standardColumn" not in locals().keys():
	standardColumn = 5
if "standardMethod" not in locals().keys():
	standardMethod = "FPKM"
if "outdir" not in locals().keys():
	outdir = "./"
if "log2" not in locals().keys():
	log2 = 1
if "fdr" not in locals().keys():
	fdr = 0.001
#print (locals().keys())
#print (locals()["__name__"])
#print (locals()["__doc__"])
#print (locals()["GeneIDColumn"])
#print ("GeneIDColumn: "+str(GeneIDColumn))
#print ("GeneExpColumn: "+str(GeneExpColumn))

def normalize(count):
	""" normalize count """
	if count > 0 :
		count = 1
	elif count < 0 :
		count = -1
	else:
		count = 0
	return count

def tree():
	""" create a tree dataSturcture, it's smarter than dict """
	return defaultdict(tree)


if __name__ == "__main__":
	samples = []
	expFiles = {}
	genes = tree()
	sampleCheck = tree()
	sampleReadsNum = defaultdict(float)

	with open (List,'r') as f:
		pattern = re.compile(r'(\S+)\s+(\S+)')
		for line in f.readlines():
			line = line.strip()
			if len(line) == 0:
				continue
			match = pattern.match(line)
			if match :
				sampleID = match.group(1)
				filePath = match.group(2)
				samples.append(sampleID)
				expFiles[sampleID] = filePath

	if len(samples) % 2 != 0:
		print ("the number of samples need be a multiple of 2")
		sys.exit(2)
	minStandardValue = 0.1
	# save expression, gene length, reads number of gene in each sample
	for s in samples:
		ff = expFiles[s]
		with open (ff,'r') as f:
			header = f.readline().strip()
			#print ("Header: " + header)
			for line in f.readlines():
				line = line.strip()
				tmp = line.split('\t')
				genes[tmp[GeneIDColumn-1]]['len'] = float(tmp[GeneLenColumn-1])
				genes[tmp[GeneIDColumn-1]][s]['reads'] = float(tmp[GeneExpColumn-1])
				genes[tmp[GeneIDColumn-1]][s][standardMethod] = float(tmp[standardColumn-1])
				if not sampleCheck[s]:
					sampleReadsNum[s] += float(tmp[GeneExpColumn-1])
				if genes[tmp[GeneIDColumn-1]][s][standardMethod] <= minStandardValue:
					minStandardValue += genes[tmp[GeneIDColumn-1]][s][standardMethod]
		sampleCheck[s] = 1

	# differential analysis

	for i in range(0,len(samples),2):
		sampleA,sampleB = samples[i],samples[i+1]
		diffs = tree()
		if sampleReadsNum[sampleA] <= 0 or sampleReadsNum[sampleB] <= 0:
			print ("total reads number of %s or %s is 0"%(sampleA,sampleB))
			sys.exit(2)
		# p-value
		PV = open("%s/%s-VS-%s.4pvalue"%(outdir,sampleA,sampleB),'w')
		for gene in genes.keys():
			if (not genes[gene][sampleA] or genes[gene][sampleA]['reads'] == 0) and (not genes[gene][sampleB] or genes[gene][sampleB]['reads'] == 0):
				continue
			if not genes[gene][sampleA]:
				genes[gene][sampleA]['reads'] == 0
				genes[gene][sampleA][standardMethod] = minStandardValue
			if not genes[gene][sampleB]:
				genes[gene][sampleB]['reads'] == 0
				genes[gene][sampleB][standardMethod] = minStandardValue
			line = "%s\t%s\t%s\t%s\t%s\n" % (gene,str(sampleReadsNum[sampleA]),str(sampleReadsNum[sampleB]),str(genes[gene][sampleA]['reads']),str(genes[gene][sampleB]['reads']))
			print (line)
			PV.writelines(line)		
		PV.close()

		os.system('/hwfssz1/ST_BIGDATA/USER/zhusitao/Bin/pythonBin/perl/GeneDiffExp/RNA_RNAdenovo_2015a_GeneDiffExp/statistic -i %s/%s-VS-%s.4pvalue -o %s/%s-VS-%s.GeneDiffExp.prepare.xls'%(outdir,sampleA,sampleB,outdir,sampleA,sampleB))
		# read perpare file
		IN = open("%s/%s-VS-%s.GeneDiffExp.prepare.xls"%(outdir,sampleA,sampleB))
		line_number = 0
		for line in IN.readlines():
			tabs = line.strip().split('\t')
			diffs[tabs[0]]["pvalue"] = float(tabs[5]) * 2
			line_number += 1
		total = line_number
		IN.close()
		# log2 fdr up/down
		sort_i = 0
		for gene,value in sorted(diffs.items(),key=lambda gene_value:gene_value[1]['pvalue'],reverse=False):
			sort_i += 1 
			diffs[gene]['fdr'] = diffs[gene]['pvalue'] * total / sort_i
			standardExpA = genes[gene][sampleA][standardMethod]
			standardExpB = genes[gene][sampleB][standardMethod]
			diffs[gene]['log2'] = log(standardExpB/standardExpA)/log(2)
			if diffs[gene]['log2'] > 0:
				diffs[gene]['ud'] = "Up"
			elif diffs[gene]['log2'] < 0:
				diffs[gene]['ud'] = "Down"
			else:
				diffs[gene]['ud'] = "-"
				
		# print resut
		OUT = open("%s/%s-VS-%s.PossionDis_Method.GeneDiffExp.xls"%(outdir,sampleA,sampleB),'w')
		FilterOut = open("%s/%s-VS-%s.PossionDis_Method.GeneDiffExpFilter.xls"%(outdir,sampleA,sampleB),'w')
		header = "GeneID\tLength\t%s-Expression\t%s-Expression\tlog2FoldChange(%s/%s)\tPvalue\tFDR\tUp/Down-Regulation(%s/%s)\n"%(sampleA,sampleB,sampleB,sampleA,sampleB,sampleA)
		OUT.writelines(header)
		FilterOut.writelines(header)
		for gene,value in sorted(diffs.items(),key=lambda gene_value:gene_value[1]['log2'],reverse=True):
			if not (fabs(diffs[gene]['log2']) >= log2 and diffs[gene]['fdr'] <= fdr):
				diffs[gene]['ud'] = "*"
			line = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(gene,genes[gene]['len'],genes[gene][sampleA][standardMethod],genes[gene][sampleB][standardMethod],diffs[gene]['log2'],diffs[gene]['pvalue'],diffs[gene]['fdr'],diffs[gene]['ud'])
			OUT.writelines(line)
			# filter
			if fabs(diffs[gene]['log2']) >= log2 and diffs[gene]['fdr'] <= fdr:
				FilterOut.writelines(line)
		OUT.close()
		FilterOut.close()
	
