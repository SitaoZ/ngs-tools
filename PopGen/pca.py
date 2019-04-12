#!/usr/bin/python 
#coding:utf-8 
import sys,re
import math
#import numpy as np
#import pandas as pd
#import matplotlib.pyplot as plt
from optparse import OptionParser
"""
Author : zhusitao 
Mail   : zhusitao@genomics.cn
Data   : 2018/1/29
Describe: This program is used to calculate population PCA and plot
	  The program was developed by python2.7 ; need change python3
"""
#usage = "usage: %prog [options] arg1 arg2"  
#parser = OptionParser(usage=usage)  
#parser.add_option("-v", "--verbose",  
#                  action="store_true", dest="verbose", default=True,  
#                  help="make lots of noise [default]")  
#parser.add_option("-q", "--quiet",  
#                  action="store_false", dest="verbose",  
#                  help="be vewwy quiet (I'm hunting wabbits)")  
#parser.add_option("-f", "--filename",  
#                  metavar="FILE", help="write output to FILE"),  
#parser.add_option("-m", "--mode",  
#                  default="intermediate",  
#              help="interaction mode: novice, intermediate, "  
#                   "or expert [default: %default]")

## step 1 snp vcf file to genotype file 
def snpToGenotype(snpvcf,genotype):
	""" Change snp vcf file to genotype file! the order is according the final snp vcf!"""
	#this hash was use for converting the genotype to a single letter according to th IUPAC
	degen={'-':'-','AA':'A','TT':'T','CC':'C','GG':'G','AG':'R','GA':'R','CT':'Y','TC':'Y','AC':'M','CA':'M','GT':'K','TG':'K','GC':'S','CG':'S','AT':'W','TA':'W','ATC':'H','ACT':'H','TAC':'H','TCA':'H','CAT':'H','CTA':'H','GTC':'B','GCT':'B','CGT':'B','CTG':'B','TCG':'B','TGC':'B','GAC':'V','GCA':'V','CGA':'V','CAG':'V','AGC':'V','ACG':'V','GAT':'D','GTA':'D','AGT':'D','ATG':'D','TGA':'D','TAG':'D','ATCG':'N','ATGC':'N','ACTG':'N','ACGT':'N','AGTC':'N','AGCT':'N','TACG':'N','TAGC':'N','TCAG':'N','TCGA':'N','TGAC':'N','TGCA':'N','CATG':'N','CAGT':'N','CTAG':'N','CTGA':'N','CGAT':'N','CGTA':'N','GATC':'N','GACT':'N','GTAC':'N','GTCA':'N','GCAT':'N','GCTA':'N'}
	number = 0
	IN = open(snpvcf,'r')
	OUT = open(genotype,'w')
	patternHeader = re.compile('^#')
	lostGenotype = re.compile('\./\.')
	realGenotype = re.compile('(\d)/(\d)')
	for line in IN.readlines():
		alt = []
		line = line.strip('\n')
		if patternHeader.search(line):
			continue
		number += 1
		if number != 1: #add the \n at the beganing of all lines except for the first line.
			OUT.writelines("\n")
		ArrayLine = line.split("\t")
		chrom = ArrayLine[0]
		posi  = ArrayLine[1]	
		ref   = ArrayLine[3]
		OUT.writelines(chrom+"\t"+posi+"\t"+ref+" ")#output the chromosome number, positions and the allele of reference
		if ',' in ArrayLine[4]:
			alt = ArrayLine[4].split(',')
		else:
			alt.append(ArrayLine[4])
		alt.insert(0,ref) # put the ref in the front of the array
		for i in range(9,len(ArrayLine)):
			match1 = lostGenotype.search(ArrayLine[i])
			if match1:
				geno = '-' # if the allele was lost in one individual, then the genotypoe will be assigned "-"
			else:
				match2 = realGenotype.search(ArrayLine[i])
				if match2:
					allele1 = alt[int(match2.group(1))]
					allele2 = alt[int(match2.group(2))]
					if ("*" in allele1) or ("*" in allele2):
						geno = "-"
					else :
						geno = allele1+allele2
			OUT.writelines(degen[geno]+" ")	
## step 2 chang genotype file into number  
def fillGenotype(ref,geno,base1,base2,i,nSNP0):
	""" To fill the heterozygous genotype! """
	pref = ref
	pgeno = geno
	major = base1
	minor = base2
	ni = i
	pnSNP = nSNP0
	pokm = 0
	psnp = 0
	if pref[0] == major:
		if pref[1] == "-":
			pref[1] = minor
			pgeno[ni-1] = 1 #het
		else:
			if minor == pref[1]:# minor allele
				pgeno[ni-1] = 1	#het
			else:
				pokm = 1
				psnp += 1
				toto = pnSNP+1
				print "Multiallele line#\t"+str(toto)+"\tat individual #\t"+str(ni)
	else:
		if pref[0] == minor:
			if pref[1] == "-": #1st snp
				pref[1] = major
				pgeno[ni-1] = 1 #het
			else:
				if major == pref[1]: # minor allele
					pgeno[ni-1] = 1 #het
				else:
					pokm = 1
					psnp += 1
					toto = pnSNP+1
					print "Multiallele line#\t"+str(toto)+"\tat individual #\t"+str(ni)
		else:
			pokm = 1
			psnp += 1
			toto = pnSNP+1
			print "Multiallele line#\t"+str(toto)+"\tat individual #\t"+str(ni)
	return pref,pgeno,pokm,psnp

## step 2 chang genotype file into number
def genotypeToNum(genotype,sampleNumbers,outputDir):
	""" The function transformates genotype file into numbers.
	(AC):M	(AG):R (AT):W (CG):S (CT):Y (GT):K """
	fileName = genotype.split("/")[-1]
	miss = []
	for i in range(int(sampleNumbers)):
		miss.append(0)
	with open(genotype,'r') as F :
		FileOut = open(outputDir+"/T-"+fileName,'w')
		FileAvg = open(outputDir+"/AVG-"+fileName,'w')
		FileSnpInd = open(outputDir+"/SNPIND-"+fileName,'w')
		geno = []
		for i in range(int(sampleNumbers)):
			geno.append(0)
		nSNP = [0,0,0,0]
		#nSNP[0]=nSNP[1]=nSNP[2]=nSNP[3]=0 # 0:total 1:missing 2:mlti 3:counted for PCA removing nonvariable
		nmiss = 0 
		for line in F.readlines():
			ArrayLine = line.strip('\n').split('\t')
			ArrayGeno = ArrayLine[2].split() # in default by blant
			ref = [0,0]
			ref[0] = ArrayGeno[0] # ref nt
			ref[1] = '-'	
			okm = 0
			i = 1
			while (i < sampleNumbers + 1) and (okm == 0):
				if ArrayGeno[i] == "N":
					ArrayGeno[i] = "-"
				if ArrayGeno[i] !=  "-":
					#homo
					if ( ArrayGeno[i]=="A") or (ArrayGeno[i]=="T") or (ArrayGeno[i]=="G") or (ArrayGeno[i]=="C"):
						if ArrayGeno[i] == ref[0]:
							geno[i-1] = 0         # homo ref
						else:
							if ref[1] == "-":
								ref[1] = ArrayGeno[i]
								geno[i-1] = 2 # homo minor
							else:
								if ArrayGeno[i] == ref[1]: # minor allele
									geno[i-1] = 2 # homo minor
								else:                 # new allele= multiple alleles: remove this SNP
									okm=1
									nSNP[2] += 1
									toto = nSNP[0]+1
									print "Multiallele line#\t"+str(toto)
					else:#hete
						if ArrayGeno[i] == "R":# AG GA
							(rec0,rec1,rec2,rec3) = fillGenotype(ref,geno,"A","G",i,nSNP[0])
							ref = rec0
							geno = rec1
							okm = rec2
							nSNP[2] += rec3	
						if ArrayGeno[i] == "S":# CG GC
							(rec0,rec1,rec2,rec3) = fillGenotype(ref,geno,"C","G",i,nSNP[0])
							ref = rec0
                                                        geno = rec1
                                                        okm = rec2
							nSNP[2] += rec3
						if ArrayGeno[i] == "W":# AT TA
							(rec0,rec1,rec2,rec3) = fillGenotype(ref,geno,"A","T",i,nSNP[0])
							ref = rec0
							geno = rec1
							okm = rec2
							nSNP[2] += rec3
						if ArrayGeno[i] == "Y":# CT TC
							(rec0,rec1,rec2,rec3) = fillGenotype(ref,geno,"C","T",i,nSNP[0])
							ref = rec0
                                                        geno = rec1
                                                        okm = rec2
                                                        nSNP[2] += rec3
						if ArrayGeno[i] == "K":# GT TG
							(rec0,rec1,rec2,rec3) = fillGenotype(ref,geno,"G","T",i,nSNP[0])
							ref = rec0
                                                        geno = rec1
                                                        okm = rec2
                                                        nSNP[2] += rec3
						if ArrayGeno[i] == "M":# AC CA
							(rec0,rec1,rec2,rec3) = fillGenotype(ref,geno,"A","C",i,nSNP[0])
                                                        ref = rec0
                                                        geno = rec1
                                                        okm = rec2
                                                        nSNP[2] += rec3
				else:
					geno[i-1]="NA"	
				i += 1
			nSNP[0] += 1
			## print out genotype if no miltiple allele
			if okm==0:
				nSNP[1] += 1
				avg=0
				ind=0
				stgeno=""
				i = 0 
				while i < sampleNumbers:
					stgeno= stgeno+str(geno[i])+"\t"
					if geno[i] != "NA":
						ind += 1
						avg += geno[i]
					else:
						miss[i] += 1
					i +=1
				if ind >0:
					avg = avg/float(ind)
				if (avg>0)and(avg<2):
					FileAvg.writelines(str(avg)+"\n")
					FileOut.writelines(stgeno+"\n")	
					nSNP[3] += 1
				else:
					Itoto = nSNP[0]
					print "Invariant site line#\t"+str(Itoto)
		print fileName+"\t"+str(nSNP[0])+"\t"+str(nSNP[1])+"\t"+str(nSNP[2])+"\t"+str(nSNP[3])+"\n"
		FileSnpInd.writelines(str(nSNP[3])+"\t"+str(sampleNumbers)+"\n")			
		FileOut.close()
		FileAvg.close()
		FileSnpInd.close()
	with open(outputDir+"/INDV-NUM-SNPs",'w') as F:
		i = 0 
		while i < sampleNumbers:
			F.writelines(str(i)+"\t"+str(miss[i])+"\n")
			i += 1

		
				
			
			
	
## step 3 calculate mean,standard deviation,covariance,covariance matrix
def calculateCovariance(InputDir,chrname,sampleNumbers):
	""" Calculate covariance matrix. dimension is individual numbers,each snp is observed value! the order is also same to the final vcf!
	    if you have same problems ,you could pay more attention to covariance matrix!"""
	snp = 0
	MtM = [[0 for i in range(sampleNumbers)]for j in range(sampleNumbers)]
	ind = sampleNumbers
	FILEIN=open(InputDir+"/T-"+chrname,'r')
	FILEAVG=open(InputDir+"/AVG-"+chrname,'r')
	FILESUM=open(InputDir+"/X-"+chrname,'w')
	
	for line in FILEIN.readlines():
		avg = FILEAVG.readline()
		avg = avg.strip()
		line = line.strip()
		avg = float(avg)
		geno = line.split('\t')
		for i in range(len(geno)):
			if geno[i] != "NA":
				geno[i] = int(geno[i])
		var = math.sqrt(avg*(1-avg/2)/2)
		# init M
		M = [0 for i in range(len(geno))]
		for j in range(len(geno)):
			if (var >0) and (geno[j] != "NA"):
				M[j] = (geno[j]-avg)/var ## varicnce
			else:
				M[j] = 0
			for k in range(j,len(geno)):
				if snp == 0:
					MtM[k][j]=MtM[j][k]=0 # init
				if (var >0) and (geno[k] != "NA"):
					M[k]= (geno[k] - avg)/var ## variance
				else:
					M[k]= 0
				if ((M[k] != "NA")and(M[j] != "NA")) and (MtM[j][k] != "NA"):
					MtM[j][k]+= M[j]*M[k] ## covariance
					MtM[k][j] = MtM[j][k] ## covariance
				else:
					MtM[k][j]=MtM[j][k] = "NA"
					print "PROBLEM: Remaining Missing Data\t"+chrname+"\t"+"ind "+j+"\t"+"SNP "+snp+"\n"
		snp += 1
		geno = []
				
	FILEIN.close()
	FILEAVG.close()
	print "IND\t"+str(ind)+"\tSNP\t"+str(snp)+"\n";
	for i in range(sampleNumbers):
		for j in range(sampleNumbers):
			if MtM[i][j] != "NA":
				MtM[i][j] = MtM[i][j]/float(snp) ## mean of covariance. MtM[i][j] is two demendion covariance.in theory should be divided (snp-1),
								 ## but here is snp
			else:
				print "PROBLEM: Remaining Missing Data\n"
			## when all loops finished finally; FILESUM get covariance matrix
			FILESUM.writelines(str(MtM[i][j])+"\t")
		FILESUM.writelines("\n")
	FILESUM.close()

def start_with(InputDir,startString):
	""" The function aims to return startswith a string files or directories in given directory! """
	import os
	if not os.path.isdir(InputDir):
		return "not true"
	else:
		d = os.listdir(InputDir)
		lst = [x for x in d]
		if len(lst)>0:
			return [x for x in lst if x.startswith(startString)]
		else:
			return "nodir"

	
## step 4: sum all chromosome covariance matrix into a total one 
def sumChromoCovarianceMatrix(InputDir,sampleNumbers):
	""" Get the final covariance matrix,put each chromosome into a total covariance matrix! """
	SNPIND_files = start_with(InputDir,"SNPIND-")
	X_files = start_with(InputDir,"X-")
	OUT = open("/hwfssz1/ST_BIGDATA/USER/zhusitao/Project/python/PCA/test"+"/Y-sum-auto",'w')
	dictZ = {}
	Mat = [[0 for i in range(sampleNumbers)] for j in range(sampleNumbers)]
	for i in SNPIND_files:
		Fi = open(InputDir+"/"+i,'r')
		ID = i.split("SNPIND-")[-1]
		for line in Fi.readlines():
			inf = line.strip().split("\t")
			dictZ[ID] = int(inf[0])
		Fi.close()
	smaple_count=SNP_total=0
	for i in X_files:
		ID = i.split("X-")[-1]
		Fi = open(InputDir+"/"+i,'r')
		if not dictZ[ID]:
			print "wrong ID"
		else:
			#print ID,dictZ[ID]
			SNP_total += dictZ[ID]
		line = 0
		for oneline in Fi.readlines():
			inf = oneline.strip().split("\t")
			if smaple_count == 0:
				smaple_count = len(inf)-1
			elif smaple_count != len(inf)-1:
				print "Bad_1 "+ID+"\tlint\t"+str(line)+"\n"
				sys.exit(1)
			for j in range(len(inf)):
				Mat[line][j] += float(inf[j]) * dictZ[ID]
			line += 1
		if smaple_count != line -1:
			print "Bad_2 "+ID+"\tline\t"+str(line)+"\n"
	covMatrix = [[0 for i in range(sampleNumbers)] for j in range(sampleNumbers)]
	for i in range(smaple_count+1):
		for j in range(smaple_count+1):
			covMatrix[i][j] = "%.15f" % (Mat[i][j]/float(SNP_total))
			OUT.writelines(str(covMatrix[i][j])+"\t")
		OUT.writelines("\n")
	OUT.close()
	return covMatrix

# step 5: calculate eigenvalue and eigenvector
def calclateEigenvalueEigenvector(covMatrix):
	"""Calculate the eigenvalue and eigenvalue matrix. Each eigenvalue corresponds to an eigenvector!!!"""
	covarianceMatrix = np.array(covMatrix,dtype=float)
	#eigenvalue = np.linalg.eigvals(covarianceMatrix) # A1 eigenvalue
	eigenvalue,eigenvector = np.linalg.eig(covarianceMatrix)   # A2 also eigenvalue, B is eigenvector
	vectorDict = {}
	for i in range(len(eigenvalue)):
		vector = []
		for j in range(len(eigenvalue)):
			vector.append(eigenvector[j][i]) # every eigenvalue coresponding to a column 
		vectorDict[eigenvalue[i]] = vector
	#for i in sorted(vectorDict.keys(),reverse = True):
	#	print i,vectorDict[i]
	dataFrame = pd.DataFrame(vectorDict).sort_index(axis=1,ascending=False)
	#print dataFrame
	return dataFrame

# step 6: pca plot
def pcaPlot(dataFrame,individual,outDir):
	""" Plot the final pca result! """
	## choose  6 eigenvalues 
	head6Eigenvalues = dataFrame.columns.tolist()[0:6]
	print "The first 6 eigenvalue: ",head6Eigenvalues
	## choose 4 eigenvectors
	head4Eigenvectors = dataFrame.columns.tolist()[0:4]
	newDataFrame = dataFrame[head4Eigenvectors]
	print newDataFrame
	individual_list = []
	with open(individual,'r') as F:
		for line in F.readlines():
			ID = line.strip('\n').split('\t')[-1]
			individual_list.append(ID)
	
	newDataFrame["individuals"] = individual_list
	# create plotDataFrame for plot in matplotlib
	plotDataFrame = pd.DataFrame()
	for i in range(1,len(newDataFrame.columns.tolist())+1):
		if i == 5:
			plotDataFrame["individuals"] = newDataFrame[newDataFrame.columns.tolist()[i-1]]
			break
		else:
			plotDataFrame["vector %d"%i] = newDataFrame[newDataFrame.columns.tolist()[i-1]]
	print plotDataFrame["vector 1"]	
	#plt.plot(plotDataFrame["vector 1"],plotDataFrame["vector 2"])
	plt.switch_backend('agg')
	plt.figure(num=1, figsize=(8, 6))
	#plt.set_title("PCA Plot")
	#plt.annotate("(3,6)", xy = (3, 6), xytext = (4, 5), arrowprops = dict(facecolor = 'black', shrink = 0.1))
	plt.xlabel("PC1")
	plt.ylabel("PC2")
	colors = np.random.rand(len(plotDataFrame["vector 1"]))
	plt.scatter(x = plotDataFrame["vector 1"],y = plotDataFrame["vector 2"],c = colors)
	#for x,y,z in zip (plotDataFrame["vector 1"],plotDataFrame["vector 2"],plotDataFrame["individuals"]):
	#	print x,y,z
	#	plt.annotate(
	#	'(%s, %s)' %(x, y),
	#	#'(%s)'%(z)
	#	xy=(x, y),
	#	xytext=(0, -10),
	#	textcoords='offset points',
	#	ha='center',
	#	va='top')
	plt.legend(plotDataFrame["individuals"],scatterpoints=1)
	plt.savefig("A.png", format='png', dpi=300)
	plt.close()

def usage():
	usage = "usage: %prog [options] args"  
	parser = OptionParser(usage)  
	parser.add_option("-f", "--file", dest="filename", help="snp vcf file from populaiton ") 
	parser.add_option("-g", "--genotype", dest="genotypePath",help="genotype file")
	parser.add_option("-n", "--numbers", dest="numbers",help="the number of individuals")
	parser.add_option("-i", "--individualTxt", dest="individualTxt",help="individual txt")
	parser.add_option("-p", "--path_num", help="the path dir of genotype numbers")
	(options, args) = parser.parse_args()
	# options is a dict; args is a list
	if not options:
		print "python %prog -h"
	else:
		return (options.filename,options.genotypePath,options.numbers,options.individualTxt)
	#if len(options.keys()) != 1:  
	#	parser.error("""Incorrect number of arguments.
	#	Add -h(--help) view help information""")  

if __name__ == '__main__':
	(snpvcf,genotype,number,individual) = usage()
	print snpvcf,genotype,number,individual
	###-------   step 0 : args for the program   -----------------------### 
	#snpvcf=sys.argv[1]
	#genotype=sys.argv[2]
	#individual = sys.argv[1]
	###------   step 1 : transformate snp to genotype   ----------------###
	#snpToGenotype(snpvcf,genotype)
	###------   step 2 : genotype to numbers   -------------------------###
	#genotypeToNum("/hwfssz1/ST_BIGDATA/USER/zhusitao/Project/python/PCA/out.genotype",16,"/hwfssz1/ST_BIGDATA/USER/zhusitao/Project/python/PCA/test")
	###------   step 3 : calculate covariance   ------------------------###
	#calculateCovariance("/hwfssz1/ST_BIGDATA/USER/zhusitao/Project/python/PCA/test","out.genotype",16)
	###------   step 4 : merge each chromosome covariance in total   ---###
	#covMatrix = sumChromoCovarianceMatrix("/hwfssz1/ST_BIGDATA/USER/zhusitao/Project/xj/reseq/result/08.PCA/PCA/data",16)	
	###------   step 5 : calculate covariance matrix   -----------------###
	#dataFrame = calclateEigenvalueEigenvector(covMatrix)
	###------   step 6 : plot pca result   -----------------------------###
	#pcaPlot(dataFrame,individual,"/hwfssz1/ST_BIGDATA/USER/zhusitao/Project/python/PCA")
