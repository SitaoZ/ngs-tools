#!/usr/bin/python 
#coding:utf-8
import sys,re,os
import getopt

"""
Author  : Zhu Siao, zhusitao@genomics.cn
Version : 1.0
Date    : 2018-2-1
Dest    : running in python3
"""

USAGE = """
        Author  : Zhu Siao, zhusitao@genomics.cn
        Version : 1.0
        Date    : 2018-2-1
        Dest    : Running in python3, For fasta file cutting.

        Usage   : python fastaDeal.py [options]

			-i , --ifile=	<inputfile>* 
			-o , --ofile=	<outputfile> 
			-f , --cutf=	<cutfiles> 
			-s , --cuts=	<cutfiles> 
			-a , --attr=	<attribute> 
			-p , --patt=	<pattern>
			-t , --subst=	<substruct>
			-u , --reform=  <reform>
			-x , --index=	<index>
			-r , (on-off key in substruct) reverse sequence
			-c , (on-off key in substruct) complement sequence
			-h , --help    

        Examples: python fastaDeal.py -i cds.fa -o ./ -f 30
		  python fastaDeal.py --ifile=cds.fa --ofile=./ --cutf=30
                  python fastaDeal.py -i cds.fa -o ./ -s 10
                  python fastaDeal.py -i cds.fa -a id:len:gc:gap > out.txt
		  python fastaDeal.py -i cds.fa -p "Pm" > out.fa
		  python fastaDeal.py -i cds.fa -t 100-200 > fragment.fa
		  python fastaDeal.py -i cds.fa -t 100-200 -r > fragment_reverse.fa
		  python fastaDeal.py -i cds.fa -t 100-200 -c > fragment_complement_reverse.fa
		  python fastaDeal.py -i cds.fa -u lower50 > reform.fa
		  python fastaDeal.py -i cds.fa -u upper > reform.fa
		  python fastaDeal.py -i cds.fa -x 1-5 > sample.fa
        """


try:
	opts,args = getopt.getopt(sys.argv[1:],"hrci:o:f:s:a:p:t:u:x:",["ifile=","ofile=","cutf=","cuts=","attr=","patt=","subst=","reform=","index="])
except getopt.GetoptError:
	print(USAGE)
	sys.exit(2)
for (opt,arg) in opts:
	if opt in ('-h', '--help'):
		print (USAGE)
		sys.exit()
	if opt in ("-i","--ifile"):
		fasta = arg
	if opt in ("-o","--ofile"):
		outputfile = arg
	if opt in ("-f","--cutf"):
		Cutf = int(arg)
	if opt in ("-s","--cuts"):
		Cuts = int(arg)
	if opt in ("-a","--attr"):
		Attr = arg
	if opt in ("-p","--patt"):
		Patt = arg
	if opt in ("-t","--subst"):
		Subst = arg
	if opt in ('-r'):
		reverse = True # on-off key
	else:
		reverse = False
	if opt in ('-c'):
		complement = True # on-off ket
	else:
		complement = False
	if opt in ("-u","-reform"):
		form = arg
	if opt in ("-x","--index"):
		index = arg
		
def delimited(file, delimiter = '\n'):
	""" a generator to return lines delimited by '>' """
	allLines = file.read()
	lines = allLines.split(delimiter)
	for line in lines:
		if line != '':
			yield line
def displaySeq(seqIn,number):
	""" printf sequence """
	seqOut = ''
	for i in range(1,len(seqIn)//number+1):
		start = (i-1)*number
		end = i*number
		seqOut += seqIn[start:end]+"\n"
	left = len(seqIn) % number
	remainder = seqIn[-left:] if left != 0 else ""
	return seqOut + remainder

def totalLength(fasta):
	""" get a fasta file total length """
	with open(fasta,'r') as F:
		seq = ''
		for line in F.readlines():
			line = line.strip()
			if line.startswith(">"):
				ID = line
			else:
				seq += line
		totalLength = len(seq)
	return totalLength
def cutFastaAdvanced(fasta,outputfile,Cutf=20):
	""" cut fasta file into Cutf parts 每个文件包含指定的碱基数目"""
	allLength = totalLength(fasta)
	file_name = os.popen('basename %s'%fasta).read().strip() 
	if os.path.exists(outputfile):
		outPutDir = outputfile+"/"+file_name+".cut"
		if os.path.exists(outPutDir):
			os.system('rm -rf %s'%outPutDir)
	else:
		outPutDir = "./"+file_name+".cut"
	os.makedirs(outPutDir)
	subLen = allLength//Cutf
	# sub file
	Cur_len = 0
	Cur_content = ''
	file_mark = 1
	with open(fasta,'r') as f:
		lines = delimited(f, '>')
		for line in lines:
			header = line.split("\n")[0]
			rawSeq = ''.join(line.split("\n")[1:])
			ruleSeq = displaySeq(rawSeq,60)
			Cur_len += len(rawSeq)
			Cur_content += ">%s\n%s\n" % (header,ruleSeq)
			if Cur_len >= subLen:
				OUT = open(outPutDir+"/"+file_name+str(file_mark),'w')
				OUT.writelines(Cur_content)
				OUT.close()
				Cur_content = ''
				Cur_len = 0
				file_mark += 1
	if len(Cur_content) != 0:
		##make the last file
		OUT = open(outPutDir+"/"+file_name+str(file_mark),'w')
		OUT.writelines(Cur_content)
		OUT.close()
def cutFasta(fasta,outputfile,Cuts=5):
	"""cut fasta into 每个文件包含指定的fasta条数"""
	file_name = os.popen('basename %s'%fasta).read().strip()
	if os.path.exists(outputfile):
		outPutDir = outputfile+"/"+file_name+".cut"
		if os.path.exists(outPutDir):
			os.system('rm -rf %s'%outPutDir)
		else:
			outPutDir = "./"+file_name+".cut"
	os.makedirs(outPutDir)
	total_num = int(os.popen("grep -c '>' %s"%(file_name)).read().strip())
	seq_num,file_num = 0,0
	if Cuts:
		seq_num = Cuts # sequence numbers in each file
		file_num = total_num//seq_num + 1 # total file numbers
	mark = int(re.sub(r'\d',"0",str(file_num)))
	mark += 1
	dir_rank = 2 if file_num > 100 else 1
	submark,num_loop = 0,0
	with open(fasta,'r') as f:
		lines = delimited(f, '>')
		for line in lines:
			header = line.split("\n")[0]
			seq = ''.join(line.split("\n")[1:])
			seq = displaySeq(seq,60)
			if dir_rank == 2 and num_loop %(seq_num* 100) ==0:
				submark += 1
				os.system('mkdir %s/%s'%(outPutDir,str(submark)))
			if num_loop % seq_num == 0:
				if dir_rank == 1:
					OUT = open(outPutDir+"/"+file_name+str(mark),'w')
				if dir_rank == 2:
					OUT = open(outPutDir+"/"+str(submark)+"/"+file_name+str(mark),'w')
				mark += 1
			OUT.writelines(">"+header+"\n"+seq+"\n")
			if num_loop % seq_num == seq_num-1:
				print ("zhusitao")
				OUT.close()
			num_loop += 1
						
##get different attributions of each sequence
##id:head:seq:len:lenwogap:gc:x:lc:uc
##lenwogap, lenght exclude gap, suppose N n as gap
##gc, includes g c G c
##x, includes x X
##lc, includes a g c t
##uc, includes A G C T
####################################################
def getAttribute(fasta,attribute):
	""" get different attributions of each sequence """
	attr = attribute.split(":")
	work = {}
	key = []
	for i in attr:
		UP = i.upper()
		work[UP] = 1
		key.append(UP)
	with open(fasta,'r') as f:
		lines = delimited(f, '>')
		for line in lines:
			header = line.split("\n")[0]
			seq = ''.join(line.split("\n")[1:])
			Tlength = len(seq)
			wogaplen = len(re.sub(r'[Nn]','',seq))
			
			result = {}
			if 'ID' in work.keys():
				result['ID'] = header
			if 'LEN' in work.keys():
				result['LEN'] = Tlength
			if 'GAP' in work.keys():
				Gaplength = Tlength - wogaplen
				result['GAP'] = Gaplength
			if 'GC' in work.keys():
				GClength = Tlength - len(re.sub(r'[GCgc]','',seq))
				result['GC'] = GClength/wogaplen
			if 'X' in work.keys():
				Xlength = Tlength -len(re.sub(r'[Xx]','',seq))
				result['X'] = Xlength/wogaplen
			if 'LC' in work.keys():
				LC = Tlength - re.sub(r'[atgc]','',seq)
				result['LC'] = LC/wogaplen
			if 'UC' in work.keys():
				UC = Tlength - re.sub(r'[ATGC]','',seq)
				result['UC'] = UC/wogaplen
			outPut = ''
			for key in result:
				outPut += str(result[key])+"\t"
			print (outPut)
def patternHeader(fasta,Patt):
	""" match the regular expression """
	with open(fasta,'r') as f:
		lines = delimited(f, '>')
		for line in lines:
			header = line.split("\n")[0]
			seq = ''.join(line.split("\n")[1:])
			if Patt:
				pat = re.compile('%s'%(Patt))
				if pat.search(header):
					print (">"+header+"\n"+displaySeq(seq,60))

def complementReverse(sequence):
	""" complement and reverse a special DNA sequence """
	sequence = sequence.upper()
	sequence = sequence.replace('A', 't')
	sequence = sequence.replace('T', 'a')
	sequence = sequence.replace('C', 'g')
	sequence = sequence.replace('G', 'c')
	sequence = sequence[::-1]
	sequence = sequence.upper()
	return sequence

							
def subStructureFragment(fasta,Subst):
	""" 截取某条染色体特定长度的序列，反向，互补"""
	with open(fasta,'r') as f:
		pattern = re.compile(r'^(\S+)')
		pattern1 = re.compile(r'(start|\d+)-(end|\d+)')
		lines = delimited(f, '>')
		for line in lines:
			header = line.split("\n")[0]
			rawSeq = ''.join(line.split("\n")[1:])
			seqLength = len(rawSeq)
			if pattern.match(header):
				seqName = pattern.match(header).group(1)
			else:
				seqName = header
			subHead = ">"+seqName
			subSeq = ''
			if pattern1.match(Subst):
				if 'start' in pattern1.match(Subst).group(1):
					subStart = pattern1.match(Subst).group(1)
				else :
					subStart = int(pattern1.match(Subst).group(1))
				if 'end' in pattern1.match(Subst).group(2):
					subEnd = pattern1.match(Subst).group(2)
				else:
					subEnd = int(pattern1.match(Subst).group(2))
				if subStart == 'start' or subStart < 1:
					subStart = 1
				if subEnd == 'end' or subEnd > seqLength:
					subEnd = seqLength
				subHead += "_"+str(subStart)+"_"+str(subEnd)
				subSeq = rawSeq[subStart:subEnd]
			if reverse:
				subSeq = subSeq[::-1]
				subHead += "_reverse"
			if complement:
				subSeq = complementReverse(subSeq)
				subHead += "_complement"
			subSeq = displaySeq(subSeq,60)
			print (subHead+"\n"+subSeq)

def reformSeq(fasta,form):
	""" all capital output or lowercase output"""
	with open(fasta,'r') as f:
		pattern = re.compile(r'([a-z]*)(\d+)')
		lines = delimited(f, '>')
		for line in lines:
			header = line.split("\n")[0]
			seq = ''.join(line.split("\n")[1:])		 
			if form == 'lower':
				seq = seq.lower()
			if form == 'upper':
				seq = seq.upper()
			match = pattern.match(form)
			if match:
				number = int(match.group(2))
				seq = displaySeq(seq,number)
				print ('>'+header+"\n"+seq)
			elif form == "pure":
				seq = seq.replace('N','')
				seq = displaySeq(seq,60)
				print ('>'+header+"\n"+seq)
			else:
				seq = displaySeq(seq,60)
				print ('>'+header+"\n"+seq)

def sampleSeq(fasta,index):
	""" sample a example fasta file in a big fasta file """
	start,end,loop = 0,0,0
	pattern1 = re.compile(r'(\d+)-(\d+)')
	pattern2 = re.compile(r'(\d+)')
	with open(fasta,'r') as f:
		lines = delimited(f, '>')
		for line in lines:
			loop += 1
			header = line.split("\n")[0]
			seq = ''.join(line.split("\n")[1:])
			match1 = pattern1.match(index)
			match2 = pattern2.match(index)
			if match1 :
				start = int(match1.group(1))
				end = int(match1.group(2))
			elif match2 :
				start = int(match2.group(1))
				end = int(match2.group(1))
			else:
				print ("please use -h or --help for more illustration!")
		
			if loop >= start and loop <= end:
				seq = displaySeq(seq,60)
				print ('>'+header+"\n"+seq)
			
#####    main execute #####
if __name__ == "__main__":
	if "Cutf" in locals().keys():
		# if Cutf defined, then execute
		cutFastaAdvanced(fasta,outputfile,Cutf)
	elif "Cuts" in locals().keys():
		cutFasta(fasta,outputfile,Cuts)
	elif "Attr" in locals().keys():
		getAttribute(fasta,Attr)
	elif "Patt" in locals().keys():
		patternHeader(fasta,Patt)
	elif "Subst" in locals().keys():
		subStructureFragment(fasta,Subst)
	elif "form" in locals().keys():
		reformSeq(fasta,form)
	elif "index" in locals().keys():
		sampleSeq(fasta,index)
