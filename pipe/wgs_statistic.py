#!/usr/bin/python
import os 
import re
from math import fabs 
import collections
from readFasta import ReadFasta
class wgs_statistic(object):
    """ this class for wgs result display including clean data, bam, snp ,indel and the variants distribution """
    def __init__(self):
        #self.fileList = bam
        #self.refer = refer
        self.samtools = "/hwfssz1/ST_BIGDATA/USER/zhusitao/Software/samtools-1.3.1/bin/samtools"
        self.outdir = os.getcwd()
    def cleanData_stat(self):
        pass
    def bam_stat(self):
        pass
    def picard_stat(self,bam):
        clean_reads,clean_bases = 0,0
        mapped_reads,mapped_bases = 0,0
        mapping_rate,uniq_rate = 0,0
        dup_reads,dup_rate = 0,0 
        mismatch_bases,mismatch_rate = 0,0
        uniq_reads,uniq_bases = 0,0
        mappedPattern = re.compile(r'XC:i:(\d+)')
        uniqPattern = re.compile(r'X0:i:(\d+)')
        mismatchPattern1 = re.compile(r'XM:i:(\d+)')
        mismatchPattern2 = re.compile(r'MD:Z:(\S+)')

        bamfh = os.popen("%s view %s"%(self.samtools,bam))
        for line in bamfh.readlines():
            tmp = line.strip().split()
            flag = int(tmp[1])
            if(flag & 0x100):
                continue 
            else:
                clean_reads += 1
                clean_bases += len(tmp[9])
                if not (flag & 0x4):
                    mapped_reads += 1
                    if mappedPattern.search(line):
                        mapped_bases += int(mappedPattern.search(line).group(1))
                    else:
                        mapped_bases += len(tmp[9])
                    if flag & 0x400 == 0:
                        if uniqPattern.search(line):
                            if int(uniqPattern.search(line).group(1)) == 1:
                                uniq_reads += 1
                            uniq_bases += len(tmp[9])
                        elif int(tmp[4]) > 0:
                            uniq_reads += 1
                            uniq_bases += len(tmp[9])
                        if mismatchPattern1.search(line):
                            mismatch_bases += int(mismatchPattern1.search(line).group(1))
                        elif mismatchPattern2.search(line):
                            mismatch_bases += self.CoutMismathNo(line)
                    else:
                        dup_reads += 1 
        mismatch_rate = mismatch_bases/mapped_bases
        mapping_rate = mapped_reads/clean_reads
        dup_rate = dup_reads/mapped_reads
        uniq_rate = uniq_reads/mapped_reads
        name = os.path.basename(bam)
        sample = name.split('.')[0]
        print ("Sample\t%s\n"%(sample))
        print ("Clean reads\t%s\n"%(clean_reads))
        print ("Clean bases (Mb)\t%.2f\n"%(clean_bases/1000000))
        print ("Mapping rate (%%)\t%.2f\n"%(100*mapping_rate))
        print ("Unique rate (%%)\t%.2f\n"%(100*uniq_rate))
        print ("Duplicate rate (%%)\t%.2f\n"%(100*dup_rate))
        print ("Mismatch rate (%%)\t%.2f\n"%(100*mismatch_rate))

    def CoutMismathNo(self,align_result):
        md_value = ''
        pattern = re.compile(r'MD:Z:(\S+)')
        if pattern.search(align_result):
            md_value =  pattern.search(align_result).group(1)
        else:
            print (align_result)
        mark = 0 
        mismatch_no = 0
        i = 0
        pattern_word = re.compile('r[A-Z]')
        pattern_int = re.compile(r'\d')
        while i < len(md_value):
            if md_value[i] == "^":
                mark = 1
            elif pattern_word.search(md_value[i]) and mark == 0:
                mismatch_no += 1
            elif pattern_int.search(md_value[i]):
                mark = 0
            i += 1
        return mismatch_no
                                

    def get_insertsize(self,bam):
       #insert = collections.defaultdict()
       insert = dict()
       insertIN = os.popen("%s view -f 0x2 -F 0x80 %s"%(self.samtools,bam))
       for line in insertIN.readlines():
           line = line.strip()
           if not line.startswith("@"):
               tmp = line.split()
               key = int(fabs(int(tmp[8])))
               if( key == 0 or key >= 700 ):
                   continue 
               else:
                   if insert.__contains__(key):
                       insert[key] += 1
                   else:
                       insert[key] = 1
       insertOUT = open("%s/insertOUT"%(self.outdir),'w')
       for key in sorted(insert.keys()):
           insertOUT.writelines(str(key)+"\t"+str(insert[key])+"\n")
       insertOUT.close() 

    def snp_stat(self,snpAnnovarResult):
        """ 1.snp count table including homozygous loci and heterozygosis loci
            2.snp position CDS,gene,mRNA sys,nonsys, bar plot
        """
        total,hom,het = 0,0,0
        with open(snpAnnovarResult,'r') as F:
            for line in F.readlines():
                lineArray = line.strip().split()
                total += 1
                if (lineArray[5] == "hom"):
                    hom += 1
                elif (lineArray[5] == "het"):
                    het += 1
            print ("SNP total")
    def indel_stat(self,indelAnnovarResult):
        """ 1. indel:insert and deletion 

        """
        pass
    def sv_stat(self):
        """ 1. structure more than 50bp
            2. breakdancer
	"""
        pass
    def sv_filter(self,inputCTX,outputCTX,Minl = 100,Maxl = 1000000,MinScore= 30,MinReads = 5):
        """ filt the result of Breakdancer(.ctx) 
            inputCTX: input the result of Breakdancer(.ctx)
            outputCTX: output the filter file(.ctx)
            Minl: the SV's Minimum size, default 100
            Maxl: the SV's Maximum size, default 1000000
            MinScore: set the minimum Score of SV, default 30
            MinReads: set the minimum number of the reads supported the SV, default 5
        """
        OUT = open(outputCTX,'w')
        IN = open(inputCTX,'r')
        for line in IN.readlines():
            line = line.strip('\n')
            if line.startswith("#"):
                OUT.writelines(line+"\n")
                continue
            Array = line.split("\t")
            chr1,pos1,chr2,pos2,svtype,size,score,num = Array[0],Array[1],Array[3],Array[4],Array[6],int(Array[7]),int(Array[8]),int(Array[9])
            if "chrUn" in chr1 or "hap" in chr1 or "random" in chr1 or "chrUn" in chr2 or "hap" in chr2 or "random" in chr2:
                continue
            if size < 0 :
                size = -size
            if ((size < Minl) or (size > Maxl) and ("DEL" in svtype or "INS" in svtype)) or ((score < MinScore) or (num < MinReads)):
                continue
            OUT.writelines(line+"\n")
        OUT.close()
        IN.close()

    def sv_anno(self,inputCTX,outputAnno):
        """ sv annotation """
        IN = open(inputCTX,'r')
        OUT = open(outputAnno,'w')
        for line in IN.readlines():
            line = line.strip('\n')
            if line.startswith("#"):
                continue
            Array = line.split("\t")
            chr1,start,chr2,end,svtype,size,num = Array[0],Array[1],Array[3],Array[4],Array[6],int(Array[7]),Array[9]
            if size < 0 :
                size = -size
            SVID = chr1+"_"+start+"_"+chr2+"_"+end+"_"+svtype
            otherpoint = chr2+":"+end
            result = "%s\t%s\t%s\t0\t0\t%s\ttype=%s\totherpoint=%s\t%d\t%s\n"%(chr1,start,start,SVID,svtype,otherpoint,size,num)
            OUT.writelines(result)
        IN.close()
        OUT.close()

    def cnv_stat(self):
        """ 1.SOAPcnv or CNVnator """
        pass





class circos():
    """ a class for snp and indel display using circos(perl package)"""
    PERL = "/hwfssz1/ST_BIGDATA/USER/zhusitao/Software/miniconda3/bin/perl"

    def __init__(self,genome):
        self.genome = genome
    def chromLen(self):
        fasta = ReadFasta(self.genome)
        fasta.readFasta()
        fastaLenDict=fasta.stat()
        for i in fastaLenDict:
            print (i,fastaLenDict[i])
        return fastaLenDict
    def get_snp_dict(self,VCF):
        snpDict = defaultdict(dict)
        with open(VCF,'r') as fh:
            for line in fh.readlines():
                line = line.strip()
                if line.startswith("#"):
                    continue
                else:
                    arr = lines.split("\t")
                    snpDict[arr[0]][arr[1]] = 0
        return snpDict
    def get_snp_density(self,win,output):
        fastaLenDict = self.chromLen()
        for chrom in fastaLenDict:
            end = fastaLenDict[chrom] - win
            i=0
            while(i<end):
                end_win = i + win
                i += win
#a = wgs_statistic()
#b = wgs_statistic("/ldfssz1/ST_BIGDATA/PMO/WORKFLOW/subworkflow_test_raw_data/wgs_test/test7/WGS/03.rmdup/ERR013153.chr11.sort.rmdup.bam")
#a.get_insertsize("/ldfssz1/ST_BIGDATA/PMO/WORKFLOW/subworkflow_test_raw_data/wgs_test/test7/WGS/02.bam/ERR013153.sort.bam")
#b.get_insertsize("/ldfssz1/ST_BIGDATA/PMO/WORKFLOW/subworkflow_test_raw_data/wgs_test/test7/WGS/03.rmdup/ERR013153.chr11.sort.rmdup.bam")
#a.picard_stat("/ldfssz1/ST_BIGDATA/PMO/WORKFLOW/subworkflow_test_raw_data/wgs_test/test7/WGS/02.bam/ERR013153.sort.bam")
#a.sv_filter("/ldfssz1/ST_BIGDATA/PMO/WORKFLOW/subworkflow_test_raw_data/wgs_test/total/WGS/10.Structure_Variation/NA12878-L3.chr6.ctx","zhusitao.filter.ctx")
#a.sv_anno("/ldfssz1/ST_BIGDATA/PMO/WORKFLOW/subworkflow_test_raw_data/wgs_test/total/WGS/10.Structure_Variation/NA12878-L3.chr6.ctx","zhusitao2.filter.ctx")
b = circos("test.fa")
b.chromLen()

