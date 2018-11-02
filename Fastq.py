#/usr/bin/python3

import os,sys
import re
import gzip
from collections import defaultdict


class fastq(object):
    """a class deal fastq file
        new fastq name: @HISEQ:310:C5MH9ANXX:1:1101:3517:2043 2:N:0:TCGGTCAC
        old fastq name: @FCD1PB1ACXX:4:1101:1799:2201#GAAGCACG/2
    """
    def __init__(self,path):
        self.path = path
        self.patternNew = re.compile(r'(?P<id>\d):N:(\d):(?P<index>\S+)',re.VERBOSE)
        self.patternOld = re.compile(r'\@\S+\#(?P<index>\S+?)\/(?P<id>\d)',re.VERBOSE)

    def fastqIN(self):
        """ read fastq """
        if self.path.endswith("gz"):
            FH = gzip.open(self.path,'r')
        else:
            FH = open(self.path, 'r')
        for line in FH:
            yield line
        FH.close()

    def qualitySystem(self):
        """ quality system from 100 lines"""
        lineCount = 0
        checklines = ''
        for line in self.fastqIN():
            lineCount += 1
            if lineCount > 100:
                break
            else:
                if lineCount % 4 == 0:
                    checklines += line.strip()

        MAX = 0
        MIN = 100
        for one in map(ord,checklines):
            if one > MAX:
                MAX = one
            if one < MIN:
                MIN = one
        if MAX <= 74 and MIN < 59:
            # sanger and Ilumina 1.8+
            return "Phred+33"
        elif MAX > 73 and MIN >=64:
            # Illumina 1.3+ and Illumina 1.5+
            return "Phred+64"
        elif MAX >73 and MIN >= 59 and MIN < 64:
            # Solexa+64
            return "Solexa+64"
        else:
            return "Unknown score encoding"



    def indexSequence(self):
        indexDict = defaultdict(int)
        indexList = []
        for line in self.fastqIN():
            matchNew = self.patternNew.search(line)
            matchOld = self.patternOld.search(line)
            if matchNew:
                indexSequence = matchNew.group("index")
                indexDict[indexSequence] += 1
            if matchOld:
                indexSequence = matchOld.group("index")
                indexDict[indexSequence] += 1
        for key in indexDict.keys():
            indexList.append(key)
        return indexList

    def pairOrSingel(self):
        count, singel, pair = 0, 0, 0
        for line in self.fastqIN():
            matchNew = self.patternNew.search(line)
            matchOld = self.patternOld.search(line)
            if matchNew:
                count += 1
                if matchNew.group("id") == "1":
                    singel += 1
                else:
                    pair += 1
            if matchOld:
                count += 1
                if matchOld.group("id") == '1':
                    singel += 1
                else:
                    pair += 1

        if count == singel:
            return "singel end"
        if count == pair :
            return "pair end"


    def to_fatsta(self,output_fasta):
        OUT=open(output_fasta,'w')
        flag = 0
        for line in self.fastqIN():
            flag += 1
            if flag == 1:
                OUT.writelines(">"+line) ## ID
            if flag == 2:
                OUT.writelines(line) ## sequence
            if flag == 4:
                flag = 0             ## break
        OUT.close()




