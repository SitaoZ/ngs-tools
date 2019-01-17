import os,sys
from itertools import groupby

"""
Author : Zhu Sitao
Date :   2018-4-14
"""


class Fasta(object):
    """ a class for Fasta file fasta operate: storage, index, item, export """
    def __init__(self,filePath):
        self.path = filePath
        self._fasta = dict()

    def __str__(self):
        """ Print the fasta dict """
        return self._fasta.__str__()

    def __len__(self):
        """ Total Length of the fasta file """
        totalLength = 0
        for ID in self._fasta:
            totalLength += len(self._fasta[ID])
        return totalLength

    def __iter__(self):
        """ Supports traversal with a for loop ,for ID loop """
        return iter(self._fasta)

    def __getitem__(self, index):
        """ Index one fasta ID and sequence"""
        outPut = self.fastaFormatOut(index, self._fasta[index], number=50)
        return outPut

    def fastaStorage(self):
        """ a function for fasta file store in dict """

        fh = open(self.path)
        faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
        for header in faiter:
            header = header.__next__()[1:].strip()
            header = header.split()[0]
            seq = "".join(s.strip() for s in faiter.__next__())
            self._fasta[header] = seq

    def fastaIndex(self):
        """ a function for fasta index """
        pass

    def fastaExport(self,filePath,number = 60):
        """ a function for fasta formatted output """
        fileOut = open(filePath,'w')
        for ID in self._fasta:
            seqIn = self._fasta[ID]
            fileOut.writelines(self.fastaFormatOut(ID, seqIn, number = number))
        fileOut.close()
    def fastaUpper(self,sequence):
        """ sequence upper output """
        return sequence.upper()

    def fastaFormatOut(self,index,sequence,number = 60):
        """ output in specical line length """
        seqOut = ''
        for i in range(1, len(sequence) // number + 1):
            start = (i - 1) * number
            end = i * number
            seqOut += sequence[start:end] + "\n"
        left = len(sequence) % number
        remainder = sequence[-left:] if left != 0 else ""
        totalOut = ">" + index + "\n" + seqOut + remainder+ "\n"
        return totalOut

    def fastaLengthStat(self):
        """ Statistic of each id length """
        for ID in self._fasta:
            print (ID, len(self._fasta[ID]))

    def fastaGC(self):
        """ Show each fasta GC and GC rate"""
        print ("ID\tGC\tGCrate(%)\tN\tNrate(%)")
        totalGC, totalN = 0, 0
        for ID in self._fasta:
            seqIn = self.fastaUpper(self._fasta[ID])  ## upper all base,especially for genome fasta file
            GC = seqIn.count("G") + seqIn.count("C")
            N = seqIn.count("N")
            totalGC += GC
            totalN += N
            GCrate = GC / len(seqIn)
            Nrate = N / len(seqIn)
            print ("%s\t%d\t%.4f\t%d\t%.4f" % (ID, GC, GCrate, N, Nrate))
        totalGCrate = totalGC / len(self)
        totalNrate = totalN / len(self)
        print ("total\t%d\t%.4f\t%d\t%.4f" % (totalGC, totalGCrate, totalN, totalNrate))

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
