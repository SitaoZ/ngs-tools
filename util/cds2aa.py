import sys
from itertools import groupby


#At current only integrate the standard translate_table


def usage():
	USAGE ="""
Description:
	cds2aa.py  --  translate nucleic acid to amino acid
	The program should be running in python3.6
Example:
	python cds2aa.py cds.fa > prot.fa
	python -check cds.fa 
		"""
	print (USAGE)

def displaySeq(seqIn,number):
	seqOut = ''
	for i in range(1,len(seqIn)/number+1):
		start = (i-1)*number
		end = i*number
		seqOut += seqIn[start:end]+"\n"
	left = len(seqIn) % number
	remainder = seqIn[-left:]
	return seqOut + remainder

def cds2aa(seq):
	STANDART = {
		'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',  # Alanine
		'TGC': 'C', 'TGT': 'C',  # Cysteine
		'GAC': 'D', 'GAT': 'D',  # Aspartic Acid
		'GAA': 'E', 'GAG': 'E',  # Glutamic Acid
		'TTC': 'F', 'TTT': 'F',  # Phenylalanine
		'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',  # Glycine
		'CAC': 'H', 'CAT': 'H',  # Histidine
		'ATA': 'I', 'ATC': 'I', 'ATT': 'I',  # Isoleucine
		'AAA': 'K', 'AAG': 'K',  # Lysine
		'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L', 'TTA': 'L', 'TTG': 'L',  # Leucine
		'ATG': 'M',  # Methionine
		'AAC': 'N', 'AAT': 'N',  # Asparagine
		'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',  # Proline
		'CAA': 'Q', 'CAG': 'Q',  # Glutamine
		'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R', 'AGA': 'R', 'AGG': 'R',  # Arginine
		'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S', 'AGC': 'S', 'AGT': 'S',  # Serine
		'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',  # Threonine
		'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',  # Valine
		'TGG': 'W',  # Tryptophan
		'TAC': 'Y', 'TAT': 'Y',  # Tyrosine
		'TAA': 'U', 'TAG': 'U', 'TGA': 'U'  # Stop
	}
	seq = seq.replace(' ','').upper()
	prot = ''
	i = 3
	if len(seq) < 3:
		print ("seq is too short")
		sys.exit(1)
	while i <= len(seq):
		codon = seq[i-3:i]
		if STANDART.__contains__(codon):
			prot += STANDART[codon]
		else:
			prot += "X"
		i += 3		
	return prot

def checkCDS(seq):
	seq = seq.replace(' ','')
	seq = seq.upper()
	start,end,mid,triple = 0,0,0,0
	length = len(seq)
	if length % 3 == 0:
		triple = 1
	if seq.startswith("ATG"):
		start = 1
	if seq.endswith("TAA") or seq.endswith("TAG") or seq.endswith("TGA"):
		end = 1
	mid = 1
	for i in range(3,length-3,3):
		codon = seq[i-3:i]
		if codon == 'TGA' or codon == 'TAG' or codon == 'TAA':
			mid = 0
			print (i)
	return start,end,mid,triple

class FastaFile:
        def __init__(self,path):
                self.path = path
                self._map = {}
                self.__fasta_iter()
        def __str__(self):
                return self._map.__str__()
        def __fasta_iter(self):
                fh = open(self.path)
                faiter =(x[1] for x in groupby(fh,lambda line : line[0] == ">")) # lambda tells to use the first item(line[0] == ">") as the grouping key
                for header in faiter:
                        header = header.next()[1:].strip()
                        header = header.split()[0]
                        seq ="".join(s.strip() for s in faiter.next())
                        self._map[header] = seq

	
	
if __name__ == "__main__":
	if len(sys.argv) == 1:
		usage()
		sys.exit(1)
	if sys.argv[1] == "-check":
		cds = sys.argv[2]
		Fasta = FastaFile(cds)._map
		for header in sorted(Fasta.keys()):
			seq = Fasta[header]
			checkCDS(seq)
	else:
		cds = sys.argv[1]
		Fasta = FastaFile(cds)._map
		for header in sorted(Fasta.keys()):
			seq = Fasta[header]	
			print (">%s [translate_table: standard]"%header)
			print (displaySeq(cds2aa(seq),50))
