import os,sys

"""
Author :Zhu Sitao
Date : 2018-4-35
Dest : for protein and protein interaction in string database
       https://string-db.org/cgi/input.pl
"""


database = "/ifs4/BC_PUB/biosoft/db/Pub/PPI/STRING/protein.sequences.v10.fa"
diamond = "/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAdenovo/RNA_RNAdenovo_2016a/PPI/../software/diamond"


def diamond(diamondPath,stringDatabase,fastaPath.M6Path):
	diamond_shell = "{diamond} blastx  --evalue  1e-5  --outfmt  6  -d {db} -a {fa} -o {out} --seg no  --threads 5  --max-target-seqs 1 --more-sensitive -b 0.5 --salltitles".format(diamond=diamondPath,db=stringDatabase,fa=fastaPath,out=outPath)
	return diamond_shell

def main(M6Path,stringDatabase,RESULT):
	gene_str = {}#string datbase relation
	result = open(RESULT,'w')
	best_gene = {}#get best result from m6
	result.write("gene1\tgene2\tprotein_cluster1\tprotein_cluster2\tscore\n")
	srtdb = {}
	protein_list = []#store unigene id which map to all database

	for line in open(protein_list_opt,'r'):
		if line.startswith('Unigene\tNr'):continue
		protein_list.append(line.strip().split('\t')[0])
	open(protein_list_opt,'r').close()
	
	#===============it had read protein_list_opt	
	candidate_dict = {}
	for line in open(diff_list_opt,'r'):
		line = line.strip().split("\t")
		if line[0] in protein_list:
			candidate_dict[line[0]] = line[1] #unigene foldchange
	#===============it had read diff_list_opt
	for line in open(M6Path,'r'):
		line = line.strip().split()
		if line[0] not in best_gene:
			best_gene[line[0]] = 1
			if line[0] in candidate_dict:#only store candidate gene
				gene_str.setdefault(line[1],[]).append(line[0])#394.NGR c10170:[unigene1,unigene2]
	os.system("date")
	for line in open(stringDatabase,'r'):
		line = line.strip().split()
		if line[0] in gene_str and line[1] in gene_str:
			for left in gene_str[line[0]]:
				for right in gene_str[line[1]]:
					result.write("%s\t%s\t%s\t%s\t%s\n"%(left,right,line[0],line[1],line[2]))
	result.close()
	os.system("date")
	
