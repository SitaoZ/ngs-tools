#coding:utf-8
import sys
import matplotlib.pyplot as plt 
from matplotlib_venn import venn2
sys.path.append("/home/zhusitao/ngs-tools/format/")
from vcf import VCF

vcf1 = sys.argv[1]
title1 = sys.argv[2]
vcf2= sys.argv[3]
title2 = sys.argv[4]

list1=[]
list2=[]

v1 = VCF(vcf1)
for record in v1.readVCF():
	chrom,pos,ref,ale = record.CHROM,record.POS,record.REF,record.ALT
	list1.append(chrom+pos+ref+ale)

v2 = VCF(vcf2)
for record in v2.readVCF():
	chrom,pos,ref,ale = record.CHROM,record.POS,record.REF,record.ALT
	list2.append(chrom+pos+ref+ale)
set1 = set(list1)
set2 = set(list2)
inter = set1&set2

inter_len = len(inter)
set1_len = len(set1)
set2_len = len(set2)

only1 = set1_len - inter_len
only2 = set2_len - inter_len

with open(title1+"_vs_"+title2+".table",'w') as F:
	F.writelines("Sample1\tSample2\tIntersection\tOnly_sample1\tOnly_sample2\n")
	F.writelines("\t".join(map(str,[set1_len,set2_len,inter_len,only1,only2]))+"\n")

plt.figure(figsize=(4, 4))
venn2(subsets=[set1,set2],set_labels=(title1, title2),set_colors=('r', 'g'))
plt.savefig("{t1}_vs_{t2}.png".format(t1=title1,t2=title2))
