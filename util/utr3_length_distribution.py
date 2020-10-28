#coding:utf-8 
import re
import sys
import numpy as np
import seaborn as sns
import pandas as pd 
from Bio.Seq import Seq
from Bio import SeqIO
from scipy import stats

transcript = "../Araport11_3_utr_20160617"
Length = dict()
for record in SeqIO.parse(transcript, 'fasta'):
    chrom = record.id
    seq = record.seq
    desc = record.description
    length = len(seq)
    Length[chrom] = length

dis = pd.DataFrame(Length.items(),columns=["transcript_id","length"])
dis.to_csv('5_id_length.csv',index=False)

ax = sns.distplot(dis["length"], rug=True, rug_kws={"color": "g"},
                  kde_kws={"color": "k", "lw": 2, "label": "3' UTR"},
                  hist_kws={"histtype": "step", "linewidth": 2,
                            "alpha": 1, "color": "g"})
ax.set(xlabel="3' UTR length",ylabel="Proportion") 
ax.set(title="The Distribution of 3' UTR of Arabidopsis thaliana")
s = pd.Series(dis['length'])
print(s.describe())
des_stat = s.describe()
print(des_stat[:8])
count,mean = map(int,des_stat[:2])
std = des_stat[2]
minc,q25,q50,q75,maxc = map(int,des_stat[3:8])
mode = stats.mode(s)[0][0]
ax.axvline(x=mode, color='r',alpha=0.3,linestyle="--")
#ax.text(800, 0.005,"count %d\n"%(count),fontsize=9)
#ax.text(2500, 0.0005,"count\nmode\nmean\nstd\nmin\nmax\n25%\n50%\n75%", fontsize=9)
#ax.text(3000, 0.0005,"{:5d}\n{:5d}\n{:5d}\n{:5.2f}\n{:5d}\n{:5d}\n{:5d}\n{:5d}\n{:5d}".format(count,mode,mean,std,minc,maxc,q25,q50,q75), fontsize=9)
p50 = np.percentile(np.array(s), 50)
p80 = np.percentile(np.array(s), 80)
p90 = np.percentile(np.array(s), 90)
print("np.percentile 50",p50)
print("np.percentile 80",p80)
print("np.percentile 90",p90)
ax = ax.get_figure()
ax.savefig('LengthDistribution3UTR.pdf')

'''
plot = sns.distplot(dis["length"]) 
# Highlight mean
plot.axvline(dis['length'].mean(), color='r')
# Highlight median
plot.axvline(dis['length'].median(), color='g')
#plot.set(xlim=(1,20)) 
plot.set(xlabel="5'UTR length",ylabel="Proportion") 
plot.set(title="The Distribution of 5' UTR of Arabidopsis thaliana")
plot_fig = plot.get_figure() 
plot_fig.savefig('LengthDistribution5UTR.png') 
'''
