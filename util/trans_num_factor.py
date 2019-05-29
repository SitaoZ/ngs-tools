#-*-coding:utf-8-*-

import sys
import os

def trans(input_file):
	DD={}
	R, S, L, FL, FR2, FR4, FR6, FR8 = [],[],[],[],[],[],[],[]
	Rt,St,Lt,FLt,FR2t,FR4t,FR6t,FR8t = [],[],[],[],[],[],[],[]
	with open(input_file,'r') as F:
		for line in F.readlines():
			if 'mi' not in line:
				continue
			r,s,l,fl,fr2,fr4,fr6,fr8 = map(float,line.strip().split()[1:])

			R.append(r)
			S.append(s)
			L.append(l)
			FL.append(fl)
			FR2.append(fr2)
			FR4.append(fr4)
			FR6.append(fr6)
			FR8.append(fr8)
	DD['R']=R
	DD['S'] = S
	DD['L'] = L
	DD['FL'] = FL
	DD['FR2'] = FR2
	DD['FR4'] = FR4
	DD['FR6'] = FR6
	DD['FR8'] = FR8

	FF ={}

	for key,value in DD.items():
		FF[key] = []
		for i in value:
			if i >100:
				FF[key].append('>100')
			elif i <= 20:
				FF[key].append('0~20')
			elif 20< i <=40:
				FF[key].append('20~40')
			elif 40 < i <= 60:
				FF[key].append('40~60')
			elif 60 < i <= 80:
				FF[key].append('60~80')
			elif 80 < i <= 100:
				FF[key].append('80~100')

	print("\t".join(['Tissue','type','count']))
	for key,value in FF.items():
		for t in ['0~20','20~40','40~60','60~80','80~100','>100']:
			print ("\t".join([key,t,str(value.count(t))]))
if __name__=="__main__":
	trans(sys.argv[1])
