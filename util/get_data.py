#-*-coding:utf-8-*-
import sys
import re,os
import glob

def get_data_from_BGI(data_path,laneid):
	fq1list = glob.glob('%s/%s_*_1.fq.gz' % (data_path, laneid))
	fq1list_sort = sorted(fq1list, key=lambda item: re.search(re.compile(r'_(\d+)_1.'), item).group(1))
	fq2list = glob.glob('%s/%s_*_2.fq.gz' % (data_path, laneid))
	fq2list_sort = sorted(fq2list, key=lambda item: re.search(re.compile(r'_(\d+)_2.'), item).group(1))
	barcode = [ os.path.basename(bar).replace('_1.fq.gz','') for bar in fq1list_sort]
	for i in range(len(barcode)):
		print (barcode[i])
		print(fq1list_sort[i])
		print(fq2list_sort[2])


if __name__=='__main__':
	get_data_from_BGI(sys.argv[1],sys.argv[2])