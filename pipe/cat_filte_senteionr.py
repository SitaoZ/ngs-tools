#coding:utf-8
import os
import re
import sys
import glob



def catfq(fqpath,laneid,outdir,percent):
	cat_dir = os.path.join(outdir,'step1_catfq')
	os.mkdir(cat_dir)
	fq1list = glob.glob('%s/%s_*_1.fq.gz'%(fqpath,laneid))
	fq1list_sort = sorted(fq1list,key=lambda item:re.search(re.compile(r'_(\d+)_1.'),item).group(1))
	fq2list = glob.glob('%s/%s_*_2.fq.gz'%(fqpath,laneid))
	fq2list_sort = sorted(fq2list,key=lambda item:re.search(re.compile(r'_(\d+)_2.'),item).group(1))
	fq1 = " ".join(fq1list_sort)
	fq2 = " ".join(fq2list_sort)
	if float(percent) < 1:
		cmd1 = 'zcat {fq1}|/hwfssz1/ST_BIGDATA/USER/zhusitao/Software/seqtk/seqtk sample -s 11 - {percent} | /hwfssz1/ST_BIGDATA/USER/zhusitao/Software/pigz-2.4/pigz -p 8 > {cat_dir}/{laneid}_1.fq.gz'.format(fq1=fq1,percent=percent,cat_dir=cat_dir,laneid=laneid)
		cmd2 = 'zcat {fq2}|/hwfssz1/ST_BIGDATA/USER/zhusitao/Software/seqtk/seqtk sample -s 11 - {percent} | /hwfssz1/ST_BIGDATA/USER/zhusitao/Software/pigz-2.4/pigz -p 8 > {cat_dir}/{laneid}_2.fq.gz'.format(fq2=fq2,percent=percent,cat_dir=cat_dir,laneid=laneid)
	else:
		cmd1 = 'zcat {fq1}|/hwfssz1/ST_BIGDATA/USER/zhusitao/Software/pigz-2.4/pigz -p 8 > {cat_dir}/{laneid}_1.fq.gz'.format(
			fq1=fq1, cat_dir=cat_dir, laneid=laneid)
		cmd2 = 'zcat {fq2}|/hwfssz1/ST_BIGDATA/USER/zhusitao/Software/pigz-2.4/pigz -p 8 > {cat_dir}/{laneid}_2.fq.gz'.format(
			fq2=fq2, cat_dir=cat_dir, laneid=laneid)

	with open(cat_dir+"/fq1.sh",'w') as F:
		F.writelines(cmd1+"\n")
	with open(cat_dir+"/fq2.sh",'w') as F:
		F.writelines(cmd2+"\n")

	final_fq1 = '{outdir}/{laneid}_1.fq.gz'.format(outdir=cat_dir,laneid=laneid)
	final_fq2 = '{outdir}/{laneid}_2.fq.gz'.format(outdir=cat_dir,laneid=laneid)
	return [final_fq1,final_fq2]

def qc(fq1,fq2,laneid,jobBin,outdir):
	qc_dir = os.path.join(outdir,'step2_qc')
	os.makedirs(qc_dir)
	os.makedirs(os.path.join(qc_dir,laneid))
	export_path = """export LIBRARY_PATH=/hwfssz1/BIGDATA_COMPUTING/gongchun/zlib/lib:$LIBRARY_PATH\nexport LD_LIBRARY_PATH=/hwfssz1/BIGDATA_COMPUTING/gongchun/zlib/lib:$LD_LIBRARY_PATH """
	cmd = """{jobbin}/soapnuke/soapnuke_v2.0/SOAPnuke2-master/SOAPnuke filter -n 0.1 -q 0.5 -T 4 -l 12 -Q 2 -G 2 -M 2 -f AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA -r AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG \\
			-1 {fq1} \\
			-2 {fq2} \\
			-C {laneid}.1.clean.fq.gz -D {laneid}.2.clean.fq.gz \\
			-o {qc_dir}/{laneid}""".format(jobbin=jobBin,fq1=fq1,fq2=fq2,laneid=laneid,qc_dir=qc_dir)
	with open(qc_dir+"/qc.sh",'w') as F:
		F.writelines(export_path+"\n")
		F.writelines(cmd+"\n")
	clean_fq1 = "{qc_dir}/{laneid}/{laneid}.1.clean.fq.gz".format(qc_dir=qc_dir,laneid=laneid)
	clean_fq2 = "{qc_dir}/{laneid}/{laneid}.2.clean.fq.gz".format(qc_dir=qc_dir, laneid=laneid)
	return [clean_fq1,clean_fq2]
def senteion(clean_fq1,clean_fq2,laneid):
	cmd = """ export SENTIEON_LICENSE=10.54.31.175:8990
export SENTIEON_HOME=/hwfssz1/BIGDATA_COMPUTING/software/tools/sentieon/sentieon-genomics-201808.05/
export SENTIEON_INSTALL_DIR=/hwfssz1/BIGDATA_COMPUTING/software/tools/sentieon/sentieon-genomics-201808.05/
export PATH=$PATH:$SENTIEON_HOME/bin

sh /nascngb/gccnt2/ST_BI/zhusitao/25.Sentieon/sentieon_quickstart_FQtoGVCF.sh out {laneid} {clean1} {clean2} {laneid}""".format(laneid=laneid,clean1=clean_fq1,clean2=clean_fq2)
	with open('work.sh','w') as F:
		F.writelines(cmd+"\n")

if __name__ == "__main__":
	fqpath = sys.argv[1]
	laneid = sys.argv[2]
	outdir = sys.argv[3]
	percent = 1
	final_fq1,final_fq2 = catfq(fqpath,laneid,outdir,percent)
	clean_fq1,clean_fq2 = qc(final_fq1,final_fq2,laneid,"/ldfssz1/ST_BIGDATA/USER/zhusitao/software",outdir)
	senteion(clean_fq1,clean_fq2,laneid)
