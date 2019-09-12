#coding:utf-8
import os
import re
import sys
import glob
from optparse import OptionParser



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
		cmd1 = 'zcat {fq1}|/hwfssz1/ST_BIGDATA/USER/zhusitao/Software/seqtk/seqtk sample -s 11 - {percent} | /hwfssz1/ST_BIGDATA/USER/zhusitao/Software/pigz-2.4/pigz -p 4 > {cat_dir}/{laneid}_1.fq.gz'.format(fq1=fq1,percent=percent,cat_dir=cat_dir,laneid=laneid)
		cmd2 = 'zcat {fq2}|/hwfssz1/ST_BIGDATA/USER/zhusitao/Software/seqtk/seqtk sample -s 11 - {percent} | /hwfssz1/ST_BIGDATA/USER/zhusitao/Software/pigz-2.4/pigz -p 4 > {cat_dir}/{laneid}_2.fq.gz'.format(fq2=fq2,percent=percent,cat_dir=cat_dir,laneid=laneid)
	else:
		cmd1 = 'zcat {fq1}|/hwfssz1/ST_BIGDATA/USER/zhusitao/Software/pigz-2.4/pigz -p 4 > {cat_dir}/{laneid}_1.fq.gz'.format(fq1=fq1, cat_dir=cat_dir, laneid=laneid)
		cmd2 = 'zcat {fq2}|/hwfssz1/ST_BIGDATA/USER/zhusitao/Software/pigz-2.4/pigz -p 4 > {cat_dir}/{laneid}_2.fq.gz'.format(fq2=fq2, cat_dir=cat_dir, laneid=laneid)
	cmd_ln = ''
	os.mkdir(os.path.join(outdir,'report','fqdata'))
	for each in fq1list_sort:
		barcode = each.split('/')[-1].replace('_1.fq.gz','')
		cmd_ln += "mkdir -p {outdir}/report/fqdata/{barcode}\nln -s {fq} {outdir}/report/fqdata/{barcode}\n".format(outdir=outdir,barcode=barcode,fq=each)
	for each in fq2list_sort:
		barcode = each.split('/')[-1].replace('_2.fq.gz','')
		cmd_ln += "mkdir -p {outdir}/report/fqdata/{barcode}\nln -s {fq} {outdir}/report/fqdata/{barcode}\n".format(outdir=outdir,barcode=barcode,fq=each)

	with open(cat_dir+"/catfq.sh",'w') as F:
		F.writelines(cmd1+"\n")
		F.writelines(cmd2+"\n")
		F.writelines(cmd_ln+'\n')

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
	cmd_stat = "perl {jobbin}/soapnuke/soapnuke_v1.0/bin/soapnuke_stat.pl {qc_dir}/{laneid}/Basic_Statistics_of_Sequencing_Quality.txt {qc_dir}/{laneid}/Statistics_of_Filtered_Reads.txt > {qc_dir}/{laneid}/{laneid}.fqstat.xls".format(jobbin=jobBin,laneid=laneid,qc_dir=qc_dir)
	cmd_ln = "mkdir -p {outdir}/report/stat\nln -s {qc_dir}/{laneid}/{laneid}.fqstat.xls {outdir}/report/stat ".format(outdir=outdir,qc_dir=qc_dir,laneid=laneid)
	with open(qc_dir+"/qc.sh",'w') as F:
		F.writelines(export_path+"\n")
		F.writelines(cmd+"\n")
		F.writelines(cmd_stat+'\n')
		F.writelines(cmd_ln+'\n')

	clean_fq1 = "{qc_dir}/{laneid}/{laneid}.1.clean.fq.gz".format(qc_dir=qc_dir,laneid=laneid)
	clean_fq2 = "{qc_dir}/{laneid}/{laneid}.2.clean.fq.gz".format(qc_dir=qc_dir, laneid=laneid)
	return [clean_fq1,clean_fq2]

def bolt(outdir,clean_fq1,clean_fq2,laneid,platform):
	samplelist = os.path.join(outdir,'sample.list')
	bolt_dir = os.path.join(outdir,'step3_bolt')
	os.makedirs(bolt_dir)
	with open(samplelist,'w') as F:
		out = ','.join([laneid,laneid,laneid,platform,clean_fq1,clean_fq2])
		F.writelines(out+"\n")

	cmd = "MegaBOLT -type allnoindex -qaligner -list {sampleList} -outputDir {bolt_dir}/{laneid}".format(sampleList=samplelist,bolt_dir=bolt_dir,laneid=laneid)
	with open(bolt_dir+"/bolt.sh",'w') as F:
		F.writelines(cmd+"\n")
	return os.path.join(bolt_dir,laneid,'output.sortdup.bqsr.bam'),os.path.join(bolt_dir,laneid,'output.sortdup.bqsr.bam.HaplotypeCaller.vcf.gz')

def bolt_full(outdir,clean_fq1,clean_fq2,laneid,platform):
	samplelist = os.path.join(outdir, 'sample.list')
	bolt_dir = os.path.join(outdir, 'step3_bolt')
	os.makedirs(bolt_dir)
	os.makedirs(os.path.join(bolt_dir,laneid))
	with open(samplelist, 'w') as F:
		out = '\t'.join([laneid, clean_fq1, clean_fq2,"AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA","AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG"])
		F.writelines(out + "\n")
	cmd = "MegaBOLT-full {sampleList} --type WGS -outdir {bolt_dir}/{laneid}".format(
		sampleList=samplelist, bolt_dir=bolt_dir, laneid=laneid)
	with open(bolt_dir + "/bolt.sh", 'w') as F:
		F.writelines(cmd + "\n")
	return os.path.join(bolt_dir,laneid),os.path.join(bolt_dir,laneid,"04.GetReport",laneid,laneid+".bam"),os.path.join(bolt_dir,laneid,"04.GetReport",laneid,laneid+".vcf.gz")

def bamStat(bam,laneid,reference,jobBin,outdir):
	bamStat_dir = os.path.join(outdir,'step4_bamStat')
	os.makedirs(bamStat_dir)
	export_path = "export PATH=/hwfssz1/ST_BIGDATA/USER/zhusitao/Software/R-3.2.0/bin/:$PATH\nexport LD_LIBRARY_PATH=/ldfssz1/ST_BIGDATA/USER/zhusitao/wdl/bolt/02.soft/lib:$LD_LIBRARY_PATH\nexport PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:$PATH"
	cmd = """
	mkdir -p {outdir}/report/stat
	{jobbin}/java/jdk1.8.0_73/bin/java -jar {jobbin}/02.soft/CollectInsertSizeMetrics.jar I={bam} O={bamStat_dir}/{laneid}.CollectInsertSizeMetrics.txt H={bamStat_dir}/{laneid}.CollectInsertSizeMetrics.pdf VALIDATION_STRINGENCY=SILENT
	/ldfssz1/ST_BIGDATA/USER/zhusitao/software/miniconda3/envs/imagemagick/bin/convert {bamStat_dir}/{laneid}.CollectInsertSizeMetrics.pdf {bamStat_dir}/{laneid}.insertsize.png
	ln -s {bamStat_dir}/{laneid}.CollectInsertSizeMetrics.txt {outdir}/report/stat
	{jobbin}/02.soft/samtools stats {bam} > {bamStat_dir}/{laneid}.samtoolsstat.xls
	ln -s {bamStat_dir}/{laneid}.samtoolsstat.xls {outdir}/report/stat
	{jobbin}/java/jdk1.8.0_73/bin/java -jar {jobbin}/02.soft/picard.jar CollectGcBiasMetrics I={bam} O={bamStat_dir}/{laneid}.gc_bias_metrics.xls CHART={bamStat_dir}/{laneid}.gc_bias_metrics.pdf S={bamStat_dir}/{laneid}.Summary.xls R={reference} VALIDATION_STRINGENCY=SILENT
	/ldfssz1/ST_BIGDATA/USER/zhusitao/software/miniconda3/envs/imagemagick/bin/convert {bamStat_dir}/{laneid}.gc_bias_metrics.pdf {bamStat_dir}/{laneid}.gcbias.png 
	ln -s {bamStat_dir}/{laneid}.gc_bias_metrics.xls {outdir}/report/stat
	ln -s {bamStat_dir}/{laneid}.gc_bias_metrics.pdf {outdir}/report/stat
	ln -s {bamStat_dir}/{laneid}.gcbias.png {outdir}/report/stat
	perl /hwfssz1/BC_COM_P0/BC_RD_P1/DNA/DNA_Human_WGS/DNA_Human_WGS_2017b_v2/processBam/processBam_v1.0/bin/depthV2.0.pl {bam} {bamStat_dir} -hg hg19 {bamStat_dir}/{laneid}.depthstat.xls
	ln -s {bamStat_dir}/{laneid}.depthstat.xls {outdir}/report/stat
	ln -s {bamStat_dir}/histPlot.pdf {outdir}/report/stat
	ln -s {bamStat_dir}/histPlot.png {outdir}/report/stat
	ln -s {bamStat_dir}/histPlot.pdf.R {outdir}/report/stat
	ln -s {bamStat_dir}/cumuPlot.png {outdir}/report/stat
	ln -s {bamStat_dir}/cumuPlot.pdf.R {outdir}/report/stat
	ln -s {bamStat_dir}/cumuPlot.pdf {outdir}/report/stat
	""".format(jobbin=jobBin,bam=bam,bamStat_dir=bamStat_dir,laneid=laneid,outdir=outdir,reference=reference)
	with open(bamStat_dir+"/bamStat.sh",'w') as F:
		F.writelines(export_path+"\n")
		F.writelines(cmd+"\n")


def bamSplit(bolt_bam,laneid,jobBin,outdir,chromlist):
	split_dir = os.path.join(outdir,'step5_split',laneid)
	os.makedirs(split_dir)
	cmd = ''
	cmd_ln='mkdir -p {outdir}/report/bam\n'.format(outdir=outdir)
	for chrom in chromlist:
		each_dir = os.path.join(split_dir,chrom)
		os.makedirs(each_dir)
		cmd+="{jobbin}/samtools/samtools-1.3.1/samtools view -b -h {bolt_bam} {chrom} > {each_dir}/{laneid}_{chrom}.bam\n".format(jobbin=jobBin,bolt_bam=bolt_bam,chrom=chrom,each_dir=each_dir,laneid=laneid)
		cmd_ln+="ln -s {each_dir}/{laneid}_{chrom}.bam {outdir}/report/bam\n".format(outdir=outdir,each_dir=each_dir,laneid=laneid,chrom=chrom)
	with open(split_dir+"/split.sh",'w') as F:
		F.writelines(cmd+'\n')
		F.writelines(cmd_ln+'\n')



def vqsr(bolt_vcf,laneid,reference,jobBin,outdir):

	# all need file
	java = os.path.join(jobBin,"java/jdk1.8.0_73/bin/java")
	gatk = os.path.join(jobBin,"gatk/3.6/GenomeAnalysisTK.jar")
	samtools = os.path.join(jobBin,"samtools/samtools-1.3.1/samtools")

	SNPparameter1 = "--maxGaussians 4 -mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 "
	SNPparameter2 = "-an QD -an MQRankSum -an ReadPosRankSum -an FS -an DP "

	Indelparameter1 = "--maxGaussians 4 -std 8.0 -mode INDEL "
	Indelparameter2 = "-an QD  -an MQRankSum -an ReadPosRankSum -an FS -an DP "

	hapmap = "/ifs4/BC_PUB/biosoft/pipeline/DNA/DNA_Human_WGS/DNA_Human_WGS_2018a/Database/hg19/gatk/hapmap_3.3.hg19.vcf.gz"
	omin = "/ifs4/BC_PUB/biosoft/pipeline/DNA/DNA_Human_WGS/DNA_Human_WGS_2018a/Database/hg19/gatk/1000G_omni2.5.hg19.vcf.gz"
	G1000 = "/ifs4/BC_PUB/biosoft/pipeline/DNA/DNA_Human_WGS/DNA_Human_WGS_2018a/Database/hg19/gatk/1000G_phase1.snps.high_confidence.hg19.vcf.gz"
	dbsnp = "/ifs4/BC_PUB/biosoft/pipeline/DNA/DNA_Human_WGS/DNA_Human_WGS_2018a/Database/hg19/gatk/dbsnp_138.hg19.vcf.gz"
	indel_gold = "/ifs4/BC_PUB/biosoft/pipeline/DNA/DNA_Human_WGS/DNA_Human_WGS_2018a/Database/hg19/gatk/Mills_and_1000G_gold_standard.indels.hg19.vcf.gz"

	vqsr_dir = os.path.join(outdir,'step6_vqsr')
	os.makedirs(vqsr_dir)
	snp_dir = os.path.join(vqsr_dir,'SNP')
	os.makedirs(snp_dir)
	indel_dir = os.path.join(vqsr_dir,'INDEL')
	os.makedirs(indel_dir)
	os.makedirs(os.path.join(outdir,'java_tmp'))

	#raw_snp_indel_vcf=os.path.join(outdir,'bolt','output.sortdup.bqsr.bam.HaplotypeCaller.vcf.gz')
	raw_snp_indel_vcf=bolt_vcf
	os.makedirs(os.path.join(snp_dir,laneid))
	snprecalFile = os.path.join(snp_dir,laneid,laneid+".snp.recalFile")
	snptranchesFile = os.path.join(snp_dir,laneid,laneid+".snp.tranchesFile")
	snprscriptFile = os.path.join(snp_dir,laneid,laneid+".snp.rscriptFile")
	snp_vqsr_vcf =os.path.join(snp_dir,laneid,laneid+".snps.VQSR.vcf")
	snp_vqsr_filter_vcf = os.path.join(snp_dir,laneid,laneid+".snps.VQSR.filter.vcf")
	snp_vqsr_filtered_vcf = os.path.join(snp_dir,laneid,laneid+".snps.VQSR.filtered.vcf")
	snp_final_vcf = os.path.join(snp_dir,laneid,laneid+".snps.final.vcf")


	os.makedirs(os.path.join(indel_dir,laneid))
	indelrecalFile = os.path.join(indel_dir,laneid,laneid+".indel.recalFile")
	indeltranchesFile = os.path.join(indel_dir,laneid,laneid+".indel.tranchesFile")
	indelrscriptFile = os.path.join(indel_dir,laneid,laneid+".indel.rscriptFile")
	indel_vqsr_vcf = os.path.join(indel_dir,laneid,laneid+".indel.VQSR.vcf")
	indel_vqsr_filter_vcf = os.path.join(indel_dir,laneid,laneid+".indel.VQSR.filter.vcf")
	indel_vqsr_filtered_vcf = os.path.join(indel_dir,laneid,laneid+".indel.VQSR.filtered.vcf")
	indel_final_vcf = os.path.join(indel_dir,laneid,laneid+".indel.final.vcf")


	cmd_snp = """{java} -Xmx6g -Djava.io.tmpdir={outdir}/java_tmp -jar {gatk} -T VariantRecalibrator -R {reference} -input {raw_snp_indel_vcf} {SNPparameter1} {SNPparameter2} -resource:hapmap,known=false,training=true,truth=true,prior=15.0 {hapmap} -resource:omni,known=false,training=true,truth=true,prior=12.0 {omin} -resource:1000G,known=false,training=true,truth=false,prior=10.0 {G1000} -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {dbsnp} -recalFile {snprecalFile} -tranchesFile {snptranchesFile} -rscriptFile {snprscriptFile}
	{java} -Xmx6g -Djava.io.tmpdir={outdir}/java_tmp -jar {gatk} -T ApplyRecalibration -R {reference} -input {raw_snp_indel_vcf} --ts_filter_level 99.5 -recalFile {snprecalFile} -tranchesFile {snptranchesFile} -mode SNP -o {snp_vqsr_vcf}
	grep -E 'PASS|#' {snp_vqsr_vcf} > {snp_vqsr_filter_vcf}
	{java} -Xmx6g -Djava.io.tmpdir={outdir}/java_tmp -jar {gatk} -T VariantFiltration -R {reference} --filterExpression "QD < 2.0 || MQ < 40.0 || ReadPosRankSum < -8.0 || FS > 60.0 || HaplotypeScore > 13.0 || MQRankSum < -12.5" --filterName LowQualFilter --missingValuesInExpressionsShouldEvaluateAsFailing --logging_level ERROR --variant {snp_vqsr_filter_vcf} -o {snp_vqsr_filtered_vcf}
	grep -E 'PASS|#' {snp_vqsr_filtered_vcf} > {snp_final_vcf}""".format(outdir=outdir,java=java,gatk=gatk,reference=reference,raw_snp_indel_vcf=raw_snp_indel_vcf,SNPparameter1=SNPparameter1,SNPparameter2=SNPparameter2,hapmap=hapmap,omin=omin,G1000=G1000,dbsnp=dbsnp,snprecalFile=snprecalFile,snptranchesFile=snptranchesFile,snprscriptFile=snprscriptFile,snp_vqsr_vcf=snp_vqsr_vcf,snp_vqsr_filter_vcf=snp_vqsr_filter_vcf,snp_vqsr_filtered_vcf=snp_vqsr_filtered_vcf,snp_final_vcf=snp_final_vcf)

	cmd_indel = """{java} -Xmx6g -Djava.io.tmpdir={outdir}/java_tmp -jar {gatk} -T VariantRecalibrator -R {reference} -input {raw_snp_indel_vcf} {Indelparameter1} {Indelparameter2} -resource:mills,known=true,training=true,truth=true,prior=12.0 {indel_gold} -recalFile {indelrecalFile} -tranchesFile {indeltranchesFile} -rscriptFile {indelrscriptFile}
	{java} -Xmx6g -Djava.io.tmpdir={outdir}/java_tmp -jar {gatk} -T ApplyRecalibration -R {reference} -input {raw_snp_indel_vcf} --ts_filter_level 99.5 -recalFile {indelrecalFile} -tranchesFile {indeltranchesFile} -mode INDEL  -o {indel_vqsr_vcf}
	grep -E 'PASS|#' {indel_vqsr_vcf} > {indel_vqsr_filter_vcf}
	{java} -Xmx6g -Djava.io.tmpdir={outdir}/java_tmp -jar {gatk} -T VariantFiltration -R {reference} --filterExpression "QD < 2.0 || ReadPosRankSum < -8.0 || SOR > 10.0 || FS > 200.0 || HaplotypeScore > 13.0 || MQRankSum < -12.5" --filterName LowQualFilter --missingValuesInExpressionsShouldEvaluateAsFailing --logging_level ERROR --variant {indel_vqsr_filter_vcf} -o {indel_vqsr_filtered_vcf}
	grep -E 'PASS|#' {indel_vqsr_filtered_vcf} > {indel_final_vcf}""".format(outdir=outdir,java=java,gatk=gatk,reference=reference,raw_snp_indel_vcf=raw_snp_indel_vcf,Indelparameter1=Indelparameter1,Indelparameter2=Indelparameter2,indel_gold=indel_gold,indelrecalFile=indelrecalFile,indeltranchesFile=indeltranchesFile,indelrscriptFile=indelrscriptFile,indel_vqsr_vcf=indel_vqsr_vcf,indel_vqsr_filter_vcf=indel_vqsr_filter_vcf,indel_vqsr_filtered_vcf=indel_vqsr_filtered_vcf,indel_final_vcf=indel_final_vcf)
				
	cmd_ln = """mkdir -p {outdir}/report/vcf
	ln -s {raw_snp_indel_vcf} {outdir}/report/vcf
	ln -s {snp_final_vcf} {outdir}/report/vcf
	ln -s {indel_final_vcf} {outdir}/report/vcf
	""".format(outdir=outdir,raw_snp_indel_vcf=raw_snp_indel_vcf,snp_final_vcf=snp_final_vcf,indel_final_vcf=indel_final_vcf)
	with open(vqsr_dir+"/vqsr.sh",'w') as F:
		F.writelines(cmd_snp+"\n")
		F.writelines(cmd_indel+"\n")
		F.writelines(cmd_ln+"\n")









def main():
	"""
	        %prog [-options] outdir

	        This program for call variant using bolt
	        e.g: %prog -1 read1.fq -2 read2.fq -r reference outdir
	    """
	parser = OptionParser()
	parser.add_option('-i', '--fqpath', help="the read1 and read2 fastq")
	parser.add_option('-s','--laneid', help="the laneid")
	parser.add_option('-c','--chip',help="the chip id ")
	parser.add_option('-r','--reference',help="the reference genome")
	parser.add_option('-j','--jobbin',help="the jobbin path")
	parser.add_option('-n', '--percent', help="the percent of library sample size")
	parser.add_option('-q', '--queue', help="the queue id")
	parser.add_option('-p', '--project', help="the project id")
	opts, args = parser.parse_args()

	if len(args) < 1:
		print('\033[0;31;40m%s\033[0m' % "\nWARM: Please give all the input parameters \n")
		sys.exit(not parser.print_help())
	elif opts.fqpath == None:
		print('\033[0;31;40m%s\033[0m' % "\nWARM: -i Please give the fastq file \n")
		sys.exit(not parser.print_help())
	elif opts.reference == None:
		print('\033[0;31;40m%s\033[0m' % "\nWARM: -r Please give the reference file \n")
		sys.exit(not parser.print_help())
	elif opts.laneid == None:
		print('\033[0;31;40m%s\033[0m' % "\nWARM: -s Please give the laneid string \n")
		sys.exit(not parser.print_help())
	elif opts.jobbin == None:
		print('\033[0;31;40m%s\033[0m' % "\nWARM: -j Please give the jobbin path \n")
		sys.exit(not parser.print_help())
	elif opts.percent == None:
		print('\033[0;31;40m%s\033[0m' % "\nWARM: -n the percent of library sample size,range 0 to 1 \n")
		sys.exit(not parser.print_help())
	elif opts.queue == None:
		print('\033[0;31;40m%s\033[0m' % "\nWARM: -q Please give the jobbin path \n")
		sys.exit(not parser.print_help())
	elif opts.project == None:
		print('\033[0;31;40m%s\033[0m' % "\nWARM: -p Please give the project string \n")
		sys.exit(not parser.print_help())

	outdir = args[0]
	if not os.path.exists(outdir):
		os.makedirs(outdir,0o755)
		os.makedirs(os.path.join(outdir, 'report'))
	fqpath = opts.fqpath
	laneid = opts.laneid
	reference = opts.reference
	jobBin = opts.jobbin
	percent = opts.percent
	#step1
	final_fq1, final_fq2 = catfq(fqpath,laneid,outdir,percent)
	#step2
	clean_fq1, clean_fq2 = qc(final_fq1,final_fq2,laneid,jobBin,outdir)
	#step3
	full_dir,bam,vcf = bolt_full(outdir,clean_fq1,clean_fq2,laneid,"COMPLETE")
	#step4
	#bamStat(bam,laneid,reference,jobBin,outdir)
	#step5
	chromlist_hg19 = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM', 'chr1_gl000191_random', 'chr1_gl000192_random', 'chr4_gl000193_random', 'chr4_gl000194_random', 'chr7_gl000195_random', 'chr8_gl000196_random', 'chr8_gl000197_random', 'chr9_gl000198_random', 'chr9_gl000199_random', 'chr9_gl000200_random', 'chr9_gl000201_random', 'chr11_gl000202_random', 'chr17_gl000203_random', 'chr17_gl000204_random', 'chr17_gl000205_random', 'chr17_gl000206_random', 'chr18_gl000207_random', 'chr19_gl000208_random', 'chr19_gl000209_random', 'chr21_gl000210_random', 'chr4_ctg9_hap1', 'chr6_apd_hap1', 'chr6_cox_hap2', 'chr6_dbb_hap3', 'chr6_mann_hap4', 'chr6_mcf_hap5', 'chr6_qbl_hap6', 'chr6_ssto_hap7', 'chr17_ctg5_hap1', 'chrUn_gl000211', 'chrUn_gl000212', 'chrUn_gl000213', 'chrUn_gl000214', 'chrUn_gl000215', 'chrUn_gl000216', 'chrUn_gl000217', 'chrUn_gl000218', 'chrUn_gl000219', 'chrUn_gl000220', 'chrUn_gl000221', 'chrUn_gl000222', 'chrUn_gl000223', 'chrUn_gl000224', 'chrUn_gl000225', 'chrUn_gl000226', 'chrUn_gl000227', 'chrUn_gl000228', 'chrUn_gl000229', 'chrUn_gl000230', 'chrUn_gl000231', 'chrUn_gl000232', 'chrUn_gl000233', 'chrUn_gl000234', 'chrUn_gl000235', 'chrUn_gl000236', 'chrUn_gl000237', 'chrUn_gl000238', 'chrUn_gl000239', 'chrUn_gl000240', 'chrUn_gl000241', 'chrUn_gl000242', 'chrUn_gl000243', 'chrUn_gl000244', 'chrUn_gl000245', 'chrUn_gl000246', 'chrUn_gl000247', 'chrUn_gl000248', 'chrUn_gl000249']
	bamSplit(bam,laneid,jobBin,outdir,chromlist=chromlist_hg19)
	#step6
	vqsr(vcf,laneid,reference,jobBin,outdir)



if __name__=="__main__":
    main()
