#coding:utf-8
import os
import sys
import glob
from optparse import OptionParser


def catfq(fqpath,sampleid,outdir):
	cat_dir = os.path.join(outdir,'catfq')
	os.mkdir(cat_dir)
	fq1list = glob.glob('%s_*_1.fq.gz'%(fqpath))
	fq2list = glob.glob('%s_*_2.fq.gz'%(fqpath))
	fq1 = " ".join(fq1list)
	fq2 = " ".join(fq2list)
	cmd1 = 'zcat {fq1}|/hwfssz1/ST_BIGDATA/USER/zhusitao/Software/seqtk/seqtk sample -s 11 - 0.999 | /hwfssz1/ST_BIGDATA/USER/zhusitao/Software/pigz-2.4/pigz -p 4 > {outdir}/{sampleid}_1.fq.gz'.format(fq1=fq1,outdir=outdir,sampleid=sampleid)
	cmd2 = 'zcat {fq2}|/hwfssz1/ST_BIGDATA/USER/zhusitao/Software/seqtk/seqtk sample -s 11 - 0.999 | /hwfssz1/ST_BIGDATA/USER/zhusitao/Software/pigz-2.4/pigz -p 4 > {outdir}/{sampleid}_1.fq.gz'.format(fq2=fq2,outdir=outdir,sampleid=sampleid)
	with open(cat_dir+"/catfq.sh",'w') as F:
		F.writelines(cmd1+"\n")
		F.writelines(cmd2+"\n")
	final_fq1 = '{outdir}/{sampleid}_1.fq.gz'.format(outdir=cat_dir,sampleid=sampleid)
	final_fq2 = '{outdir}/{sampleid}_2.fq.gz'.format(outdir=cat_dir,sampleid=sampleid)
	return [final_fq1,final_fq2]

def qc(fq1,fq2,sampleid,jobBin,outdir):
	qc_dir = os.path.join(outdir,'qc')
	os.mkdir(qc_dir)
	export_path = """export LIBRARY_PATH=/hwfssz1/BIGDATA_COMPUTING/gongchun/zlib/lib:$LIBRARY_PATH\nexport LD_LIBRARY_PATH=/hwfssz1/BIGDATA_COMPUTING/gongchun/zlib/lib:$LD_LIBRARY_PATH """
	cmd += """{jobbin}/soapnuke/soapnuke_v2.0/SOAPnuke2-master/SOAPnuke filter -n 0.1 -q 0.5 -T 4 -l 12 -Q 2 -G 2 -M 2 -f AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA -r AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG 
			-1 {fq1} 
			-2 {fq2} 
			-o {qc_dir}/{sampleid} 
			-C {sampleid}.1.clean.fq 
			-D {sampleid}.2.clean.fq """.format(jobbin=jobBin,fq1=fq1,fq2=fq2,sampleid=sampleid,qc_dir=qc_dir)
	with open(qc_dir+"/qc.sh",'w') as F:
		F.writelines(export_path+"\n")
		F.writelines(cmd+"\n")
	clean_fq1 = "{qc_dir}/{sampleid}/{sampleid}.1.clean.fq".format(qc_dir=qc_dir,sampleid=sampleid)
	clean_fq2 = "{qc_dir}/{sampleid}/{sampleid}.2.clean.fq".format(qc_dir=qc_dir, sampleid=sampleid)
	return [clean_fq1,clean_fq2]

def bolt(outdir,clean_fq1,clean_fq2,sampleid,platform):
	samplelist = os.path.join(outdir,'sample.list')
	bolt_dir = os.path.join(outdir,'bolt')
	os.makedirs(bolt_dir)
	with open(samplelist,'w') as F:
		out = ','.join([sampleid,sampleid,sampleid,platform,clean_fq1,clean_fq2])
		F.writelines(outdir+"\n")

	cmd = "MegaBOLT -type allnoindex -qaligner -list {sampleList} -outputDir {bolt_dir}/{sampleid}".format(sampleList=samplelist,bolt_dir=bolt_dir,sampleid=sampleid)
	with open(bolt_dir+"/bolt.sh",'w') as F:
		F.writelines(cmd+"\n")
	return os.path.join(bolt_dir,sampleid,'output.sortdup.bqsr.bam'),os.path.join(bolt_dir,sampleid,'output.sortdup.bqsr.bam.HaplotypeCaller.vcf.gz')

def bamStat(bam,sampleid,reference,jobBin,outdir):
	bamStat_dir = os.path.join(outdir,'bamStat')
	os.makedirs(stat_dir)
	export_path = "export PATH=/hwfssz1/ST_BIGDATA/USER/zhusitao/Software/R-3.2.0/bin/:$PATH\nexport LD_LIBRARY_PATH=/ldfssz1/ST_BIGDATA/USER/zhusitao/wdl/bolt/02.soft/lib:$LD_LIBRARY_PATH\nexport PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:$PATH"
	cmd = """
	{jobbin}/java/jdk1.8.0_73/bin/java -jar {jobbin}/02.soft/CollectInsertSizeMetrics.jar I={bam} O={bamStat_dir}/{sampleid}.CollectInsertSizeMetrics.txt H={bamStat_dir}/{sampleid}.CollectInsertSizeMetrics.pdf VALIDATION_STRINGENCY=SILENT
	/ldfssz1/ST_BIGDATA/USER/zhusitao/software/miniconda3/envs/imagemagick/bin/convert {bamStat_dir}/{sampleid}.CollectInsertSizeMetrics.pdf {bamStat_dir}/{sampleid}.insertsize.png
	ln -s {bamStat_dir}/{sampleid}.CollectInsertSizeMetrics.txt {outdir}/report/stat
	{jobbin}/02.soft/samtools stats {bam} > {bamStat_dir}/{sampleid}.samtoolsstat.xls
	{jobbin}/java/jdk1.8.0_73/bin/java -jar {jobbin}/02.soft/picard.jar CollectGcBiasMetrics I={bam} O={bamStat_dir}/{sampleid}.gc_bias_metrics.xls CHART={bamStat_dir}/{sampleid}.gc_bias_metrics.pdf S={bamStat_dir}/{sampleid}.Summary.xls R={reference} VALIDATION_STRINGENCY=SILENT
	/ldfssz1/ST_BIGDATA/USER/zhusitao/software/miniconda3/envs/imagemagick/bin/convert {bamStat_dir}/{sampleid}.gc_bias_metrics.pdf {bamStat_dir}/{sampleid}.gcbias.png 
	ln -s {bamStat_dir}/{sampleid}.gc_bias_metrics.xls {outdir}/report/stat
	perl {jobbin}/02.soft/depthV2.0.pl {bam} {bamStat_dir} -b {jobbin}/02.soft/nonN_region > {bamStat_dir}/{sampleid}.depthstat.xls
	ln -s {bamStat_dir}/{sampleid}.depthstat.xls {outdir}/report/stat
	""".format(jobbin=jobBin,bam=bam,bamStat_dir=bamStat_dir,sampleid=sampleid,outdir=outdir,reference=reference)
	with open(bamStat_dir+"/bamStat.sh",'w') as F:
		F.writelines(export_path+"\n")
		F.writelines(cmd+"\n")


def bamSplit(bolt_bam,sampleid,jobBin,outdir,chromlist):
	split_dir = os.path.join(outdir,'split',sampleid)
	os.makedirs(split_dir)
	cmd = ''
	for chrom in chromlist:
		each_dir = os.path.join(split_dir,chrom)
		os.makedirs(each_dir)
		cmd+="{jobbin}/samtools/samtools-1.3.1/samtools view -b -h {bolt_bam} {chrom} > {each_dir}/${sampleid}.${chrom}.sort.bam\n".format(jobbin=jobBin,bolt_bam=bolt_bam,chrom=chrom,each_dir=each_dir,sampleid=sampleid)
	with(split_dir+"/split.sh",'w') as F:
		F.writelines(cmd+'\n')


def vqsr(bolt_vcf,sampleid,jobBin,outdir):

	# all need file
	java = os.path.join(jobBin,"/java/jdk1.8.0_73/bin/java")
	gatk = os.path.join(jobBin,"/gatk/3.6/GenomeAnalysisTK.jar")
	samtools = os.path.join(jobBin,"/samtools/samtools-1.3.1/samtools")

	SNPparameter1 = "--maxGaussians 4 -mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 "
	SNPparameter2 = "-an QD -an MQRankSum -an ReadPosRankSum -an FS -an DP "

	Indelparameter1 = "--maxGaussians 4 -std 8.0 -mode INDEL "
	Indelparameter2 = "-an QD  -an MQRankSum -an ReadPosRankSum -an FS -an DP "

	hapmap = "/ifs4/BC_PUB/biosoft/pipeline/DNA/DNA_Human_WGS/DNA_Human_WGS_2018a/Database/hg19/gatk/hapmap_3.3.hg19.vcf.gz"
	omin = "/ifs4/BC_PUB/biosoft/pipeline/DNA/DNA_Human_WGS/DNA_Human_WGS_2018a/Database/hg19/gatk/1000G_omni2.5.hg19.vcf.gz"
	G1000 = "/ifs4/BC_PUB/biosoft/pipeline/DNA/DNA_Human_WGS/DNA_Human_WGS_2018a/Database/hg19/gatk/1000G_phase1.snps.high_confidence.hg19.vcf.gz"
	dbsnp = "/ifs4/BC_PUB/biosoft/pipeline/DNA/DNA_Human_WGS/DNA_Human_WGS_2018a/Database/hg19/gatk/dbsnp_138.hg19.vcf.gz"
	indel_gold = "/ifs4/BC_PUB/biosoft/pipeline/DNA/DNA_Human_WGS/DNA_Human_WGS_2018a/Database/hg19/gatk/Mills_and_1000G_gold_standard.indels.hg19.vcf.gz"

	vqsr_dir = os.path.join(outdir,'vqsr')
	os.makedirs(vqsr_dir)
	snp_dir = os.path.join(vqsr_dir,'SNP')
	os.makedirs(snp_dir)
	indel_dir = os.path.join(vqsr_dir,'INDEL')
	os.makedirs(indel_dir)

	raw_snp_indel_vcf=os.path.join(outdir,bolt,'output.sortdup.bqsr.bam.HaplotypeCaller.vcf.gz')
	os.makedirs(snp_dir,sampleid)
	snprecalFile = os.path.join(snp_dir,sampleid,sampleid+".snp.recalFile")
	snptranchesFile = os.path.join(snp_dir,sampleid,sampleid+".snp.tranchesFile")
	snprscriptFile = os.path.join(snp_dir,sampleid,sampleid+".snp.rscriptFile")
	snp_vqsr_vcf =os.path.join(snp_dir,sampleid,sampleid+".snps.VQSR.vcf")
	snp_vqsr_filter_vcf = os.path.join(snp_dir,sampleid,sampleid+".snps.VQSR.filter.vcf")
	snp_vqsr_filtered_vcf = os.path.join(snp_dir,sampleid,sampleid+".snps.VQSR.filtered.vcf")
	snp_final_vcf = os.path.join(snp_dir,sampleid,sampleid+".snps.final.vcf")


	os.makedirs(indel_dir, sampleid)
	indelrecalFile = os.path.join(indel_dir,sampleid,sampleid+".indel.recalFile")
	indeltranchesFile = os.path.join(indel_dir,sampleid,sampleid+".indel.tranchesFile")
	indelrscriptFile = os.path.join(indel_dir,sampleid,sampleid+".indel.rscriptFile")
	indel_vqsr_vcf = os.path.join(indel_dir,sampleid,sampleid+".indel.VQSR.vcf")
	indel_vqsr_filter_vcf = os.path.join(indel_dir,sampleid,sampleid+".indel.VQSR.filter.vcf")
	indel_vqsr_filtered_vcf = os.path.join(indel_dir,sampleid,sampleid+".indel.VQSR.filtered.vcf")
	indel_final_vcf = os.path.join(indel_dir,sampleid,sampleid+".indel.final.vcf")


	cmd_snp = """{java} -Xmx6g -Djava.io.tmpdir=${outdir}/java_tmp -jar {gatk} -T VariantRecalibrator -R {reference} -input {raw_snp_indel_vcf} {SNPparameter1} {SNPparameter2} -resource:hapmap,known=false,training=true,truth=true,prior=15.0 {hapmap} -resource:omni,known=false,training=true,truth=true,prior=12.0 {omin} -resource:1000G,known=false,training=true,truth=false,prior=10.0 {G1000} -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {dbsnp} -recalFile {snprecalFile} -tranchesFile {snptranchesFile} -rscriptFile {snprscriptFile}
	{java} -Xmx6g -Djava.io.tmpdir=${outdir}/java_tmp -jar {gatk} -T ApplyRecalibration -R {reference} -input {raw_snp_indel_vcf} --ts_filter_level 99.5 -recalFile {snprecalFile} -tranchesFile {snptranchesFile} -mode SNP -o {snp_vqsr_vcf}"
	"grep -E 'PASS|#' {snp_vqsr_vcf} > {snp_vqsr_filter_vcf}"
	{java} -Xmx6g -Djava.io.tmpdir={outdir}/java_tmp -jar {gatk} -T VariantFiltration -R {reference} --filterExpression "QD < 2.0 || MQ < 40.0 || ReadPosRankSum < -8.0 || FS > 60.0 || HaplotypeScore > 13.0 || MQRankSum < -12.5" --filterName LowQualFilter --missingValuesInExpressionsShouldEvaluateAsFailing --logging_level ERROR --variant {snp_vqsr_filter_vcf} -o {snp_vqsr_filtered_vcf}
	grep -E 'PASS|#' {snp_vqsr_filtered_vcf} > {snp_final_vcf}""".format(outdir=outdir,gatk=gatk,reference=reference,raw_snp_indel_vcf=raw_snp_indel_vcf,SNPparameter1=SNPparameter1,SNPparameter2=SNPparameter2,hapmap=hapmap,omin=omin,G1000=G1000,dbsnp=dbsnp,snprecalFile=snprecalFile,snptranchesFile=snptranchesFile,snprscriptFile=snprscriptFile,snp_vqsr_vcf=snp_vqsr_vcf,snp_vqsr_filter_vcf=snp_vqsr_filter_vcf,snp_vqsr_filtered_vcf=snp_vqsr_filtered_vcf,snp_final_vcf=snp_final_vcf)

	cmd_indel = """{java} -Xmx6g -Djava.io.tmpdir={outdir}/java_tmp -jar {gatk} -T VariantRecalibrator -R {reference} -input {raw_snp_indel_vcf} {Indelparameter1} {Indelparameter2} -resource:mills,known=true,training=true,truth=true,prior=12.0 {indel_gold} -recalFile {indelrecalFile} -tranchesFile {indeltranchesFile} -rscriptFile {indelrscriptFile}
	{java} -Xmx6g -Djava.io.tmpdir={outdir}/java_tmp -jar {gatk} -T ApplyRecalibration -R {reference} -input {raw_snp_indel_vcf} --ts_filter_level 99.5 -recalFile {indelrecalFile} -tranchesFile {indeltranchesFile} -mode INDEL  -o {indel_vqsr_vcf}
	grep -E 'PASS|#' {indel_vqsr_vcf} > {indel_vqsr_filter_vcf}
	{java} -Xmx6g -Djava.io.tmpdir={outdir}/java_tmp -jar {gatk} -T VariantFiltration -R {reference} --filterExpression "QD < 2.0 || ReadPosRankSum < -8.0 || SOR > 10.0 || FS > 200.0 || HaplotypeScore > 13.0 || MQRankSum < -12.5" --filterName LowQualFilter --missingValuesInExpressionsShouldEvaluateAsFailing --logging_level ERROR --variant {indel_vqsr_filter_vcf} -o {indel_vqsr_filtered_vcf}
	grep -E 'PASS|#' {indel_vqsr_filtered_vcf} > {indel_final_vcf}""".format(outdir=outdir,gatk=gatk,reference=reference,raw_snp_indel_vcf=raw_snp_indel_vcf,Indelparameter1=Indelparameter1,Indelparameter2=Indelparameter2,indel_gold=indel_gold,indelrecalFile=indelrecalFile,indeltranchesFile=indeltranchesFile,indelrscriptFile=indelrscriptFile,indel_vqsr_vcf=indel_vqsr_vcf,indel_vqsr_filter_vcf=indel_vqsr_filter_vcf,indel_vqsr_filtered_vcf=indel_vqsr_filtered_vcf,indel_final_vcf=indel_final_vcf)
				
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
	        e.g: %prog -1 read1.fq -2 read2.fq -r reference -o outdir
	    """
	parser = OptionParser()
	parser.add_option('-i', '--fqpath', help="the read1 and read2 fastq")
	parser.add_option('-s','--sampleid', help="the sampleid")
	parser.add_option('-r','--reference',help="the reference genome")
	parser.add_option('-j','--jobbin',help="the jobbin path")
	parser.add_option('-q', '--queue', help="the queue id")
	parser.add_option('-p', '--project', help="the project id")
	opts, args = parser.parse_args()
	print("opts:", opts)
	print('args:', args)
	if len(args) < 1:
		print('\033[0;31;40m%s\033[0m' % "\nWARM: Please give all the input parameters \n")
		sys.exit(not parser.print_help())
	elif opts.fqpath == None:
		print('\033[0;31;40m%s\033[0m' % "\nWARM: Please give the fastq file \n")
		sys.exit(not parser.print_help())
	elif opts.reference == None:
		print('\033[0;31;40m%s\033[0m' % "\nWARM: Please give the reference file \n")
		sys.exit(not parser.print_help())
	elif opts.sampleid == None:
		print('\033[0;31;40m%s\033[0m' % "\nWARM: Please give the sampleid string \n")
		sys.exit(not parser.print_help())
	elif opts.jobbin == None:
		print('\033[0;31;40m%s\033[0m' % "\nWARM: Please give the jobbin path \n")
		sys.exit(not parser.print_help())
	elif opts.queue == None:
		print('\033[0;31;40m%s\033[0m' % "\nWARM: Please give the jobbin path \n")
		sys.exit(not parser.print_help())
	elif opts.project == None:
		print('\033[0;31;40m%s\033[0m' % "\nWARM: Please give the project string \n")
		sys.exit(not parser.print_help())

	outdir = args[0]
	if not os.path.exists(outdir):
		os.makedirs(outdir,0o755)
	fqpath = opts.fqpath
	sampleid = opts.sampleid
	reference = opts.reference
	jobBin = opts.jobbin
	final_fq1, final_fq2 = catfq(fqpath,sampleid,outdir)
	clean_fq1, clean_fq2 = qc(final_fq1,final_fq2,sampleid,jobBin,outdir)
	bam,vcf = bolt(outdir,clean_fq1,clean_fq2,sampleid,"COMPLETE")
	bamStat(bam,sampleid,reference,jobBin,outdir)
	bamSplit(bam,sampleid,jobBin,outdir,chromlist=['chr1','chr2','chr3'])
	vqsr(vcf,sampleid,jobBin,outdir)



if __name__=="__main__":
    main()