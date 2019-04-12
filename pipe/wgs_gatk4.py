import os
import re
import sys
from optparse import OptionParser


def make_dir(outdir,step):
    """ make dir for all step"""
    if not os.path.exists(outdir):
        sys.exit("%s not exist"%(outdir))
    else:
        os.makedirs(os.path.join(outdir,step),0o755)
    return os.path.join(outdir,step)

def write_shell(cmd,step,shell_dir):
    """ write the cmd to shell """
    F = open("{outdir}/{action}.sh".format(outdir=shell_dir,action=step),'w')
    F.write(cmd+"\n")
    F.close()
    return "{outdir}/{action}.sh".format(outdir=shell_dir,action=step)


def filter(sampleid,fq1,fq2,outdir):
    """filter fastq using soapnuke2 """
    soapnuke2 = "/ldfssz1/ST_BIGDATA/USER/zhusitao/software/soapnuke/soapnuke_v2.0/SOAPnuke2-master/SOAPnuke"
    filter_outdir = make_dir(outdir, 'filter')
    cmd = """{soapnuke} filter -n 0.1 -q 0.5 -T 4 -l 12 -Q 2 -G 2 -M 2 \\
            -f AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA -r AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG \\
            -1 {fq1} -2 {fq2} -C {clean_fq1} -D {clean_fq2} -o {filterDir}""".format(soapnuke=soapnuke2,
                                                                                  fq1=fq1,fq2=fq2,
                                                                                  clean_fq1=sampleid+".1.clean.fq",
                                                                                  clean_fq2=sampleid+".2.clean.fq",
                                                                                  filterDir=filter_outdir)
    cmd_shell = write_shell(cmd,'filter',outdir)
    return cmd_shell

def minimap2(clean_fq1,clean_fq2,reference,outdir):
    """ alignment read to genome using minimap2 """
    minimap2 = "/ldfssz1/ST_BIGDATA/USER/xujunhao/software/minimap2/minimap2"
    align_dir = make_dir(outdir,'minimap2')
    cmd = """{minimap2} -ax sr -t 4 {ref} {clean_fq1} {clean_fq2} -o {alignDir}""".format(minimap2=minimap2,
                                                                                        ref=reference,
                                                                                        clean_fq1=clean_fq1,
                                                                                        clean_fq2=clean_fq2,
                                                                                        alignDir=align_dir)
    cmd_shell = write_shell(cmd,'alignment',outdir)
    return cmd_shell

def rmdup(input_bam,outdir):
    """ rmdup using latest picards """
    java="/hwfssz1/ST_BIGDATA/USER/liuxing2/software/source/jdk1.8.0_131/bin/java"
    picards = "/ldfssz1/ST_BIGDATA/USER/zhusitao/software/picards/picard-2.18.29.jar"
    markdup_bam = "zhusitao"
    matrix = "zhusitao"
    cmd = """{java} -jar {picards} MarkDuplicates I={input_bam} O={output} M={matrix_txt}""".format(java=java,
                                                                                                    picards=picards,
                                                                                                    input_bam=input_bam,
                                                                                                    output=markdup_bam,
                                                                                                    matrix_txt=matrix)

    cmd_shell = write_shell(cmd,'markdup',outdir)

def bamStat(sortbam,reference):
    """ stat bam file samtools and picards """
    java="/hwfssz1/ST_BIGDATA/USER/liuxing2/software/source/jdk1.8.0_131/bin/java"
    picards = "/ldfssz1/ST_BIGDATA/USER/zhusitao/software/picards/picard-2.18.29.jar"
    stat_dir = make_dir(outdir,'bamStat')
    cmd_path = "export PATH=/hwfssz1/ST_BIGDATA/USER/zhusitao/Software/R-3.2.0/bin/:$PATH\n"
    # insert szie
    cmd_insert = "{java} -jar {picards} CollectInsertSizeMetrics I={bam} O={stat_dir}/CollectInsertSizeMetrics_txt H={stat_dir}/CollectInsertSizeMetrics_pdf VALIDATION_STRINGENCY=SILENT".format(java=java,picards=picards,bam=sortbam,stat_dir=stat_dir)
    # gc Bias
    cmd_gcBias = "{java} -jar {picards} CollectGcBiasMetrics  I={bam} O={stat_dir}/gc_bias_metrics_xls CHART={stat_dir}/gc_bias_metrics_pdf S={stat_dir}/Summary_xls R={reference} VALIDATION_STRINGENCY=SILENT\n".format(java=java,picards=picards,bam=sortbam,stat_dir=stat_dir,reference=reference)



def gatk4():
    """ get the variant file using gatk4 """
    gatk="/ldfssz1/ST_BIGDATA/USER/zhusitao/software/gatk/4.0/gatk"

def main():
    """
        %prog [-options]

        This program for call variant using gatk4
        e.g: %prog -1 read1.fq -2 read2.fq -r reference -o outdir
    """
    parser = OptionParser()
    parser.add_option('-1', '--fastq1', help="the read1 fastq")
    parser.add_option('-2', '--fastq2', help="the read2 fastq")
    parser.add_option('-r', '--reference', help="the reference")
    parser.add_option('-o', '--outDir', help="the out put dir ")
    opts, args = parser.parse_args()
    if len(args) < 5:
        print('\033[0;31;40m%s\033[0m' % "\nWARM: Please give all the input parameters \n")
        sys.exit(not parser.print_help())
    elif opts.reference == None:
        print('\033[0;31;40m%s\033[0m' % "\nWARM: Please give the reference file \n")
        sys.exit(not parser.print_help())

    referencefile = args[0]


if __name__=="__main__":
    main()