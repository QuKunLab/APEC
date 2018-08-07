#!/usr/bin/python
import os
import sys
from optparse import OptionParser
import subroutines
from multiprocessing import Pool
#
#
opts = OptionParser()
usage = "Map sequencing data\nusage: %prog -d data -w work --index bowtie2-index --picard picard.jar --tss TSS.txt --np 4"
opts = OptionParser(usage=usage, version="%prog 2.1")
opts.add_option("-d", help="data folder, the folder including all trimmed fastq files. "
                          +"These should be pair-end data, and named as: "
                          +"type1-001_1.trim.fastq, type1-001_1.trim.fastq, etc. "
                          +"{type1, type2, ......} can be composed by any character except (_) and (-). "
                          +"No need to organize them by yourself if you have run s01_trimming to get trimmed fastq files.")
opts.add_option("-w", help="work folder, defined by user.")
opts.add_option("--index", default="../reference/hg19",
                help="The reference file path in bowtie2/indexes folder, default=../reference/hg19")
opts.add_option("--picard", default="../reference/picard.jar",
                help="The picard.jar file path, default=../reference/picard.jar")
opts.add_option("--tss", default="../reference/hg19_refseq_genes_TSS.txt",
                help="The TSS file path, can be downloaded from RefSeq, default=../reference/hg19_refseq_genes_TSS.txt")
opts.add_option("--np", default=1, help="Number of CPUs used for mapping, default=1")
options, arguments = opts.parse_args()
#
#
#
def mapping(par):
    work_dir, cell, options, input1, input2, chr_list = par[0], par[1], par[2], par[3], par[4], par[5]
    sam = work_dir + cell + '.sam'
    bam = work_dir + cell + '.bam'
    log = work_dir + cell + '.map.log'
    sorted_bam = work_dir + cell + '.sorted.bam'
    filtered_bam = work_dir + cell + '.filtered.bam'
    marked_bam = work_dir + cell + '.marked.bam'
    removed_duplicate = work_dir + cell + '.dups.log'
    quality_state = work_dir + cell + '.stats.log'
    hist_log = work_dir + cell + '.hist.log'
    hist_pdf = work_dir + cell + '.hist.pdf'
    refSeqTSS = work_dir + cell + '.RefSeqTSS'
    print cell, work_dir
#
#    os.popen('bowtie2 -X2000 -p ' + options.np + ' --rg-id ' + cell + ' -x ' + options.index
#        + ' -1 ' + input1 + ' -2 ' + input2 + ' -S ' + sam + ' 2> ' + log)
    os.popen('bowtie2 -X2000 -p 1 --rg-id ' + cell + ' -x ' + options.index
        + ' -1 ' + input1 + ' -2 ' + input2 + ' -S ' + sam + ' 2> ' + log)
    os.popen('samtools view -bS ' + sam + ' -o ' + bam)
    os.popen('rm ' + sam)
    os.popen('java -XX:+UseSerialGC -Xmx1g -jar ' + options.picard + ' SortSam SO=coordinate VALIDATION_STRINGENCY=SILENT I='
        + bam + ' O=' + sorted_bam)
    os.popen('samtools index ' + sorted_bam)
    os.popen('samtools view -b -q 30 ' + sorted_bam + ' -o ' + filtered_bam + ' ' + chr_list)
    os.popen('java -XX:+UseSerialGC -Xmx1g -jar ' + options.picard + ' MarkDuplicates INPUT=' + filtered_bam +' OUTPUT='
        + marked_bam + ' METRICS_FILE=' + removed_duplicate
        + ' REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT')
    os.popen('samtools index ' + marked_bam)
    os.popen('echo -e "Chromosome\tLength\tProperPairs\tBadPairs:Raw" >> ' + quality_state)
    os.popen('samtools idxstats ' + sorted_bam + ' >> ' + quality_state)
#    os.popen('echo -e "Chromosome\tLength\tProperPairs\tBadPairs:Filtered" >> ' + quality_state)
#    os.popen('samtools idxstats ' + marked_bam + ' >> ' + quality_state)
    os.popen('java -XX:+UseSerialGC -Xmx1g -jar ' + options.picard + ' CollectInsertSizeMetrics VALIDATION_STRINGENCY=SILENT I='
        + marked_bam + ' O=' + hist_log + ' H=' + hist_pdf + ' W=1000')
    subroutines.draw_TSS_insert(options.tss, marked_bam, refSeqTSS)
    return
#
#
chr_list = 'chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY'
#
fastq_list = [x for x in os.listdir(options.d) if x[-11:]=='.trim.fastq']
fastq_list.sort()
cell_list = [x.split('_')[0] for x in fastq_list]
cell_list = list(set(cell_list))
cell_list.sort()
if not os.path.exists(options.w): os.popen('mkdir ' + options.w)
#
parameters = []
for cell in cell_list:
    work_dir = options.w + '/' + cell + '/'
    if os.path.exists(work_dir): os.popen('rm -rf' + work_dir)
    os.popen('mkdir ' + work_dir)
    input1, input2 = options.d+'/'+cell+'_1.trim.fastq', options.d+'/'+cell+'_2.trim.fastq'
    par = [work_dir, cell, options, input1, input2, chr_list]
    parameters.append(par)
#    mapping(par)
#
pool = Pool(int(options.np))
pool.map(mapping, parameters)
pool.close()
pool.join()
#
#
#
