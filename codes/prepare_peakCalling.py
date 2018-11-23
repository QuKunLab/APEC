#!/usr/bin/python
#
import warnings
warnings.filterwarnings("ignore")
#
import os
import numpy
import sys
from optparse import OptionParser
import subroutines
#
#
opts = OptionParser()
usage = "Call peaks\nusage: %prog -s project --blist blacklist.bed --fa genome_chr.fa --tss tssFile --logq 3"
opts = OptionParser(usage=usage, version="%prog 1.0")
opts.add_option("-s", help="The project folder.")
opts.add_option("--picard", default="../reference/picard.jar",
                help="The picard.jar file path, default=../reference/picard.jar")
opts.add_option("--blist", default='../reference/hg19_blacklist.JDB.bed', 
                help="Blacklist.bed, default=../reference/hg19_blacklist.JDB.bed")
opts.add_option("--fa", default='../reference/hg19_chr.fa', 
                help="Genome_chr.fa, default=../reference/hg19_chr.fa")
opts.add_option('--tss', default='../reference/hg19_refseq_genes_TSS.txt',
                help='TSS file, default=../reference/hg19_refseq_genes_TSS.txt')
opts.add_option('--ref', default='hg19', help='Name of genome reference, default=hg19')
opts.add_option('--logq', default='3',
                help='Threshold of -log(p-value) for top peaks, default=3.')
options, arguments = opts.parse_args()
#
workspace_folder = options.s + '/work/'
peak_folder = options.s + '/peak/'
genome_fasta = options.fa
tssFile = options.tss
os.popen('mkdir ' + peak_folder)
#
#
print '!!!!!!  merge all marked bam files  !!!!!!'
bam_folder = [x for x in os.listdir(workspace_folder)]
bam_folder.sort()
print 'cells number:', len(bam_folder)
marked_bam = []
merged_bam = peak_folder + 'mergeAll.bam'
for folder in bam_folder:
    path = workspace_folder + folder + '/'
    if len(folder.split('.'))<=1:
        marked_bam.extend([path + x for x in os.listdir(path) if x[-10:]=='marked.bam'])
if len(marked_bam)<=1000:
    marked_bam = ' '.join(marked_bam)
    os.popen('samtools merge -f ' + merged_bam + ' ' + marked_bam)
else:
    n_batch = len(marked_bam)//1000 + 1
    temps = []
    for i_batch in range(0, n_batch):
        temp_bam = peak_folder+'temp_'+str(i_batch)+'.bam'
        temps.append(temp_bam)
        start, end = i_batch*1000, min((i_batch+1)*1000, len(marked_bam))
        marked = ' '.join(marked_bam[start:end])
        os.popen('samtools merge -f ' + temp_bam + ' ' + marked)
        os.popen('samtools index ' + temp_bam)
    all_temp = ' '.join(temps)
    os.popen('samtools merge -f ' + merged_bam + ' ' + all_temp)
#        
os.popen('samtools index ' + merged_bam)
print '!!!!!!  merge done  !!!!!!'
print
#
hist_log = peak_folder + 'mergeAll.hist.log'
hist_pdf = peak_folder + 'mergeAll.hist.pdf'
os.popen('java -XX:+UseSerialGC -Xmx1g -jar '+options.picard+' CollectInsertSizeMetrics VALIDATION_STRINGENCY=SILENT I='
    + merged_bam + ' O=' + hist_log + ' H=' + hist_pdf + ' W=1000')
#
refSeqTSS = peak_folder + 'mergeAll.RefSeqTSS'
subroutines.draw_TSS_insert(tssFile, merged_bam, refSeqTSS)
#
print '!!!!!!  call peak by macs2  !!!!!!'
peak_file = peak_folder + 'peaks'
os.popen('macs2 callpeak --nomodel -t ' + merged_bam + ' -n ' 
    + peak_file + ' --nolambda --keep-dup all --call-summits')
print '!!!!!!  call peak done  !!!!!!'
print
#
summit = peak_folder + 'peaks_summits.bed'
filtered_peak = peak_folder + 'filtered_peaks.bed'
if options.blist:
    print '!!!!!!  filter peaks  !!!!!!'
    os.popen('bedtools intersect -v -a ' + summit + ' -b ' + options.blist
        + " | sort -k5 -nr > " + filtered_peak)
    print '!!!!!!  filter peaks done  !!!!!!'
    print
else:
    os.popen('sort -k5 -nr ' + summit + ' > ' + filtered_peak)

print '!!!!!!  get top N peaks by q-value  !!!!!!'
fold_rank = numpy.loadtxt(filtered_peak, 'str', delimiter='\t')
fold_rank[:, 1] = numpy.array(map(int, fold_rank[:, 1])) - 250
fold_rank[:, 2] = numpy.array(map(int, fold_rank[:, 2])) + 250
toppeaks = peak_folder + 'top_peaks.bed'
with open(toppeaks, 'w') as output:
    for peak in fold_rank:
        if float(peak[-1])>=float(options.logq):
            print >> output, peak[0]+'\t'+peak[1]+'\t'+peak[2]
#numpy.savetxt(toppeaks, fold_rank[:50000, :3], fmt='%s', delimiter='\t')
print '!!!!!!  get top peaks done  !!!!!!'
print
#
#
print '!!!!!!  annotate peaks  !!!!!!'
annotate_output = peak_folder + 'annotate_output.bed'
annotate_peak = peak_folder + 'annotate_peak.bed'
os.popen("annotatePeaks.pl " + toppeaks + " " + options.ref + " -size given > " + annotate_output)
temp01_file = peak_folder + 'temp01.bed'
with open(annotate_output) as annotate_file, open(temp01_file, 'w') as temp01:
    for i, line in enumerate(annotate_file):
        if i>0:
            words = line.split('\t')
            words = ['NA' if x=='' else x for x in words]
            leave = words[1:7] + words[9:13] + [words[15]]
            print >> temp01, '\t'.join(leave)
os.popen('bedtools sort -i ' + temp01_file + ' > ' + annotate_peak)
print '!!!!!!  annotate done  !!!!!!'
print
#
#
print '!!!!!!  get transposase bias by GC content  !!!!!!'
trans_bias = peak_folder + 'transposase_bias.bed'
temp02_file = peak_folder + 'temp02.bed'
temp03_file = peak_folder + 'temp03.bed'
with open(annotate_peak) as annotate_file, open(temp02_file, 'w') as temp02:
    for i, line in enumerate(annotate_file):
        words = line.split('\t')
        leave = words[0:3]
        print >> temp02, '\t'.join(leave)
# use GC contents to estimate the transposase bias
os.popen('bedtools nuc -fi ' + genome_fasta + ' -bed ' + temp02_file + ' > ' + temp03_file)
with open(temp03_file) as temp03, open(trans_bias, 'w') as bias:
    for i, line in enumerate(temp03):
        if i>0:
            words = line.split('\t')
            leave = words[0:3] + [words[4]]
            print >> bias, '\t'.join(leave)
print '!!!!!!  get bias done  !!!!!!'
#
#

