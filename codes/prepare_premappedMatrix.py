#!/usr/bin/python
import warnings
warnings.filterwarnings("ignore")
import os
import numpy
import sys
from optparse import OptionParser
import subroutines
#
#
opts = OptionParser()
usage = "Organize premapped data\nusage: %prog -s source_folder --ref ref --fa chr.fa --bg bg.txt --meme motifs.meme --np 4"
opts = OptionParser(usage=usage, version="%prog 2.1")
opts.add_option("-s", help="Source folder")
opts.add_option("--ref", default='hg19', help="Genome reference, default=hg19")
opts.add_option("--fa", default='../reference/hg19_chr.fa', help="Genome fasta file, default=../reference/hg19_chr.fa")
opts.add_option("--bg", default='../reference/tier1_markov1.norc.txt',
                help="Background file, default=../reference/tier1_markov1.norc.txt")
opts.add_option("--pvalue", default=0.00005, help="P-value threshold for FIMO, default=0.00005")
opts.add_option("--meme", default='../reference/JASPAR2018_CORE_vertebrates_redundant_pfms_meme.txt',
                help="Motif file, default=../reference/JASPAR2018_CORE_vertebrates_redundant_pfms_meme.txt")
opts.add_option("--np", default=1, help="Number of CPU cores to use, default=1")
options, arguments = opts.parse_args()
#
#
top_peaks = options.s + '/peak/top_peaks.bed'
print '!!!!!!  annotate peaks  !!!!!!'
annotate_output = options.s + '/peak/annotate_output.bed'
annotate_peak = options.s + '/peak/annotate_peak.bed'
#
os.popen("annotatePeaks.pl " + top_peaks + " " + options.ref + " -size given > " + annotate_output)
temp01_file = options.s + '/peak/temp01.bed'
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
print '!!!!!!  get transposase bias by GC content  !!!!!!'
trans_bias = options.s + '/peak/transposase_bias_filtered.bed'
temp02_file = options.s + '/peak/temp02.bed'
temp03_file = options.s + '/peak/temp03.bed'
with open(annotate_peak) as annotate_file, open(temp02_file, 'w') as temp02:
    for i, line in enumerate(annotate_file):
        words = line.split('\t')
        leave = words[0:3]
        print >> temp02, '\t'.join(leave)
# use GC contents to estimate the transposase bias
os.popen('bedtools nuc -fi ' + options.fa + ' -bed ' + temp02_file + ' > ' + temp03_file)
with open(temp03_file) as temp03, open(trans_bias, 'w') as bias:
    for i, line in enumerate(temp03):
        if i>0:
            words = line.split('\t')
            leave = words[0:3] + [words[4]]
            print >> bias, '\t'.join(leave)
print '!!!!!!  get bias done  !!!!!!'
#
#
#
motif_folder = options.s + '/matrix/motif'
peaks_file = options.s + '/peak/annotate_peak.bed'
if os.path.exists(motif_folder): os.popen('rm -rf ' + motif_folder)
os.popen('mkdir ' + motif_folder)
# run FIMO for motif-site searching
motifFasta = options.s + '/matrix/motif.fasta'
os.popen('bedtools getfasta -fi ' + options.fa + ' -bed ' + peaks_file + ' -fo ' + motifFasta)
subroutines.batch_fimo(options.bg, options.pvalue, options.meme, motifFasta, motif_folder, int(options.np))
# motif annotation
TFmatrix_file = options.s + '/matrix/motif_filtered.csv'
subroutines.score_peaks(peaks_file, motif_folder, TFmatrix_file)
#
#
#
#
