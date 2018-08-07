#!/usr/bin/python
#
import os
from optparse import OptionParser
import subroutines
#
opts = OptionParser()
usage = "Peak alignment\nusage: %prog -w work -p peaks -m matrix --cellinfo cell_info.csv " \
    + "--fa chr.fa --bg bg --meme motif.meme --np 4"
opts = OptionParser(usage=usage, version="%prog 2.1")
opts.add_option("-w", help="work folder")
opts.add_option("-p", help="peak folder")
opts.add_option("-m", help="matrix folder")
opts.add_option("--cellinfo", help="cell_info.csv file with 2 columns (name, notes). "
                          +"It's usually automatically built by trimming step, and placed in data folder.")
opts.add_option("--fa", default='../reference/hg19_chr.fa', help="Genome fasta file, default=../reference/hg19_chr.fa")
opts.add_option("--bg", default='../reference/tier1_markov1.norc.txt',
                help="Background file, default=../reference/tier1_markov1.norc.txt")
opts.add_option("--pvalue", default=0.00005, help="P-value threshold for FIMO, default=0.00005")
opts.add_option("--meme", default='../reference/JASPAR2018_Homo_sapiens.meme',
                help="Motif file, default=../reference/JASPAR2018_Homo_sapiens.meme")
opts.add_option("--np", default=1, help="Number of CPU cores used for motif searching, default=1")
options, arguments = opts.parse_args()
#
matrix_folder = options.m
peaks_folder = options.p
motif_folder = matrix_folder + '/motif'
peaks_file = peaks_folder + '/annotate_peak.bed'
#
if not os.path.exists(matrix_folder): os.popen('mkdir ' + matrix_folder)
if os.path.exists(motif_folder): os.popen('rm -rf ' + motif_folder)
os.popen('mkdir ' + motif_folder)
#
#### run FIMO for motif-site searching
motifFasta = matrix_folder + '/motif.fasta'
os.popen('bedtools getfasta -fi ' + options.fa + ' -bed ' + peaks_file + ' -fo ' + motifFasta)
subroutines.batch_fimo(options.bg, options.pvalue, options.meme, motifFasta, motif_folder, int(options.np))
#
#### motif annotation
TFmatrix_file = matrix_folder + '/motif_TF.csv'
subroutines.score_peaks(peaks_file, motif_folder, TFmatrix_file)
#
#### count reads for peaks
bam_file = peaks_folder + "/mergeAll.bam"
reads_matrix = matrix_folder + "/reads.csv"
subroutines.counts_per_peak(bam_file, peaks_file, reads_matrix)
#
#
subroutines.QC_table(options.cellinfo, options.w, matrix_folder)
#
#
