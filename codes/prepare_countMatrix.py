#!/usr/bin/python
import warnings
warnings.filterwarnings("ignore")
import os
from optparse import OptionParser
import subroutines
#
opts = OptionParser()
usage = "Align reads to build matrix\nusage: %prog -s source_folder --fa chr.fa --bg bg --meme motif.meme --np 4"
opts = OptionParser(usage=usage, version="%prog 2.1")
opts.add_option("-s", help="Source folder.")
opts.add_option("--fa", default='../reference/hg19_chr.fa', help="Genome fasta file, default=../reference/hg19_chr.fa")
opts.add_option("--bg", default='../reference/tier1_markov1.norc.txt',
                help="Background file, default=../reference/tier1_markov1.norc.txt")
opts.add_option("--pvalue", default=0.00005, help="P-value threshold for FIMO, default=0.00005")
opts.add_option("--meme", default='../reference/JASPAR2018_CORE_vertebrates_redundant_pfms_meme.txt',
                help="Motif file, default=../reference/JASPAR2018_CORE_vertebrates_redundant_pfms_meme.txt")
opts.add_option("--np", default=1, help="Number of CPU cores used for motif searching, default=1")
options, arguments = opts.parse_args()
#
matrix_folder = options.s + '/matrix/'
peaks_folder = options.s + '/peak/'
work_folder = options.s + '/work/'
motif_folder = matrix_folder + '/motif'
peaks_file = peaks_folder + '/annotate_peak.bed'
cell_info = options.s + '/data/cell_info.csv'
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
matrix_df = subroutines.counts_per_peak(bam_file, peaks_file, reads_matrix)
matrix_df.to_csv(reads_matrix, sep=',')
#
#
subroutines.QC_table(cell_info, work_folder, matrix_folder)
#
#
