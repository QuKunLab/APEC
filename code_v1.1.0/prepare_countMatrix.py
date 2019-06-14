#!/usr/bin/python
import warnings
warnings.filterwarnings("ignore")
import os,pandas
from optparse import OptionParser
import subroutines
import scipy.io,scipy.sparse
#
opts = OptionParser()
usage = "Align reads to build matrix\nusage: %prog -s project --fa chr.fa --bg bg --meme motif.meme --np 4"
opts = OptionParser(usage=usage, version="%prog 1.0")
opts.add_option("-s", help="The project folder.")
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
peaks_file = peaks_folder + '/top_peaks.bed'
cell_info = options.s + '/data/cell_info.csv'
motifFasta = matrix_folder + '/motif.fasta'
TFmatrix_file = matrix_folder + '/motif_TF.csv'
bam_file = peaks_folder + "/mergeAll.bam"
reads_matrix = matrix_folder + "/reads.mtx"
#
if not os.path.exists(matrix_folder): os.popen('mkdir ' + matrix_folder)
#if os.path.exists(motif_folder): os.popen('rm -rf ' + motif_folder)
#os.popen('mkdir ' + motif_folder)
#
#### run FIMO for motif-site searching
#os.popen('bedtools getfasta -fi ' + options.fa + ' -bed ' + peaks_file + ' -fo ' + motifFasta)
#subroutines.batch_fimo(options.bg, options.pvalue, options.meme, motifFasta, motif_folder, int(options.np))
#
#### motif annotation
#subroutines.score_peaks(peaks_file, motif_folder, TFmatrix_file)
#
#### count reads for peaks
cell_info_df = pandas.read_csv(cell_info, sep='\t', index_col=0)
matrix = subroutines.counts_per_peak(bam_file, peaks_file, reads_matrix, cell_info_df)
matrix = scipy.sparse.coo_matrix(matrix.T)
scipy.io.mmwrite(reads_matrix, matrix)
#
#
subroutines.QC_table(cell_info, work_folder, matrix_folder)
#
#
