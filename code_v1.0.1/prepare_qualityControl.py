#!/usr/bin/python
import warnings
warnings.filterwarnings("ignore")
import os
import numpy
import scipy.io
import pandas
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from optparse import OptionParser
#
#
opts = OptionParser()
usage = "Build QC table\nusage: %prog -s project --frag 0.2 --lib 2000"
opts = OptionParser(usage=usage, version="%prog 1.0")
opts.add_option("-s", help="The project folder.")
opts.add_option("--pfrag", default=0.2, help="Threshold for percentage of fragments in peaks, "
                +"default=0.2, decrease it for more filtered-samples, increase it for better quality")
opts.add_option("--lib", default=2000, help="Threshold for fragment number, default=2000, "
                +"decrease it for more filtered-samples, increase it for better quality")
options, arguments = opts.parse_args()
#
if not os.path.exists(options.s+'/figure'): os.popen('mkdir ' + options.s+'/figure')
frag_thresh = float(options.pfrag)
lib_size_thresh = int(options.lib)
#
reads_df = pandas.read_csv(options.s+'/matrix/reads.csv', sep=',', index_col=0,
                   engine='c', na_filter=False, low_memory=False)
reads = reads_df.values
cell_info = pandas.read_csv(options.s+'/matrix/cell_info.merged.csv', sep='\t', index_col=0)
cells_name = reads_df.index.values
cell_info = cell_info.loc[cells_name]
libSize = cell_info['final_reads'].values
double = 1
readsInPeaks = reads.sum(axis=1).astype(float)*double/cell_info['final_reads'].values
cell_quality = numpy.vstack((libSize, readsInPeaks))
cell_quality_df = pandas.DataFrame(cell_quality.T, index=cell_info.index.values, columns=['lib_size', 'frag_in_peak'])
cell_quality_df.to_csv(options.s+'/matrix/cell_quality.csv', sep='\t')
#
plt.scatter(libSize+1, readsInPeaks, s=10)
plt.vlines(lib_size_thresh, -0.05, 1.05, linestyles='dashed')
plt.hlines(frag_thresh, 1e2, 1e7, linestyles='dashed')
plt.xscale('log')
plt.xlabel('final mapped reads')
plt.ylabel('fragments in peaks(%)')
plt.xlim(1e2, 1e7)
plt.ylim(-0.05, 1.05)
plt.savefig(options.s + '/figure/cell_quality.pdf')
#
with open(options.s+'/matrix/filtered_cells.txt', 'w') as output:
    for ic,cell in enumerate(cell_info.index.values):
        cell_info = cell_quality_df.loc[cell].values
        if (cell_info[0]>lib_size_thresh) & (cell_info[1]>frag_thresh):
            print >> output, ic, cell 
#
drop_peaks = [peak for peak in reads_df.columns.values if len(numpy.where(reads_df[peak].values>0)[0])<3]
reads_df = reads_df.drop(drop_peaks, axis=1)
motif_df = pandas.read_csv(options.s+'/matrix/motif_TF.csv', sep=',', index_col=0)
motif_df = motif_df.loc[reads_df.columns.values]
peak_num = [int(x[4:]) for x in reads_df.columns.values]
bias = open(options.s+'/peak/transposase_bias.bed').readlines()
with open(options.s+'/peak/transposase_bias_filtered.bed', 'w') as output:
    for ipeak in peak_num:
        print >> output, bias[ipeak][:-1]
cell_names = [x.split()[1] for x in open(options.s+'/matrix/filtered_cells.txt').readlines()]
reads_df = reads_df.loc[cell_names]
reads_df.to_csv(options.s+'/matrix/filtered_reads.csv', sep=',')
motif_df.to_csv(options.s+'/matrix/motif_filtered.csv', sep=',')
print reads_df.shape, motif_df.shape
#
#
