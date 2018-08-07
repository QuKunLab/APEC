#!/usr/bin/python
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
usage = "QC table\nusage: %prog -m matrix -p peak -f figure --frag 0.2 --lib 2000"
opts = OptionParser(usage=usage, version="%prog 2.1")
opts.add_option("-m", help="matrix folder")
opts.add_option("-p", help="peak folder")
opts.add_option("-f", help="figure folder")
opts.add_option("--frag", default=0.2, help="Threshold for fragment ratio in peaks, "
                +"default=0.2, decrease it for more filtered-samples, increase it for better quality")
opts.add_option("--lib", default=2000, help="Threshold for final count, default=2000, "
                +"decrease it for more filtered-samples, increase it for better quality")
options, arguments = opts.parse_args()
#
if not os.path.exists(options.f): os.popen('mkdir ' + options.f)
frag_thresh = float(options.frag)
lib_size_thresh = int(options.lib)
#
reads_df = pandas.DataFrame.from_csv(options.m+'/reads.csv', sep=',')
reads = reads_df.values
cells = pandas.DataFrame.from_csv(options.m+'/cell_info.merged.csv', sep='\t')
cells_name = reads_df.index.values
cell_info = cells.loc[cells_name]
libSize = cell_info['final_reads'].values
double = 1
readsInPeaks = reads.sum(axis=1).astype(float)*double/cell_info['final_reads'].values
cell_quality = numpy.vstack((libSize, readsInPeaks))
cell_quality_df = pandas.DataFrame(cell_quality.T, index=cell_info.index.values, columns=['lib_size', 'frag_in_peak'])
cell_quality_df.to_csv(options.m + '/cell_quality.csv', sep='\t')
#
plt.scatter(libSize+1, readsInPeaks, s=10)
plt.vlines(lib_size_thresh, -0.05, 1.05, linestyles='dashed')
plt.hlines(frag_thresh, 1e2, 1e7, linestyles='dashed')
plt.xscale('log')
plt.xlabel('final mapped reads')
plt.ylabel('fragments in peaks(%)')
plt.xlim(1e2, 1e7)
plt.ylim(-0.05, 1.05)
plt.savefig(options.f + '/cell_quality.pdf')
#
with open(options.m+'/filtered_cells.txt', 'w') as output:
    for ic,cell in enumerate(cell_info.index.values):
        cell_info = cell_quality_df.loc[cell].values
        if (cell_info[0]>lib_size_thresh) & (cell_info[1]>frag_thresh):
            print >> output, ic, cell 
#
drop_peaks = [peak for peak in reads_df.columns.values if len(numpy.where(reads_df[peak].values>0)[0])<3]
reads_df = reads_df.drop(drop_peaks, axis=1)
motif_df = pandas.DataFrame.from_csv(options.m+'/motif_TF.csv', sep=',')
motif_df = motif_df.loc[reads_df.columns.values]
peak_num = [int(x[4:]) for x in reads_df.columns.values]
bias = open(options.p+'/transposase_bias.bed').readlines()
with open(options.p+'/transposase_bias_filtered.bed', 'w') as output:
    for ipeak in peak_num:
        print >> output, bias[ipeak][:-1]
cell_names = [x.split()[1] for x in open(options.m+'/filtered_cells.txt').readlines()]
reads_df = reads_df.loc[cell_names]
reads_df.to_csv(options.m+'/filtered_reads.csv', sep=',')
motif_df.to_csv(options.m+'/motif_filtered.csv', sep=',')
print reads_df.shape, motif_df.shape
#
#
