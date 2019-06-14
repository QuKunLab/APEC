#!/usr/bin/python
import warnings
warnings.filterwarnings("ignore")
import os,numpy,scipy.sparse,scipy.io,pandas,sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from optparse import OptionParser
#
#
opts = OptionParser()
usage = "Build QC table\nusage: %prog -s project --pfrag 0.2 --lib 2000"
opts = OptionParser(usage=usage, version="%prog 1.0.6")
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
reads = scipy.sparse.csr_matrix(scipy.io.mmread(options.s+'/matrix/reads.mtx')).T
cell_info = pandas.read_csv(options.s+'/matrix/cell_info.merged.csv', sep='\t', index_col=0)
libSize = cell_info['final_reads'].values
readsInPeaks = reads.sum(axis=1).A[:,0].astype(float)/cell_info['final_reads'].values
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
filtered_quality_df = cell_quality_df.loc[cell_quality_df['lib_size']>lib_size_thresh]
filtered_quality_df = filtered_quality_df.loc[filtered_quality_df['frag_in_peak']>frag_thresh]
filtered_cell_names = filtered_quality_df.index.values
filtered_cell_index = [list(cell_quality_df.index.values).index(x) for x in filtered_cell_names]
filtered_cells_df = cell_info.loc[filtered_cell_names, 'notes']
filtered_cells_df = pandas.DataFrame(filtered_cells_df.values, index=filtered_cell_names, columns=['notes'])
filtered_cells_df.to_csv(options.s+'/matrix/filtered_cells.csv', sep='\t')
#
nonzero_per_peak = numpy.array([len(numpy.where(x.A>0)[0]) for x in reads.T])
filtered_peak_index = numpy.where(nonzero_per_peak>=3)[0]
with open(options.s+'/peak/top_peaks.bed') as in_file, \
     open(options.s+'/peak/top_filtered_peaks.bed', 'w') as out_file:
    for iline,line in enumerate(in_file):
        if iline in filtered_peak_index:
            out_file.write(line)
#
reads = reads[filtered_cell_index, :]
reads = reads[:, filtered_peak_index]
scipy.io.mmwrite(options.s+'/matrix/filtered_reads.mtx', scipy.sparse.coo_matrix(reads.T))
#
#
