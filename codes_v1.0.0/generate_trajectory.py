#!/usr/bin/python
import warnings
warnings.filterwarnings("ignore")
#
import numpy
import pandas
from optparse import OptionParser
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn
import os
from sklearn.decomposition import PCA
#
opts = OptionParser()
usage = "Build monocle trajectory\nusage: %prog -s project --npc 5 --cfile cluster.csv"
opts = OptionParser(usage=usage, version="%prog 1.0")
opts.add_option("-s", help="The project folder.")
opts.add_option("--npc", default=5, help="Number of principle components used for pseudo-time trajectory, defaul=5")
opts.add_option("--cfile", help='Cell-types file, e.g. cell_info.csv or cluster.csv')
opts.add_option("--dim", default=3, help="Plot 2D or 3D trajectory, default=3")
opts.add_option("--angle", default='30,30', help='Angles to rotate the 3D trajectory, default=30,30')
options, arguments = opts.parse_args()
#
#
def trajectory(options):
    reads_csv = options.s+'/matrix/Accesson_reads.csv'
    celltype_csv = options.cfile
    reads = pandas.read_csv(reads_csv, sep=',', index_col=0, engine='c', na_filter=False, low_memory=False)
#
    matrix = reads.values
    normal = numpy.array([x/x.sum()*1000000 for x in matrix])
    reads = pandas.DataFrame(normal, index=reads.index.values, columns=reads.columns.values)
#
    npc = int(options.npc)
    pca_result = PCA(n_components=npc, svd_solver='full').fit_transform(reads)
    reads = pandas.DataFrame(pca_result[:, :npc], columns=['pc'+str(x) for x in xrange(npc)], index=reads.index.values)
    print reads.shape
    if 'cluster' in options.cfile:
        celltype_df = pandas.read_csv(celltype_csv, sep='\t', index_col=0)
    else:
        cellinfo = pandas.read_csv(options.cfile, sep='\t', index_col=0)
        celltype_df = pandas.DataFrame(cellinfo['notes'].values, index=cellinfo.index.values, columns=['cluster'])
        celltype_df = celltype_df.ix[reads.index.values]
    input_csv = options.s+'/matrix/monocle_reads.csv'
    cells_csv = options.s+'/matrix/monocle_cells.tsv'
    peaks_csv = options.s+'/matrix/monocle_peaks.tsv'
    trajectory_csv = options.s+'/result/monocle_trajectory.csv'
    reduced_csv = options.s+'/result/monocle_reduced_dimension.csv'
    peaks_df = pandas.DataFrame(reads.columns.values, index=reads.columns.values, columns=['gene_short_name'])
    peaks_df.to_csv(peaks_csv, sep='\t') 
    celltype_df.to_csv(cells_csv, sep='\t')
    reads = reads.ix[celltype_df.index.values, peaks_df.index.values]
    print reads.shape, peaks_df.shape, celltype_df.shape
    reads.T.to_csv(input_csv, sep=',')
#
    os.system('Rscript run_monocle.R '+input_csv+' '+cells_csv+' '+peaks_csv+' '+reduced_csv+' '+trajectory_csv)
    return
#
#
def plot_traj(reduced_df, cells_df, out_fig, options):
    size = 20
    wordsize = 20
    colors = numpy.array(['pink', 'red', '#377eb8', 'green', 'skyblue', 'lightgreen', 'gold',
                      '#ff7f00', '#000066', '#ff3399', '#a65628', '#984ea3', '#999999',
                      '#e41a1c', '#dede00', 'b', 'g', 'r', 'c', 'm', 'y', 'k'])
    fig = plt.figure(1, figsize=(10,10))
    if int(options.dim)==3:
        ax = fig.add_subplot(111, projection='3d')
        cell_types = list(set(list(cells_df['cluster'].values)))
        cell_types.sort()
        for itype,ctype in enumerate(cell_types):
            cluster = cells_df.loc[cells_df['cluster']==ctype]
            cluster = cluster.index.values
            components = reduced_df[cluster].values
            ax.scatter(components[0], components[1], components[2], c=colors[itype], edgecolors='none', 
                       s=size, label=cell_types[itype])
        ax.set_zlabel('Component 3')
        beta1, beta2 = int(options.angle.split(',')[0]), int(options.angle.split(',')[1])
        ax.view_init(beta1, beta2)
        ax.set_zticklabels([])
        ax.set_zlim(reduced_df.iloc[2].min()-0.02, reduced_df.iloc[2].max()+0.02)
    elif int(options.dim)==2:
        ax = fig.add_subplot(111)
        cell_types = list(set(list(cells_df['cluster'].values)))
        cell_types.sort()
        for itype,ctype in enumerate(cell_types):
            cluster = cells_df.loc[cells_df['cluster']==ctype]
            cluster = cluster.index.values
            components = reduced_df[cluster].values
            ax.scatter(components[0], components[1], c=colors[itype], edgecolors='none', 
                       s=size, label=cell_types[itype])
    ax.legend(fontsize=wordsize, bbox_to_anchor=(1.0, 1.0))
    ax.set_xlabel('Component 1')
    ax.set_ylabel('Component 2')
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_xlim(reduced_df.iloc[0].min()-0.02, reduced_df.iloc[0].max()+0.02)
    ax.set_ylim(reduced_df.iloc[1].min()-0.02, reduced_df.iloc[1].max()+0.02)
    fig.savefig(out_fig, bbox_inches='tight')
    return
#
trajectory(options)
reduced_df = pandas.read_csv(options.s+'/result/monocle_reduced_dimension.csv', sep=',', index_col=0)
cells_df = pandas.read_csv(options.s+'/matrix/monocle_cells.tsv', sep='\t', index_col=0)
plot_traj(reduced_df, cells_df, options.s+'/figure/pseudotime_trajectory.pdf', options)
#
