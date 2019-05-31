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
import sys
import scipy.stats
#
#
opts = OptionParser()
usage = "Render marker on tSNE or trajectory plot\nusage: %prog -s project --cfile TSNE_by_Accesson.csv --type motif --name RELA"
opts = OptionParser(usage=usage, version="%prog 1.0")
opts.add_option("-s", help="The project folder.")
opts.add_option("--cfile", help="TSNE_by_Accesson.csv or monocle_reduced_dimension.csv file to use for rendering")
opts.add_option("--type", default='motif', help="type of marker to plot, can be motif, gene, or accesson")
opts.add_option("--name", help="name of marker to plot")
opts.add_option("--angle", default='30,30', help='Angles to rotate the 3D trajectory, default=30,30')
opts.add_option("--sharp", default='0', help='Cutoff range for deviation or expression, default=0, i.e. no sharpening')
options, arguments = opts.parse_args()
#
#
def draw_marker(options):
    if 'monocle_reduced' in options.cfile:
        tsne_df = pandas.read_csv(options.cfile, sep=',', index_col=0).T
    else:
        tsne_df = pandas.read_csv(options.cfile, sep='\t', index_col=0)
    tsne, cells = tsne_df.values, tsne_df.index.values
    if options.type=='motif':
        reads_df = pandas.read_csv(options.s+'/result/deviation_chromVAR.csv', sep=',', index_col=0,
                   engine='c', na_filter=False, low_memory=False).T
        motifs = []
        for tf in reads_df.columns.values:
            names = tf.split('-')[:-1]
            if options.name in names: motifs.append(tf)
        print(motifs)
        if len(motifs)==0:
            print("No corresponding marker!")
            sys.exit()
        else:
            reads = reads_df.ix[cells, motifs].values.sum(axis=1)
    elif options.type=='gene':
        reads_df = pandas.read_csv(options.s+'/matrix/genes_scored_by_peaks.csv', sep=',', index_col=0,
                   engine='c', na_filter=False, low_memory=False).T
        gene = list(set([options.name]).intersection(set(reads_df.columns.values)))
        if len(gene)==0:
            print("No corresponding marker!")
            sys.exit()
        else:
            reads = reads_df.loc[cells, gene].values.sum(axis=1)
    elif options.type=='accesson':
        reads_df = pandas.read_csv(options.s+'/matrix/Accesson_reads.csv', sep=',', index_col=0,
                   engine='c', na_filter=False, low_memory=False)
        normal = numpy.array([x/x.sum()*10000 for x in reads_df.values])
        reads_df = pandas.DataFrame(normal, index=reads_df.index, columns=reads_df.columns)
        reads = reads_df.loc[cells, options.name].values
    order = numpy.argsort(reads)
    if str(options.sharp)!='0':
        up, down = int(options.sharp.split(',')[-1]), int(options.sharp.split(',')[0])
        reads = numpy.clip(reads, down, up)
#
    if 'monocle_reduced' in options.cfile:
        clist = ['blue', 'silver', 'red']
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list('mylist', clist, N=256)
        fig1 = plt.figure(1, figsize=(12,10))
        ax = fig1.add_subplot(111, projection='3d')
        im = ax.scatter(tsne[order,0], tsne[order,1], tsne[order,2], cmap=cmap, c=reads[order], edgecolors='none', s=10)
        beta1, beta2 = int(options.angle.split(',')[0]), int(options.angle.split(',')[1])
        ax.view_init(beta1, beta2)
        cbar = plt.colorbar(im, shrink=0.15, ticks=[reads.min(), reads.max()], aspect=8)
        cbar.ax.set_yticklabels([round(reads.min(),2), round(reads.max(),2)])
        width, height, rad = tsne[:,0].max()-tsne[:,0].min(), tsne[:,1].max()-tsne[:,1].min(), tsne[:,2].max()-tsne[:,2].min()
        ax.set_xlim((tsne[:,0].min()-0.01*width, tsne[:,0].max()+0.01*width))
        ax.set_ylim((tsne[:,1].min()-0.01*height, tsne[:,1].max()+0.01*height))
        ax.set_zlim((tsne[:,2].min()-0.01*rad, tsne[:,2].max()+0.01*rad))
        ax.set_zticks([])
    else:
        fig1 = plt.figure(1, figsize=(6,5))
        ax = fig1.add_subplot(111)
        im = ax.scatter(tsne[order,0], tsne[order,1],  cmap='Spectral_r', c=reads[order], edgecolors='none', s=20)
        cbar = plt.colorbar(im, shrink=0.15, ticks=[reads.min(), reads.max()], aspect=8)
        cbar.ax.set_yticklabels([round(reads.min(),2), round(reads.max(),2)])
#        width, height = tsne[:,0].max()-tsne[:,0].min(), tsne[:,1].max()-tsne[:,1].min()
#        ax.set_xlim((tsne[:,0].min()-0.01*width, tsne[:,0].max()+0.01*width))
#        ax.set_ylim((tsne[:,1].min()-0.01*height, tsne[:,1].max()+0.01*height))
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(options.name)
    fig1.savefig(options.s+'/figure/'+options.type+'_'+options.name+'_on_'+options.cfile.split('/')[-1].split('.')[0]+'.pdf',
                 bbox_inches='tight')
    return
#
#
draw_marker(options)
#
#
#
#
