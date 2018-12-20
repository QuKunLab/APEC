#!/usr/bin/python
import warnings
warnings.filterwarnings("ignore")
#
import numpy
import subroutines
import scipy.io
import Bias_corrected_deviation
import pandas
import sys
import os
from optparse import OptionParser
from multiprocessing import Pool
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn
#
#
global TFmotif, samples, reads, expected, TFnames, cell_names, GC_bias
#
opts = OptionParser()
usage = "Cluster by ChromVAR\nusage: %prog -s project --np 4 --nc 0"
opts = OptionParser(usage=usage, version="%prog 1.0")
opts.add_option("-s", help="The project folder.")
opts.add_option("--ns", default=50, help="Number of permuted samplings, default=50")
opts.add_option("--np", default=1, help="CPU cores used for samplings, default=1")
opts.add_option("--nc", default=0, 
                help="Number of groups for cell to be clustered, default=0, i.e. predicted by Louvain algorithm")
opts.add_option("--npc", default=40, help="Number of PCs used for clustering, default=40")
opts.add_option("--space", default='pca', help="transform space used for clustering, can be pca or tsne, default=pca")
opts.add_option("--hc", default='yes', help="Run hierarchical clustering or not, would be very slow for more than 5000 cells."
                + " if hc=no, only KNN clustering will be applied. default=yes.")
opts.add_option("--tsneLR", default=100, help="Learning rate of tSNE analysis, default=100")
opts.add_option("--tsneRS", default=1, help="Random state of tSNE analysis, default=1")
opts.add_option("--tsneIter", default=1000, help="Muximum number of Iterations for tSNE analysis, default=1000")
options, arguments = opts.parse_args()
#
#
if not os.path.exists(options.s+'/result'): os.popen('mkdir ' + options.s+'/result')
if not os.path.exists(options.s+'/figure'): os.popen('mkdir ' + options.s+'/figure')
#
#
def initiation(options):
    reads_df = pandas.read_csv(options.s+'/matrix/filtered_reads.csv', sep=',', index_col=0,
                   engine='c', na_filter=False, low_memory=False)
    reads_df = reads_df + 0.0001
    reads, cell_names = reads_df.values+0, reads_df.index.values
    TFmotif_df = pandas.read_csv(options.s+'/matrix/motif_filtered.csv', sep=',', index_col=0)
    TFmotif_origin = TFmotif_df.values.T
    TFnames = TFmotif_df.columns.values
    print 'read-counts matrix:', reads.shape
    print 'TFmotif matrix:', TFmotif_origin.shape
    GC_bias = numpy.array([float(x.split()[3]) for x in open(options.s+'/peak/transposase_bias_filtered.bed').readlines()])
    TFmotif = numpy.asarray([x for x in TFmotif_origin if x.sum() > 0])
    TFmotif[numpy.where(TFmotif > 0)] = 1
    TFnames = [x for i,x in enumerate(TFnames) if TFmotif_origin[i, :].sum() > 0]
    return TFmotif, reads, TFnames, cell_names, GC_bias
#
#
def permuted_sampling(options):
    peak_reads = numpy.log10(reads.sum(axis=0)+1.0)
    ngrid, std = 50, 1
    print GC_bias.shape, peak_reads.shape
    samples = Bias_corrected_deviation.batch_sampling(GC_bias, peak_reads, ngrid, std, int(options.ns), int(options.np))
    print 'permuted sampling done!'
    return samples
#
#
def raw_deviation(options):
    expected = Bias_corrected_deviation.expected_matrix(reads, TFmotif, GC_bias)
    raw_dev = Bias_corrected_deviation.raw_deviation(TFmotif, reads, expected)
    numpy.savetxt(options.s+'/result/raw_deviation.txt', raw_dev)
    print 'raw deviation done!'
    return expected, raw_dev
#
#
def background_deviation(iIter):
    numpy.random.seed(12345+iIter)
    (nCell, nPeak) = reads.shape
    (nTF, nPeak) = TFmotif.shape
    background_dev = numpy.zeros((nTF, nCell))
    B_matrix = numpy.zeros((nPeak, nPeak))
    for iPeak in range(0, nPeak):
        B_matrix[iPeak, int(samples[iIter, iPeak])] = 1
    background_dev = Bias_corrected_deviation.deviation(TFmotif, B_matrix, reads, expected)
    print 'background deviation for sample '+str(iIter+1)+' done!'
    return background_dev
#
#
def corrected_deviation(options):
    kIterations = numpy.arange(0, int(options.ns), 1, dtype=int)
    pool = Pool(int(options.np))
    bg_dev = pool.map(background_deviation, kIterations)
    pool.close()
    pool.join()
    bg_dev = numpy.array(bg_dev)
    print 'background deviations done!'
    bg_dev_mean = bg_dev.mean(axis=0)
    bg_dev_std = bg_dev.std(axis=0)
    raw_dev = numpy.loadtxt(options.s+'/result/raw_deviation.txt')
    corrected_dev = (raw_dev - bg_dev_mean) / bg_dev_std
    dev_df = pandas.DataFrame(corrected_dev, index=TFnames, columns=cell_names)
    dev_df.to_csv(options.s+'/result/deviation_chromVAR.csv', sep=',')
    return
#
#
def cell_cluster(options):
    n_clust = int(options.nc)
    dev_df = pandas.read_csv(options.s+'/result/deviation_chromVAR.csv', sep=',', index_col=0,
                   engine='c', na_filter=False, low_memory=False)
    corr = dev_df.corr()
    if n_clust==0:
        n_clust = subroutines.predict_cluster(corr)
        print "predicted number of clusters: ", n_clust
#
    pca_result, tsne_result = subroutines.PCA_tSNE(options, corr, "PCA_by_chromVAR.pdf", "TSNE_by_chromVAR.pdf")
    tsne_df = pandas.DataFrame(tsne_result, index=corr.index.values, columns=['TSNE1', 'TSNE2'])
    tsne_df.to_csv(options.s+'/result/TSNE_by_chromVAR.csv', sep='\t')
    pca_df = pandas.DataFrame(pca_result, index=corr.index.values, columns=['PC'+str(i) for i in range(0, int(options.npc))])
    if options.space=='tsne':
        matrix_df = tsne_df
    else:
        matrix_df = pca_df
    subroutines.plot_knn_cluster(options, matrix_df, n_clust, tsne_result, "KNN_cluster_by_chromVAR.pdf", 
                            "KNN_cluster_by_chromVAR.csv")
#
    if options.hc=='yes':
        subroutines.hierarchy_cluster(options, corr, n_clust, "cell_cell_correlation_by_chromVAR.png",
                            "Hierarchical_cluster_by_chromVAR.csv")
        KNN_df = pandas.read_csv(options.s+"/result/KNN_cluster_by_chromVAR.csv", sep='\t', index_col=0)
        subroutines.heatmap_compare(options, corr, KNN_df, "HC_KNN_compare_by_chromVAR.png")
    return
#
#
TFmotif, reads, TFnames, cell_names, GC_bias = initiation(options)
samples = permuted_sampling(options)
expected, raw_dev = raw_deviation(options)
corrected_deviation(options)
cell_cluster(options)
#
#
#
#
#
