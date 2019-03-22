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
from sklearn.manifold import TSNE
from sklearn import cluster
from sklearn.neighbors import kneighbors_graph
import time
#
#
global TFmotif, samples, reads, expected, TFnames, cell_names, GC_bias
#
opts = OptionParser()
usage = "Cluster by ChromVAR\nusage: %prog -s project --np 8"
opts = OptionParser(usage=usage, version="%prog 1.0.5")
opts.add_option("-s", help="The project folder.")
opts.add_option("--format", default='csv', help="=csv or mtx; read fragment count matrix from csv or mtx file, default=csv.")
opts.add_option("--ns", default=50, help="Number of permuted samplings, default=50")
opts.add_option("--np", default=1, help="CPU cores used for samplings, default=1")
opts.add_option("--nc", default=0, help="Number of cell clusters, default=0, i.e. predicted by Louvain algorithm")
opts.add_option("--hc", default='no', help="Run hierarchical clustering or not, default=no.")
opts.add_option("--rs", default=1, help="Random state of tSNE analysis, default=1")
options, arguments = opts.parse_args()
#
#
if not os.path.exists(options.s+'/result'): os.popen('mkdir ' + options.s+'/result')
if not os.path.exists(options.s+'/figure'): os.popen('mkdir ' + options.s+'/figure')
#
#
def initiation(options):
    if options.format=='mtx':
        reads = scipy.io.mmread(options.s+'/matrix/filtered_reads.mtx')
        reads = scipy.sparse.csr_matrix(reads)
        cells = pandas.read_csv(options.s+'/matrix/filtered_cells.csv', sep='\t', index_col=0,
                   engine='c', na_filter=False, low_memory=False)
        cells_names = cells.index.values
    else:
        reads_df = pandas.read_csv(options.s+'/matrix/filtered_reads.csv', sep=',', index_col=0,
                   engine='c', na_filter=False, low_memory=False)
        reads, cell_names = reads_df.values + 0.0001, reads_df.index.values
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
    reads_df = pandas.read_csv(options.s+'/result/deviation_chromVAR.csv', sep=',', index_col=0,
                   engine='c', na_filter=False, low_memory=False)
    matrix = reads_df.T
    connect = kneighbors_graph(matrix, n_neighbors=20, include_self=False)
    connectivity = 0.5*(connect + connect.T)
    if int(options.nc)==0:
        n_clust, clusters = subroutines.predict_cluster(matrix, connectivity.todense())
        print "predicted number of cell-clusters: ", n_clust
        clusters.to_csv(options.s+'/result/louvain_cluster_by_chromVAR.csv', sep='\t')
        tsne_result = TSNE(n_components=2, random_state=int(options.rs)).fit_transform(matrix.values)
        subroutines.plot_cluster(options, clusters, n_clust, tsne_result, 'louvain_cluster_by_chromVAR.pdf')
    else:
        n_clust = int(options.nc)
        clusters = subroutines.knn_cluster(options, matrix, n_clust, connectivity, "KNN_cluster_by_chromVAR.csv")
        tsne_result = TSNE(n_components=2, random_state=int(options.rs)).fit_transform(matrix.values)
        subroutines.plot_cluster(options, clusters, n_clust, tsne_result, 'KNN_cluster_by_chromVAR.pdf')
#
    subroutines.plot_tSNE(options, matrix, tsne_result, "TSNE_by_chromVAR.pdf")
    tsne_df = pandas.DataFrame(tsne_result, index=matrix.index, columns=['TSNE1', 'TSNE2'])
    tsne_df.to_csv(options.s+'/result/TSNE_by_chromVAR.csv', sep='\t')
#
    if options.hc=='yes':
        subroutines.hierarchy_cluster(options, matrix, n_clust, "cell_cell_correlation_by_chromVAR.png",
                                  "Hierarchical_cluster_by_chromVAR.csv")
    return
#
#
#t1=time.time()
TFmotif, reads, TFnames, cell_names, GC_bias = initiation(options)
samples = permuted_sampling(options)
expected, raw_dev = raw_deviation(options)
corrected_deviation(options)
cell_cluster(options)
#print time.time()-t1
#
#
#
#
#
