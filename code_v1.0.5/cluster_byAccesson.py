#!/usr/bin/python
import warnings
warnings.filterwarnings("ignore")
#
import os
import numpy
import pandas
from optparse import OptionParser
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn
import scipy.spatial.distance
import scipy.io
from scipy import sparse
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn import cluster
from sklearn.neighbors import kneighbors_graph
from sklearn.cluster import MiniBatchKMeans
from sklearn.feature_extraction.text import TfidfTransformer
from sklearn.decomposition.truncated_svd import TruncatedSVD
from sklearn.preprocessing import MinMaxScaler
from scipy import stats
import subroutines
import random
import time
import copy
#
#
opts = OptionParser()
usage = "Clustering by accesson\nusage: %prog -s project"
opts = OptionParser(usage=usage, version="%prog 1.0.5")
opts.add_option("-s", help="The project folder.")
opts.add_option("--format", default='csv', help="=csv or mtx; read fragment count matrix from csv or mtx file, default=csv.")
opts.add_option("--nc", default=0, help="Number of cell clusters, default=0, i.e. predicted by Louvain algorithm")
opts.add_option("--ngroup", default=600, help="Number of accessons, default=600")
opts.add_option("--npc", default=40, help="Number of PCs used for peak grouping, default=40")
opts.add_option("--ncell", default=10000, help="Maximum number of cells used for peak grouping, default=10000")
opts.add_option("--norm", default='zscore', help="=zscore or probability; set normalization method, default=zscore.")
opts.add_option("--filter", default='yes', help="=yes or no; filter accessons or not, default=yes.")
opts.add_option("--rs", default=0, help="Random state of tSNE analysis, default=0")
opts.add_option("--wt", default=0.2, help="Weight coefficient of tSNE analysis, default=0.2")
opts.add_option("--hc", default='no', help="Run hierarchical clustering or not, default=no.")
options, arguments = opts.parse_args()
#
#
def build_accesson(options):
    ngroups, ncell_cut = int(options.ngroup), int(options.ncell)
    if options.format=='mtx':
        reads = scipy.io.mmread(options.s+'/matrix/filtered_reads.mtx')
        reads = scipy.sparse.csr_matrix(reads)
        cells = pandas.read_csv(options.s+'/matrix/filtered_cells.csv', sep='\t', index_col=0,
                   engine='c', na_filter=False, low_memory=False)
        cells = cells.index.values
        peaks = ['peak'+str(x) for x in range(0, reads.shape[1])]
        if len(cells)>ncell_cut:
            normal = copy.deepcopy(reads)
            normal.data = numpy.ones(len(normal.data))
            index_sampled = random.sample(numpy.arange(0, len(cells), 1, dtype=int), ncell_cut)
            normal2 = normal[index_sampled, :]
        else:
            reads = numpy.array(reads.todense())
            normal = numpy.array([x*10000.0/x.sum() for x in reads])
            normal = numpy.log2(normal+1)
            normal2 = normal
    else:
        reads_df = pandas.read_csv(options.s+'/matrix/filtered_reads.csv', sep=',', index_col=0,
                   engine='c', na_filter=False, low_memory=False)
        cells, peaks = reads_df.index.values, reads_df.columns.values
        reads = reads_df.values
        normal = numpy.array([x*10000.0/x.sum() for x in reads])
        normal = numpy.log2(normal+1)
        normal2 = normal
        normal_df = pandas.DataFrame(normal2, index=reads_df.index, columns=reads_df.columns)
        normal_df.to_csv(options.s+'/matrix/normal_reads.csv', sep=',')
    npc = min(int(options.npc), normal2.shape[0], normal2.shape[1])
    if len(cells)>ncell_cut:
        pca_result = TruncatedSVD(n_components=npc).fit_transform(normal2.T)
    else:
        pca_result = PCA(n_components=npc, svd_solver='full').fit_transform(normal2.T)
    connectivity = kneighbors_graph(pca_result, n_neighbors=10, include_self=False)
    connectivity = 0.5*(connectivity + connectivity.T)
    ward_linkage = cluster.AgglomerativeClustering(n_clusters=ngroups, linkage='ward', connectivity=connectivity)
    ward_linkage.fit(pca_result)
    y_predict = ward_linkage.labels_.astype(numpy.int)
    peak_labels_df = pandas.DataFrame(y_predict, index=peaks, columns=['group'])
    peak_labels_df.to_csv(options.s+'/matrix/Accesson_peaks.csv', sep='\t')
    groups = list(set(y_predict))
    coAccess_matrix = numpy.array([normal[:, numpy.where(y_predict==x)[0]].sum(axis=1) for x in groups]).T
    if (options.format=='mtx')&(len(cells)>ncell_cut):
        coAccess_df = pandas.DataFrame(coAccess_matrix[0], index=cells, columns=groups)
    else:
        coAccess_df = pandas.DataFrame(coAccess_matrix, index=cells, columns=groups)
    coAccess_df.to_csv(options.s+'/matrix/Accesson_reads.csv', sep=',')
    return
#
#
def filtering(raw_matrix):
    matrix = numpy.array([x/x.sum()*1000000 for x in raw_matrix]) 
    gini = lambda x : 0.5 * numpy.abs(numpy.subtract.outer(x, x)).mean() / numpy.mean(x)
    if raw_matrix.shape[0]<=10000:
        coef_gini = numpy.array([gini(x) for x in matrix.T])
    else:
        coef_gini = numpy.array([len(numpy.where(x>0)[0]) for x in matrix.T])
    means = matrix.mean(axis=0)
    mean_cutoff = (means.mean() - numpy.std(means))*2
    kk, bb = numpy.polyfit(numpy.log10(means), coef_gini, 1)
    yy = coef_gini + numpy.log10(means) * (-kk)
    index = numpy.where((yy>bb)&(means>mean_cutoff))[0]
    matrix_filtered = matrix[:, index]
    return matrix_filtered
#
#
def weighted_tsne(matrix, clusters, options):
    cell_types = list(set(clusters['cluster'].values))
    adjusted = copy.deepcopy(matrix.values)
    for ctype in cell_types:
        cluster_cells = numpy.where(clusters==ctype)[0]
        weight = adjusted[cluster_cells, :].mean(axis=0) * float(options.wt)
        adjusted[cluster_cells, :] = numpy.array([x+weight for x in adjusted[cluster_cells, :]])
    tsne_result = TSNE(n_components=2, random_state=int(options.rs)).fit_transform(adjusted)
    return tsne_result
#
#
def cluster_plot(options):
    if not os.path.exists(options.s+'/result'): os.popen('mkdir '+options.s+'/result')
    if not os.path.exists(options.s+'/figure'): os.popen('mkdir '+options.s+'/figure')
    reads_df = pandas.read_csv(options.s+'/matrix/Accesson_reads.csv', sep=',', index_col=0,
                   engine='c', na_filter=False, low_memory=False)
    if options.filter=='yes':
        normal = filtering(reads_df.values)
    else:
        normal = numpy.array([x/x.sum()*1000000 for x in reads_df.values])
    if options.norm=='zscore':
        normal = stats.zscore(normal, axis=1)
    elif options.norm=='probability':
        normal = numpy.array([x/x.sum() for x in normal])
    else:
        print "Error: --norm should be zscore or probability"
    matrix = pandas.DataFrame(normal, index=reads_df.index, columns=['pc'+str(x+1) for x in xrange(len(normal.T))])
    connect = kneighbors_graph(matrix, n_neighbors=20, include_self=False)
    connectivity = 0.5*(connect + connect.T)
    if int(options.nc)==0:
        n_clust, clusters = subroutines.predict_cluster(matrix, connectivity.todense())
        print "predicted number of cell-clusters: ", n_clust
        clusters.to_csv(options.s+'/result/louvain_cluster_by_Accesson.csv', sep='\t')
        tsne_result = weighted_tsne(matrix, clusters, options)
        subroutines.plot_cluster(options, clusters, n_clust, tsne_result, 'louvain_cluster_by_Accesson.pdf')
    else:
        n_clust = int(options.nc)
        clusters = subroutines.knn_cluster(options, matrix, n_clust, connectivity, "KNN_cluster_by_Accesson.csv")
        tsne_result = weighted_tsne(matrix, clusters, options)
        subroutines.plot_cluster(options, clusters, n_clust, tsne_result, 'KNN_cluster_by_Accesson.pdf')
#
    subroutines.plot_tSNE(options, reads_df, tsne_result, "TSNE_by_Accesson.pdf")
    tsne_df = pandas.DataFrame(tsne_result, index=reads_df.index, columns=['TSNE1', 'TSNE2'])
    tsne_df.to_csv(options.s+'/result/TSNE_by_Accesson.csv', sep='\t')
#   
    if options.hc=='yes':
        clipped = numpy.clip(normal, -1, 6)
        corr = numpy.corrcoef(clipped)
        corr_df = pandas.DataFrame(corr, index=reads_df.index, columns=reads_df.index)
        subroutines.hierarchy_cluster(options, corr_df, n_clust, "cell_cell_correlation_by_Accesson.png",
                                  "Hierarchical_cluster_by_Accesson.csv")
    return
#
#
#t1=time.time()
#
build_accesson(options)
cluster_plot(options)
#
#t2=time.time()
#print t2-t1
#
#
