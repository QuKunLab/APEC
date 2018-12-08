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
from sklearn.cluster import MiniBatchKMeans
import scipy.spatial.distance
from scipy import sparse
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.decomposition import PCA
from sklearn.decomposition import TruncatedSVD
from sklearn.manifold import TSNE
from sklearn import cluster
from sklearn.neighbors import kneighbors_graph
from scipy import stats
import subroutines
#
#
opts = OptionParser()
usage = "Clustering by accesson\nusage: %prog -s project --nc 0"
opts = OptionParser(usage=usage, version="%prog 1.0")
opts.add_option("-s", help="The project folder.")
opts.add_option("--nc", default=0, 
                help="Number of cell-groups to be clustered, default=0, which means it will be predicted by Louvain algorithm")
opts.add_option("--npeak", default=5, help="Minimum number of peaks for each accesson, default=5, decrease it if too few accessons")
opts.add_option("--ngroup", default=600, help="Number of accessons, default=600")
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
def normalize(options):
    reads = pandas.read_csv(options.s+'/matrix/filtered_reads.csv', sep=',', index_col=0,
                   engine='c', na_filter=False, low_memory=False)
    normal = []
    for cell in reads.index.values:
        reads_per_cell = reads.loc[cell].values
        reads_per_cell = reads_per_cell * 10000 / reads_per_cell.sum()
        normal.append(numpy.log2(reads_per_cell + 1))
    normal = numpy.array(normal)
    normal_df = pandas.DataFrame(normal, index=reads.index.values, columns=reads.columns.values)
    normal_df.to_csv(options.s+'/matrix/normal_reads.csv', sep=',')
    return
#
#
def lonely(peak_labels_df, cutoff):
    peak_labels = peak_labels_df.values
    cluster_index, lonely_index = {}, {}
    for i in range(0, peak_labels.max()):
        index = numpy.where(peak_labels==i)[0]
        if len(index)>=cutoff:
            cluster_index[i] = peak_labels_df.index.values[index]
        elif len(index)>0:
            lonely_index[i] = peak_labels_df.index.values[index]
    lonely_peaks = [y for x in lonely_index.values() for y in x]
    return lonely_peaks, cluster_index
#
#
def build_accesson(options):
    normal_df = pandas.read_csv(options.s+'/matrix/normal_reads.csv', sep=',', index_col=0,
                   engine='c', na_filter=False, low_memory=False)
    normal = normal_df.values
    npc = min(int(options.npc), len(normal[:,0]), len(normal[0,:]))
    ngroups = int(options.ngroup)
    pca_result = PCA(n_components=npc, svd_solver='full').fit_transform(normal.T)
    connectivity = kneighbors_graph(pca_result, n_neighbors=10, include_self=False)
    connectivity = 0.5*(connectivity + connectivity.T)
    ward_linkage = cluster.AgglomerativeClustering(n_clusters=ngroups, linkage='ward', connectivity=connectivity)
    ward_linkage.fit(pca_result)
    y_predict = ward_linkage.labels_.astype(numpy.int)
    peak_labels_df = pandas.DataFrame(y_predict, index=normal_df.columns.values, columns=['group'])
    peak_labels_df.to_csv(options.s+'/matrix/Accesson_peaks.csv', sep='\t')
    lonely_peaks, cluster_index = lonely(peak_labels_df, int(options.npeak))
    print 'number of lonely peaks:', len(lonely_peaks)
    print 'number of valid accessons:', len(cluster_index)
    coAccess_matrix, groups, STDs = [], [], []
    for key in cluster_index.keys():
        peaks = cluster_index[key]
        average = normal_df[peaks].values.sum(axis=1) #/ len(peaks)
        coAccess_matrix.append(average)
        groups.append(key)
    coAccess_matrix, groups = numpy.array(coAccess_matrix).T, numpy.array(groups)
    coAccess_df = pandas.DataFrame(coAccess_matrix, index=normal_df.index.values, columns=groups)
    coAccess_df.to_csv(options.s+'/matrix/Accesson_reads.csv', sep=',')
    return
#
#
def predict_nc(corr, options):
    n_clust = int(options.nc)
    if n_clust==0:
        n_clust = subroutines.predict_cluster(corr)
        print "predicted number of cell-clusters: ", n_clust
        options.nc = n_clust
    return n_clust
#
#
def cluster_plot(options):
    if not os.path.exists(options.s+'/result'): os.popen('mkdir '+options.s+'/result')
    if not os.path.exists(options.s+'/figure'): os.popen('mkdir '+options.s+'/figure')
    reads_df = pandas.read_csv(options.s+'/matrix/Accesson_reads.csv', sep=',', index_col=0,
                   engine='c', na_filter=False, low_memory=False)
#
    normal = numpy.array([x/x.sum()*1000000 for x in reads_df.values])
    reads_df = pandas.DataFrame(normal, index=reads_df.index, columns=reads_df.columns)
    npc = min(200, len(normal[:,0]), len(normal[0,:]))
    pca_result = PCA(n_components=npc, svd_solver='full').fit_transform(reads_df)
    corr = pandas.DataFrame(numpy.corrcoef(pca_result), index=reads_df.index, columns=reads_df.index)
#
    n_clust = predict_nc(corr, options)
#
    pca_result, tsne_result = subroutines.PCA_tSNE(options, corr, "PCA_by_Accesson.pdf", "TSNE_by_Accesson.pdf")
    tsne_df = pandas.DataFrame(tsne_result, index=corr.index.values, columns=['TSNE1', 'TSNE2'])
    tsne_df.to_csv(options.s+'/result/TSNE_by_Accesson.csv', sep='\t')
    pca_df = pandas.DataFrame(pca_result, index=corr.index.values, columns=['PC'+str(i) for i in range(0, int(options.npc))])
    if options.space=='tsne':
        matrix_df = tsne_df
    else:
        matrix_df = pca_df
    subroutines.plot_knn_cluster(options, matrix_df, n_clust, tsne_result, "KNN_cluster_by_Accesson.pdf",
                                 "KNN_cluster_by_Accesson.csv")
#
    if options.hc=='yes':
        subroutines.hierarchy_cluster(options, corr, n_clust, "cell_cell_correlation_by_Accesson.png",
                                  "Hierarchical_cluster_by_Accesson.csv")
        KNN_df = pandas.read_csv(options.s+"/result/KNN_cluster_by_Accesson.csv", sep='\t', index_col=0)
        subroutines.heatmap_compare(options, corr, KNN_df, "HC_KNN_compare_by_Accesson.png")
    return
#
#
#
#
normalize(options)
build_accesson(options)
cluster_plot(options)
#
#
#
