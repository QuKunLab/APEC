#!/usr/bin/python
import warnings
warnings.filterwarnings("ignore")
#
import os,numpy,pandas,random,time,copy
from optparse import OptionParser
import scipy.spatial.distance
import scipy.io,scipy.sparse,scipy.stats
from sklearn.decomposition import PCA
from sklearn import cluster
import sklearn.utils.sparsefuncs
from sklearn.neighbors import kneighbors_graph
from sklearn.decomposition.truncated_svd import TruncatedSVD
from MulticoreTSNE import MulticoreTSNE as McTSNE
import subroutines
#
#
opts = OptionParser()
usage = "Clustering by accesson\nusage: %prog -s project"
opts = OptionParser(usage=usage, version="%prog 1.0.6")
opts.add_option("-s", help="The project folder.")
opts.add_option("--nc", default=0, help="Number of cell clusters, default=0, i.e. predicted by Louvain algorithm")
opts.add_option("--ngroup", default=600, help="Number of accessons, default=600")
opts.add_option("--npc", default=40, help="Number of PCs used for peak grouping, default=40")
opts.add_option("--ncell", default=10000, help="Maximum cell-number for dataset that uses exact PCA, default=10000. "
                +"Dataset with more cells will use truncated SVD reduction.")
opts.add_option("--norm", default='zscore', help="=zscore or probability; set normalization method, default=zscore.")
opts.add_option("--rs", default=0, help="Random state of tSNE analysis, default=0")
opts.add_option("--wt", default=0.2, help="Weight coefficient of tSNE analysis, default=0.2")
opts.add_option("--hc", default='no', help="Run hierarchical clustering or not, default=no.")
options, arguments = opts.parse_args()
#
#
def build_accesson(options):
    ngroups, ncell_cut = int(options.ngroup), int(options.ncell)
    reads = scipy.io.mmread(options.s+'/matrix/filtered_reads.mtx')
    reads = scipy.sparse.csr_matrix(reads).T*1.0
    cells = pandas.read_csv(options.s+'/matrix/filtered_cells.csv', sep='\t', index_col=0,
                            engine='c', na_filter=False, low_memory=False)
    cells = cells.index.values
    peaks = ['peak'+str(x) for x in range(0, reads.shape[0])]
    scale = numpy.array(10000.0 / reads.sum(axis=0))[0]
    sklearn.utils.sparsefuncs.inplace_column_scale(reads, scale)
    reads.data = numpy.log2(reads.data+1)
    npc = min(int(options.npc), reads.shape[0], reads.shape[1])
    if len(cells)>ncell_cut:
        pca_result = TruncatedSVD(n_components=npc, algorithm='arpack', random_state=0).fit_transform(reads)
    else:
        pca_result = PCA(n_components=npc, svd_solver='arpack', random_state=0).fit_transform(reads.A)
    ward_linkage = cluster.AgglomerativeClustering(n_clusters=ngroups, linkage='ward')
    y_predict = ward_linkage.fit_predict(pca_result)
    peak_labels_df = pandas.DataFrame(y_predict, index=peaks, columns=['group'])
    peak_labels_df.to_csv(options.s+'/matrix/Accesson_peaks.csv', sep='\t')
    groups = list(set(y_predict))
    coAccess_matrix = numpy.array([reads[numpy.where(y_predict==x)[0],:].sum(axis=0) for x in groups])
    coAccess_matrix = coAccess_matrix[:,0,:].T
    coAccess_df = pandas.DataFrame(coAccess_matrix, index=cells, columns=groups)
    coAccess_df.to_csv(options.s+'/matrix/Accesson_reads.csv', sep=',')
    return
#
#
def weighted_tsne(matrix, clusters, options):
    cell_types = list(set(clusters['cluster'].values))
    adjusted = copy.deepcopy(matrix.values)
    for ctype in cell_types:
        cluster_cells = numpy.where(clusters==ctype)[0]
        weight = adjusted[cluster_cells, :].mean(axis=0) * float(options.wt)
        adjusted[cluster_cells, :] = numpy.array([x+weight for x in adjusted[cluster_cells, :]])
    tsne_result = McTSNE(n_components=2, random_state=int(options.rs)).fit_transform(adjusted)
    return tsne_result
#
#
def cluster_plot(options):
    if not os.path.exists(options.s+'/result'): os.popen('mkdir '+options.s+'/result')
    if not os.path.exists(options.s+'/figure'): os.popen('mkdir '+options.s+'/figure')
    reads_df = pandas.read_csv(options.s+'/matrix/Accesson_reads.csv', sep=',', index_col=0,
                   engine='c', na_filter=False, low_memory=False)
    normal = numpy.array([x/x.sum() for x in reads_df.values])
    if options.norm=='zscore':
        normal = scipy.stats.zscore(normal, axis=1)
    matrix = pandas.DataFrame(normal, index=reads_df.index, columns=['pc'+str(x+1) for x in xrange(len(normal.T))])
    connect = kneighbors_graph(matrix, n_neighbors=20, include_self=False)
    connectivity = 0.5*(connect + connect.T)
    if int(options.nc)==0:
        n_clust, clusters = subroutines.predict_cluster(matrix, connectivity.todense())
        print("predicted number of cell-clusters: ", n_clust)
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
build_accesson(options)
cluster_plot(options)
#
#
