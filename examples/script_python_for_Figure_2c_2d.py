import os,sys,numpy
import pandas
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn
import scipy.stats
#
#
def get_matrix(markers, gene_csv, cluster_csv):
    gene_df = pandas.read_csv(gene_csv, sep=',', index_col=0, 
                              engine='c', na_filter=False, low_memory=False).T
    cluster_df = pandas.read_csv(cluster_csv, sep='\t', index_col=0)
    clusters = list(set(cluster_df['cluster'].values))
    matrix, reads = [], []
    for cluster in clusters:
        cells = cluster_df.loc[cluster_df['cluster']==cluster].index.values
        expr = gene_df.loc[markers, cells].values.mean(axis=1)
        matrix.append(expr)
        reads.append(gene_df[cells].values.mean(axis=1))
    matrix = numpy.array(matrix)
    matrix_df = pandas.DataFrame(matrix, index=['c_'+str(x) for x in clusters], columns=markers)
    matrix_df.T.to_csv('cluster_vs_genes.csv', sep=',')
    return
#
#
def plot_heatmap(matrix_csv, cluster_order, gene_order, out_fig):
    matrix_df = pandas.read_csv(matrix_csv, sep=',', index_col=0)
    zscore = scipy.stats.zscore(matrix_df.values, axis=1)
    zscore = scipy.stats.zscore(zscore, axis=0)
    z_df = pandas.DataFrame(zscore, index=matrix_df.index, columns=matrix_df.columns)
    xx = [0, 520, 1040, 1260, 1380]
    yy = [0, 520, 1040, 1460, 1880]
    xclust = [500, 500, 200, 100, 100]
    ygene = [500, 500, 400, 400, 300]
    fig0 = plt.figure(figsize=(10,17))
    x_clust = 0
    for ic,cluster in enumerate(cluster_order):
        y_gene = 0
        for ig,gene in enumerate(gene_order):
            ax = plt.subplot2grid((2200, 1500),
                 (yy[ig], xx[ic]), rowspan=ygene[ig], colspan=xclust[ic])
            y_gene += len(gene)
            im = ax.imshow(z_df.loc[gene,cluster], cmap='bwr', aspect='auto',
                           interpolation='none', vmax=2.0, vmin=-2)
            if y_gene==len(z_df.index.values):
                ax.set_xticks(numpy.arange(0, len(cluster)))
                ax.set_xticklabels(cluster, fontsize=15)
            else:
                ax.set_xticklabels([])
            if ic==0: 
                ax.set_yticks(numpy.arange(0, len(gene)))
                ax.set_yticklabels(gene, fontsize=15)
            else:
                ax.set_yticklabels([])
        x_clust += len(cluster)
    plt.savefig(out_fig, bbox_inches='tight')
    plt.close()
    return
#
#
def pearson(csv, clusts, outfig):
    matrix_df = pandas.read_csv(csv, sep=',', index_col=0)
    matrix_df = matrix_df[clusts]
    seaborn.set_context('poster')
    corr = matrix_df.corr()
    seaborn.clustermap(corr, method='ward', cmap='YlOrRd')
    plt.savefig(outfig, bbox_inches='tight')
    plt.close()
    return
#
#
#### Please run script_for_project02.py first !!!!
#### Please change "./project02/" to the path that you placed project02 !!!!  
#
#
gene_score_file = './project02/matrix/genes_scored_by_TSS_peaks.csv'
#
if not os.path.exists(gene_score_file):
    print('Error !!!!')
    print('Please run script_for_project02.py first !!!!')
    print('Please change "./project02/" to the path that you placed project02 !!!!')
    sys.exit()
#
#
#### These are marker genes reported by Preissl et al.
markers = ['Neurod1', 'Neurod2', 'Neurod6', 'Tbr1', 'Slc17a7',
           'Gad1', 'Gad2', 'Slc32a1', 'Dlx1', 'Dlx5',
           'Bcan', 'Aldh1l1', 'Slc1a2', 'Slc1a3',
           'Mobp', 'Mag', 'Plp1', 'Mog',
           'C1qb', 'Ctss', 'Spi1']
#
#
get_matrix(markers, gene_score_file, 'project02-result/cluster_by_APEC.csv')
#
c_order = [['c_6', 'c_5', 'c_2', 'c_3', 'c_8'], 
           ['c_0', 'c_10', 'c_1', 'c_12', 'c_11'],
           ['c_4', 'c_13'], ['c_9'], ['c_7']]
g_order = [['Neurod1', 'Neurod2', 'Neurod6', 'Tbr1', 'Slc17a7'],
           ['Gad1', 'Gad2', 'Slc32a1', 'Dlx1', 'Dlx5'], 
           ['Bcan', 'Aldh1l1', 'Slc1a2', 'Slc1a3'],
           ['Mobp', 'Mag', 'Plp1', 'Mog'],
           ['C1qb', 'Ctss', 'Spi1']]
plot_heatmap('cluster_vs_genes.csv', c_order, g_order, 'cluster_vs_gene.png')
#
clusts = ['c_'+str(x) for x in range(0,14)]
pearson('cluster_vs_genes.csv', clusts, 'cluster_pearson_corr.png')
#
#