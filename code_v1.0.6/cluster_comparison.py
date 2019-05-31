import numpy, pandas
from optparse import OptionParser
import scipy.special
import sklearn.metrics
from sklearn.metrics.cluster import contingency_matrix
#
#
pts = OptionParser()
usage = "Compare clustering method\nusage: %prog --c1 cell_info.csv --c2 KNN_cluster.csv"
opts = OptionParser(usage=usage, version="%prog 1.0.6")
opts.add_option("--c1", help="filtered_cells.csv in <matrix> folder, or cluster_by_XXX.csv in <result> folder")
opts.add_option("--c2", help="cell cluster file different with c1")
options, arguments = opts.parse_args()
#
#
cluster_df1 = pandas.read_csv(options.c1, sep='\t', index_col=0)
cluster_df2 = pandas.read_csv(options.c2, sep='\t', index_col=0)
if 'notes' in cluster_df1.columns.values: cluster_df1['cluster']=cluster_df1['notes']
if 'notes' in cluster_df2.columns.values: cluster_df2['cluster']=cluster_df2['notes']
clusters1 = cluster_df1['cluster'].values
clusters2 = cluster_df2['cluster'].values
ari = sklearn.metrics.adjusted_rand_score(clusters1, clusters2)
nmi = sklearn.metrics.mutual_info_score(clusters1, clusters2)
ami = sklearn.metrics.adjusted_mutual_info_score(clusters1, clusters2)
cont_mat = contingency_matrix(clusters1, clusters2)
print('ARI=', ari)
print('NMI=', nmi)
print('AMI=', ami)
print('contingency matrix:')
print(cont_mat)
#
#
#
#
#
