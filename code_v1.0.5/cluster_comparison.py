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
import scipy.special
#
#
pts = OptionParser()
usage = "Compare clustering method\nusage: %prog --c1 cell_info.csv --c2 KNN_cluster.csv"
opts = OptionParser(usage=usage, version="%prog 1.0.5")
opts.add_option("--c1", help="cell_info.csv in <data> folder, or cluster_by_XXX.csv in <result> folder")
opts.add_option("--c2", help="cell cluster file different with c1")
options, arguments = opts.parse_args()
#
#
cluster_df1 = pandas.read_csv(options.c1, sep='\t', index_col=0)
cluster_df2 = pandas.read_csv(options.c2, sep='\t', index_col=0)
if 'notes' in cluster_df1.columns.values: cluster_df1['cluster']=cluster_df1['notes']
if 'notes' in cluster_df2.columns.values: cluster_df2['cluster']=cluster_df2['notes']
clusters1 = list(set(cluster_df1['cluster'].values))
clusters2 = list(set(cluster_df2['cluster'].values))
contingency = numpy.zeros((len(clusters1), len(clusters2)), dtype=int)
for i,clust_i in enumerate(clusters1):
    index_i = numpy.where(cluster_df1['cluster'].values==clust_i)[0]
    cells_i = cluster_df1.index.values[index_i]
    for j,clust_j in enumerate(clusters2):
        index_j = numpy.where(cluster_df2['cluster'].values==clust_j)[0]
        cells_j = cluster_df2.index.values[index_j]
        overlap = list(set(cells_i).intersection(set(cells_j)))
        contingency[i,j] = len(overlap)
contingency_df = pandas.DataFrame(contingency, index=clusters1, columns=clusters2)
print contingency_df
#
sum_ai = 0
for ai in contingency.sum(axis=1):
    sum_ai += scipy.special.binom(ai, 2)
sum_bj = 0
for bj in contingency.sum(axis=0):
    sum_bj += scipy.special.binom(bj, 2)
sum_nij = 0
for row in contingency:
    for nij in row:
        sum_nij += scipy.special.binom(nij, 2)
n_binom = scipy.special.binom(contingency.sum(), 2)
ari = (sum_nij - sum_ai*sum_bj/n_binom) / (0.5*(sum_ai+sum_bj)
      - sum_ai*sum_bj/n_binom)
print 'ARI=', ari
#
#
#
#
#
