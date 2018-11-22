#!/usr/bin/python
import warnings
warnings.filterwarnings("ignore")
#
import numpy
import pandas
from optparse import OptionParser
from scipy import stats
import os
import sys
from multiprocessing import Pool
import subroutines
#
#
opts = OptionParser()
usage = "Enriched motifs/genes/peaks of a cluster (batch)\nusage: %prog -s source_folder --cfile cluster.csv --cluster 1 --vs 2,3"
opts = OptionParser(usage=usage, version="%prog 2.1")
opts.add_option("-s", help="Source folder.")
opts.add_option("--cfile", help="cluster.csv file of a clustering method, e.g. KNN_cluster_by_Accesson.csv in result folder")
opts.add_option("--cluster", help="The cluster for specific markers analysis, can be {0, 1, ..., nCluster}, or a batch of clusters like 0,2,3")
opts.add_option("--vs", default='all', help="vs which clusters to search specific markers for target clusters, e.g. 1,4,2, default=all")
opts.add_option("--pvalue", default=0.001, help='P-value threshold for specific markers, default=0.001')
opts.add_option("--fold", default=2, help='Fold change cutoff of specific markers, default=2')
opts.add_option("--motif", default='no', help='Whether to search differential motifs for target cluster, default=no.')
opts.add_option("--gene", default='no', help='Whether to search differential genes for target cluster, default=no.')
options, arguments = opts.parse_args()
#
#
def group_cells(options):
    deviation = pandas.read_csv(options.s+'/result/deviation_chromVAR.csv', sep=',', index_col=0,
                   engine='c', na_filter=False, low_memory=False)
    cluster_df = pandas.read_csv(options.cfile, sep='\t', index_col=0)
    kCluster = options.cluster.split(',')
    if options.vs!='all': vsCluster = options.vs.split(',')
    if 'cluster' not in cluster_df.columns.values:
        cluster_df['cluster'] = cluster_df['notes']
    else:
        kCluster = map(int, kCluster)
        if options.vs!='all': vsCluster = map(int, vsCluster)
    if options.vs=='all': vsCluster = list(set(cluster_df['cluster'].values)-set(kCluster))
    cluster_df = cluster_df.loc[deviation.columns.values]
    cell_inCluster = cluster_df.loc[cluster_df['cluster'].isin(kCluster)].index.values
    cell_outCluster = cluster_df.loc[cluster_df['cluster'].isin(vsCluster)].index.values
    print len(cell_inCluster), len(cell_outCluster)
    return cell_inCluster, cell_outCluster
#
#
def get_diff_motifs(options, subname, cell_inCluster, cell_outCluster):
    deviation = pandas.read_csv(options.s+'/result/deviation_chromVAR.csv', sep=',', index_col=0,
                   engine='c', na_filter=False, low_memory=False)
    dev_in, dev_out = deviation[cell_inCluster].values, deviation[cell_outCluster].values
    mean_in, mean_out = dev_in.mean(axis=1), dev_out.mean(axis=1)
    delta = mean_in - mean_out
    ttest, pvalues = stats.ttest_ind(dev_in.T, dev_out.T, equal_var=False)
    matrix = numpy.array([mean_in, mean_out, delta, pvalues]).T
    columns = ['dev_inCluster', 'dev_outCluster', 'd_dev', 'P-value']
    compare_df = pandas.DataFrame(matrix, index=deviation.index, columns=columns)
    compare_df = compare_df.loc[compare_df['d_dev']>=1]
    compare_df = compare_df.loc[compare_df['P-value']<=float(options.pvalue)]
    compare_df = compare_df.sort_values(by=['P-value'])
    compare_df.to_csv(options.s+'/result/motifs_of_cluster_'+subname+'.csv', sep='\t')
    return
#
#
def get_diff_genes(options, subname, cell_inCluster, cell_outCluster):
    expr = pandas.read_csv(options.s+'/matrix/genes_scored_by_peaks.csv', sep=',', index_col=0,
                   engine='c', na_filter=False, low_memory=False)
    dev_in, dev_out = expr[cell_inCluster].values, expr[cell_outCluster].values
    mean_in, mean_out = dev_in.mean(axis=1), dev_out.mean(axis=1)
    delta = (mean_in + 1e-4) / (mean_out + 1e-4)
    ttest, pvalues = stats.ttest_ind(dev_in.T, dev_out.T, equal_var=False)
    matrix = numpy.array([mean_in, mean_out, delta, pvalues]).T
    columns = ['expr_inCluster', 'expr_outCluster', 'fold', 'P-value']
    compare_df = pandas.DataFrame(matrix, index=expr.index, columns=columns)
    compare_df = compare_df.loc[compare_df['fold']>=float(options.fold)]
    compare_df = compare_df.loc[compare_df['P-value']<=float(options.pvalue)]
    compare_df = compare_df.sort_values(by=['P-value'])
    compare_df.to_csv(options.s+'/result/genes_of_cluster_'+subname+'.csv', sep='\t')
    return
#
#
subname = '_'.join(options.cluster.split(','))
if options.vs!='all':
    subname += '_VS_' + '_'.join(options.vs.split(','))
subroutines.specific_peak(options, subname)
cells_in, cells_out = group_cells(options)
if options.motif=='yes':
    get_diff_motifs(options, subname, cells_in, cells_out)
if options.gene=='yes':
    get_diff_genes(options, subname, cells_in, cells_out)
#
#
#
#
