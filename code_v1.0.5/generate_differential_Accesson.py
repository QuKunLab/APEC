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
usage = "Enriched accessons of a cluster (batch) \nusage: %prog -s project --cfile cluster.csv --cluster 1 --vs 2,3"
opts = OptionParser(usage=usage, version="%prog 1.0")
opts.add_option("-s", help="The project folder.")
opts.add_option("--cfile", help="cluster.csv file of a clustering method, e.g. louvain_cluster_by_Accesson.csv in result folder")
opts.add_option("--cluster", help="the cluster(s) for specific-TF analysis, can be {0, 1, ..., nCluster}, or a batch of clusters like 0,2,3")
opts.add_option("--vs", default='all', help="vs which cluster(s) to search specific TF for target cluster(s), e.g. 1,4,2, default=all")
opts.add_option("--pvalue", default=0.001, help='P-value threshold for specific peaks, default=0.001')
opts.add_option("--fold", default=2, help='Fold change cutoff of specific peaks, default=2')
options, arguments = opts.parse_args()
#
#
subname = '_'.join(options.cluster.split(','))
if options.vs!='all':
    subname += '_VS_' + '_'.join(options.vs.split(','))
subroutines.specific_accesson(options, subname)
#
#