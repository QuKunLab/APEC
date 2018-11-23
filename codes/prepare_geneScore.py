#!/usr/bin/python
import warnings
warnings.filterwarnings("ignore")
import os
import numpy
import pandas
import sys
from optparse import OptionParser
#
#
opts = OptionParser()
usage = "Evaluate gene score by TSS peaks\nusage: %prog -s project --gtf hg19.gtf --distal 20000"
opts = OptionParser(usage=usage, version="%prog 1.0")
opts.add_option("-s", help="The project folder.")
opts.add_option("--gtf", default='../reference/hg19_RefSeq_genes.gtf', 
                help="gtf file for genome, default=../reference/hg19_RefSeq_genes.gtf")
opts.add_option("--distal", default=20000, 
                help="distal region around TSS for peak searching, default=20000")
options, arguments = opts.parse_args()
#
#
def get_tss_region(options):
    mm10_df = pandas.read_csv(options.gtf, sep='\t', index_col=0)
    genes = list(set(mm10_df['name2']))
    genes.sort()
    mm10_df.index = mm10_df['name']
    names, tss = [], []
    for symbol in genes:
        sub_df = mm10_df.loc[mm10_df['name2']==symbol]
        if len(sub_df.index.values)>=1:
            chrom = list(set(sub_df['chrom'].values))
            strand = list(set(sub_df['strand'].values))
            if len(chrom)==1:
                if strand[0]=='+':
                    starts = list(set(map(str, sub_df['txStart'].values)))
                    start = ','.join(starts)
                elif strand[0]=='-':
                    starts = list(set(map(str, sub_df['txEnd'].values)))
                    start = ','.join(starts)
                names.append(symbol)
                tss.append([chrom[0], start])
    tss = numpy.array(tss)
    tss_df = pandas.DataFrame(tss, index=names, columns=['chrom', 'tss'])
    tss_df.to_csv(options.s+'/peak/genes_tss_region.csv', sep='\t')
    return
#
#
def get_tss_peaks(options):
    peaks = [[x.split()[0], (int(x.split()[1])+int(x.split()[2]))/2]
             for x in open(options.s+'/peak/annotate_peak.bed').readlines()]
    peaks_df = pandas.DataFrame(peaks, index=['peak'+str(x) for x in xrange(len(peaks))], 
                                columns=['chrom', 'center'])
    tss_df = pandas.read_csv(options.s+'/peak/genes_tss_region.csv', sep='\t', index_col=0)
    for gene in tss_df.index.values:
        chrom, tsses = tss_df.ix[gene, 'chrom'], tss_df.ix[gene, 'tss']
        tsses = map(int, tsses.split(','))
        chr_peaks = peaks_df.loc[peaks_df['chrom']==chrom]
        proxim_peaks, distal_peaks = [], []
        for tss in tsses:
            peaks1 = chr_peaks.loc[abs(chr_peaks['center']-tss)<=2000].index.values
            peaks2 = chr_peaks.loc[abs(chr_peaks['center']-tss)<=int(options.distal)].index.values
            proxim_peaks.extend(peaks1)
            distal_peaks.extend(peaks2)
        proxim_peaks = list(set(proxim_peaks))
        distal_peaks = list(set(distal_peaks)-set(proxim_peaks))
        if len(proxim_peaks)==0: proxim_peaks = ['NONE'] 
        if len(distal_peaks)==0: distal_peaks = ['NONE']
        proxim_peaks = ';'.join(proxim_peaks)
        tss_df.ix[gene, 'proximal'] = proxim_peaks
        distal_peaks = ';'.join(distal_peaks)
        tss_df.ix[gene, 'distal'] = distal_peaks
    tss_df.to_csv(options.s+'/peak/genes_tss_peaks.csv', sep='\t')
    return
#
#
def get_score_from_peaks(options):
    tss_df = pandas.read_csv(options.s+'/peak/genes_tss_peaks.csv', sep='\t', index_col=0)
    reads_df = pandas.read_csv(options.s+'/matrix/filtered_reads.csv', sep=',', index_col=0)
    all_peaks = reads_df.columns.values
    genes, score = [], []
    for igene,gene in enumerate(tss_df.index.values):
        distal = tss_df.loc[gene, 'distal'].split(';')
        proximal = tss_df.loc[gene, 'proximal'].split(';')
        distal = list(set(distal).union(set(proximal)))
        distal = list(set(distal).intersection(set(all_peaks)))
        if len(distal)>0:
            signal = reads_df[distal].values.mean(axis=1)
            genes.append(gene)
            score.append(signal)
    score = numpy.array(score)
    score_df = pandas.DataFrame(score, index=genes, columns=reads_df.index)
    score_per_cell = score.sum(axis=0)
    R_wave = [numpy.log(x*10000.0/score_per_cell[i]+1) for i,x in enumerate(score.T)]
    R_wave = numpy.array(R_wave)
    normal_df = pandas.DataFrame(R_wave.T, index=genes, columns=reads_df.index)
    normal_df.to_csv(options.s+'/matrix/genes_scored_by_peaks.csv', sep=',')
    return
#
#
get_tss_region(options)
get_tss_peaks(options)
get_score_from_peaks(options)
#
#
