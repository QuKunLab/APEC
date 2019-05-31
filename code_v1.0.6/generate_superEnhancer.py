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
import sys
import scipy.stats
#
#
opts = OptionParser()
usage = "Search super enhancers\nusage: %prog -s project"
opts = OptionParser(usage=usage, version="%prog 1.0")
opts.add_option("-s", help="The project folder.")
options, arguments = opts.parse_args()
#
def search_supper(options):
    supper_range = 1e6
    peaks = numpy.array([x.split()[:3] for x in open(options.s+'/peak/top_peaks.bed').readlines()])
    chroms = list(set(peaks[:, 0]))
    chroms.sort()
    accesson_df = pandas.DataFrame.from_csv(options.s+'/matrix/Accesson_peaks.csv', sep='\t')
    accessons = list(set(accesson_df['group'].values))
    all_suppers, all_locate, all_base = [], [], []
    for access in accessons:
        peaks_group = accesson_df.loc[accesson_df['group']==access].index.values
        peaks_index = [int(x[4:]) for x in peaks_group]
        peaks_info = peaks[peaks_index, :]
        peaks_dict = {x:[] for x in chroms}
        peaks_label = {x:[] for x in chroms}
        for ip,pp in enumerate(peaks_info):
            peaks_dict[pp[0]].append((int(pp[1])+int(pp[2]))/2)
            peaks_label[pp[0]].append(peaks_index[ip])
        peaks_label = {x:numpy.array(peaks_label[x]) for x in peaks_label.keys()}
        suppers, locate, base = [], [], []
        for chrom in peaks_dict.keys():
            if len(peaks_dict[chrom])>=2:
                position = numpy.array(peaks_dict[chrom])
                supper = [numpy.where(abs(position-x)<=supper_range)[0] for x in position if 
                      len(numpy.where(abs(position-x)<=supper_range)[0])>1]
                supper = ['-'.join(map(str, x)) for x in supper]
                supper = list(set(supper))
                supper = [list(map(int, x.split('-'))) for x in supper]
                supper_peaks = numpy.array(['-'.join(map(str, peaks_label[chrom][x])) for x in supper])
                if len(supper_peaks)>0:
                    for ii,ss in enumerate(supper_peaks):
                        peaks_in = map(int, ss.split('-'))
                        start = peaks[peaks_in, 1:].astype(int).min()
                        end = peaks[peaks_in, 1:].astype(int).max()
                        delta = numpy.array([abs(peaks_in[i+1]-x) for i,x in enumerate(peaks_in[:-1])])
                        close = numpy.where(delta<=2)[0]
                        percent = len(close)/float(len(delta))
                        if (len(delta)>=2) & (percent>0.5) :
                            suppers.append(ss)
                            locate.append(access)
                            base.append(chrom+':'+str(start)+'-'+str(end))
        all_suppers.extend(suppers)
        all_locate.extend(locate)
        all_base.extend(base)
    supper_df = pandas.DataFrame(numpy.array([all_locate, all_base]).T, index=all_suppers, columns=['peaks', 'position'])
    supper_df.to_csv(options.s+'/result/potential_super_enhancer.csv', sep='\t')
    return
#
#
#
#
search_supper(options)
#
#
#
#
