#
import os
import numpy
import pysam
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn
from multiprocessing import Pool
import random
import pandas
import copy
import scipy.io
from sklearn.manifold import TSNE
from sklearn.neighbors import kneighbors_graph
from sklearn.decomposition import PCA
from sklearn import cluster
import scipy.cluster.hierarchy
import scipy.stats
import networkx
import community
import numba
#
global colors
colors = numpy.array(['pink', 'red', '#377eb8', 'green', 'skyblue', 'lightgreen', 'gold',
                      '#ff7f00', '#000066', '#ff3399', '#a65628', '#984ea3', '#999999',
                      '#e41a1c', '#dede00', 'b', 'g', 'r', 'c', 'm', 'y', 'k'])
#
@numba.jit()
def counts_per_peak(bam_file, peak_bed, out_file):
    peak_info = numpy.loadtxt(peak_bed, 'str', delimiter="\t")
    bamfile = pysam.Samfile(bam_file, "rb")
    outData = {}
    for ipeak in range(0,len(peak_info)):
        start, end = int(peak_info[ipeak][1]), int(peak_info[ipeak][2])
        rList = []
        for read in bamfile.fetch(peak_info[ipeak][0].tolist(), max(0,start-2000), end+2000):
            if read.mapq<30:  sys.exit('not prefiltered correctly')
            if  read.qname in rList:  continue
            if read.is_reverse: 
                position = read.pos + read.alen
            else:
                position = read.pos
            if (position>start) & (position<end):
                barcode = dict(read.tags)['RG']
                if not outData.has_key(barcode): outData[barcode] = numpy.zeros(len(peak_info))
                outData[barcode][ipeak] += 1
                rList.append(read.qname)
        if ipeak%5000==0: 
           print 'Peak:\t'+str(ipeak)
    headers = outData.keys()
    headers.sort()
    outFinal = []
    for barcode in headers:
        outFinal.append(outData[barcode])
    outFinal = numpy.array(outFinal)
    peaks=['peak'+str(i) for i in range(0,len(outFinal[0, :]))]
    matrix_df = pandas.DataFrame(outFinal, index=headers, columns=peaks)
    print matrix_df.shape
    return matrix_df
#
#
def motif_search(info):
    motif, bgFile, threshold, motifFile, motifFasta, outFolder = info[0], info[1], info[2], info[3], info[4], info[5]
    motif_name = motif.split('-')[-1]
    print motif_name
    fimoFile = outFolder + '/' + motif +'.fimo'
    os.popen('fimo --bgfile ' + bgFile + ' --text --thresh ' + threshold + ' --motif ' + motif_name 
        + ' --no-qvalue --verbosity 1 ' + motifFile + ' ' + motifFasta + ' > ' + fimoFile)
    bedFile = outFolder + '/' + motif +'.bed'
    with open(fimoFile) as fimo, open(bedFile, 'w') as bed:
        for line in fimo:
            if line[0]=='#':
                continue
            else:
                words = line.split('\t')
                chrom = words[1].split(':')[0]
                start = int(words[1].split(':')[1].split('-')[0]) + int(words[2])
                end = int(words[1].split(':')[1].split('-')[0]) + int(words[3])
                strand, score, pvalue, name = words[4], words[5], words[6], words[0]
                newLine = chrom+'\t'+str(start)+'\t'+str(end)+'\t'+strand+'\t'+score+'\t'+pvalue+'\t'+name
                print >> bed, newLine
    os.popen('gzip ' + bedFile)
    os.popen('rm ' + fimoFile)
    return
#
#
def batch_fimo(backgroudFile, threshold, motifFile, motifFasta, outFolder, n_processor):
    motifs = []
    with open(motifFile) as mfile:
        for line in mfile:
            words = line.split(' ')
            if words[0] == 'MOTIF':
                info = words[2][:-1]
                info = info.replace('::', '-')
                info = info.replace(' ', '-')
                info = info.replace(':', '-')
                motifs.append(info+'-'+words[1])
    nMotif = len(motifs)
    info = numpy.vstack((motifs, [backgroudFile]*nMotif, [threshold]*nMotif, [motifFile]*nMotif, 
                         [motifFasta]*nMotif, [outFolder]*nMotif)).T
    pool = Pool(n_processor)
    pool.map(motif_search, info)
    pool.close()
    pool.join()
    return
#
#
def assign_TSS_insert(position, matrix, start, end, t, i, rows, tss_info):
    if (float(position)>=start) & (float(position)<end-1) & (t<rows):
        base = position - start
        if len(tss_info[0]) == 3:
            matrix[t, base] += 1
        elif tss_info[i,3] == "-":     
            matrix[t, len(matrix[0])-base-1] += 1    
        else:
            matrix[t, base] += 1
    return matrix
#
#
def TSS_insert_matrix(rows, cols, bam_file, tss_info, half_width):
    matrix = numpy.zeros((rows,cols))
    bam = pysam.Samfile(bam_file, "rb")
    for i in range(0, len(tss_info)):
        center = int(tss_info[i,1])+(int(tss_info[i,2])-int(tss_info[i,1]))/2
        start, end = center-int(half_width), center+int(half_width)
        for read in bam.fetch(str(tss_info[i,0]), max(0,start-2000), end+2000):
            if read.mapq<30:
                continue
            if read.is_reverse:
                continue
            else:
                left, real_len = read.pos+4, abs(read.tlen)-9
                right = left + real_len
                matrix = assign_TSS_insert(left, matrix, start, end, real_len, i, rows, tss_info)
                matrix = assign_TSS_insert(right, matrix, start, end, real_len, i, rows, tss_info)
    return matrix
#
def draw_TSS_insert(tss_file, bam_file, out_file):
    tss_info = numpy.loadtxt(tss_file,'str')
    half_width = 2000
    rows, cols = 1000, half_width*2
    matrix = TSS_insert_matrix(rows, cols, bam_file, tss_info, half_width)
    matrix = matrix.sum(axis=0)
    numpy.savetxt(out_file, matrix, delimiter='\t',fmt='%s')
    fig = plt.figure(figsize=(8.0, 5.0))
    plt.plot(matrix/numpy.mean(matrix[1:200]), 'k.')
    plt.plot(numpy.convolve(matrix,numpy.ones(20),'same')/20/numpy.mean(matrix[1:200]), 'r')
    plt.xlabel('Position relative to center')
    plt.ylabel('Insertions')
    fig.savefig(out_file+'.pdf')
    plt.close(fig)
    return
#
#
def score_peaks(peaks_file, motif_folder, out_file):
    peaks_info = numpy.loadtxt(peaks_file,'str',delimiter="\t")
    files = [motif_folder+'/'+x for x in os.listdir(motif_folder)]
    files.sort()
    outData = numpy.zeros([len(peaks_info),len(files)])
    headers = []
    for i in range(0,len(files)):
        file = files[i]
        fName = file.split('/')[-1].split('.bed')[0].split('.narrowPeak')[0]
        chipData = numpy.loadtxt(file, 'str') 
        headers.append(fName)
        if len(chipData)>0:
            chip = {}
            for line in chipData:
                chrom, start, end = line[0], int(line[1]), int(line[2])
                if not chip.has_key(chrom): 
                    chip[chrom] = [[start, end]]
                else:
                    chip[chrom].append([start, end])
            for j in range(0,len(peaks_info)):
                peakChr, peakStart, peakEnd = peaks_info[j,0], int(peaks_info[j,1]), int(peaks_info[j,2]) 
                try:
                    for site in chip[peakChr]:
                        if (site[0]>=peakStart) & (site[1]<=peakEnd):
                            outData[j,i]+=1
                            break
                except:
                    continue
    TFmotif_df = pandas.DataFrame(outData, index=['peak'+str(i) for i in xrange(len(peaks_info))], columns=headers)
    TFmotif_df.to_csv(out_file, sep=',')
    return
#
#
def QC_table(cell_info_file, work_folder, out_folder):
    cell_info = pandas.read_csv(cell_info_file, sep='\t', index_col=0)
    out_info = copy.deepcopy(cell_info)
    folders = [work_folder+'/'+x for x in os.listdir(work_folder)]
    folders.sort()
    for fold in folders:
        cell_name = fold.split('/')[-1]
        for ffile in os.listdir(fold):
            if ffile[-9:]=='.dups.log':  cell_info.ix[cell_name, 'dups_file'] = fold+'/'+ffile
            elif ffile[-10:]=='.stats.log':  cell_info.ix[cell_name, 'stat_file'] = fold+'/'+ffile
            elif ffile[-10:]=='.RefSeqTSS':  cell_info.ix[cell_name, 'tss_file'] = fold+'/'+ffile
            elif ffile[-8:]=='.map.log':  cell_info.ix[cell_name, 'align_file'] = fold+'/'+ffile
#
    for cell_name in cell_info.index.values:
#        print cell_name, cell_info.ix[cell_name,'align_file']
        alignFile = numpy.loadtxt(cell_info.ix[cell_name,'align_file'],dtype='str',delimiter=',')
        for k in range(0,len(alignFile)):
            if 'overall alignment rate' in alignFile[k]: id_align_rate = k
            if 'reads; of these:' in alignFile[k]: id_reads_number = k
        out_info.ix[cell_name, 'align_rate'] = float(alignFile[id_align_rate].split('%')[0])/100.
        out_info.ix[cell_name, 'all_reads'] = int(alignFile[id_reads_number].split(' reads;')[0])*2 
        statData = pandas.read_csv(cell_info.ix[cell_name,'stat_file'], sep='\t', index_col=0)
        out_info.ix[cell_name, 'chrM_reads'] = statData.ix['chrM', 'ProperPairs']
        out_info.ix[cell_name, 'mapped_reads'] = statData['ProperPairs'].values.sum()
        out_info.ix[cell_name, 'chrM_rate'] = statData.ix['chrM', 'ProperPairs'] / float(statData['ProperPairs'].values.sum())
        dupFile = numpy.loadtxt(cell_info.ix[cell_name,'dups_file'], dtype='str', delimiter=',')
        dupLine = dupFile[1].split('\t')[1:]
        for iword,word in enumerate(dupLine):
            if word=='': dupLine[iword]='0'
        try:
            dupValue = numpy.array(dupLine,float)
        except:
            dupValue = numpy.zeros(len(dupLine))
        out_info.ix[cell_name, 'filtered_reads'] = int(dupValue[0]+dupValue[1]*2)
        out_info.ix[cell_name, 'duplicate_rate'] = round(dupValue[-2], 4)
        out_info.ix[cell_name, 'final_reads'] = int((dupValue[1]-dupValue[5])*2+dupValue[0]-dupValue[4])
        out_info.ix[cell_name, 'tss_reads'] = numpy.loadtxt(cell_info.ix[cell_name,'tss_file']).sum()
    out_file = out_folder + '/cell_info.merged.csv'
    out_info.to_csv(out_file, sep='\t')
    return
#
#
#
#
def predict_cluster(reads_df):
    connect = kneighbors_graph(reads_df, 20, include_self=False)
    connectivity = 0.5*(connect + connect.T).todense()
    graph = networkx.from_numpy_matrix(connectivity)
    partition = community.best_partition(graph)
    return len(list(set(partition.values())))
#
#
#
#
def hierarchy_cluster(options, matrix_df, n_clust, outfig0, outcsv0):
    cellinfo_df = pandas.read_csv(options.s+'/data/cell_info.csv', sep='\t', index_col=0)
    cell_type = cellinfo_df.ix[matrix_df.index.values, 'notes']
    cTypes = list(set(cell_type))
    cTypes.sort()
    cTypeIndex = [numpy.where(cell_type==x) for x in cTypes]
#
    lut = dict(zip(cTypes, colors[:len(cTypes)]))
    row_colors = cell_type.map(lut)
    cmap = seaborn.diverging_palette(h_neg=210, h_pos=350, s=90, l=30, as_cmap=True)
    seaborn.set(font_scale=1.2)
    fig0 = seaborn.clustermap(matrix_df, method='ward', metric='euclidean', cmap=cmap, 
           row_colors=row_colors, figsize=(20,20))
    plt.setp(fig0.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    plt.setp(fig0.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)
    fig0.savefig(options.s + '/figure/' + outfig0, bbox_inches='tight')
    z_linkage = scipy.cluster.hierarchy.linkage(matrix_df, method='ward', metric='euclidean')
    f_cluster = scipy.cluster.hierarchy.fcluster(z_linkage, n_clust, criterion='maxclust')
    clust_df = pandas.DataFrame(f_cluster, index=matrix_df.index.values, columns=['cluster'])
    clust_df.to_csv(options.s+'/result/'+outcsv0, sep='\t')
    return
#
#
#
#
def PCA_tSNE(options, matrix_df, outfig1, outfig2):
    cellinfo_df = pandas.read_csv(options.s+'/data/cell_info.csv', sep='\t', index_col=0)
    cell_type = cellinfo_df.ix[matrix_df.index.values, 'notes']
    cTypes = list(set(cell_type))
    cTypes.sort()
    cTypeIndex = [numpy.where(cell_type==x) for x in cTypes]
#
    npc = min(int(options.npc), len(matrix_df.values[:,0]), len(matrix_df.values[0,:]))
    pca = PCA(n_components=npc, svd_solver='full')
    pca_result = pca.fit_transform(matrix_df)
    seaborn.set(font_scale=2)
    fig1, axes = plt.subplots(1, figsize=(15,15))
    for ict,ct in enumerate(cTypes):
        axes.scatter(pca_result[cTypeIndex[ict], 0], pca_result[cTypeIndex[ict], 1], c=colors[ict], label=ct, s=50)
    axes.legend(bbox_to_anchor=(1, 1), loc=2, borderaxespad=0.)
    fig1.savefig(options.s+'/figure/'+outfig1, bbox_inches='tight')
#
    tsne_result = TSNE(n_components=2, learning_rate=int(options.tsneLR), random_state=int(options.tsneRS),
                       n_iter=int(options.tsneIter)).fit_transform(pca_result)
    fig2, axes = plt.subplots(1, figsize=(15,15))
    for ict,ct in enumerate(cTypes):
        axes.scatter(tsne_result[cTypeIndex[ict], 0], tsne_result[cTypeIndex[ict], 1], c=colors[ict], label=ct, s=50)
    axes.legend(bbox_to_anchor=(1, 1), loc=2, borderaxespad=0.)
    fig2.savefig(options.s+'/figure/'+outfig2, bbox_inches='tight')
    return pca_result, tsne_result
#
#
#
#
def plot_knn_cluster(options, matrix_df, n_clust, tsne_result, outfig, outcsv):
    connectivity = kneighbors_graph(matrix_df, n_neighbors=20, include_self=False)
    connectivity = 0.5*(connectivity + connectivity.T)
    ward_linkage = cluster.AgglomerativeClustering(n_clusters=n_clust, linkage='ward', connectivity=connectivity)
    ward_linkage.fit(matrix_df)
    y_predict = ward_linkage.labels_.astype(numpy.int)
    seaborn.set(font_scale=2)
    fig3, axes = plt.subplots(1, figsize=(15,15))
    for i in range(0,n_clust):
        index = numpy.where(y_predict==i)[0]
        axes.scatter(tsne_result[index,0], tsne_result[index,1], color=colors[i], label=str(i), s=50)
    axes.legend(bbox_to_anchor=(1, 1), loc=2, borderaxespad=0.)
    fig3.savefig(options.s+'/figure/'+outfig, bbox_inches='tight')
    predict_df = pandas.DataFrame(y_predict, index=matrix_df.index.values, columns=['cluster'])
    predict_df.to_csv(options.s+'/result/'+outcsv, sep='\t')
    return
#
#
#
#
def plot_cluster(options, cluster_df, n_clust, tsne_result, outfig):
    seaborn.set(font_scale=2)
    fig3, axes = plt.subplots(1, figsize=(15,15))
    for i in range(0,n_clust):
        index = numpy.where(cluster_df.values==i)[0]
        axes.scatter(tsne_result[index,0], tsne_result[index,1], color=colors[i], label=str(i), s=50)
    axes.legend(bbox_to_anchor=(1, 1), loc=2, borderaxespad=0.)
    fig3.savefig(options.s+'/figure/'+outfig, bbox_inches='tight')
    return
#
#
#
#
def heatmap_compare(options, matrix_df, clust_df, outfig0):
    cell_type = clust_df.ix[matrix_df.index.values, 'cluster']
    cTypes = list(set(cell_type))
    cTypes.sort()
    cTypeIndex = [numpy.where(cell_type==x) for x in cTypes]
#
    lut = dict(zip(cTypes, colors[:len(cTypes)]))
    row_colors = cell_type.map(lut)
    cmap = seaborn.diverging_palette(h_neg=210, h_pos=350, s=90, l=30, as_cmap=True)
    seaborn.set(font_scale=1.2)
    fig0 = seaborn.clustermap(matrix_df, method='ward', metric='euclidean', cmap=cmap,
           row_colors=row_colors, figsize=(20,20))
    plt.setp(fig0.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    plt.setp(fig0.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)
    fig0.savefig(options.s + '/figure/' + outfig0, bbox_inches='tight')
    return
#
#
#
#
def specific_accesson(options, subname):
    cluster_df = pandas.read_csv(options.cfile, sep='\t', index_col=0)
    kCluster = options.cluster.split(',')
    if options.vs!='all': vsCluster = options.vs.split(',')
    if 'cluster' not in cluster_df.columns.values:
        cluster_df['cluster'] = cluster_df['notes']
    else:
        kCluster = map(int, kCluster)
        if options.vs!='all': vsCluster = map(int, vsCluster)
    if options.vs=='all': vsCluster = list(set(cluster_df['cluster'].values)-set(kCluster))
    accesson_df = pandas.read_csv(options.s + '/matrix/Accesson_reads.csv', sep=',', index_col=0,
                   engine='c', na_filter=False, low_memory=False)
    normal = numpy.array([x/x.sum()*1000000 for x in accesson_df.values])
    accesson_df = pandas.DataFrame(normal, index=accesson_df.index, columns=accesson_df.columns)
    cluster_df = cluster_df.loc[accesson_df.index.values]
    peaks_df = pandas.read_csv(options.s + '/matrix/Accesson_peaks.csv', sep='\t', index_col=0,
                   engine='c', na_filter=False, low_memory=False)
#    peaks_bed = numpy.array(open(options.s+'/peak/annotate_peak.bed').readlines())
    peaks_bed = numpy.array(open(options.s+'/peak/top_peaks.bed').readlines())
    cell_inCluster = cluster_df.loc[cluster_df['cluster'].isin(kCluster)].index.values
    cell_outCluster = cluster_df.loc[cluster_df['cluster'].isin(vsCluster)].index.values
    print len(cell_inCluster), len(cell_outCluster)
    read_inCluster = accesson_df.loc[cell_inCluster].values
    read_outCluster = accesson_df.loc[cell_outCluster].values
    mean_in, mean_out = read_inCluster.mean(axis=0), read_outCluster.mean(axis=0)
    ttest, pvalues = scipy.stats.ttest_ind(read_inCluster, read_outCluster, equal_var=False)
    fold = (mean_in + 1e-4) / (mean_out + 1e-4)
    matrix = numpy.vstack((mean_in, mean_out, fold, pvalues)).T
    columns = ['mean_inCluster','mean_outCluster', 'fold', 'p-value']
    matrix_df = pandas.DataFrame(matrix, index=accesson_df.columns, columns=columns)
    matrix_df = matrix_df.loc[matrix_df['p-value'] <= float(options.pvalue)]
    matrix_df = matrix_df.loc[matrix_df['fold'] >= float(options.fold)]
    matrix_df = matrix_df.sort_values(by=['p-value'])
    for acc in matrix_df.index.values:
        peaks = peaks_df.loc[peaks_df['group']==int(acc)].index.values
        matrix_df.loc[acc, 'N_peaks'] = len(peaks)
        peaks_index = [int(x[4:]) for x in peaks]
        if len(peaks)>=5:
            with open(options.s+'/result/accesson_'+str(acc)+'_peaks.bed', 'w') as output1:
                for line in peaks_bed[peaks_index]:
                    words = line.split()
                    print >> output1, '\t'.join(words)
    matrix_df.to_csv(options.s+'/result/Accessons_of_Cluster_'+subname+'.csv', sep='\t')
    return
#
#
def specific_peak(options, subname):
    cluster_df = pandas.read_csv(options.cfile, sep='\t', index_col=0)
    kCluster = options.cluster.split(',')
    if options.vs!='all': vsCluster = options.vs.split(',')
    if 'cluster' not in cluster_df.columns.values:
        cluster_df['cluster'] = cluster_df['notes']
    else:
        kCluster = map(int, kCluster)
        if options.vs!='all': vsCluster = map(int, vsCluster)
    if options.vs=='all': vsCluster = list(set(cluster_df['cluster'].values)-set(kCluster))
    reads_df = pandas.read_csv(options.s+'/matrix/filtered_reads.csv', sep=',', index_col=0,
                   engine='c', na_filter=False, low_memory=False)
#    peaks_bed = open(options.s+'/peak/annotate_peak.bed').readlines()
    peaks_bed = open(options.s+'/peak/top_peaks.bed').readlines()
    cluster_df = cluster_df.loc[reads_df.index.values]
    cell_inCluster = cluster_df.loc[cluster_df['cluster'].isin(kCluster)].index.values
    cell_outCluster = cluster_df.loc[cluster_df['cluster'].isin(vsCluster)].index.values
    print len(cell_inCluster), len(cell_outCluster)
#
    read_inCluster = reads_df.loc[cell_inCluster].values
    read_outCluster = reads_df.loc[cell_outCluster].values
    mean_in, mean_out = read_inCluster.mean(axis=0), read_outCluster.mean(axis=0)
    ttest, pvalues = scipy.stats.ttest_ind(read_inCluster, read_outCluster, equal_var=False)
    fold = (mean_in + 1e-4) / (mean_out + 1e-4)
    matrix = numpy.vstack((mean_in, mean_out, fold, pvalues)).T
    columns = ['mean_inCluster','mean_outCluster', 'fold', 'p-value']
    matrix_df = pandas.DataFrame(matrix, index=reads_df.columns, columns=columns)
    matrix_df = matrix_df.loc[matrix_df['p-value'] <= float(options.pvalue)]
    matrix_df = matrix_df.loc[matrix_df['fold'] >= float(options.fold)]
    matrix_df = matrix_df.sort_values(by=['p-value'])
#
    peaks_specific = [peaks_bed[int(peak[4:])] for peak in matrix_df.index.values]
    peaks_position = [x.split()[0]+':'+x.split()[1]+'-'+x.split()[2] for x in peaks_specific]
    matrix_df['position'] = peaks_position
    matrix_df.to_csv(options.s+'/result/Peaks_of_Cluster_'+subname+'.csv', sep='\t')
    with open(options.s+'/result/Peaks_of_Cluster_'+subname+'.bed', 'w') as output:
        for line in peaks_specific:
            center = (int(line.split()[1]) + int(line.split()[2])) // 2
            print >> output, line[:-1]
    return
#
#
#
#
#
#
