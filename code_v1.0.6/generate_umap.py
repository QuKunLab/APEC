import sys,getopt,numpy,pandas,umap,scipy.stats
from sklearn.decomposition import PCA
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#
#
def get_parameters(argv):
    project, cellinfo, rand_stat, norm_method = '', '', 0, 'zscore'
    help_info = ['Plot UMAP figure \n',
                 'python generate_umap.py -s <project path> \n',
                 '\t\t\t -c <cluster.csv or filtered_cell.csv file> \n',
                 '\t\t\t -r <seed of random_state for UMAP, default=0 \n',
                 '\t\t\t -n <normalization method, zscore or probability, default=zscore> \n']
    help_info = ''.join(help_info)
    try:
        opts, args = getopt.getopt(argv,"hs:c:r:n:")
    except getopt.GetoptError:
        print('Incorrect input parameters!')
        print(help_info)
        sys.exit(2)
    for opt,arg in opts:
        if opt=='-h':
            print(help_info)
            sys.exit()
        elif opt=='-s':
            project = arg
        elif opt=='-c':
            cellinfo = arg
        elif opt=='-r':
            rand_stat = int(arg)
        elif opt=='-n':
            norm_method = arg
    return project, cellinfo, rand_stat, norm_method
#
#
def run_umap(project, cellinfo, rand_stat=0, norm_method='zscore'):
    mat_df = pandas.read_csv(project+'/matrix/Accesson_reads.csv', sep=',', index_col=0,
                             engine='c', na_filter=False, low_memory=False)
    if norm_method=='zscore':
        matrix = scipy.stats.zscore(mat_df.values, axis=1)
    elif norm_method=='probability':
        matrix = numpy.array([x*10000.0/x.sum() for x in mat_df.values])
    else:
        print('-n should be zscore or probability')
        sys.exit()
    umap_result = umap.UMAP(n_components=2, random_state=rand_stat).fit_transform(matrix)
    cellinfo_df = pandas.read_csv(cellinfo, sep='\t', index_col=0, engine='c', na_filter=False, low_memory=False)
    umap_df = pandas.DataFrame(umap_result, index=cellinfo_df.index, columns=['UMAP1', 'UMAP2'])
    umap_df.to_csv(project+'/result/umap.csv', sep='\t')
#
    if 'notes' in cellinfo_df.columns.values: cellinfo_df['cluster'] = cellinfo_df['notes']
    cTypes = list(set(cellinfo_df['cluster'].values))
    cTypes.sort()
    cTypeIndex = [numpy.where(cellinfo_df['cluster'].values==x) for x in cTypes]
    colors = matplotlib.colors.CSS4_COLORS.keys()
    fig2, axes = plt.subplots(1, figsize=(15,15))
    for ict,ct in enumerate(cTypes):
        axes.scatter(umap_result[cTypeIndex[ict], 0], umap_result[cTypeIndex[ict], 1], c=colors[ict], label=ct, s=50)
    axes.legend(bbox_to_anchor=(1, 1), loc=2, borderaxespad=0.)
    fig2.savefig(project+'/figure/umap.pdf', bbox_inches='tight')
    return
#
#
if __name__=='__main__':
    project, cellinfo, rand_stat, norm_method = get_parameters(sys.argv[1:])
    run_umap(project, cellinfo, rand_stat, norm_method)
#
#
#
