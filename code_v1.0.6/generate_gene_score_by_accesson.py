import os,numpy,pandas,time,numba,scipy.stats
from optparse import OptionParser
#
#
opts = OptionParser()
usage = "Generate gene score from accesson\nusage: %prog -s project --genome hg19"
opts = OptionParser(usage=usage, version="%prog 1.0")
opts.add_option("-s", help="The project folder.")
opts.add_option("--genome", default='hg19', help="Genome reference for gene annotation, can be hg19 or mm10, default=hg19")
opts.add_option("--width", default=1000000, help="Width of genome window for fisher exact test, default=1000000")
opts.add_option("--pvalue", default=0.01, help='P-value threshold for fisher exact test, default=0.01')
options, arguments = opts.parse_args()
#
#
def annotate_peak(project, genome='hg19'):
    peaks_bed = project + '/peak/top_filtered_peaks.bed'
    annotated = project + '/peak/top_annotated_peaks.bed'
    os.popen('annotatePeaks.pl '+peaks_bed+'  '+genome+' -strand both -size given > '+annotated)
    return
#
#
def annotate_accessons(project, width=1e6, pvalue=0.01):
    peak_bed = project + '/peak/top_filtered_peaks.bed'
    accesson_csv = project + '/matrix/Accesson_peaks.csv'
    annotated_bed = project + '/peak/top_annotated_peaks.bed'
    peak_list = [[x.split()[0], (int(x.split()[1])+int(x.split()[2]))//2] for x in open(peak_bed).readlines()]
    peak_index = ['peak'+str(x) for x in range(0, len(peak_list))]
    peak_df = pandas.DataFrame(peak_list, index=peak_index, columns=['chromosome', 'base'])
    annotated_df = pandas.read_csv(annotated_bed, sep='\t', index_col=0,
                                   engine='c', na_filter=False, low_memory=False)
    annotated_df.index = ['peak'+x.split('-')[-1] for x in annotated_df.index.values]
    accesson_df = pandas.read_csv(accesson_csv, sep='\t', index_col=0,
                                  engine='c', na_filter=False, low_memory=False)
    accessons = list(set(accesson_df['group'].values))
    accessons.sort()
    accesson_annotate = pandas.DataFrame(columns=['genes', '-log10(P-value)'])
    gene_annotate = pandas.DataFrame(columns=['accessons', '-log10(P-value)'])
    for acc in accessons:
        accPeak_df = peak_df.loc[accesson_df.loc[accesson_df['group']==acc].index]
        genes = {}
        for peak in accPeak_df.index.values:
            chrom, base = accPeak_df.loc[peak, 'chromosome'], accPeak_df.loc[peak, 'base']
            base_up, base_down = base - width, base + width
            sameChrom_in_acc = accPeak_df.loc[accPeak_df['chromosome']==chrom]
            sameChrom_overall = peak_df.loc[peak_df['chromosome']==chrom]
            sameRegion_in_acc = numpy.where(abs(sameChrom_in_acc['base']-base)<=width)[0]
            sameRegion_overall = numpy.where(abs(sameChrom_overall['base']-base)<=width)[0]
            matrix = numpy.array([[len(sameRegion_in_acc), len(accPeak_df)-len(sameRegion_in_acc)],
                                  [len(sameRegion_overall), len(peak_df)-len(sameRegion_overall)]])
            odd, p_value = scipy.stats.fisher_exact(matrix)
            if p_value<=pvalue:
                gene_symbol = annotated_df.loc[peak, 'Gene Name']
                log10_P = -numpy.log10(p_value)
                if gene_symbol not in genes.keys():
                    genes[gene_symbol] = log10_P
                elif genes[gene_symbol] < log10_P:
                    genes[gene_symbol] = log10_P
        for gene in genes.keys():
            if gene not in gene_annotate.index.values:
                gene_annotate.loc[gene, 'accessons'] = str(acc)
                gene_annotate.loc[gene, '-log10(P-value)'] = str(genes[gene])
            else:
                gene_annotate.loc[gene, 'accessons'] += ';' + str(acc)
                gene_annotate.loc[gene, '-log10(P-value)'] += ';' + str(genes[gene])
        accesson_annotate.loc[acc, 'genes'] = ';'.join(genes.keys())
        accesson_annotate.loc[acc, '-log10(P-value)'] = ';'.join(list(map(str, genes.values())))
    accesson_annotate.to_csv(project+'/matrix/Accesson_annotated.csv', sep='\t')
    gene_annotate.to_csv(project+'/matrix/gene_annotated.csv', sep='\t')
    return
#
#
def get_gene_score(project):
    gene_annotate = pandas.read_csv(project+'/matrix/gene_annotated.csv', sep='\t', index_col=0)
    accesson_matrix = pandas.read_csv(project+'/matrix/Accesson_reads.csv', sep=',', index_col=0,
                                      engine='c', na_filter=False, low_memory=False)
    genes, matrix = [], []
    for gene in gene_annotate.index.values:
        accessons = gene_annotate.loc[gene, 'accessons'].split(';')
        weight = list(map(float, gene_annotate.loc[gene, '-log10(P-value)'].split(';')))
        if len(accessons)>0:
            sub_matrix = accesson_matrix[accessons].values
            expression = numpy.average(sub_matrix, axis=1, weights=weight).T
            matrix.append(expression)
            genes.append(gene)
    matrix = numpy.array(matrix).T
    expression_df = pandas.DataFrame(matrix, index=accesson_matrix.index, columns=genes)
    expression_df.to_csv(project+'/matrix/gene_score.csv', sep=',')
    return
#
#
t1 = time.time()
annotate_peak(options.s, genome=options.genome)
annotate_accessons(options.s, int(options.width), float(options.pvalue))
get_gene_score(options.s)
print(time.time()-t1)
#
#
