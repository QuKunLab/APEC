from APEC import clustering, plot, generate
from optparse import OptionParser
opts = OptionParser()
usage = "Clustering by accesson\nusage: %prog -p project03"
opts = OptionParser(usage=usage, version="%prog 1.0")
opts.add_option("-p", help="The project03 folder.")
options, arguments = opts.parse_args()
#
#
#### Cluster cells, plot tSNE map, and calculate ARI; cost ~5 minutes.
#### The output file project02/figure/TSNE_by_APEC_with_notes_label.pdf is Fig. S3b of the APEC paper.
#
clustering.build_accesson(options.p, ngroup=700)
clustering.cluster_byAccesson(options.p, norm='zscore')
plot.plot_tsne(options.p, wt=1)
clustering.cluster_comparison(options.p+'/matrix/filtered_cells.csv',
                              options.p+'/result/cluster_by_APEC.csv')
#
