from APEC import clustering, plot, generate
from optparse import OptionParser
opts = OptionParser()
usage = "Clustering by accesson\nusage: %prog -p project01 -r reference"
opts = OptionParser(usage=usage, version="%prog 1.0")
opts.add_option("-p", help="The project01 folder.")
opts.add_option("-r", help="The reference folder")
options, arguments = opts.parse_args()
#
#
#### Cluster cells, plot tSNE map, and calculate ARI; cost ~10 minutes.
#### The output file project01/figure/TSNE_by_APEC_with_notes_label.pdf is Fig. 1b of the APEC paper.
#
clustering.build_accesson(options.p, ngroup=600)
clustering.cluster_byAccesson(options.p, norm='probability')
plot.plot_tsne(options.p)
clustering.cluster_comparison(options.p+'/matrix/filtered_cells.csv',
                              options.p+'/result/cluster_by_APEC.csv',
                              exclude='UNK')
#
#
#### Caculate motif enrichment and plot motifs on tSNE map; cost <1 hour on 8-core CPU.
#### The output files project01/figure/motif_XXX_on_tsne_by_APEC.pdf are Fig. 1d of the APEC paper.
#
# generate.motif_matrix(options.p, genome_fa=options.r+'/hg19_chr.fa',
#                       background=options.r+'/tier1_markov1.norc.txt',
#                       meme=options.r+'/JASPAR2018_CORE_vertebrates_redundant_pfms_meme.txt',
#                       np=8)
# clustering.cluster_byMotif(options.p, np=8)
#
plot.plot_feature(options.p, space='tsne', feature='motif', name='Erg', clip=[0,7])
plot.plot_feature(options.p, space='tsne', feature='motif', name='CTCF', clip=[0,5])
plot.plot_feature(options.p, space='tsne', feature='motif', name='GATA1', clip=[0,25])
#
#
#### Construct trajectory and plot motifs on trajectory; cost ~1 minute.
#### The output file project01/figure/pseudotime_trajectory_with_notes_label.pdf is Fig. 3c of the APEC paper.
#### The output files project01/figure/motif_XXX_on_trajectory_by_APEC.pdf are Fig. 3d of the APEC paper.
#
generate.monocle_trajectory(options.p, around=[6, 25])
plot.plot_trajectory(options.p, angles=[60,90])
#
plot.plot_feature(options.p, space='trajectory', feature='motif', name='Hoxa9', angles=[60,90], clip=[-6,9])
plot.plot_feature(options.p, space='trajectory', feature='motif', name='GATA1', angles=[60,90], clip=[-10,10])
plot.plot_feature(options.p, space='trajectory', feature='motif', name='CEBPB', angles=[60,90], clip=[-10,10])
plot.plot_feature(options.p, space='trajectory', feature='motif', name='TCF4', angles=[60,90], clip=[-5,10])
#
