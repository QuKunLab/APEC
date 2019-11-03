from APEC import clustering, plot, generate
#
#
#### Cluster cells, plot tSNE map, and calculate ARI; cost ~5 minutes.
#### The output file project02/figure/TSNE_by_APEC_with_notes_label.pdf is Fig. S3b of the APEC paper.
#
clustering.build_accesson('project03/', ngroup=700)
clustering.cluster_byAccesson('project03/', norm='zscore')
plot.plot_tsne('project03/', wt=1)
clustering.cluster_comparison('project03/matrix/filtered_cells.csv',
                              'project03/result/cluster_by_APEC.csv')
#
