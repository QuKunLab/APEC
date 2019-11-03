from APEC import clustering, plot, generate
#
#
#### Cluster cells, plot tSNE map and cell-cell correlation; cost ~90 minutes.
#### The output file project02/figure/TSNE_by_APEC_with_cluster_label.pdf is Fig. 2a of the APEC paper.
#### The output file project02/figure/cell_cell_correlation_by_APEC_with_cluster_label.png is Fig. 2b of the APEC paper
#
clustering.build_accesson('project02/', ngroup=600)
clustering.cluster_byAccesson('project02/', norm='probability')
plot.plot_tsne('project02/', rs=3)
plot.correlation('project02/', cell_label='cluster', clip=[0,1])
#
#
#### Generate gene scores from ATAC signals; cost ~30 minutes.
#
generate.gene_score('project02/', genome_gtf='../reference/mm10_RefSeq_genes.gtf')
#
#
#### Generate differential accessons between sub-clusters of excitatory neurons and inhibitory neurons;
#### cost ~1 minute.
#
log2Fold = 0.848
generate.get_nearby_genes('project02/')
generate.differential_feature('project02/', feature='accesson', target='0', vs='1,10,11,12', log2_fold=log2Fold)
generate.differential_feature('project02/', feature='accesson', target='1', vs='0,10,11,12', log2_fold=log2Fold)
generate.differential_feature('project02/', feature='accesson', target='10', vs='0,1,11,12', log2_fold=log2Fold)
generate.differential_feature('project02/', feature='accesson', target='11', vs='0,1,10,12', log2_fold=log2Fold)
generate.differential_feature('project02/', feature='accesson', target='12', vs='0,1,10,11', log2_fold=log2Fold)
generate.differential_feature('project02/', feature='accesson', target='2', vs='3,5,6,8', log2_fold=0.678)
generate.differential_feature('project02/', feature='accesson', target='3', vs='2,5,6,8', log2_fold=log2Fold)
generate.differential_feature('project02/', feature='accesson', target='5', vs='2,3,6,8', log2_fold=log2Fold)
generate.differential_feature('project02/', feature='accesson', target='6', vs='2,3,5,8', log2_fold=log2Fold)
generate.differential_feature('project02/', feature='accesson', target='8', vs='2,3,5,6', log2_fold=log2Fold)
#
#
#### Caculate motif enrichment; cost ~4 hours on 8-core CPU.
#### Note: ../reference/ should be the reference folder built by the user.
#
generate.motif_matrix('project02/', genome_fa='../reference/mm10_chr.fa',
                      background='../reference/tier1_markov1.norc.txt',
                      meme='../reference/JASPAR2018_CORE_vertebrates_redundant_pfms_meme.txt',
                      np=8)
clustering.cluster_byMotif('project02/', np=8)
#
#
#### Plot motifs on the tSNE map; cost ~1 minutes.
#### The output files project02/figure/motif_XXX_on_tsne_by_APEC.pdf are Fig. 2e of the APEC paper.
#
plot.plot_feature('project02/', space='tsne', feature='motif', name='NEUROD1', clip=[-5,10])
plot.plot_feature('project02/', space='tsne', feature='motif', name='OLIG2', clip=[-5,8])
plot.plot_feature('project02/', space='tsne', feature='motif', name='MEF2C', clip=[-5,8])
plot.plot_feature('project02/', space='tsne', feature='motif', name='MEIS2', clip=[-4,6])
plot.plot_feature('project02/', space='tsne', feature='motif', name='Dlx2', clip=[-5,11])
plot.plot_feature('project02/', space='tsne', feature='motif', name='NOTO', clip=[-6,10])
plot.plot_feature('project02/', space='tsne', feature='motif', name='Sox2', clip=[-22,27])
plot.plot_feature('project02/', space='tsne', feature='motif', name='ETS1', clip=[-4,8])
#
#