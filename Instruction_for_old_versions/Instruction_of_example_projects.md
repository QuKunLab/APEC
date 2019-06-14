Please visit the following website to download example projects:

http://galaxy.ustc.edu.cn:30803/APEC/

(1) **project01** is an example of running APEC from the raw fastq files. Users only need to copy the raw fastq files in **data** folder to start the analysis. This project contains single cell samples of human's leukemic stem cells, leukemic blast cells, LMPP, monocytes and HL60 cells, from "Schep, A.N., Wu, B., Buenrostro, J.D. & Greenleaf, W.J. chromVAR: inferring transcription-factor-associated accessibility from single-cell epigenomic data. Nat Methods 14, 975-978 (2017)".

How to run:

    #### get fragment count matrix from raw fastq files, take 10~20 hours on 8-core/32GB computer ####

    bash APEC_prepare_steps.sh -s $project01 -g hg19 -n 10 -l 8 -p 0.05 -f 800

    #### cluster cells by APEC algorithm, take <5 minutes on one CPU-core of computer ####

    python cluster_byAccesson.py -s $project01 --ngroup 720

    #### cluster cells by chromVAR algorithm, take ~30 minutes on 8-core computer ####

    python cluster_byMotif.py -s $project01 --np 8

    #### get differential peaks/genes/motifs for cell cluster 1, take several minutes ####

    python generate_gene_score_by_accesson.py -s $project01 --genome hg19
    python generate_differential_markers.py -s $project01 --cfile $project01/result/louvain_cluster_by_Accesson.csv \
                                            --cluster 1 --vs all --motif yes --gene yes

    #### plot enrichment of motif RUNX1 on tSNE diagram, take several minutes ####

    python generate_markers_on_plots.py -s $project01 --cfile $project01/result/TSNE_by_Accesson.csv \
                                        --type motif --name RUNX1

    #### search for potential super enhancer, take several minutes ####

    python generate_superEnhancer.py -s $project01


(2) **project02** is an example of running APEC from the user's own fragment count matrix. Users need to copy "filtered_cells.csv" and "filtered_reads.mtx" in **matrix** folder, and "top_filtered_peaks.bed" in **peak** folder to start the clustering. This project contains single cell samples from the forebrain of adult mice, from "Preissl, S. et al. Single-nucleus analysis of accessible chromatin in developing mouse forebrain reveals cell-type-specific transcriptional regulation. Nat Neurosci 21, 432-439 (2018)".

How to run:

    #### prepare dataset from user's own fragment count matrix，take 1~2 hours on on 8-core/32GB computer ####

    python prepare_premappedMatrix.py -s $project02 --ref mm10 --fa ../reference/mm10_chr.fa --np 8

    #### cluster cells by APEC algorithm, take ~40 minutes on one CPU-core of computer ####

    python cluster_byAccesson.py -s $project02 --ngroup 700 --hc yes

    #### cluster cells by chromVAR algorithm, take ~8 hours on 8-core computer ####

    python cluster_byMotif.py -s $project02 --np 8 --hc no

    #### plot motifs on tSNE diagram of APEC. It takes ~1 minute for each motif ####

    python generate_markers_on_plots.py -s ../results/project02 --cfile ../results/project02/result/TSNE_by_Accesson.csv \
                                        --type motif --name NEUROD1 --sharp -5,10
    python generate_markers_on_plots.py -s ../results/project02 --cfile ../results/project02/result/TSNE_by_Accesson.csv \
                                        --type motif --name OLIG2 --sharp -5,8
    python generate_markers_on_plots.py -s ../results/project02 --cfile ../results/project02/result/TSNE_by_Accesson.csv \
                                        --type motif --name MEF2C --sharp -5,8
    python generate_markers_on_plots.py -s ../results/project02 --cfile ../results/project02/result/TSNE_by_Accesson.csv \
                                        --type motif --name MEIS2 --sharp -4,6
    python generate_markers_on_plots.py -s ../results/project02 --cfile ../results/project02/result/TSNE_by_Accesson.csv \
                                        --type motif --name Dlx2 --sharp -5,11
    python generate_markers_on_plots.py -s ../results/project02 --cfile ../results/project02/result/TSNE_by_Accesson.csv \
                                        --type motif --name NOTO --sharp -6,10
    python generate_markers_on_plots.py -s ../results/project02 --cfile ../results/project02/result/TSNE_by_Accesson.csv \
                                        --type motif --name Sox2 --sharp -22,27
    python generate_markers_on_plots.py -s ../results/project02 --cfile ../results/project02/result/TSNE_by_Accesson.csv \
                                        --type motif --name ETS1 --sharp -4,8


(3) **project03** contains single cell samples of hematopoietic stem cell differentiation (including HSC, MPP, CMP, GMP, MEP, LMPP, CLP, and pDC cells), from "Buenrostro, J.D. et al. Integrated Single-Cell Analysis Maps the Continuous Regulatory Landscape of Human Hematopoietic Differentiation. Cell 173, 1535-1548 e1516 (2018)".

How to run:

    #### 1. Cell clustering by APEC algorithm. It takes ~10 minutes.

    python cluster_byAccesson.py -s $project03 --norm probability

    #### 2. Generate pseudo-time trajectory by monocle. It takes ~4 minutes.

    python generate_trajectory.py -s $project03 --cfile $project03/matrix/filtered_cells.csv --npc 5 --angle 39,94

    #### 3 Plot genes on tSNE. It takes several minutes.

    python generate_gene_score_by_accesson.py -s $project01 --genome hg19
    python generate_markers_on_plots.py -s $project03 --cfile $project03/result/TSNE_by_APEC.csv \
                                        --type gene --name FOXO1

    #### 4. Plot motifs on trajectory. It takes ~1 minute for each motif.

    python generate_markers_on_plots.py -s $project03 --cfile $project03/result/monocle_reduced_dimension.csv \
                                        --type motif --name Hoxa9 --angle 39,94 --sharp -6,9
    python generate_markers_on_plots.py -s $project03 --cfile $project03/result/monocle_reduced_dimension.csv \
                                        --type motif --name GATA1 --angle 39,94 --sharp -10,10
    python generate_markers_on_plots.py -s $project03 --cfile $project03/result/monocle_reduced_dimension.csv \
                                        --type motif --name CEBPB --angle 39,94 --sharp -10,10
    python generate_markers_on_plots.py -s $project03 --cfile $project03/result/monocle_reduced_dimension.csv \
                                        --type motif --name TCF4 --angle 39,94 --sharp -5,10

**NOTE**: Monocle is very sensitive to the input matrix. For example, when we use different versions of scipy (such as 1.0.0 and 1.2.1), the PCA matrices generated by sklearn will be slightly different (difference < 1.0e-7), and even this small difference will produce completely different shapes of trajectories by monocle. Different versions of many other related libraries/packages in python and R environment will also result in different shapes of trajectories. However, APEC can still accurately profile the developmental and differentiation pathway, even if the shapes of trajectories are different. Here we recommend that users install scipy=1.0.0 to generate the trajectory more similar to the Figure in our paper.
