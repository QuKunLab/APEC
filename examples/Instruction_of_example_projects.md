Please visit the following website to download example projects:

http://galaxy.ustc.edu.cn:30803/APEC/

(1) **project01** is an example of running APEC from the raw fastq files. Users only need to copy the raw fastq files in **data** folder to start the analysis. This project contains single cell samples of human's leukemic stem cells, leukemic blast cells, LMPP, monocytes and HL60 cells, from "Schep, A.N., Wu, B., Buenrostro, J.D. & Greenleaf, W.J. chromVAR: inferring transcription-factor-associated accessibility from single-cell epigenomic data. Nat Methods 14, 975-978 (2017)".

How to run:

    #### get fragment count matrix from raw fastq files, take 10~20 hours on 8-core/32GB computer ####
    
    bash APEC_prepare_steps.sh -s $project01 -g hg19 -n 10 -l 8 -p 0.05 -f 800
    
    #### cluster cells by APEC algorithm, take <5 minutes on one CPU-core of computer ####
    
    python cluster_byAccesson.py -s $project01
    
    #### cluster cells by chromVAR algorithm, take ~30 minutes on 8-core computer ####
    
    python cluster_byMotif.py -s $project01 --np 8
    
    #### get differential peaks/genes/motifs for cell cluster 1, take several minutes ####
    
    python generate_differential_markers.py -s $project01 --cfile $project01/result/KNN_cluster_by_Accesson.csv \
                                            --cluster 1 --vs all --motif yes --gene yes
    
    #### plot enrichment of motif RUNX1 on tSNE diagram, take several minutes ####
    
    python generate_markers_on_plots.py -s $project01 --cfile $project01/result/TSNE_by_Accesson.csv \
                                        --type motif --name RUNX1
    
    #### search for potential super enhancer, take several minutes ####
    
    python generate_superEnhancer.py -s $project01


(2) **project02** is an example of running APEC from the user's own fragment count matrix. Users need to copy "cell_info.csv" in **data** folder, "filtered_reads.csv" in **matrix** folder, and "filtered_reads.csv" in **peak** folder to start the clustering. This project contains single cell samples from the forebrain of adult mice, from "Preissl, S. et al. Single-nucleus analysis of accessible chromatin in developing mouse forebrain reveals cell-type-specific transcriptional regulation. Nat Neurosci 21, 432-439 (2018)".

How to run:

    #### prepare dataset from user's own fragment count matrix，take 1~2 hours on on 8-core/32GB computer ####

    python prepare_premappedMatrix.py -s $project02 --ref mm10 --fa ../reference/mm10_chr.fa --np 8

    #### cluster cells by APEC algorithm, take ~40 minutes on one CPU-core of computer ####
    
    python cluster_byAccesson.py -s $project02 --nc 10 --hc no
    
    #### cluster cells by chromVAR algorithm, take ~8 hours on 8-core computer ####
    
    python cluster_byMotif.py -s $project02 --np 8 --nc 10 --hc no

