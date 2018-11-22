# User Guide for APEC

(Accessibility Pattern based Epigenomic Clustering)

## Introduction

APEC can perform fine cell type clustering on single cell chromatin accessibility data from scATAC-seq, snATAC-seq, sciATAC-seq or any other relevant experiment. It can also be used to control data quality, map fragment count matrices, search for important differential motifs/genes for each cell cluster, find super enhancers, and construct pseudo-time trajectory (by calling Monocle).

**If users want to process the raw fastq data from scATAC-seq experiment, please run APEC from section 2 “Fragment count matrix”. If users have their own fragment count matrix, where each element is the number of fragments per-cell-per-peak, please run APEC from section 3 “Clustering”.**

## 1. Requirements and installation

### 1.1	Requirements

APEC requires users to use Linux system, as well as Python (version 2.7.5+) and R (version 3.4+) environment. Users also need to install the following packages:

(1) Python packages and libraries: 

    numpy, scipy, pandas, sklearn, multiprocessing, numba, pysam,
    matplotlib, seaborn, networkx, python-louvain, python-Levenshtein
    
    all upon python packages can be installed by: 
    pip install package_name

(2) R packages and libraries: 

    Monocle: http://cole-trapnell-lab.github.io/monocle-release/

(3) Other necessary software:

    Bowtie2: https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.9/
    Samtools: https://github.com/samtools/samtools
    Picard: https://github.com/broadinstitute/picard/releases/tag/2.18.14
    Bedtools: http://bedtools.readthedocs.io/en/latest/content/installation.html
    HOMER: http://homer.ucsd.edu/homer/download.html
    Macs2: https://github.com/taoliu/MACS.git
    Meme: http://meme-suite.org/doc/download.html?man_type=web
  
### 1.2	Installation
Users simply completes the APEC installation by copying the APEC folder to any path on the computer (i.e. $APEC). There are two subfolders in APEC: a **codes** folder, which contains all APEC programs for data processing; a **reference** folder, which contains all necessary index and reference files for the hg19 and mm10 genomes. So users can run APEC program directly in $APEC/codes/, or put this path in the system environment to use it elsewhere. The **reference** folder is required for APEC and should be placed in the same path with the **codes** folder. It contains the following files:
 
    hg19_refseq_genes_TSS.txt, hg19_RefSeq_genes.gtf, hg19_blacklist.JDB.bed,
    hg19_chr.fa, hg19_chr.fa.fai, hg19.chrom.sizes,
    hg19.1.bt2, hg19.2.bt2, hg19.3.bt2, hg19.4.bt2,
    mm10_refseq_genes_TSS.txt, mm10_RefSeq_genes.gtf, mm10_blacklist.BIN.bed,
    mm10_chr.fa, mm10_chr.fa.fai, mm10.chrom.sizes,
    mm10.1.bt2, mm10.2.bt2, mm10.3.bt2, mm10.4.bt2,
    JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt, tier1_markov1.norc.txt

## 2. Fragment count matrix

### 2.1	Arrangement of raw data

Users need to build a source folder (i.e. $source), which contains a **data** folder, then copy all raw sequencing fastq files into the $source/**data**/ folder. All these pair-end fastq files should be named as:
 
    type1-001_1.fastq, type1-001_2.fastq, type1-002_1.fastq, type1-002_2.fastq, ……;
    type2-001_1.fastq, type2-001_2.fastq, type2-002_1.fastq, type2-002_2.fastq, ……;
    ……

where "\_1" and "\_2" indicate forward and backward reads for pair-end sequencing. {type1, type2, ...} can be cell-types or batches of samples, such as {GM, K562, ...}, or {batch1, batch2, ...}, or any other words without underline “_” or dash “-”.
The **work**, **matrix**, **peak**, **result** and **figure** folders will be automatically built by subsequent steps, and placed in $source folder.
 
### 2.2	Easy-run of matrix preparation

Users can use the script APEC_prepare_steps.sh to finish the process from raw data to fragment count matrix.  This script includes steps of “trimming”, “mapping”, “peak calling”, “aligning read counts matrix”, “quality contral”, “estimating gene score”.

Example:

    bash APEC_prepare_steps.sh -s $source -g hg19 -n 4 -l 3 -t 0.2 -f 2000

Input parameters:

    -s: source path, which should contain data folder before running APEC.
    -g: hg19 or mm10.
    -n: Number of CPU cores.
    -l: Threshold for the –log(Q-value) of peaks, used to filter peaks.
    -p: Threshold of the percentage of fragments in peaks, used to filter cells.
    -f: Threshold of the fragment number of each cell, used to filter cells.

Output files:

The script ***APEC_prepare_steps.sh*** will generate **work**, **peak**, **matrix**, and **figure** folders with many output files. Here, we only introduce files that are useful to users.

(1) In **data** folder:
 
    cell_info.csv: Two-column list of cell information.

(2) In **work** folder:

For each cell, the mapping step can generate a subfolder (with cell name) in the **work** folder. There are several useful files in each subfolder:

    cell_name.hist.pdf: A histogram of fragment length distribution of each cell.
    cell_name.RefSeqTSS.pdf: Insert enrichment around TSS regions of each cell.

(3) In **peak** folder:

    mergeAll.hist.pdf: A histogram of fragment length distribution of all cells.
    mergeAll.RefSeqTSS.pdf: Insert enrichment around TSS regions of all cells.
    top_peaks.bed: Top peaks obtained by Q-value filtering.
    annotate_peak.bed: Annotation information of peaks. 
    genes_scored_by_peaks.csv: Gene scores evaluated by TSS peaks.

(4) In **matrix** folder:

    reads.csv: Fragment count matrix.
    cell_info.merged.csv: Data quality report of each cell.
    filtered_reads.csv: Filtered fragment count matrix.

(5) In **figure** folder:

    cell_quality.pdf: A scatter plot of the fragment number and the percentage of fragments in peaks.

## 3.	Clustering

### 3.1	If users want to apply their own matrix

**Skip this step** if users have run APEC_prepare_steps.sh in Section 2 “Fragment count matrix” and generated fragment count matrix from the raw scATAC-seq data.

If users have their own fragment count matrix, please build data, work, matrix, peak, result and figure folders in $source path, and place “cell_info.csv” file in data folder, “top_peaks.bed” in peak folder, “filtered_reads.csv” in matrix folder. Then users need to run script prepare_premappedMatrix.py before clustering and further analysis.

Example:

    python prepare_premappedMatrix.py -s $source --ref hg19 --fa $APEC/reference/hg19_chr.fa --np 4

Input parameters:

    -s: source path that contains data, work, matrix, peak, result and figure folders.
    --ref: hg19 or mm10.
    --fa: hg19_chr.fa or mm10_chr.fa in reference folder.
    --np: Number of CPU cores.

Details about initial files:

    cell_info.csv: Two-column (separated by tabs) list of cell information, such as:
                        	name    notes
                        	CD4-001 CD4
                        	CD4-002 CD4
                        	CD8-001 CD8
                        	CD8-002 CD8
    top_peaks.bed: Three-column list of peaks (500 BP accessible regions), such as:
                        	chr1    3094060 3094560
                        	chr1    3094768 3095268
                        	chr1    3113480 3113980
                        	chr1    3119987 3120487
    filtered_reads.csv: Fragment count matrix, where each row is a cell and each column represents peak. 
                        The name of cells should be same with “cell_info.csv”, and the order of peaks 
                        should be same with “top_peaks.bed”. All numbers are separated by commas, for example:
				,peak1,peak2,peak3,peak4,peak5,peak6
			        CD4-001,1,0,2,0,0,1
			        CD4-002,0,1,0,0,2,0
			        CD8-001,0,0,0,1,1,0
			        CD8-002,2,1,0,0,0,1

### 3.2 Clustering based on accessons

To clustering cells by accessons, users can run the script cluster_byAccesson.py on the fragment count matrix.

Example:

    Python cluster_byAccesson.py -s source --nc 0 --space pca

Input parameters:

    -s: source path that contains data, work, matrix, peak, result and figure folders.
    --nc: Number of cell clusters. If nc=0, it will be predicted by Louvain algorithm. 
    --space: Transfromed space used for clustering, can be “pca” or “tsne”.

Output files important to users:

(1) In **result** folder:

    TSNE_by_Accesson.csv: Top 2 components of tSNE transformation.
    KNN_cluster_by_Accesson.csv: Result of K-nearest-neighbor (KNN) clustering.
    Hierarchical_cluster_by_Accesson.csv: Result of Hierarchical clustering.

(2) In **figure** folder:

    PCA_by_Accesson.pdf: PCA diagram of cells, colored by the “notes” column in the ~/data/cell_info.csv file.
    TSNE_by_Accesson.pdf: tSNE diagram of cells, colored by the “notes” column in the ~/data/cell_info.csv file.
    KNN_cluster_by_Accesson.pdf: KNN clustering result of cells, plotted on tSNE map.
    cell_cell_correlation_by_Accesson.png: Hierarchical clustering heat map of cell-cell correlation matrix.
                                           Colors of sidebar are defined by the “notes” column in 
                                           the ~/data/cell_info.csv file.
    HC_KNN_compare_by_Accesson.png: Hierarchical clustering heat map of cell-cell correlation matrix.
                                    Colors of sidebar are defined by KNN clustering result.

### 3.3 Clustering based on motifs (as chromVAR)

The script cluster_byMotif.py provides a python version of motif-based cell clustering, with the same algorithm as chromVAR.

Example:

    python cluster_byMotif.py -s $source --ns 8 --np 8 --nc 0 --space pca

Input parameters:

    -s: source path that contains data, work, matrix, peak, result and figure folders.
    --ns: Number of permuted samplings to calculate the bias-corrected deviations of motifs.
          Its default value is 50.
    --np: Number of CPU cores.
    --nc: Number of cell clusters. If nc=0, it will be predicted by Louvain algorithm. 
    --space: Transfromed space used for clustering, can be “pca” or “tsne”.

Output files description:

(1) In result folder:

    deviation_chromVAR.csv: The bias-corrected deviation matrix of motifs
    TSNE_by_chromVAR.csv, KNN_cluster_by_chromVAR.csv, and Hierarchical_cluster _by_chromVAR.csv 
    are similar to the results in Section 3.2, except that they are generated by motif-based 
    clustering (same as chromVAR).

(2) In figure folder:

    PCA_by_chromVAR.pdf, TSNE_by_chromVAR.pdf, KNN_cluster_by_chromVAR.pdf, cell_cell_correlation
    _by_chromVAR.png, and HC_KNN_compare_by_chromVAR.png are similar to the results in Section 3.2,
    except that they are generated by motif-based clustering (same as chromVAR).


