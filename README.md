# User Guide for APEC (v1.1.0)

(Accessibility Pattern based Epigenomic Clustering)

APEC can perform fine cell type clustering on single cell chromatin accessibility data from scATAC-seq, snATAC-seq, sciATAC-seq or any other relevant experiment. It can also be used to evaluate gene expression from relevant accesson, search for differential motifs/genes for each cell cluster, find super enhancers, and construct pseudo-time trajectory (by calling Monocle). **If users have already obtained the fragment-count-per-peak matrix from other mapping pipelines (such as CellRanger), please run APEC from the first Section "Run APEC from fragment count matrix". If users have only the raw fastq files, please jump to the second Section "Get fragment count matrix from raw data".**

## Run AEPC from fragment count matrix

### 1. Requirements and installation

#### 1.1 Requirements

APEC requires Linux system (CentOS 7.3+ or Ubuntu 16.04+), as well as Python (2.7.15+ or 3.6.8+). If users want to build pseudotime trajectory with APEC, please install R (3.5.1) environment and monocle (2.10.0). Also, the following software are required for APEC:

    Bedtools: http://bedtools.readthedocs.io/en/latest/content/installation.html
    Meme 4.11.2: http://meme-suite.org/doc/download.html?man_type=web
    Homer: http://homer.ucsd.edu/homer/

**notes: Users need to download genome reference for Homer by "perl /path-to-homer/configureHomer.pl -install hg19" and "perl /path-to-homer/configureHomer.pl -install mm10".**

The files in **reference** folder are required for APEC. **But we didn't upload reference files to GitHub since they are too big. Users can download all reference files from http://galaxy.ustc.edu.cn:30803/APEC/**. The **reference** folder should contains the following files:

    hg19_RefSeq_genes.gtf, hg19_chr.fa, hg19_chr.fa.fai,
    mm10_RefSeq_genes.gtf, mm10_chr.fa, mm10_chr.fa.fai,
    JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt, tier1_markov1.norc.txt

#### 1.2 Install and import APEC

Users can install APEC by:

    pip install APEC==1.1.0.7

Due to the compatibility problem (especially for rpy2), we don't recommend conda environment. Users can use **pyenv** to build a sub environment for APEC. If users want to call Paga (instead of monocle) to construct pseudotime trajectory, please use APEC in Python3 environment and install the following packages:

    pip install scanpy anndata

In Ipython, Jupyter-notebook or a python script, users can import packages of APEC by:

    from APEC import clustering,plot,generate

Users can inquire the manual for each function of APEC by using "help()" in Ipython or Jupyter, for example:

    help(clustering.cluster_byAccesson)

### 2. Input data

Users need to prepare a project folder (termed '$project'), which contains **matrix**, **peak**, **result** and **figure** folders. Please place "filtered_cells.csv" and "filtered_reads.mtx" in **matrix** folder, "top_filtered_peaks.bed" in **peak** folder. Here is the instruction for three input files:

    filtered_cells.csv: Two-column (separated by tabs) list of cell information ('name' and 'notes'):
                        The 'name' column stores cell names (or barcodes); the 'notes' column can be cell-type,
                        development stage, batch index or any other cell information, such as:
                        	name    notes
                        	CD4-001 CD4
                        	CD4-002 CD4
                        	CD8-001 CD8
                        	CD8-002 CD8
    top_filtered_peaks.bed: Three-column list of peaks, which is a standard bed format file.
                            It is similar to the "peaks.bed" file in the CellRanger output of a 10X scATAC-seq dataset.
    filtered_reads.mtx: Fragment count matrix in mtx format, where each row is a peak and each column represents a cell.
                        It is similar to the "matrix.mtx" file in the CellRanger output of a 10X scATAC-seq dataset.
                        The order of cells should be the same with "filtered_cells.csv", and the order of peaks should
                        be the same with "top_filtered_peaks.bed".

### 3. Functions of APEC (step by step)

#### 3.1 Clustering by APEC

Use the following codes to cluster cells by APEC algorithm:

    clustering.build_accesson('$project', ngroup=600)
    clustering.cluster_byAccesson('$project', norm='zscore')

input parameters:

    ngroup:   Number of accessons, default=600.
    nc:       Number of cell clusters, set it to 0 if users want to predict cluster number by Louvain algorithm, default=0.
    norm:     Normalization method for accesson matrix, can be 'zscore' or 'probability', default='zscore'.
    filter:   Filter high dispersion accessons or not, can be 'yes' or 'no', default='yes'.

output files:

    $project/matrix/Accesson_peaks.csv
    $project/matrix/Accesson_reads.csv
    $project/result/louvain_cluster_by_APEC.csv

Then users can plot tSNE, UMAP or corrlation heatmap for cells:

    plot.plot_tsne('$project', cell_label='notes')
    plot.plot_tsne('$project', cell_label='cluster')
    plot.plot_umap('$project', cell_label='notes')
    plot.correlation('$project', cell_label='notes')

input parameters:

    matrix_type:    Type of input matrix, can be 'APEC' or 'chromVAR', default='APEC'.
                    If matrix_type='APEC', it will use accesson matrix yielded by clustering.cluster_byAccesson();
                    if matrix_type='chromVAR', it will use deviation matrix yielded by clustering.cluster_byMotif().
    cell_label:     Color labels for cells, can be 'notes' or 'cluster', default='notes'.
                    If cell_label='cluster', it will use clustering result of clustering.cluster_byXXX().
    cluster:        Clustering algorithm used in clustering.cluster_byXXX(), default='louvain'.

output files:

    $project/result/TSNE_by_APEC.csv
    $project/figure/TSNE_by_APEC_with_notes_label.pdf
    $project/figure/TSNE_by_APEC_with_cluster_label.pdf
    $project/result/UMAP_by_APEC.csv
    $project/figure/UMAP_by_APEC_with_notes_label.pdf
    $project/figure/cell_cell_correlation_by_APEC_with_notes_label.png

<img src="images/TSNE_by_APEC_with_notes_label.jpg" width="500">

_TSNE_by_APEC_with_notes_label.pdf_

<img src="images/cell_cell_correlation_by_APEC_with_notes_label.jpg" width="500">

_cell_cell_correlation_by_APEC_with_notes_label.png_

#### 3.2 Clustering by chromVAR (optional, required for motif analysis)

Use the following codes to cluster cells by chromVAR algorithm:

    generate.motif_matrix('$project', genome_fa='$reference/hg19_chr.fa',
                          background='$reference/tier1_markov1.norc.txt',
                          meme='$reference/JASPAR2018_CORE_vertebrates_redundant_pfms_meme.txt',
                          np=4)
    clustering.cluster_byMotif('$project', np=4)

input parameters:

    genome_fa:   Path to hg19_chr.fa or mm10_chr.fa in $reference folder.
    background:  Path to tier1_markov1.norc.txt in $reference folder.
    meme:        Path to JASPAR2018_CORE_vertebrates_redundant_pfms_meme.txt in $reference folder.
    np:          Number of CPU cores used for parallel calculation, default=4.
    nc:          Number of cell clusters, set it to 0 if users want to predict cluster number by Louvain algorithm, default=0.
    ns:          Number of permuted sampling, default=50.

output files:

    $project/result/deviation_chromVAR.csv
    $project/result/louvain_cluster_by_chromVAR.csv

#### 3.3 Evaluate ARI, NMI and AMI for clustering result

If users have the real cell type in the 'notes' column of '$project/matrix/filtered_cells.csv', please use the following code to calculate ARI, NMI and AMI to estimate the accuracy of the clustering algorithm.

    clustering.cluster_comparison('$project/matrix/filtered_cells.csv',
                                  '$project/result/louvain_cluster_by_APEC.csv')

The output ARI, NMI and AMI values will present on the screen directly. Please make sure filtered_cells.csv contains the FACS label for each cell. For some datasets, such as Hematopoietic cells, the user should ignore all "unknown" cells before the calculation of ARI.

#### 3.4 Generate pseudotime trajectory

By default, APEC adapts monocle to generate pseudotime trajectory from accesson matrix:

    generate.monocle_trajectory('$project', npc=5)
    plot.plot_trajectory('$project', cell_label='notes', angles=[30,30])

input parameters:

    npc:            Number of principal components used to build trajectory, default=5.
    cell_label:     Color labels for cells, can be 'notes' or 'cluster', default='notes'.
    cluster:        Clustering algorithm used in clustering.cluster_byXXX(), default='louvain'.
    angles:         Rotation angles for 3D trajectory, e.g. [100,20], default=[30,30].

output files:

    $project/result/monocle_trajectory.csv
    $project/result/monocle_reduced_dimension.csv
    $project/figure/pseudotime_trajectory_with_notes_label.pdf

In Python3 environment, users can also use Paga to generate trajectory:

    generate.paga_trajectory('$project', cell_label='notes', npc=5)

output files:

    $project/figure/paga_skeleton_with_notes_label.pdf
    $project/figure/paga_trajectory_with_notes_label.pdf

<img src="images/pseudotime_trajectory_with_notes_label.jpg" width="500">

_pseudotime_trajectory_with_notes_label.pdf_

#### 3.5 Generate gene expression

    generate.gene_expression('$project', genome='hg19', width=1000000, pvalue=0.01)

input parameters:

    genome:      Genome reference for Homer, can be "hg19" or "mm10", default="hg19".
    width:       Width of Genome region for fisher exact test, default=1000000.
    pvalue:      P-value threshold for fisher exact test, default=0.01.

output files:

    $project/matrix/Accesson_annotated.csv
    $project/matrix/gene_annotated.csv
    $project/matrix/gene_expression.csv

Users can also score genes by the peaks around their TSS regions:

    generate.gene_score('$project', genome_gtf='hg19_RefSeq_genes.gtf', distal=20000)

output file:

    $project/matrix/genes_scored_by_TSS_peaks.csv

#### 3.6 Generate differential feature for a cell cluster

Get differential accessons:

    generate.nearby_genes('$project', genome_gtf='hg19_RefSeq_genes.gtf', distal=20000)   # optional step
    generate.differential_feature('$project', feature='accesson', target='0', vs='all')

Get differential motifs/genes:

    generate.differential_feature('$project', feature='motif', target='0', vs='all')
    generate.differential_feature('$project', feature='gene', target='0', vs='all')

input parameters:

    feature:     Type of feature, can be 'accesson' or 'motif' or 'gene', default='accesson'.
                 If feature='accesson', run clustering.cluster_byAccesson() first;
                 if feature='motif', run clustering.cluster_byMotif() first;
                 if feature='gene', run generate.gene_expression() first.
    cell_label:  Cell labels used for differential analysis, can be 'notes' or 'cluster', default='cluster'.
    cluster:     Clustering algorithm used in clustering.cluster_byXXX(), default='louvain'.
    target:      The target cluster that users search for differential features, default='1'.
                 If cell_label='cluster', target is one element in the 'cluster' column of XXX_cluster_by_XXX.csv file;
                 if cell_label='notes', target is one element in the 'notes' column of filtered_cells.csv file.
    vs:          Versus which clusters, can be '2,3,4' or 'all', default='all' (means all the rest clusters).
    pvalue:      P-value for student-t test, default=0.01.
    log2_fold:   Cutoff for log2(fold_change), default=1.

output file:

    $project/result/differential_accesson_of_cluster_X_vs_XXX.csv
    $project/result/differential_motif_of_cluster_X_vs_XXX.csv
    $project/result/differential_gene_of_cluster_X_vs_XXX.csv

#### 3.7 Plot motif/gene on tSNE/trajectory diagram

    plot.plot_feature('$project', space='tsne', feature='gene', name='FOXO1')
    plot.plot_feature('$project', space='trajectory', feature='motif', name='GATA1')

input parameters:

    space:          In which space we draw the feature, can be 'tsne' or 'umap' or 'trajectory', default='tsne'.
                    If space='tsne', run plot.plot_tsne() first;
                    if space='umap', run plot.plot_umap() first;
                    if space='trajectory', run generate.monocle_trajectory() first.
    feature:        Type of the feature, can be 'accesson' or 'motif' or 'gene', default='accesson'.
                    If feature='accesson', run clustering.cluster_byAccesson() first;
                    if feature='motif', run clustering.cluster_byMotif() first;
                    if feature='gene', run generate.gene_expression() first.
    matrix_type:    Type of input matrix for tSNE/UMAP, can be 'APEC' or 'chromVAR', default='APEC'.
                    If matrix_type='APEC', it will use tSNE/UMAP result of APEC;
                    if matrix_type='chromVAR', it will use tSNE/UMAP result of chromVAR.
    name:           Name of the feature.
                    If feature='accesson', name=accesson number, i.e. '1';
                    if feature='motif', name=motif symbol, i.e. 'GATA1';
                    if feature='gene', name=gene symbol, i.e. 'CD36'.
    clip:           Clip range for the input matrix, can be [min, max] or 'none', default='none'.
    angles:         Rotation angles for 3D trajectory, e.g. [100,20], default=[30,30].

output files:

    $project/figure/gene_FOXO1_on_tsne_by_APEC.pdf
    $project/figure/motif_GATA1_on_trajectory_by_APEC.pdf

<img src="images/gene_FOXO1_on_tsne_by_APEC.jpg" width="500">

_gene_FOXO1_on_tsne_by_APEC.pdf_

<img src="images/motif_GATA1_on_trajectory_by_APEC.jpg" width="500">

_motif_GATA1_on_trajectory_by_APEC.pdf_

**Notes: Plotting feature on tSNE diagram requires the running of plot.plot_tsne() beforehand (see 3.1), and plotting feature on trajectory requires the running of generate.monocle_trajectory() beforehand (see 3.4).**

#### 3.8 Generate potential super enhancer

    generate.search_super_enhancer('$project', super_range=1000000)

input parameter:

    super_range:    Genome range to search for super enhancer, default=1000000.

output file:

    $project/result/potential_super_enhancer.csv


## Get fragment count matrix from raw data
## (this part is only available on GitHub:https://github.com/QuKunLab/APEC)

### 1. Requirements and installation

All of the following software needs to be placed in the global environment of the Linux system to ensure that they can be called in any path/folder. Picard is also required, but we have placed it into $APEC/reference folder, and users don't need to install it. We recommend that users adopt the latest version of these software, except Meme (version 4.11.2).

    Bowtie2: https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.9/
    Samtools: https://github.com/samtools/samtools
    Bedtools: http://bedtools.readthedocs.io/en/latest/content/installation.html
    Macs2: https://github.com/taoliu/MACS.git
    Meme 4.11.2: http://meme-suite.org/doc/download.html?man_type=web
    pysam for python: set up by "pip install pysam"
    Levenshtein for python: set up by "pip install python-Levenshtein"

### 1.2	Installation

Users can simply install this part by copying the **code_v1.1.0** folder and **reference** folder into a same path. Users **must** run ***APEC_prepare_steps.sh*** directly in code_v1.1.0/, since each program will invoke the reference files automatically. The **reference** folder is required, **but we didn't upload reference files to GitHub since they are too big. Users can download all reference files from http://galaxy.ustc.edu.cn:30803/APEC/**. The **reference** folder should contains the following files:

    hg19_refseq_genes_TSS.txt, hg19_RefSeq_genes.gtf, hg19_blacklist.JDB.bed,
    hg19_chr.fa, hg19_chr.fa.fai, hg19.chrom.sizes,
    hg19.1.bt2, hg19.2.bt2, hg19.3.bt2, hg19.4.bt2, hg19.rev.1.bt2, hg19.rev.2.bt2,
    mm10_refseq_genes_TSS.txt, mm10_RefSeq_genes.gtf, mm10_blacklist.BIN.bed,
    mm10_chr.fa, mm10_chr.fa.fai, mm10.chrom.sizes,
    mm10.1.bt2, mm10.2.bt2, mm10.3.bt2, mm10.4.bt2, mm10.rev.1.bt2, mm10.rev.2.bt2,
    JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt, tier1_markov1.norc.txt, picard.jar

## 2. Fragment count matrix

### 2.1	Arrangement of raw data

The **raw_data** folder should contain all raw sequencing fastq files into the. All these pair-end fastq files should be named as:

    type1-001_1.fastq, type1-001_2.fastq, type1-002_1.fastq, type1-002_2.fastq, ……;
    type2-001_1.fastq, type2-001_2.fastq, type2-002_1.fastq, type2-002_2.fastq, ……;
    ……

where "\_1" and "\_2" indicate forward and backward reads for pair-end sequencing. {type1, type2, ...} can be cell-types or batches of samples, such as {GM, K562, ...}, or {batch1, batch2, ...}, or any other words without underline "\_" or dash "-".
Users need to build a **project** folder to store the result. The **work**, **matrix**, **peak** and **figure** folders will be automatically built by subsequent steps, and placed in **project** folder.

### 2.2	Easy-run of matrix preparation

Users can use the script ***APEC_prepare_steps.sh*** to finish the process from raw data to fragment count matrix.  This script includes steps of "trimming", "mapping", "peak calling", "aligning read counts matrix", and "quality contral". Running this step on our example project (i.e. project01 with 672 cells) will take 10~20 hours on an 8-core 32 GB computer, since the sequence mapping step is the slowest step.

Example:

    bash APEC_prepare_steps.sh -r $raw_data -s $project -g hg19 -n 4 -l 3 -p 0.2 -f 2000

Input parameters:

    -r: The raw_data folder
    -s: The project folder.
    -g: "hg19" or "mm10".
    -n: Number of CPU cores.
    -l: Threshold for the –log(Q-value) of peaks, used to filter peaks.
    -p: Threshold of the percentage of fragments in peaks, used to filter cells.
    -f: Threshold of the fragment number of each cell, used to filter cells.

Output files:

The script ***APEC_prepare_steps.sh*** will generate **work**, **peak**, **matrix**, and **figure** folders with many output files. Here, we only introduce files that are useful to users. For our example projects, all of these results can be reproduced on a general computer system.

(1) In **work** folder:

For each cell, the mapping step can generate a subfolder (with cell name) in the **work** folder. There are several useful files in each subfolder:

    cell_name.hist.pdf: A histogram of fragment length distribution of each cell.
    cell_name.RefSeqTSS.pdf: Insert enrichment around TSS regions of each cell.

(2) In **peak** folder:

    mergeAll.hist.pdf: A histogram of fragment length distribution of all cells.
    mergeAll.RefSeqTSS.pdf: Insert enrichment around TSS regions of all cells.
    top_filtered_peaks.bed: Filtered top peaks, ranked by Q-value.

(3) In **matrix** folder:

    reads.csv: Fragment count matrix.
    cell_info.merged.csv: Data quality report of each cell.
    filtered_cells.csv: Filtered cells information in csv format.
    filtered_reads.mtx: Filtered fragment count matrix in mtx format.

(4) In **figure** folder:

    cell_quality.pdf: A scatter plot of the fragment number and the percentage of fragments in peaks.
