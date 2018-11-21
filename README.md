# User Guide for APEC

(Accessibility Pattern based Epigenomic Clustering)

## Introduction

APEC can perform fine cell type clustering on single cell chromatin accessibility data from scATAC-seq, snATAC-seq, sciATAC-seq or any other relevant experiment. It can also be used to control data quality, map fragment count matrices, search for important differential motifs/genes for each cell cluster, find super enhancers, and construct pseudo-time trajectory (by calling Monocle).

**If users want to process the raw fastq data from scATAC-seq experiment, please run APEC from section 2 “Fragment count matrix”. If users have obtained the fragment count matrix, where each element is the number of fragments per-cell-per-peak, please run APEC from section 3 “Clustering”.**

## 1.	Requirements and installation

### 1.1	Requirements

APEC requires users to use Linux system, as well as Python (version 2.7.5+) and R (version 3.4+) environment. Users also need to install the following packages:

(1)	Python packages and libraries: 

    numpy, scipy, pandas, sklearn, multiprocessing, numba, pysam,
    matplotlib, seaborn, networkx, python-louvain, python-Levenshtein
    
    all upon python packages can be installed by: 
    pip install package_name

(2)	R packages and libraries: 

    Monocle: http://cole-trapnell-lab.github.io/monocle-release/

(3)	Other necessary software:

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

## 2.	Fragment count matrix

### 2.1	Arrangement of raw data

Users need to build a source folder (i.e. $source), which contains a **data** folder, then copy all raw sequencing fastq files into the <$source/data/> folder. All these pair-end fastq files should be named as:
 
    type1-001_1.fastq, type1-001_2.fastq, type1-002_1.fastq, type1-002_2.fastq, ……;
    type2-001_1.fastq, type2-001_2.fastq, type2-002_1.fastq, type2-002_2.fastq, ……;
    ……

where "\_1" and "\_2" indicate forward and backward reads for pair-end sequencing. {type1, type2, ...} can be cell-types or batches of samples, such as {GM, K562, ...}, or {batch1, batch2, ...}, or any other words without underline “_” or dash “-”.
The **work**, **matrix**, **peak**, **result** and **figure** folders will be automatically built by subsequent steps, and placed in $source folder.
 
### 2.2	Easy-run of matrix preparation

Users can use the script APEC_prepare_steps.sh to finish the process from raw data to fragment count matrix.  This script includes steps of “trimming”, “mapping”, “peak calling”, “aligning read counts matrix”, “quality contral”, “estimating gene score”.

Usage: 

    bash APEC_prepare_steps.sh -s source_folder -g genome_index
                               -n nCPUs -l logq -t tssfrag -f frag

Example:

    bash APEC_prepare_steps.sh -s /home/user/test -g hg19 -n 4 
                               -l 3 -t 0.2 -f 2000
Input parameters:

    -s: source path, which should contain <data> folder before running APEC.
    -g: hg19 or mm10.
    -n: Number of CPU cores.
    -l: Threshold for the –log(Q-value) of peaks, used to filter peaks.
    -p: Threshold of the percentage of fragments in peaks, used to filter cells.
    -f: Threshold of the fragment number of each cell, used to filter cells.
Output files:

The script ***APEC_prepare_steps.sh*** will generate **work**, **peak**, **matrix**, and **figure** folders with many output files. Below are the details about these output files.

(1) Files in **data** folder:
 
    cell_name.trim.fastq: Trimmed data from raw fastq file. 
    cell_info.csv: A list of cell information containing two columns, for example::
                name    notes
                CD4-001   CD4
                CD4-002   CD4
                CD8-001   CD8
                CD8-002   CD8

(2) Files in **work** folder:

For each cell, the mapping step can generate a subfolder (with cell name) in the **work** folder. There are several useful files in each subfolder:

    cell_name.hist.pdf: A histogram of fragment length distribution of each cell.
    cell_name.RefSeqTSS.pdf: Insert enrichment around TSS regions of each cell.

(3) Files in **peak** folder:

    mergeAll.hist.pdf: A histogram of fragment length distribution of all cells.
    mergeAll.RefSeqTSS.pdf: Insert enrichment around TSS regions of all cells.
    top_peaks.bed: Top peaks obtained by Q-value filtering.
    annotate_peak.bed: Annotation information of peaks. 
    genes_scored_by_peaks.csv: Gene scores evaluated by TSS peaks.

(4) Files in **matrix** folder:

    reads.csv: Fragment count matrix.
    cell_info.merged.csv: Data quality report of each cell.
    filtered_reads.csv: Filtered fragment count matrix.

(5) Files in **figure** folder:

    cell_quality.pdf: A scatter plot of the fragment number and the percentage of fragments in peaks.



