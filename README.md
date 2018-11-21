# User Guide for APEC

(Accessibility Pattern based Epigenomic Clustering)

## Introduction

APEC can perform fine cell type clustering on single cell chromatin accessibility data from scATAC-seq, snATAC-seq, sciATAC-seq or any other relevant experiment. It can also be used to control data quality, map fragment count matrices, search for important differential motifs/genes for each cell cluster, find super enhancers, and construct pseudo-time trajectory (by calling Monocle).

The primary algorithm of APEC is based on grouping accessible regions throuhg their accessibility patterns prior to cell clustering. We define a group of accessible regions as an **accesson**, and apply the accesson count matrix for further analysis.

**If the user wants to process the raw fastq data from scATAC-seq experiment, please run APEC from section 2 “Fragment count matrix”. If the user has obtained the fragment count matrix, where each element is the number of fragments per-cell-per-peak, please run APEC from section 3 “Clustering”.**

## 1.	Requirements and installation

### 1.1	Requirements

APEC requires users to use Linux system, as well as Python (version 2.7.5+) and R (version 3.4+) environment. Users also need to install the following packages:

(1)	Python packages and libraries: 

numpy, scipy, pandas, sklearn, multiprocessing, numba, pysam, matplotlib, seaborn, networkx, python-louvain, python-Levenshtein

all upon python packages can be installed by:
**pip install package_name**

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
Users simply completes the APEC installation by copying the APEC folder to any path on the computer (i.e. $PATH). There are two subfolders in APEC: a <codes> folder, which contains all APEC programs for data processing; a <reference> folder, which contains all necessary index and reference files for the hg19 and mm10 genomes. So users can run APEC program directly in $PATH/APEC/codes/, or put this path in the system environment to use it elsewhere. The <reference> folder is required for APEC and should be placed in the same path with the <codes> folder. It contains the following files:
 
hg19_refseq_genes_TSS.txt, hg19_RefSeq_genes.gtf, hg19_blacklist.JDB.bed, 

hg19_chr.fa, hg19_chr.fa.fai, hg19.chrom.sizes,

hg19.1.bt2, hg19.2.bt2, hg19.3.bt2, hg19.4.bt2,

mm10_refseq_genes_TSS.txt, mm10_RefSeq_genes.gtf, mm10_blacklist.BIN.bed, 

mm10_chr.fa, mm10_chr.fa.fai, mm10.chrom.sizes, 

mm10.1.bt2, mm10.2.bt2, mm10.3.bt2, mm10.4.bt2,

JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt, tier1_markov1.norc.txt


