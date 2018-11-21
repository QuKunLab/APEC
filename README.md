User Guide for APEC

(Accessibility Pattern based Epigenomic Clustering)
======
Introduction
------

APEC can perform fine cell type clustering on single cell chromatin accessibility data from scATAC-seq, snATAC-seq, sciATAC-seq or any other relevant experiment. It can also be used to control data quality, map fragment count matrices, search for important differential motifs/genes for each cell cluster, find super enhancers, and construct pseudo-time trajectory (by calling Monocle).

The primary algorithm of APEC is based on grouping accessible regions throuhg their accessibility patterns prior to cell clustering. We define a group of accessible regions as an accesson, and apply the accesson count matrix for further analysis.

If the user wants to process the raw fastq data from scATAC-seq experiment, please run APEC from section 2 “Fragment count matrix”. If the user has obtained the fragment count matrix, where each element is the number of fragments per-cell-per-peak, please run APEC from section 3 “Clustering”.

