**Updates on 2019-05-28**

    Fit the problem that cluster_byMotif.py cannot be used on the count matrix of mtx format (i.e. filtered_reads.mtx).
    Also, in cluster_byMotif.py we adopted MulticoreTSNE instead of sklearn.manifold.TSNE

**Updates on 2019-05-11**
    
    Update generate_differential_markers.py:
        If '--motif' set to 'no', 'projec/result/deviation_chromVAR.csv' is not required.

**Updates on 2019-04-27**

    We corrected the error on line 388 of subroutines.py

**Updates on 2019-04-25**

    We corrected the error on line 123 of cluster_byAccesson.py.

**Updates on 2019-04-21**

    There is a problem with the TSNE package in scikit-learn when analyzing a large number of cells.
    Now We use mulitcore-tsne in cluster_byAccesson.py to analyze datasets containing more than 10,000 cells.