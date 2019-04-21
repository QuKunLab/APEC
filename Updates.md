**Updates on 2019-04-21**

    There is a problem with the TSNE package in scikit-learn when analyzing a large number of cells.
    Now We use mulitcore-tsne in cluster_byAccesson.py to analyze datasets containing more than 10,000 cells.