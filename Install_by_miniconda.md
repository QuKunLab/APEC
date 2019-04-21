
#### download miniconda

    wget https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86_64.sh

#### install miniconda

    bash Miniconda2-latest-Linux-x86_64.sh

#### set channels, all channels must be arranged in the following order

    conda config --add channels bioconda
    conda config --add channels r
    conda config --add channels defaults
    conda config --add channels conda-forge

#### set a new environment for APEC, i.e. apec_env

    conda create -n apec_env python=2.7

#### install python packages

    conda install -n apec_env bowtie2
    conda install -n apec_env samtools
    conda install -n apec_env bedtools
    conda install -n apec_env macs2
    conda install -n apec_env meme=4.11.2
    conda install -n apec_env ucsc-bedgraphtobigwig

    conda install -n apec_env numpy
    conda install -n apec_env scipy=1.0.0
    conda install -n apec_env pandas
    conda install -n apec_env -c auto multiprocessing
    conda install -n apec_env numba
    conda install -n apec_env pysam
    conda install -n apec_env matplotlib
    conda install -n apec_env seaborn
    conda install -n apec_env setuptools
    conda install -n apec_env networkx
    conda install -n apec_env python-louvain=0.11
    conda install -n apec_env python-Levenshtein
    conda install -n apec_env scikit-learn=0.20.0
    conda install -n apec_env multicore-tsne

    conda install -n apec_env bioconductor-monocle=2.4.0 (It's not recommended to install monocle with conda)
    conda install -n apec_env libiconv   (required if you use conda to install monocle)
    conda install -n apec_env r-cluster   (required if you use conda to install monocle)
    conda install -n apec_env r-stringr=1.2.0   (required if you use conda to install monocle)

**Note**: We found some problems for the R enviroment installed by conda, so we recommend that users do not use conda to install the R environment and Monocle.

**Note**: We also found some problems when installing non-python software via conda/bioconda, so we recommend that users do not use conda/bioconda to install bowtie2, samtools, bedtools, macs2, meme (4.11.2) and bedgraghtobigwig.

#### activate apec_env

    conda activate apec_env
