
#### download miniconda

    wget https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86_64.sh

#### install miniconda

    bash Miniconda2-latest-Linux-x86_64.sh

#### set channels

    conda config --add channels conda-forge
    conda config --add channels defaults
    conda config --add channels r
    conda config --add channels bioconda


#### set a new environment for APEC, i.e. apec_env

    conda create -n apec_env python=2.7

#### install python packages

    conda install -n apec_env bowtie2
    conda install -n apec_env samtools
    conda install -n apec_env bedtools
    conda install -n apec_env homer
    conda install -n apec_env macs2
    conda install -n apec_env meme=4.11.2
    conda install -n apec_env ucsc-bedgraphtobigwig

    conda install -n apec_env numpy
    conda install -n apec_env scipy
    conda install -n apec_env pandas
    conda install -n apec_env -c auto multiprocessing
    conda install -n apec_env numba
    conda install -n apec_env pysam
    conda install -n apec_env matplotlib
    conda install -n apec_env seaborn
    conda install -n apec_env setuptools
    conda install -n apec_env networkx
    conda install -n apec_env python-louvain
    conda install -n apec_env python-Levenshtein
    conda install -n apec_env scikit-learn
    
    conda install -n apec_env bioconductor-monocle

#### activate apec_env

    source activate apec_env

#### download Homer packages

    perl ~/miniconda2/envs/apec_env/share/homer-X.X.X/configureHomer.pl -install hg19
    perl ~/miniconda2/envs/apec_env/share/homer-X.X.X/configureHomer.pl -install mm10
    
    