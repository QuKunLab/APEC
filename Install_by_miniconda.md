
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
    conda install -n apec_env homer  (not required for code_v1.0.1)
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
    
    conda install -n apec_env bioconductor-monocle (not recommended, please install R and Monocle independently)

**Note**: We found some problems for the R enviroment installed by conda, so we recommend that users do not use conda to install the R environment and Monocle.

**Note**: We also found some problems when installing non-python software via conda/bioconda, so we recommend that users do not use conda/bioconda to install bowtie2, samtools, bedtools, macs2, meme (4.11.2) and bedgraghtobigwig.

#### activate apec_env

    conda activate apec_env

#### download Homer packages (not required for code_v1.0.1)

    perl ~/miniconda2/envs/apec_env/share/homer-X.X.X/configureHomer.pl -install hg19
    perl ~/miniconda2/envs/apec_env/share/homer-X.X.X/configureHomer.pl -install mm10
    
    
