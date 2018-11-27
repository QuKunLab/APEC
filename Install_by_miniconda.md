
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

#### activate apec_env

source activate apec_env

#### install python packages

conda install -n apec_env numpy
conda install -n apec_env scipy
conda install -n apec_env pandas
conda install -n apec_env scikit-learny
conda install -n apec_env -c auto multiprocessing
conda install -n apec_env numba
conda install -n apec_env pysam
conda install -n apec_env matplotlib
conda install -n apec_env seaborn
conda install -n apec_env networkx
conda install -n apec_env python-louvain
conda install -n apec_env python-Levenshtein