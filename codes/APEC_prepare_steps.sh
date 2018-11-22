#!/usr/bin/bash
#
#### input parameters defined by users #############################################
#
ARGS=`getopt -o hs:g:n:l:t:f: -l help,source:,genome:,np:,logq:,pfrag:,frag: -- "$@"`
if [ $? != 0 ] ; then echo "Terminating..." >&2 ; exit 1 ; fi
eval set -- "$ARGS"
while true ; do
    case "$1" in
        -h|--help)
           echo "
bash APEC_prepare_steps.sh -s source_folder -g genome_index -n nCPUs -l logq -t pfrag -f frag
     -s/--source:   source path, which should contain <data> folder before running APEC.
     -g/--genome:   hg19 or mm10.
     -n/--np:       Number of CPU cores.
     -l/--logq:     Threshold for the -log(Q-value) of peaks, used to filter peaks.
     -p/--pfrag:    Threshold of the percentage of fragments in peaks, used to filter cells.
     -f/--frag:     Threshold of the fragment number of each cell, used to filter cells."
           exit 1 ;;
        -s|--source) source="$2" ; shift 2;;
        -g|--genome) genome="$2" ; shift 2;;
        -n|--np) np="$2" ; shift 2;;    
        -l|--logq) logq="$2" ; shift 2;;
        -t|--pfrag) pfrag="$2" ; shift 2;;
        -f|--frag) frag="$2" ; shift 2;;
        --) shift; break ;;
        *) echo "unknown parameter: {$1}" ; exit 1 ;;
    esac
done
#
picard=../reference/picard.jar
ref=$genome
fa="../reference/"$genome"_chr.fa"
index="../reference/"$genome
tss="../reference/"$genome"_refseq_genes_TSS.txt"
if [[ $genome == "hg19" ]]; then
    blist=../reference/hg19_blacklist.JDB.bed
elif [[ $genome == "mm10" ]]; then
    blist=../reference/mm10_blacklist.BIN.bed
fi
gtf="../reference/"$genome"_RefSeq_genes.gtf"
np=$np
logq=$logq
pfrag=$pfrag
frag=$frag
#
#
#
#### processes to prepare raw data ###########
#
python prepare_trimming.py -s $source --np $np
#
python prepare_mapping.py -s $source --index $index --picard $picard --tss $tss --np $np
#
python prepare_peakCalling.py -s $source --blist $blist --fa $fa --tss $tss --ref $ref --logq $logq
#
python prepare_countMatrix.py -s $source --fa $fa --np $np
#
python prepare_qualityControl.py -s $source --pfrag $pfrag --lib $frag
#
python prepare_geneScore.py -s $source --gtf $gtf
#
#
#
#
