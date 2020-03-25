#!/usr/bin/bash
#
#### input parameters defined by users #############################################
#
ARGS=`getopt -o hr:s:g:n:l:p:f: -l help,raw:,project:,genome:,np:,logq:,pfrag:,frag: -- "$@"`
if [ $? != 0 ] ; then echo "Terminating..." >&2 ; exit 1 ; fi
eval set -- "$ARGS"
while true ; do
    case "$1" in
        -h|--help)
           echo "
bash APEC_prepare_steps.sh -r raw_data -s project -g genome_index -n nCPUs -l logq -p pfrag -f frag
     -r/--raw:      The raw_data path, where the pair-end sequencing data locate.
     -s/--project:  The project path.
     -g/--genome:   hg19 or mm10.
     -n/--np:       Number of CPU cores.
     -l/--logq:     Threshold for the -log(Q-value) of peaks, used to filter peaks.
     -p/--pfrag:    Threshold of the percentage of fragments in peaks, used to filter cells.
     -f/--frag:     Threshold of the fragment number of each cell, used to filter cells."
           exit 1 ;;
        -r|--raw) raw_data="$2" ; shift 2;;
        -s|--project) project="$2" ; shift 2;;
        -g|--genome) genome="$2" ; shift 2;;
        -n|--np) np="$2" ; shift 2;;
        -l|--logq) logq="$2" ; shift 2;;
        -p|--pfrag) pfrag="$2" ; shift 2;;
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
#
#
#
#### processes to prepare raw data ###########
#
python prepare_trimming.py -r $raw_data -s $project --np $np
#
python prepare_mapping.py -s $project --index $index --picard $picard --tss $tss --np $np
#
python prepare_peakCalling.py -s $project --blist $blist --fa $fa --tss $tss --ref $ref --logq $logq
#
python prepare_countMatrix.py -s $project --fa $fa --np $np
#
python prepare_qualityControl.py -s $project --pfrag $pfrag --lib $frag
#
#
#
#
