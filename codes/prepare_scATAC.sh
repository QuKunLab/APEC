#!/usr/bin/bash
#
#
#
#
#### input parameters defined by users ##########
#
data=../Example23_LY/data
work=../Example23_LY/work
peak=../Example23_LY/peak
matrix=../Example23_LY/matrix
result=../Example23_LY/result
figure=../Example23_LY/figure
np=40
picard=../reference/picard.jar
ref=mm10
fa=../reference/mm10_chr.fa
index=../reference/mm10
tss=../reference/mm10_refseq_genes_TSS.txt
blist=../reference/mm10_blacklist.BIN.bed
meme=../reference/JASPAR2018_Mus_musculus.meme
bgfile=../reference/tier1_markov1.norc.txt
TSSfragment=0.1
readsPerCell=2000
#
#
#
#
python prepare_trimming.py -d $data --np $np
#
python prepare_mapping.py -d $data -w $work --index $index --picard $picard --tss $tss --np $np
#
python prepare_peakCalling.py -w $work -p $peak --blist $blist --fa $fa --tss $tss --ref $ref
#
python prepare_countMatrix.py -w $work -p $peak -m $matrix --cellinfo $data/cell_info.csv --fa $fa --bg $bgfile --meme $meme --np $np
#
python prepare_qualityControl.py -m $matrix -p $peak -f $figure --frag $TSSfragment --lib $readsPerCell
#
#
