#
library('monocle')
args <- commandArgs(T)
#
####### build CellDataSet type object ##########
#
expr_matrix <- read.csv(args[1], header=TRUE, row.names=1, sep=',', check.names=FALSE)
cells <- read.delim(args[2], row.names=1)
genes <- read.delim(args[3], row.names=1)
pd <- new("AnnotatedDataFrame", data=cells)
fd <- new("AnnotatedDataFrame", data=genes)
HSMM <- newCellDataSet(as.matrix(expr_matrix), phenoData=pd, featureData=fd, expressionFamily=negbinomial.size())
#HSMM <- detectGenes(HSMM, min_expr=0.001)
#expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= 4))
#HSMM <- HSMM[expressed_genes, ]
HSMM <- estimateSizeFactors(HSMM)
#
#
####### reduce dimension ##########
#
HSMM <- reduceDimension(HSMM, max_components=3, norm_method='none', num_dim=20, reduction_method='DDRTree', verbose=T)
HSMM <- orderCells(HSMM)
write.csv(HSMM@reducedDimS, args[4])
write.csv(HSMM@phenoData@data, args[5])
#
#
#
