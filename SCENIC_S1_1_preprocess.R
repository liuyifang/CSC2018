# source("~/profile.R")
source("/BIGDATA1/pku_hkdeng_1/R/R-3.4.2/profile.R")

library(viridis)
library(Biobase)
library(RcisTarget.mm9.motifDatabases.20k)

setwd("SCENIC_Mouse") # Or in the first chunk if running a notebook

load("data/esetMouse.RData")

exprMat <- exprs(esetMouse)
dim(exprMat)
exprMat[1:3, 1:3]

cellInfo <- pData(esetMouse)[colnames(exprMat), "Sample", drop=F]

Batch_vec <- c("Ctrl", "Drug")
Batch_cols <- viridis(2, option = "viridis")

# Color to assign to the variables (same format as for NMF::aheatmap)
colVars <- list(level1class = setNames(Batch_cols, Batch_vec))
save(colVars, file="data/colVars.RData")
# plot.new(); legend(0,1, fill=colVars$level1class, legend=names(colVars$level1class))

# Get genes in databases:
data(mm9_500bpUpstream_motifRanking) # or 10kbp, they should have the same genes
genesInDatabase <- mm9_500bpUpstream_motifRanking@rankings$rn

# Get TFS in databases:
data(mm9_direct_motifAnnotation)
allTFs <- mm9_direct_motifAnnotation$allTFs

nCellsPerGene <- apply(exprMat, 1, function(x) sum(x>0))
nCountsPerGene <- apply(exprMat, 1, sum)

summary(nCellsPerGene)
summary(nCountsPerGene)

max(exprMat)

sum(exprMat>0) / sum(exprMat==0)

minReads <- ncol(exprMat)*.01
length(nCountsPerGene)
genesLeft_minReads <- names(nCountsPerGene)[which(nCountsPerGene > minReads)]
length(genesLeft_minReads)

minSamples <- ncol(exprMat)*.01
nCellsPerGene2 <- nCellsPerGene[genesLeft_minReads]
genesLeft_minCells <- names(nCellsPerGene2)[which(nCellsPerGene2 > minSamples)]
length(genesLeft_minCells)

genesLeft_minCells_inDatabases <- genesLeft_minCells[which(genesLeft_minCells %in% genesInDatabase)]
length(genesLeft_minCells_inDatabases)

exprMatrix_filtered <- exprMat[genesLeft_minCells_inDatabases, ]
dim(exprMatrix_filtered)
exprMatrix_filtered[1:3, 1:3]
save(exprMatrix_filtered, file="int/1.1_exprMatrix_filtered.RData")

# # Check whether any relevant gene / potential gene of interest is missing:
# interestingGenes <- c("Thy1","Sox17","Sall4","Pou5f1","Nanog","Zscan4c","Prrx1","Twist2","Zeb2","Fbn1","Gata4","Gata6","Foxa2","Epcam","Zscan4d","Zscan4f","Gm4340","Tcstv1","Esrrb","Dppa2","Pecam1","Fgf4","Sox2","Dppa3","Dppa4","Zfp42")
# interestingGenes[which(!interestingGenes %in% rownames(exprMatrix_filtered))]

rm(exprMat)

inputTFs <- allTFs[allTFs %in% rownames(exprMatrix_filtered)]
save(inputTFs, file="int/1.2_inputTFs.RData")

c(allTFs=length(allTFs), inputTFs=length(inputTFs))

date()
sessionInfo()
