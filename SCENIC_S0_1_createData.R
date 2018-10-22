library(cellrangerRkit)

pipestance_path <- "rawData"
gbm <- load_cellranger_matrix(pipestance_path)
exprSM <- exprs(gbm)
pData_AGG <- pData(gbm)
fData_AGG <- fData(gbm)

exprMatrix <- as.matrix(exprSM)
rm(exprSM)
dim(exprMatrix)
exprMatrix[1:3,1:3]

geneNames <- as.character(fData_AGG$symbol)
rownames(exprMatrix) <- geneNames
exprMatrix[1:3,1:3]
exprMatrix <- exprMatrix[unique(rownames(exprMatrix)),] # Remove duplicated rows
dim(exprMatrix)

save(exprMatrix, fData_AGG, pData_AGG, file = "data/S0_createData.RData")

date()
sessionInfo()
