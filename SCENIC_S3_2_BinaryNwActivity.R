options(stringsAsFactors=FALSE)
today <- format(Sys.time(), "%y%m%d")

library(data.table)
library(Biobase)
library(AUCell)
library(SCENIC)
suppressWarnings(library(NMF, verbose=FALSE, warn.conflicts=FALSE, quietly=TRUE))

setwd("SCENIC_Mouse")

## ----loadAssignmentB-----------------------------------------------------
load("int/3.4_AUCellThresholds.RData")

## ----createBinaryMatrix--------------------------------------------------
# Get cells assigned to each regulon
regulonsCells <- lapply(cells_AUCellThresholds, function(x) x$assignment)
length(regulonsCells)

# Conver to matrix (regulons with zero assigned cells are lost)
regulonActivity <- reshape2::melt(regulonsCells)
binaryRegulonActivity <- t(table(regulonActivity[,1], regulonActivity[,2]))
class(binaryRegulonActivity) <- "matrix"
save(binaryRegulonActivity, file="int/3.6_binaryRegulonActivity.RData")

dim(binaryRegulonActivity)
binaryRegulonActivity[1:10,1:3]

## ------------------------------------------------------------------------
binaryRegulonActivity_nonDupl <- binaryRegulonActivity[which(rownames(binaryRegulonActivity) %in% onlyNonDirectExtended(rownames(binaryRegulonActivity))),]
save(binaryRegulonActivity_nonDupl, file="int/3.7_binaryRegulonActivity_nonDupl.RData")

## ------------------------------------------------------------------------
cbind(nCellsOn=sort(rowSums(binaryRegulonActivity), decreasing=TRUE)[1:15])
summary(rowSums(binaryRegulonActivity))
nCellsOn <- sort(rowSums(binaryRegulonActivity), decreasing=TRUE)
nCellsOn <- as.data.frame(nCellsOn)
colnames(nCellsOn) <- "nCellsOn"
fname <- paste0("regulon_activity_percentage_", today, ".RData")
save(nCellsOn, file = fname)

# ## ----boxplots, fig.height=4, fig.width=8---------------------------------
# par(mfrow=c(1,2))
# boxplot(rowSums(binaryRegulonActivity_nonDupl), main="nCells per regulon",
#         sub='number of cells \nthat have the regulon active',
#         col="darkolivegreen1", border="#001100", lwd=2, frame=FALSE)
# boxplot(colSums(binaryRegulonActivity_nonDupl), main="nRegulons per Cell",
#         sub='number of regulons \nactive per cell',
#         col="darkolivegreen1", border="#001100", lwd=2, frame=FALSE)

## ----loadColors----------------------------------------------------------
load("data/colVars.RData")
load("data/esetMouse.RData")
cellInfo <- pData(esetMouse)[, c("Sample", "Cell"), drop=F]
minCells <- ncol(esetMouse) * .01

## ----regulonCorrelation--------------------------------------------------
load("int/3.6_binaryRegulonActivity.RData")
load("int/3.7_binaryRegulonActivity_nonDupl.RData")

x <- as.data.frame(rowSums(binaryRegulonActivity_nonDupl) /  ncol(esetMouse))
x$regulon <- row.names(x)

regulonSelection <- list()

# All regulons.
regulonSelection[["All regulons \n (including duplicated regulons)"]] <- rownames(binaryRegulonActivity)

# Active in > 1% cells
regMinCells <- names(which(rowSums(binaryRegulonActivity_nonDupl) > minCells))
regulonSelection[["Regulons active in more than 1% of cells"]] <- regMinCells

# Correlation across regulons (based on binary cell activity)
reguCor <- cor(t(binaryRegulonActivity_nonDupl[regMinCells,]))
reguCor[which(is.na(reguCor))] <- 0
diag(reguCor) <- 0

# Regulons that co-ocurr in similar cells. If a regulon is relevant by itself it will not be shown, also check the regulons ignored.
corrRegs <- names(which(rowSums(abs(reguCor) > 0.30) > 0))
regulonSelection[["Regulons with any other regulon correlated\n with abs(cor)>0.30 \n(and active in at least 1% of cells)"]]  <- corrRegs

missingRegs <- rownames(binaryRegulonActivity_nonDupl)[which(!rownames(binaryRegulonActivity_nonDupl) %in% corrRegs)]
regulonSelection[["Regulons no other regulons correlated\n with abs(cor)>0.30 \n or active in fewer than 1% of cells"]]  <- missingRegs

save(regulonSelection,file="int/3.8_regulonSelections.RData")

## Set regulon order (for plotting)
binaryRegulonOrder <- hclust(as.dist(1-reguCor[corrRegs,corrRegs]))
binaryRegulonOrder <- binaryRegulonOrder$labels[binaryRegulonOrder$order]
save(binaryRegulonOrder,file="int/3.9_binaryRegulonOrder.RData")

## ----heatmapPlot, eval=TRUE----------------------------------------------
for(i in seq_len(length(regulonSelection)))
{
    selRegs <- names(regulonSelection)[i]
    if(length(regulonSelection[[selRegs]])>1)
    {
        binaryMat <- binaryRegulonActivity[regulonSelection[[selRegs]],,drop=FALSE]
        NMF::aheatmap(binaryMat, scale="none", revC=TRUE, main=selRegs,
                    annCol=cellInfo[colnames(binaryMat),, drop=FALSE],
                    annColor=colVars,
                    color = c("white", "black"),
                    filename=paste0("output/Step3_3.3_binaryRegulonActivity_Heatmap_",i,".pdf"))
    }
}

## ----htmlPreview, echo=FALSE, fig.height=7, fig.width=7, eval=TRUE-------
# selRegs <- names(regulonSelection)[3]
# binaryMat <- binaryRegulonActivity[regulonSelection[[selRegs]],,drop=FALSE]
# NMF::aheatmap(binaryMat, scale="none", revC=TRUE, main=selRegs,
#                 annCol=cellInfo[colnames(binaryMat),, drop=FALSE],
#                 annColor=colVars,
#                 color = c("white", "black"))

fname <- paste0("output/regulon1_", today, ".RData")
selRegs <- names(regulonSelection)[1]
binaryMat <- binaryRegulonActivity[regulonSelection[[selRegs]],,drop=FALSE]
save(binaryMat, file = fname)

fname <- paste0("output/regulon2_", today, ".RData")
selRegs <- names(regulonSelection)[2]
binaryMat <- binaryRegulonActivity[regulonSelection[[selRegs]],,drop=FALSE]
save(binaryMat, file = fname)

fname <- paste0("output/regulon3_", today, ".RData")
selRegs <- names(regulonSelection)[3]
binaryMat <- binaryRegulonActivity[regulonSelection[[selRegs]],,drop=FALSE]
save(binaryMat, file = fname)

fname <- paste0("output/regulon4_", today, ".RData")
selRegs <- names(regulonSelection)[4]
binaryMat <- binaryRegulonActivity[regulonSelection[[selRegs]],,drop=FALSE]
save(binaryMat, file = fname)

## ----sessionInfo---------------------------------------------------------
date()
sessionInfo()
