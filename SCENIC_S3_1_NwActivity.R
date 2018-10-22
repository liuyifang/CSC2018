# source("~/profile.R")
source("/BIGDATA1/pku_hkdeng_1/R/R-3.4.2/profile.R")

library(data.table)
library(Biobase)
library(AUCell)
library(SCENIC)
suppressWarnings(library(NMF, verbose=FALSE, warn.conflicts=FALSE, quietly=TRUE))

setwd("SCENIC_Mouse")

## ----loadEset------------------------------------------------------------
load("data/esetMouse.RData")
exprMat <- exprs(esetMouse)
dim(exprMat)

# load("data/colVars.RData")
# cellInfo <- pData(esetMouse)[,names(colVars), drop=F]

## ----loadRegulons--------------------------------------------------------
load("int/2.6_regulons_asGeneSet.RData")
regulons <- regulons[order(lengths(regulons), decreasing=TRUE)]
regulons <- regulons[lengths(regulons)>=10]

# Add the TF & rename
regulons <- setNames(lapply(names(regulons), function(tf) sort(unique(c(gsub("_extended", "", tf), regulons[[tf]])))), names(regulons))
names(regulons) <- paste(names(regulons), " (",lengths(regulons), "g)", sep="")
save(regulons, file="int/3.0_regulons_forAUCell.RData")
length(regulons)
cbind(names(regulons)[1:10])

## ----runAUCell, eval=TRUE------------------------------------------------
# 1. Create rankings
# increase Plots window
aucellRankings <- AUCell.buildRankings(exprMat, nCores=8 * parallel::detectCores(logical = FALSE) - 2, plotStats=TRUE)
# abline(v=aucellRankings@nGenesDetected["1%"], col="skyblue3", lwd=5, lty=3)
save(aucellRankings, file="int/3.1_aucellRankings.RData")

# 2. Calculate AUC
regulonAUC <- AUCell.calcAUC(regulons, aucellRankings, aucMaxRank=aucellRankings@nGenesDetected["1%"], nCores=8 * parallel::detectCores(logical = FALSE) - 2)

## ----sortAUC-------------------------------------------------------------
# Order the modules by similarity, for easier exploration in the upcoming steps & save
variableRegulons <- names(which(apply(getAuc(regulonAUC), 1, sd) > 0))
reguDist <-as.dist(1-cor(t(getAuc(regulonAUC)[variableRegulons,]), method="spear"))
reguClust <- hclust(reguDist, method="ward.D2")
regulonClusters <- setNames(dynamicTreeCut::cutreeDynamic(reguClust, distM=as.matrix(reguDist), verbose = FALSE), reguClust$labels)
regulonOrder <- reguClust$labels[reguClust$order]
regulonOrder <- regulonOrder[order(regulonClusters[regulonOrder], decreasing = TRUE)]
regulonAUC@matrix <- regulonAUC@matrix[regulonOrder,]
save(regulonAUC, file="int/3.2_regulonAUC.RData")

load("int/3.0_regulons_forAUCell.RData")
load("int/3.1_aucellRankings.RData")
load("int/3.2_regulonAUC.RData")

## ----tSNE_AUC, eval=TRUE-------------------------------------------------
# (It is recommended to try different perplexity values)
regulonAUC_subset <- subset(regulonAUC, onlyNonDirectExtended(rownames(regulonAUC)))

# PCA-based t-SNE
set.seed(123)
tsneAUC <- Rtsne::Rtsne(t(getAuc(regulonAUC_subset)), initial_dims=10, perplexity=10)
rownames(tsneAUC$Y) <- colnames(regulonAUC_subset)
colnames(tsneAUC$Y) <- c("tsne1", "tsne2")
save(tsneAUC, file="int/3.3_tsneRegulonAUC_PCA.RData")

# Alternative: Distance-based t-SNE:
corDist <- as.dist(1-cor(getAuc(regulonAUC_subset)))
set.seed(123)
tsneAUC <- Rtsne::Rtsne(corDist, is_distance=TRUE, perplexity=10)
rownames(tsneAUC$Y) <- labels(corDist)
colnames(tsneAUC$Y) <- c("tsne1", "tsne2")
save(tsneAUC, file="int/3.3_tsneRegulonAUC_Dist.RData")

## ----tSNE_plot, fig.height=5, fig.width=10-------------------------------
load("int/3.3_tsneRegulonAUC_PCA.RData")
tSNE <- tsneAUC$Y
# par(mfrow=c(1,2))

# tSNE2 <- read.csv("../tSNE2.csv", row.names = 1)
# tSNE2 <- as.matrix(tSNE2)

# # Number of genes detected:
# nGenesPerCell <- apply(exprMat, 2, function(x) sum(x>0))
# colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))
# cellColorNgenes <- setNames(adjustcolor(colorPal(10), alpha=.8)[as.numeric(cut(nGenesPerCell,breaks=10, right=F,include.lowest=T))], names(nGenesPerCell))

# plot(tSNE, col=cellColorNgenes[rownames(tSNE)], pch=16, main="nGenes", sub="t-SNE on regulon activity (AUC)")

# # Other known properties:
# for(varName in names(colVars))
# {
#   cellColor <- setNames(colVars[[varName]][cellInfo[,varName]], rownames(cellInfo))
#   plot(tSNE, col=cellColor[rownames(tSNE)], pch=16, main=varName, sub="t-SNE on regulon activity (AUC)")
# }

## ----histogramsTsne, eval=TRUE-------------------------------------------
# Cairo::CairoPDF("output/Step3_RegulonActivity_AUCtSNE.pdf", width=20, height=5)
# par(mfrow=c(1,4))

# # tSNE (colored by number of genes detected per cell)
# plot(tSNE, col=cellColorNgenes[rownames(tSNE)], pch=16, main="nGenes", sub="t-SNE on regulon activity (AUC)")
# plot(tSNE, col=cellColor[rownames(tSNE)], pch=16, main=varName, sub="t-SNE on regulon activity (AUC)")
# plot.new(); plot.new()

# Plot module activity, thresholds & assignment:
cells_AUCellThresholds <- plot_aucTsne(tSNE=tSNE, exprMat=exprMat, regulonAUC=regulonAUC, alphaOff=0.1)
# dev.off()
save(cells_AUCellThresholds, file="int/3.4_AUCellThresholds.RData")

## ----htmlPreview, echo=FALSE, fig.height=6, fig.width=7, eval=TRUE-------
x <- as.data.frame(rownames(regulonAUC))
regOrder <- c("Sox2 (142g)", "Nanog (200g)", "Pou5f1 (32g)")
# tiff(file=paste0("SIII_tSNE2_", today,".tiff"), height = 6, width = 8, units = 'in', res = 300, compression = 'none')
# par(mfrow=c(3,4))
# tmp <- plot_aucTsne(tSNE=tSNE2, exprMat=exprMat, regulonAUC=regulonAUC[regOrder,], alphaOff=0.1, cex=.8)
# dev.off()

## ----thresholds2edit-----------------------------------------------------
load("int/3.4_AUCellThresholds.RData")

# Get cells assigned to each regulon
regulonsCells <- lapply(cells_AUCellThresholds, function(x) x$assignment)

### Save threshold info as text (e.g. to edit/modify...)
trhAssignment <- sapply(cells_AUCellThresholds, function(x) unname(x$aucThr$selected))
commentsThresholds <- sapply(cells_AUCellThresholds, function(x) unname(x$aucThr$comment))

table2edit <- cbind(regulon=names(trhAssignment),
                    threshold=trhAssignment,
                    nCellsAssigned=lengths(regulonsCells)[names(trhAssignment)],
                    AUCellComment=commentsThresholds,
                    nGenes=gsub("[\\(g\\)]", "", regmatches(names(cells_AUCellThresholds), gregexpr("\\(.*?\\)", names(cells_AUCellThresholds)))),
                    clusteringOrder=1:length(trhAssignment),
                    clusterGroup=regulonClusters[names(trhAssignment)],
                    onlyNonDirectExtended=(names(trhAssignment) %in% onlyNonDirectExtended(names(trhAssignment))),
                    personalNotes="")
write.table(table2edit, file="int/3.5_1_AUCellThresholds.txt", row.names=F, quote=F, sep="\t")

## ----sessionInfo---------------------------------------------------------
date()
sessionInfo()


