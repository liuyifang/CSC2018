source("/BIGDATA1/pku_hkdeng_1/R/R-3.4.2/profile.R")

## ----setup, echo=FALSE, message=FALSE, warning=FALSE---------------------
# Suppress loading messages when building the HTML
suppressPackageStartupMessages({
    library(Biobase)
    library(data.table)
    library(reshape2)
    library(GENIE3)
})

# Do not convert strings to factors (IMPORTANT! Specially if reading-in GENIE3 text output)
options(stringsAsFactors=FALSE)

# To build a personalized report, update this working directory:
setwd("SCENIC_Mouse") # Or in the first chunk if running a notebook

## ----loadGenie3, eval=TRUE-----------------------------------------------
library(GENIE3)
# Convert the weight matrix into links:
load("int/1.3_GENIE3_weightMatrix.RData")
linkList <- getLinkList(weightMatrix, threshold=0.001) # (slighly faster)
# linkList <- getLinkList(weightMatrix)
colnames(linkList) <- c("TF", "Target", "weight")
# order by weight
linkList <- linkList[order(linkList[,"weight"], decreasing=TRUE),]
save(linkList, file="int/1.5_GENIE3_linkList.RData")

## ----checkLinkList-------------------------------------------------------
load("int/1.5_GENIE3_linkList.RData")
dim(linkList)
head(linkList)

## ----weightStats---------------------------------------------------------
quantile(linkList$weight, probs=c(0.75, 0.90))
# plot(linkList$weight[1:1000000], type="l", ylim=c(0, max(linkList$weight)), main="Weight of the links",
#      ylab="Weight", xlab="Links sorted decreasingly")
# abline(h=0.001, col="blue") # Threshold
sum(linkList$weight>0.001)/nrow(linkList)

## ----filterLinkList------------------------------------------------------
linkList_001 <- linkList[which(linkList[,"weight"]>0.001),]
# Number of links over the threshold:
nrow(linkList_001)

## ----splitLinkList, eval=TRUE--------------------------------------------
tfModules <- list()

linkList_001$TF <- as.character(linkList_001$TF)
linkList_001$Target <- as.character(linkList_001$Target)

#### Create TF-modules:
# 1: Weight > 0.001 (filtered in previous step)
tfModules[["w001"]] <- split(linkList_001$Target, factor(linkList_001$TF))

# 2: Weight > 0.005
llminW <- linkList_001[which(linkList_001[,"weight"]>0.005),]
tfModules[["w005"]] <- split(llminW$Target, factor(llminW$TF))

# 3: Top 50 targets for each TF
# ("w001" should be ordered decreasingly by weight)
tfModules[["top50"]] <- lapply(tfModules[["w001"]], function(x) x[1:(min(length(x), 50))])

# 4-6: Top regulators per target
# (linkList_001 should be ordered by weight!)
linkList_001_byTarget <- split(linkList_001, factor(linkList_001$Target))
save(linkList_001_byTarget, file="int/1.5_linkList_001_byTarget.RData")

nTopTfs <- c(5, 10, 50)
nTopTfs <- setNames(nTopTfs, paste("top", nTopTfs, "perTarget", sep=""))

library(reshape2); library(data.table)
topTFsperTarget <- lapply(linkList_001_byTarget, function(llt) {
   nTFs <- nTopTfs[which(nTopTfs <= nrow(llt))]
   melt(lapply(nTFs, function(x) llt[1:x,"TF"]))
})
topTFsperTarget <- topTFsperTarget[which(!sapply(sapply(topTFsperTarget, nrow), is.null))]
topTFsperTarget.asDf <-  data.frame(rbindlist(topTFsperTarget, idcol=TRUE))
head(topTFsperTarget.asDf)
colnames(topTFsperTarget.asDf) <- c("Target", "TF", "method")

# Merge the all the gene-sets:
tfModules.melted <- melt(tfModules)
colnames(tfModules.melted) <- c("Target", "TF", "method")
tfModules <- rbind(tfModules.melted, topTFsperTarget.asDf)

save(tfModules, file="int/1.6_tfModules.RData")

## ----showTfModules-------------------------------------------------------
load("int/1.6_tfModules.RData")
# Basic counts:
rbind(nGeneSets=nrow(tfModules),
      nTFs=length(unique(tfModules$TF)),
      nTargets=length(unique(tfModules$Target)))

## ----addCorr, eval=TRUE--------------------------------------------------
load("int/1.4_corrMat.RData")
# Keep only correlation between TFs and potential targets
tfs <- unique(tfModules$TF)
corrMat <- corrMat[tfs,]

# Split TF modules according to correlation
tfModules_byTF <- split(tfModules, factor(tfModules$TF))
tfModules_withCorr_byTF <- lapply(tfModules_byTF, function(tfGeneSets)
{
    tf <- unique(tfGeneSets$TF)
    targets <- tfGeneSets$Target
    cbind(tfGeneSets, corr=c(as.numeric(corrMat[tf,targets] > 0.03) - as.numeric(corrMat[tf,targets] < -0.03)))
})
tfModules_withCorr <- data.frame(rbindlist(tfModules_withCorr_byTF))
save(tfModules_withCorr, file="int/1.7_tfModules_withCorr.RData")

## ----showTfModules_withCorr----------------------------------------------
load("int/1.7_tfModules_withCorr.RData")
head(tfModules_withCorr)
dim(tfModules_withCorr)

## ----sessionInfo---------------------------------------------------------
date()
sessionInfo()

