# source("~/profile.R")
source("/BIGDATA1/pku_hkdeng_1/R/R-3.4.2/profile.R")

setwd("SCENIC_Mouse") # Or in the first chunk if running a notebook

load("int/1.1_exprMatrix_filtered.RData")
corrMat <- cor(t(exprMatrix_filtered), method="spearman")
save(corrMat, file="int/1.4_corrMat.RData")

date()
sessionInfo()
