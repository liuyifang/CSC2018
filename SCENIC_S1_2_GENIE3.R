source("/BIGDATA1/pku_hkdeng_1/R/R-3.4.2/profile.R")

library(foreach)
library(GENIE3)

setwd("SCENIC_Mouse")

load("int/1.1_exprMatrix_filtered.RData")
exprMatrix_filtered <- log2(exprMatrix_filtered+1)
load("int/1.2_inputTFs.RData")

set.seed(1011)
weightMatrix <- GENIE3(exprMatrix_filtered, regulators = inputTFs, nCores = 8 * parallel::detectCores(logical = FALSE) - 2)
save(weightMatrix, file = paste0("int/1.3_GENIE3_weightMatrix.RData"))

date()
sessionInfo()
