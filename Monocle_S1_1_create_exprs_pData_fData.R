update.packages(ask = FALSE)
source("http://bioconductor.org/biocLite.R")
biocLite()

library(cellrangerRkit)

pipestance_path <- "/Users/m/tmp/data/_AGG11_170901/AGG11_mapped"
gbm <- load_cellranger_matrix(pipestance_path)
exprs_gbm <- exprs(gbm)
save(exprs_gbm, file = "exprs_gbm.Robj")
x <- pData(gbm)
write.csv(x, file = "pData_raw.csv", quote = F)
y <- fData(gbm)
write.csv(y, file = "fData_raw.csv", quote = F)

sessionInfo()
