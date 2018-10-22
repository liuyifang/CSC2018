.libPaths("/lustre1/lch3000_pkuhpc/liuyf/R/R-3.4.1/library")

library(Biobase)
library(reshape2)
library(ggplot2)
library(igraph)
library(monocle)

load("exprs_gbm.Robj")
pd <- read.csv("pData_to_monocle.csv", row.names = 1)
fd <- read.csv("fData_to_monocle.csv", row.names = 1)

AGG <- newCellDataSet(exprs_gbm,
                      phenoData = new("AnnotatedDataFrame", pd),
                      featureData = new("AnnotatedDataFrame", fd),
                      lowerDetectionLimit=0.5,
                      expressionFamily=negbinomial.size())

## ----show_mRNA_totals, eval = TRUE, fig.width = 4, fig.height = 2, fig.align="center"----
pData(AGG)$Total_mRNAs <- Matrix::colSums(exprs(AGG))
upper_bound <- 10^(mean(log10(pData(AGG)$Total_mRNAs)) + 2*sd(log10(pData(AGG)$Total_mRNAs)))
lower_bound <- 10^(mean(log10(pData(AGG)$Total_mRNAs)) - 2*sd(log10(pData(AGG)$Total_mRNAs)))
upper_bound
lower_bound

# check
length(unique(pData(AGG)$Batch)) # 13 batch (include ES)
length(pData(AGG)$Barcode) # total cell is 39342

sum(pData(AGG)$Batch == "MEF")
sum(pData(AGG)$Batch == "MEF" &
      pData(AGG)$Total_mRNAs <= lower_bound)
sum(pData(AGG)$Batch == "MEF" &
      pData(AGG)$Total_mRNAs >= upper_bound)

sum(pData(AGG)$Batch == "XEN")
sum(pData(AGG)$Batch == "XEN" &
      pData(AGG)$Total_mRNAs <= lower_bound)
sum(pData(AGG)$Batch == "XEN" &
      pData(AGG)$Total_mRNAs >= upper_bound)

sum(pData(AGG)$Batch == "ES")
sum(pData(AGG)$Batch == "ES" &
      pData(AGG)$Total_mRNAs <= lower_bound)
sum(pData(AGG)$Batch == "ES" &
      pData(AGG)$Total_mRNAs >= upper_bound)

sum(pData(AGG)$Batch == "SI5D")
sum(pData(AGG)$Batch == "SI5D" &
      pData(AGG)$Total_mRNAs <= lower_bound)
sum(pData(AGG)$Batch == "SI5D" &
      pData(AGG)$Total_mRNAs >= upper_bound)

sum(pData(AGG)$Batch == "SI12D")
sum(pData(AGG)$Batch == "SI12D" &
      pData(AGG)$Total_mRNAs <= lower_bound)
sum(pData(AGG)$Batch == "SI12D" &
      pData(AGG)$Total_mRNAs >= upper_bound)

sum(pData(AGG)$Batch == "SII8D")
sum(pData(AGG)$Batch == "SII8D" &
      pData(AGG)$Total_mRNAs <= lower_bound)
sum(pData(AGG)$Batch == "SII8D" &
      pData(AGG)$Total_mRNAs >= upper_bound)

sum(pData(AGG)$Batch == "SII12D")
sum(pData(AGG)$Batch == "SII12D" &
      pData(AGG)$Total_mRNAs <= lower_bound)
sum(pData(AGG)$Batch == "SII12D" &
      pData(AGG)$Total_mRNAs >= upper_bound)

sum(pData(AGG)$Batch == "SIII3D")
sum(pData(AGG)$Batch == "SIII3D" &
      pData(AGG)$Total_mRNAs <= lower_bound)
sum(pData(AGG)$Batch == "SIII3D" &
      pData(AGG)$Total_mRNAs >= upper_bound)

sum(pData(AGG)$Batch == "SIII6D")
sum(pData(AGG)$Batch == "SIII6D" &
      pData(AGG)$Total_mRNAs <= lower_bound)
sum(pData(AGG)$Batch == "SIII6D" &
      pData(AGG)$Total_mRNAs >= upper_bound)

sum(pData(AGG)$Batch == "SIII8D")
sum(pData(AGG)$Batch == "SIII8D" &
      pData(AGG)$Total_mRNAs <= lower_bound)
sum(pData(AGG)$Batch == "SIII8D" &
      pData(AGG)$Total_mRNAs >= upper_bound)

sum(pData(AGG)$Batch == "SIII10D")
sum(pData(AGG)$Batch == "SIII10D" &
      pData(AGG)$Total_mRNAs <= lower_bound)
sum(pData(AGG)$Batch == "SIII10D" &
      pData(AGG)$Total_mRNAs >= upper_bound)

sum(pData(AGG)$Batch == "SIII15D")
sum(pData(AGG)$Batch == "SIII15D" &
      pData(AGG)$Total_mRNAs <= lower_bound)
sum(pData(AGG)$Batch == "SIII15D" &
      pData(AGG)$Total_mRNAs >= upper_bound)

sum(pData(AGG)$Batch == "SIII21D")
sum(pData(AGG)$Batch == "SIII21D" &
      pData(AGG)$Total_mRNAs <= lower_bound)
sum(pData(AGG)$Batch == "SIII21D" &
      pData(AGG)$Total_mRNAs >= upper_bound)

# filter out lower, upper & ES
AGG <- AGG[,pData(AGG)$Total_mRNAs > lower_bound &
             pData(AGG)$Total_mRNAs < upper_bound]
AGG <- AGG[,pData(AGG)$Batch != "ES"]
pData(AGG)$Batch <- factor(pData(AGG)$Batch,
                           levels = c("MEF", "SI5D", "SI12D", "XEN",
                                      "SII8D", "SII12D",
                                      "SIII3D", "SIII6D", "SIII8D", "SIII10D", "SIII15D", "SIII21D"))

# check
unique(pData(AGG)$Batch)
length(unique(pData(AGG)$Batch))
sum(pData(AGG)$Batch %in% c("ES"))
# all batch
sum(pData(AGG)$Batch %in% c("MEF","XEN","ES","SI5D","SI12D","SII8D","SII12D","SIII3D","SIII6D","SIII8D","SIII10D","SIII15D","SIII21D"))
# without ES
sum(pData(AGG)$Batch %in% c("MEF","XEN","SI5D","SI12D","SII8D","SII12D","SIII3D","SIII6D","SIII8D","SIII10D","SIII15D","SIII21D"))

monocle_raw <- exprs(AGG)
monocle_raw <- as.matrix(monocle_raw)
dim(monocle_raw)

# Normalize
x <- monocle_raw
monocle_norm <- sweep(x,2,colSums(x),'/')*median(colSums(x))

# Convert to data frame
monocle_raw <- as.data.frame(monocle_raw)
monocle_norm <- as.data.frame(monocle_norm)

## ----estimate_size_and_dispersion, eval=TRUE-----------------------------
AGG <- estimateSizeFactors(AGG)
AGG <- estimateDispersions(AGG, cores = 6)
AGG <- detectGenes(AGG, min_expr = 0.1)

monocle_pData <- pData(AGG)
monocle_fData <- fData(AGG)

write.csv(monocle_pData, file = "monocle_pData.csv", quote = F)
write.csv(monocle_fData, file = "monocle_fData.csv", quote = F)
save(monocle_raw, monocle_norm, file = "monocle_raw_norm.Robj")
save(AGG, file = "AGG_estimate_size_and_dispersion.Robj")
sessionInfo()
