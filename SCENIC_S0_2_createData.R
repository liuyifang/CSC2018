source("/BIGDATA1/pku_hkdeng_1/R/R-3.4.2/profile.R")

library(Biobase)

# load("data/S0_createData.RData")
load("/BIGDATA1/pku_hkdeng_1/CapitalBio/data/S0_createData.RData")

clusters <- read.csv("clusters_anno.csv", row.names = 1)
# cellLabels <- pData_AGG[, "Batch", drop = FALSE]
cellLabels <- clusters[, "Cell", drop = FALSE]
colnames(cellLabels) <- "level1class"

# filter by Cell
subGroup <- c("EC",
              "Fibro",
              "Myloid")
cellLabelsSubGroup <- subset(clusters, Cell %in% subGroup)

# # filter by Barcode
# subGroup <- read.csv("successfulFailedTerminal.csv")
# barcodeSubGroup <- subGroup$Barcode
# cellLabelsSubGroup <- cellLabels[barcodeSubGroup, , drop = FALSE]

barcodeSubGroup <- row.names(cellLabelsSubGroup)
exprMatrixSubGroup <- exprMatrix[, barcodeSubGroup]
dim(cellLabelsSubGroup)
dim(exprMatrixSubGroup)
cellLabelsSubGroup[1:3, ]
exprMatrixSubGroup[1:3, 1:3]

dir.create("SCENIC_Mouse")
setwd("SCENIC_Mouse") # Or in the first chunk if running a notebook
dir.create("int")
dir.create("output")
dir.create("data")

esetMouse <- new("ExpressionSet",
                      exprs = exprMatrixSubGroup,
                      phenoData = new("AnnotatedDataFrame",
                                    data = data.frame(cellLabelsSubGroup[colnames(exprMatrixSubGroup),, drop = FALSE])))
save(esetMouse, file="data/esetMouse.RData")

date()
sessionInfo()
