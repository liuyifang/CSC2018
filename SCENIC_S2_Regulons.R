source("/BIGDATA1/pku_hkdeng_1/R/R-3.4.2/profile.R")

library(data.table)
library(Biobase)
library(AUCell)
library(SCENIC)
suppressWarnings(library(NMF, verbose=FALSE, warn.conflicts=FALSE, quietly=TRUE))

# To build a personalized report, update this working directory:
setwd("SCENIC_Mouse") # Or in the first chunk if running a notebook

## ----loadTFmodules, eval=TRUE--------------------------------------------
load("int/1.7_tfModules_withCorr.RData")

## ----chooseOrg-----------------------------------------------------------
library(RcisTarget.mm9.motifDatabases.20k)

# Get genes in databases:
data(mm9_500bpUpstream_motifRanking) # or 10kbp, they should have the same genes
genesInDatabase <- mm9_500bpUpstream_motifRanking@rankings$rn

# Get TFS in databases:
data(mm9_direct_motifAnnotation)
allTFs <- mm9_direct_motifAnnotation$allTFs

# Motif rankings (genes x motifs)
data(mm9_10kbpAroundTss_motifRanking)
motifRankings <- list()
motifRankings[["500bp"]] <- mm9_500bpUpstream_motifRanking
motifRankings[["10kbp"]] <- mm9_10kbpAroundTss_motifRanking

# Motif annotation (TFs)
direct_motifAnnotation <- mm9_direct_motifAnnotation
data(mm9_inferred_motifAnnotation) # optional
inferred_motifAnnotation <- mm9_inferred_motifAnnotation

# Remove genes missing from RcisTarget databases
#  (In case the input matrix wasn't already filtered)
tfModules_withCorr <- tfModules_withCorr[which(as.character(tfModules_withCorr$TF) %in% allTFs),]
geneInDb <- tfModules_withCorr$Target %in% motifRankings[["500bp"]]@rankings$rn
# Genes in co-expression modules not available in RcisTargetDatabases:
missingGenes <- sort(unique(tfModules_withCorr[which(!geneInDb),"Target"]))
missingGenes
tfModules_withCorr <- tfModules_withCorr[which(geneInDb),]

# Targets with positive correlation
tfModules_Selected <- tfModules_withCorr[which(tfModules_withCorr$corr==1),]

# Add a column with the geneSet name (TF_method)
tfModules_Selected <- cbind(tfModules_Selected, geneSetName=paste(tfModules_Selected$TF, tfModules_Selected$method, sep="_"))
head(tfModules_Selected)

# Split into tfModules (TF-modules, with several methods)
tfModules <- split(tfModules_Selected$Target, tfModules_Selected$geneSetName)

# Keep gene sets with at least 20 genes
tfModules <- tfModules[which(lengths(tfModules)>=20)]

# Add TF to the gene set (used in the following steps, careful if editing)
tfModules <- setNames(lapply(names(tfModules), function(gsn) {
    tf <- strsplit(gsn, "_")[[1]][1]
    unique(c(tf, tfModules[[gsn]]))
    }), names(tfModules))
save(tfModules, file="int/2.1_tfModules_forMotifEnrichmet.RData")

## ----statsTFmodules------------------------------------------------------
load("int/2.1_tfModules_forMotifEnrichmet.RData")
tfModulesSummary <- t(sapply(strsplit(names(tfModules), "_"), function(x) x[1:2]))
sort(table(tfModulesSummary[,2]))

## ----RcisTarget, eval=TRUE-----------------------------------------------
library(RcisTarget)
################################################################
# 1. Calculate motif enrichment for each TF-module

### 1.1 Calculate enrichment
motifs_AUC <- lapply(motifRankings, function(ranking) calcAUC(tfModules, ranking, aucMaxRank=0.01*nrow(ranking@rankings), nCores=8 * parallel::detectCores(logical = FALSE) - 2, verbose=FALSE))
save(motifs_AUC, file="int/2.2_motifs_AUC.RData") # renamed from: 2.2_motifs_AUC_500bp_10kbp.RData
# load(file="int/2.2_motifs_AUC.RData")

### 1.2 Conver to table, filter by NES & add the TFs to which the motif is annotated
# (For each database...)
motifEnrichment <- lapply(motifs_AUC, function(aucOutput)
{
  # Extract the TF of the gene-set name (i.e. MITF_w001):
  tf <- sapply(setNames(strsplit(rownames(aucOutput), "_"), rownames(aucOutput)), function(x) x[[1]])

  # Calculate NES and add motif annotation (provide tf in 'highlightTFs'):
  addMotifAnnotation(aucOutput, highlightTFs=tf, nesThreshold=3.0, digits=3,
                  motifAnnot_direct=direct_motifAnnotation,
                  motifAnnot_inferred=inferred_motifAnnotation)
})

# Merge both tables, adding a column that contains the 'motifDb'
motifEnrichment <- do.call(rbind, lapply(names(motifEnrichment), function(dbName){
  cbind(motifDb=dbName, motifEnrichment[[dbName]])
}))
save(motifEnrichment, file="int/2.3_motifEnrichment.RData")
cat("Number of motifs in the initial enrichment: ", nrow(motifEnrichment))

### 1.3 Keep only the motifs annotated to the initial TF
motifEnrichment_selfMotifs <- motifEnrichment[which(motifEnrichment$TFinDB != ""),, drop=FALSE]
save(motifEnrichment_selfMotifs, file="int/2.4_motifEnrichment_selfMotifs.RData")
cat("Number of motifs annotated to the initial TF: ", nrow(motifEnrichment_selfMotifs))
rm(motifEnrichment)

################################################################
# 2. Prune targets
library(data.table)

motifEnrichment_selfMotifs_wGenes <- lapply(names(motifRankings), function(motifDbName){
  addSignificantGenes(resultsTable=motifEnrichment_selfMotifs[motifDb==motifDbName],
                      geneSets=tfModules,
                      rankings=motifRankings[[motifDbName]],
                      maxRank=5000, method="aprox", nCores=8 * parallel::detectCores(logical = FALSE) - 2)
  })

motifEnrichment_selfMotifs_wGenes <- rbindlist(motifEnrichment_selfMotifs_wGenes)
save(motifEnrichment_selfMotifs_wGenes, file="int/2.5_motifEnrichment_selfMotifs_wGenes.RData")

# Save as text:
write.table(motifEnrichment_selfMotifs_wGenes, file="output/Step2_MotifEnrichment.tsv",
            sep="\t", quote=FALSE, row.names=FALSE)

## ----showSelfMotifs------------------------------------------------------
load("int/2.5_motifEnrichment_selfMotifs_wGenes.RData")
dim(motifEnrichment_selfMotifs_wGenes)
motifEnrichment_selfMotifs_wGenes[order(NES,decreasing=TRUE)][1:5,-"enrichedGenes", with=F]

## ----regulonTargetsInfo, eval=TRUE---------------------------------------
motifEnrichment.asIncidList <- apply(motifEnrichment_selfMotifs_wGenes, 1, function(oneMotifRow) {
  genes <- strsplit(oneMotifRow["enrichedGenes"], ";")[[1]]
  oneMotifRow <- data.frame(rbind(oneMotifRow), stringsAsFactors=FALSE)
  data.frame(oneMotifRow[rep(1, length(genes)),c("NES", "motif", "highlightedTFs", "TFinDB")], genes, stringsAsFactors = FALSE)
})
motifEnrichment.asIncidList <- rbindlist(motifEnrichment.asIncidList)
colnames(motifEnrichment.asIncidList) <- c("NES", "motif", "TF", "annot", "gene")
motifEnrichment.asIncidList <- data.frame(motifEnrichment.asIncidList, stringsAsFactors = FALSE)

# Get targets for each TF, but keep info about best motif/enrichment
# (directly annotated motifs are considered better)
regulonTargetsInfo <- lapply(split(motifEnrichment.asIncidList, motifEnrichment.asIncidList$TF), function(tfTargets){
  # print(unique(tfTargets$TF))
  tfTable <- as.data.frame(do.call(rbind, lapply(split(tfTargets, tfTargets$gene), function(enrOneGene){
    directAnnot <- "**" %in% enrOneGene$annot
    enrOneGeneByAnnot <- enrOneGene
    if(directAnnot) enrOneGeneByAnnot <- enrOneGeneByAnnot[which(enrOneGene$annot == "**"),]
    bestMotif <- which.max(enrOneGeneByAnnot$NES)

    cbind(TF=unique(enrOneGene$TF), gene=unique(enrOneGene$gene), nMotifs=nrow(enrOneGene),
          bestMotif=as.character(enrOneGeneByAnnot[bestMotif,"motif"]), NES=as.numeric(enrOneGeneByAnnot[bestMotif,"NES"]),
          directAnnot=directAnnot)
  })), stringsAsFactors=FALSE)
  tfTable[order(tfTable$NES, decreasing = TRUE),]
})
regulonTargetsInfo <- rbindlist(regulonTargetsInfo)
colnames(regulonTargetsInfo) <- c("TF", "gene", "nMotifs", "bestMotif", "NES", "directAnnot")

# Optional: Add Genie3 score
load("int/1.5_GENIE3_linkList.RData")
linkList <- linkList[which(linkList$weight>=0.001),]
rownames(linkList) <- paste(linkList$TF, linkList$Target,sep="__")
regulonTargetsInfo <- cbind(regulonTargetsInfo, Genie3Weight=linkList[paste(regulonTargetsInfo$TF, regulonTargetsInfo$gene,sep="__"),"weight"])

save(regulonTargetsInfo, file="int/2.6_regulonTargetsInfo.RData")
write.table(regulonTargetsInfo, file="output/Step2_regulonTargetsInfo.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

## ----tfRegulons, eval=TRUE-----------------------------------------------
regulonTargetsInfo_splitByAnnot <- split(regulonTargetsInfo, regulonTargetsInfo$directAnnot)
regulons <- sapply(split(regulonTargetsInfo_splitByAnnot[["TRUE"]], regulonTargetsInfo_splitByAnnot[["TRUE"]][,"TF"]), function(x) sort(as.character(unlist(x[,"gene"]))))
regulons_extended <- sapply(split(regulonTargetsInfo_splitByAnnot[["FALSE"]],regulonTargetsInfo_splitByAnnot[["FALSE"]][,"TF"]), function(x) unname(x[,"gene"]))
regulons_extended <- sapply(names(regulons_extended), function(tf) sort(unique(c(regulons[[tf]], regulons_extended[[tf]]))))
names(regulons_extended) <- paste(names(regulons_extended), "_extended", sep="")
regulons <- c(regulons, regulons_extended)
save(regulons, file="int/2.6_regulons_asGeneSet.RData")

## ----tfRegulonInfo-------------------------------------------------------
load("int/2.6_regulons_asGeneSet.RData")
# Number of regulons and summary of sizes:
length(regulons)
summary(lengths(regulons))

## ----tfRegulons_faster, eval=FALSE---------------------------------------
## selfMotifs_byTF <- split(motifEnrichment_selfMotifs_wGenes, motifEnrichment_selfMotifs_wGenes$highlightedTFs)
## regulons <- lapply(selfMotifs_byTF,
##                             function(x) unique(unlist(strsplit(x$enrichedGenes, ";"))))
## save(regulons, file="int/2.6_regulons_B_asGeneSet.RData")

## ----incidMats-----------------------------------------------------------
incidList <- melt(regulons)
incidMat <- table(incidList[,2], incidList[,1])
save(incidMat, file="int/2.6_regulons_asIncidMat.RData")
dim(incidMat)

## ----selfRegTfs----------------------------------------------------------
table(sapply(names(regulons), function(x) x %in% regulons[[x]]))

# ## ----exampleTfMotifs-----------------------------------------------------
# selTF <- "Sox2"
# subsetTable <- motifEnrichment_selfMotifs_wGenes[highlightedTFs %in% selTF][order(NES,decreasing=TRUE)][,-"enrichedGenes", with=F]
#
# subsetTable <- addLogo(subsetTable)
# library(DT)
# datatable(subsetTable, escape=FALSE, filter="top", options=list(pageLength=5))
#
# ## ----showOneTfEnrichment-------------------------------------------------
# geneSetName <- "Sox2_top50"
# motifDbName <- "10kbp"
# selectedMotifs <- subsetTable[geneSet==geneSetName & motifDb==motifDbName, motif]
# selectedMotifs <- selectedMotifs[1:3]

## ----signifGenesPlot-----------------------------------------------------
# pdf("int/2.8_RCC_selectedMotifs.pdf")
# par(mfrow=c(2,2))
# signifGenes_SelectedMotifs <- getSignificantGenes(tfModules[[geneSetName]],
#                                         motifRankings[[motifDbName]],
#                                         signifRankingNames=selectedMotifs,
#                                         plotCurve=TRUE, maxRank=5000, nCores=8 * parallel::detectCores(logical = FALSE) - 2,
#                                         genesFormat="geneList", method="aprox")
# dev.off()

# # Motif & number of genes:
# cbind(lengths(signifGenes_SelectedMotifs$enrichedGenes))

## ----sessionInfo---------------------------------------------------------
date()
sessionInfo()

