monocle_tSNE <- read.csv("t-SNE.csv", row.names = 1)
rownames(monocle_tSNE) <- monocle_tSNE$Barcode
m_tSNE <- monocle_tSNE[,c("data_dim_1","data_dim_2","Batch")]

is_MEF_stage <- m_tSNE[m_tSNE$Batch %in% "MEF",]
tiff(file = paste0("is_MEF_stage","",".tiff"), width = 7, height = 7, units = 'in', res = 300, compression = 'none')
par(mar = c(3, 3, 3, 3))
plot(m_tSNE[,c("data_dim_1","data_dim_2")], pch = 20, cex = 0.3, col = "#F2F2F2", frame.plot = F, xaxt = 'n', yaxt = 'n', ann = FALSE)
for ( k in 1:length(is_MEF_stage$Batch) ){
  points(is_MEF_stage[k,1], is_MEF_stage[k,2], col = "#440154FF", pch = 20, cex = 0.3)
}
dev.off()

is_SI5D_stage <- m_tSNE[m_tSNE$Batch %in% "SI5D",]
tiff(file = paste0("is_SI5D_stage","",".tiff"), width = 7, height = 7, units = 'in', res = 300, compression = 'none')
par(mar = c(3, 3, 3, 3))
plot(m_tSNE[,c("data_dim_1","data_dim_2")], pch = 20, cex = 0.3, col = "#F2F2F2", frame.plot = F, xaxt = 'n', yaxt = 'n', ann = FALSE)
for ( k in 1:length(is_SI5D_stage$Batch) ){
  points(is_SI5D_stage[k,1], is_SI5D_stage[k,2], col = "#482173FF", pch = 20, cex = 0.3)
}
dev.off()

is_SI12D_stage <- m_tSNE[m_tSNE$Batch %in% "SI12D",]
tiff(file = paste0("is_SI12D_stage","",".tiff"), width = 7, height = 7, units = 'in', res = 300, compression = 'none')
par(mar = c(3, 3, 3, 3))
plot(m_tSNE[,c("data_dim_1","data_dim_2")], pch = 20, cex = 0.3, col = "#F2F2F2", frame.plot = F, xaxt = 'n', yaxt = 'n', ann = FALSE)
for ( k in 1:length(is_SI12D_stage$Batch) ){
  points(is_SI12D_stage[k,1], is_SI12D_stage[k,2], col = "#433E85FF", pch = 20, cex = 0.3)
}
dev.off()

is_XEN_stage <- m_tSNE[m_tSNE$Batch %in% "XEN",]
tiff(file = paste0("is_XEN_stage","",".tiff"), width = 7, height = 7, units = 'in', res = 300, compression = 'none')
par(mar = c(3, 3, 3, 3))
plot(m_tSNE[,c("data_dim_1","data_dim_2")], pch = 20, cex = 0.3, col = "#F2F2F2", frame.plot = F, xaxt = 'n', yaxt = 'n', ann = FALSE)
for ( k in 1:length(is_XEN_stage$Batch) ){
  points(is_XEN_stage[k,1], is_XEN_stage[k,2], col = "#38598CFF", pch = 20, cex = 0.3)
}
dev.off()

is_SII8D_stage <- m_tSNE[m_tSNE$Batch %in% "SII8D",]
tiff(file = paste0("is_SII8D_stage","",".tiff"), width = 7, height = 7, units = 'in', res = 300, compression = 'none')
par(mar = c(3, 3, 3, 3))
plot(m_tSNE[,c("data_dim_1","data_dim_2")], pch = 20, cex = 0.3, col = "#F2F2F2", frame.plot = F, xaxt = 'n', yaxt = 'n', ann = FALSE)
for ( k in 1:length(is_SII8D_stage$Batch) ){
  points(is_SII8D_stage[k,1], is_SII8D_stage[k,2], col = "#2D708EFF", pch = 20, cex = 0.3)
}
dev.off()

is_SII12D_stage <- m_tSNE[m_tSNE$Batch %in% "SII12D",]
tiff(file = paste0("is_SII12D_stage","",".tiff"), width = 7, height = 7, units = 'in', res = 300, compression = 'none')
par(mar = c(3, 3, 3, 3))
plot(m_tSNE[,c("data_dim_1","data_dim_2")], pch = 20, cex = 0.3, col = "#F2F2F2", frame.plot = F, xaxt = 'n', yaxt = 'n', ann = FALSE)
for ( k in 1:length(is_SII12D_stage$Batch) ){
  points(is_SII12D_stage[k,1], is_SII12D_stage[k,2], col = "#25858EFF", pch = 20, cex = 0.3)
}
dev.off()

is_SIII3D_stage <- m_tSNE[m_tSNE$Batch %in% "SIII3D",]
tiff(file = paste0("is_SIII3D_stage","",".tiff"), width = 7, height = 7, units = 'in', res = 300, compression = 'none')
par(mar = c(3, 3, 3, 3))
plot(m_tSNE[,c("data_dim_1","data_dim_2")], pch = 20, cex = 0.3, col = "#F2F2F2", frame.plot = F, xaxt = 'n', yaxt = 'n', ann = FALSE)
for ( k in 1:length(is_SIII3D_stage$Batch) ){
  points(is_SIII3D_stage[k,1], is_SIII3D_stage[k,2], col = "#1E9B8AFF", pch = 20, cex = 0.3)
}
dev.off()

is_SIII6D_stage <- m_tSNE[m_tSNE$Batch %in% "SIII6D",]
tiff(file = paste0("is_SIII6D_stage","",".tiff"), width = 7, height = 7, units = 'in', res = 300, compression = 'none')
par(mar = c(3, 3, 3, 3))
plot(m_tSNE[,c("data_dim_1","data_dim_2")], pch = 20, cex = 0.3, col = "#F2F2F2", frame.plot = F, xaxt = 'n', yaxt = 'n', ann = FALSE)
for ( k in 1:length(is_SIII6D_stage$Batch) ){
  points(is_SIII6D_stage[k,1], is_SIII6D_stage[k,2], col = "#2BB07FFF", pch = 20, cex = 0.3)
}
dev.off()

is_SIII8D_stage <- m_tSNE[m_tSNE$Batch %in% "SIII8D",]
tiff(file = paste0("is_SIII8D_stage","",".tiff"), width = 7, height = 7, units = 'in', res = 300, compression = 'none')
par(mar = c(3, 3, 3, 3))
plot(m_tSNE[,c("data_dim_1","data_dim_2")], pch = 20, cex = 0.3, col = "#F2F2F2", frame.plot = F, xaxt = 'n', yaxt = 'n', ann = FALSE)
for ( k in 1:length(is_SIII8D_stage$Batch) ){
  points(is_SIII8D_stage[k,1], is_SIII8D_stage[k,2], col = "#51C56AFF", pch = 20, cex = 0.3)
}
dev.off()

is_SIII10D_stage <- m_tSNE[m_tSNE$Batch %in% "SIII10D",]
tiff(file = paste0("is_SIII10D_stage","",".tiff"), width = 7, height = 7, units = 'in', res = 300, compression = 'none')
par(mar = c(3, 3, 3, 3))
plot(m_tSNE[,c("data_dim_1","data_dim_2")], pch = 20, cex = 0.3, col = "#F2F2F2", frame.plot = F, xaxt = 'n', yaxt = 'n', ann = FALSE)
for ( k in 1:length(is_SIII10D_stage$Batch) ){
  points(is_SIII10D_stage[k,1], is_SIII10D_stage[k,2], col = "#85D54AFF", pch = 20, cex = 0.3)
}
dev.off()

is_SIII15D_stage <- m_tSNE[m_tSNE$Batch %in% "SIII15D",]
tiff(file = paste0("is_SIII15D_stage","",".tiff"), width = 7, height = 7, units = 'in', res = 300, compression = 'none')
par(mar = c(3, 3, 3, 3))
plot(m_tSNE[,c("data_dim_1","data_dim_2")], pch = 20, cex = 0.3, col = "#F2F2F2", frame.plot = F, xaxt = 'n', yaxt = 'n', ann = FALSE)
for ( k in 1:length(is_SIII15D_stage$Batch) ){
  points(is_SIII15D_stage[k,1], is_SIII15D_stage[k,2], col = "#C2DF23FF", pch = 20, cex = 0.3)
}
dev.off()

is_SIII21D_stage <- m_tSNE[m_tSNE$Batch %in% "SIII21D",]
tiff(file = paste0("is_SIII21D_stage","",".tiff"), width = 7, height = 7, units = 'in', res = 300, compression = 'none')
par(mar = c(3, 3, 3, 3))
plot(m_tSNE[,c("data_dim_1","data_dim_2")], pch = 20, cex = 0.3, col = "#F2F2F2", frame.plot = F, xaxt = 'n', yaxt = 'n', ann = FALSE)
for ( k in 1:length(is_SIII21D_stage$Batch) ){
  points(is_SIII21D_stage[k,1], is_SIII21D_stage[k,2], col = "#FDE725FF", pch = 20, cex = 0.3)
}
dev.off()

G1_MEFs <- read.csv("G1_MEFs.csv")
row.names(G1_MEFs) <- G1_MEFs$Barcode
G1_MEFs <- G1_MEFs[,c("Data.Dim.1","Data.Dim.2")]
colnames(G1_MEFs) <- c("data_dim_1","data_dim_2")

feeders <- read.csv("feeders.csv")
row.names(feeders) <- feeders$Barcode
feeders <- feeders[,c("Data.Dim.1","Data.Dim.2")]
colnames(feeders) <- c("data_dim_1","data_dim_2")

feeders_II_III <- read.csv("rectangle.csv")
feeders_II_III <- subset(feeders_II_III, Batch %in% c("SIII8D","SIII3D","SII12D","SII8D","SIII10D","SIII6D","SIII15D","SIII21D"))
feeders_II_III <- feeders_II_III[,c("Data.Dim.1","Data.Dim.2")]
colnames(feeders_II_III) <- c("data_dim_1","data_dim_2")

G2_RPG_I <- read.csv("G2_RPG_I.csv")
row.names(G2_RPG_I) <- G2_RPG_I$Barcode
G2_RPG_I <- G2_RPG_I[,c("Data.Dim.1","Data.Dim.2")]
colnames(G2_RPG_I) <- c("data_dim_1","data_dim_2")

RPG_I_all <- read.csv("RPG_I_all.csv")
row.names(RPG_I_all) <- RPG_I_all$Barcode
RPG_I_all <- RPG_I_all[,c("Data.Dim.1","Data.Dim.2")]
colnames(RPG_I_all) <- c("data_dim_1","data_dim_2")

G3_XENs <- read.csv("G3_XENs.csv")
row.names(G3_XENs) <- G3_XENs$Barcode
G3_XENs <- G3_XENs[,c("Data.Dim.1","Data.Dim.2")]
colnames(G3_XENs) <- c("data_dim_1","data_dim_2")

G4_CiPSCs <- read.csv("G4_CiPSCs.csv")
row.names(G4_CiPSCs) <- G4_CiPSCs$Barcode
G4_CiPSCs <- G4_CiPSCs[,c("Data.Dim.1","Data.Dim.2")]
colnames(G4_CiPSCs) <- c("data_dim_1","data_dim_2")

G5_Ci2Cs <- read.csv("G5_Ci2Cs.csv")
row.names(G5_Ci2Cs) <- G5_Ci2Cs$Barcode
G5_Ci2Cs <- G5_Ci2Cs[,c("Data.Dim.1","Data.Dim.2")]
colnames(G5_Ci2Cs) <- c("data_dim_1","data_dim_2")

G6_RPG_II_III <- read.csv("G6_RPG_II_III.csv")
row.names(G6_RPG_II_III) <- G6_RPG_II_III$Barcode
G6_RPG_II_III <- G6_RPG_II_III[,c("Data.Dim.1","Data.Dim.2")]
colnames(G6_RPG_II_III) <- c("data_dim_1","data_dim_2")

dim(G1_MEFs)
dim(G2_RPG_I)
dim(G3_XENs)
dim(G4_CiPSCs)
dim(G5_Ci2Cs)
dim(G6_RPG_II_III)

tiff(file = paste0("groups","_v5",".tiff"), width = 7, height = 7, units = 'in', res = 300, compression = 'none')
par(mar = c(3, 3, 3, 3))
plot(m_tSNE[,c("data_dim_1","data_dim_2")], pch = 20, cex = 0.3, col = "#F2F2F2", frame.plot = F, xaxt = 'n', yaxt = 'n', ann = FALSE)
for ( k in 1:length(row.names(G1_MEFs)) ){
  points(G1_MEFs[k,"data_dim_1"],G1_MEFs[k,"data_dim_2"],col="#0365C0",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(feeders)) ){
  points(feeders[k,"data_dim_1"],feeders[k,"data_dim_2"],col="#8B572A",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(G2_RPG_I)) ){
  points(G2_RPG_I[k,"data_dim_1"],G2_RPG_I[k,"data_dim_2"],col="#01882B",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(RPG_I_all)) ){
  points(RPG_I_all[k,"data_dim_1"],RPG_I_all[k,"data_dim_2"],col="#01882B",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(feeders_II_III)) ){
  points(feeders_II_III[k,"data_dim_1"],feeders_II_III[k,"data_dim_2"],col="#8B572A",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(G3_XENs)) ){
  points(G3_XENs[k,"data_dim_1"],G3_XENs[k,"data_dim_2"],col="#DCBD23",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(G4_CiPSCs)) ){
  points(G4_CiPSCs[k,"data_dim_1"],G4_CiPSCs[k,"data_dim_2"],col="#C82506",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(G6_RPG_II_III)) ){
  points(G6_RPG_II_III[k,"data_dim_1"],G6_RPG_II_III[k,"data_dim_2"],col="#DF6A0F",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(G5_Ci2Cs)) ){
  points(G5_Ci2Cs[k,"data_dim_1"],G5_Ci2Cs[k,"data_dim_2"],col="#773F9B",pch=20,cex=0.3)
}
dev.off()

G1_MEFs_in_MEF_stage <- G1_MEFs[row.names(G1_MEFs) %in% row.names(is_MEF_stage),]
G2_RPG_I_in_MEF_stage <- G2_RPG_I[row.names(G2_RPG_I) %in% row.names(is_MEF_stage),]
G3_XENs_in_MEF_stage <- G3_XENs[row.names(G3_XENs) %in% row.names(is_MEF_stage),]
G4_CiPSCs_in_MEF_stage <- G4_CiPSCs[row.names(G4_CiPSCs) %in% row.names(is_MEF_stage),]
G5_Ci2Cs_in_MEF_stage <- G5_Ci2Cs[row.names(G5_Ci2Cs) %in% row.names(is_MEF_stage),]
G6_RPG_II_III_in_MEF_stage <- G6_RPG_II_III[row.names(G6_RPG_II_III) %in% row.names(is_MEF_stage),]

length(row.names(G1_MEFs_in_MEF_stage))
length(row.names(G2_RPG_I_in_MEF_stage))
length(row.names(G3_XENs_in_MEF_stage))
length(row.names(G4_CiPSCs_in_MEF_stage))
length(row.names(G5_Ci2Cs_in_MEF_stage))
length(row.names(G6_RPG_II_III_in_MEF_stage))

feeders_in_MEF_stage <- feeders[row.names(feeders) %in% row.names(is_MEF_stage),]
RPG_I_all_in_MEF_stage <- RPG_I_all[row.names(RPG_I_all) %in% row.names(is_MEF_stage),]
feeders_II_III_in_MEF_stage <- feeders_II_III[row.names(feeders_II_III) %in% row.names(is_MEF_stage),]

length(row.names(feeders_in_MEF_stage))
length(row.names(RPG_I_all_in_MEF_stage))

tiff(file = paste0("groups_in_","MEF","_stage.tiff"), width = 7, height = 7, units = 'in', res = 300, compression = 'none')
par(mar = c(3, 3, 3, 3))
plot(m_tSNE[,c("data_dim_1","data_dim_2")], pch = 20, cex = 0.3, col = "#F2F2F2", frame.plot = F, xaxt = 'n', yaxt = 'n', ann = FALSE)
for ( k in 1:length(row.names(is_MEF_stage)) ){
  points(is_MEF_stage[k,1], is_MEF_stage[k,2], col = "#C82406", pch = 20, cex = 0.3)
}
for ( k in 1:length(row.names(G1_MEFs_in_MEF_stage)) ){
  points(G1_MEFs_in_MEF_stage[k,"data_dim_1"],G1_MEFs_in_MEF_stage[k,"data_dim_2"],col="#0365C0",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(feeders_in_MEF_stage)) ){
  points(feeders_in_MEF_stage[k,"data_dim_1"],feeders_in_MEF_stage[k,"data_dim_2"],col="#8B572A",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(RPG_I_all_in_MEF_stage)) ){
  points(RPG_I_all_in_MEF_stage[k,"data_dim_1"],RPG_I_all_in_MEF_stage[k,"data_dim_2"],col="#01882B",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(feeders_II_III_in_MEF_stage)) ){
  points(feeders_II_III_in_MEF_stage[k,"data_dim_1"],feeders_II_III_in_MEF_stage[k,"data_dim_2"],col="#8B572A",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(G3_XENs_in_MEF_stage)) ){
  points(G3_XENs_in_MEF_stage[k,"data_dim_1"],G3_XENs_in_MEF_stage[k,"data_dim_2"],col="#DCBD23",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(G4_CiPSCs_in_MEF_stage)) ){
  points(G4_CiPSCs_in_MEF_stage[k,"data_dim_1"],G4_CiPSCs_in_MEF_stage[k,"data_dim_2"],col="#C82506",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(G6_RPG_II_III_in_MEF_stage)) ){
  points(G6_RPG_II_III_in_MEF_stage[k,"data_dim_1"],G6_RPG_II_III_in_MEF_stage[k,"data_dim_2"],col="#DF6A0F",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(G5_Ci2Cs_in_MEF_stage)) ){
  points(G5_Ci2Cs_in_MEF_stage[k,"data_dim_1"],G5_Ci2Cs_in_MEF_stage[k,"data_dim_2"],col="#773F9B",pch=20,cex=0.3)
}
dev.off()

# SI5D_stage
G1_MEFs_in_SI5D_stage <- G1_MEFs[row.names(G1_MEFs) %in% row.names(is_SI5D_stage),]
G2_RPG_I_in_SI5D_stage <- G2_RPG_I[row.names(G2_RPG_I) %in% row.names(is_SI5D_stage),]
G3_XENs_in_SI5D_stage <- G3_XENs[row.names(G3_XENs) %in% row.names(is_SI5D_stage),]
G4_CiPSCs_in_SI5D_stage <- G4_CiPSCs[row.names(G4_CiPSCs) %in% row.names(is_SI5D_stage),]
G5_Ci2Cs_in_SI5D_stage <- G5_Ci2Cs[row.names(G5_Ci2Cs) %in% row.names(is_SI5D_stage),]
G6_RPG_II_III_in_SI5D_stage <- G6_RPG_II_III[row.names(G6_RPG_II_III) %in% row.names(is_SI5D_stage),]

length(row.names(G1_MEFs_in_SI5D_stage))
length(row.names(G2_RPG_I_in_SI5D_stage))
length(row.names(G3_XENs_in_SI5D_stage))
length(row.names(G4_CiPSCs_in_SI5D_stage))
length(row.names(G5_Ci2Cs_in_SI5D_stage))
length(row.names(G6_RPG_II_III_in_SI5D_stage))

feeders_in_SI5D_stage <- feeders[row.names(feeders) %in% row.names(is_SI5D_stage),]
RPG_I_all_in_SI5D_stage <- RPG_I_all[row.names(RPG_I_all) %in% row.names(is_SI5D_stage),]
feeders_II_III_in_SI5D_stage <- feeders_II_III[row.names(feeders_II_III) %in% row.names(is_SI5D_stage),]

length(row.names(feeders_in_SI5D_stage))
length(row.names(RPG_I_all_in_SI5D_stage))

tiff(file = paste0("groups_in_","SI5D_stage.tiff"), width = 7, height = 7, units = 'in', res = 300, compression = 'none')
par(mar = c(3, 3, 3, 3))
plot(m_tSNE[,c("data_dim_1","data_dim_2")], pch = 20, cex = 0.3, col = "#F2F2F2", frame.plot = F, xaxt = 'n', yaxt = 'n', ann = FALSE)
for ( k in 1:length(row.names(is_SI5D_stage)) ){
  points(is_SI5D_stage[k,1], is_SI5D_stage[k,2], col = "#C82406", pch = 20, cex = 0.3)
}
for ( k in 1:length(row.names(G1_MEFs_in_SI5D_stage)) ){
  points(G1_MEFs_in_SI5D_stage[k,"data_dim_1"],G1_MEFs_in_SI5D_stage[k,"data_dim_2"],col="#0365C0",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(feeders_in_SI5D_stage)) ){
  points(feeders_in_SI5D_stage[k,"data_dim_1"],feeders_in_SI5D_stage[k,"data_dim_2"],col="#8B572A",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(RPG_I_all_in_SI5D_stage)) ){
  points(RPG_I_all_in_SI5D_stage[k,"data_dim_1"],RPG_I_all_in_SI5D_stage[k,"data_dim_2"],col="#01882B",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(feeders_II_III_in_SI5D_stage)) ){
  points(feeders_II_III_in_SI5D_stage[k,"data_dim_1"],feeders_II_III_in_SI5D_stage[k,"data_dim_2"],col="#8B572A",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(G3_XENs_in_SI5D_stage)) ){
  points(G3_XENs_in_SI5D_stage[k,"data_dim_1"],G3_XENs_in_SI5D_stage[k,"data_dim_2"],col="#DCBD23",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(G4_CiPSCs_in_SI5D_stage)) ){
  points(G4_CiPSCs_in_SI5D_stage[k,"data_dim_1"],G4_CiPSCs_in_SI5D_stage[k,"data_dim_2"],col="#C82506",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(G6_RPG_II_III_in_SI5D_stage)) ){
  points(G6_RPG_II_III_in_SI5D_stage[k,"data_dim_1"],G6_RPG_II_III_in_SI5D_stage[k,"data_dim_2"],col="#DF6A0F",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(G5_Ci2Cs_in_SI5D_stage)) ){
  points(G5_Ci2Cs_in_SI5D_stage[k,"data_dim_1"],G5_Ci2Cs_in_SI5D_stage[k,"data_dim_2"],col="#773F9B",pch=20,cex=0.3)
}
dev.off()

# SI12D_stage
G1_MEFs_in_SI12D_stage <- G1_MEFs[row.names(G1_MEFs) %in% row.names(is_SI12D_stage),]
G2_RPG_I_in_SI12D_stage <- G2_RPG_I[row.names(G2_RPG_I) %in% row.names(is_SI12D_stage),]
G3_XENs_in_SI12D_stage <- G3_XENs[row.names(G3_XENs) %in% row.names(is_SI12D_stage),]
G4_CiPSCs_in_SI12D_stage <- G4_CiPSCs[row.names(G4_CiPSCs) %in% row.names(is_SI12D_stage),]
G5_Ci2Cs_in_SI12D_stage <- G5_Ci2Cs[row.names(G5_Ci2Cs) %in% row.names(is_SI12D_stage),]
G6_RPG_II_III_in_SI12D_stage <- G6_RPG_II_III[row.names(G6_RPG_II_III) %in% row.names(is_SI12D_stage),]

length(row.names(G1_MEFs_in_SI12D_stage))
length(row.names(G2_RPG_I_in_SI12D_stage))
length(row.names(G3_XENs_in_SI12D_stage))
length(row.names(G4_CiPSCs_in_SI12D_stage))
length(row.names(G5_Ci2Cs_in_SI12D_stage))
length(row.names(G6_RPG_II_III_in_SI12D_stage))

feeders_in_SI12D_stage <- feeders[row.names(feeders) %in% row.names(is_SI12D_stage),]
RPG_I_all_in_SI12D_stage <- RPG_I_all[row.names(RPG_I_all) %in% row.names(is_SI12D_stage),]
feeders_II_III_in_SI12D_stage <- feeders_II_III[row.names(feeders_II_III) %in% row.names(is_SI12D_stage),]

length(row.names(feeders_in_SI12D_stage))
length(row.names(RPG_I_all_in_SI12D_stage))

tiff(file = paste0("groups_in_","SI12D_stage.tiff"), width = 7, height = 7, units = 'in', res = 300, compression = 'none')
par(mar = c(3, 3, 3, 3))
plot(m_tSNE[,c("data_dim_1","data_dim_2")], pch = 20, cex = 0.3, col = "#F2F2F2", frame.plot = F, xaxt = 'n', yaxt = 'n', ann = FALSE)
for ( k in 1:length(row.names(is_SI12D_stage)) ){
  points(is_SI12D_stage[k,1], is_SI12D_stage[k,2], col = "#C82406", pch = 20, cex = 0.3)
}
for ( k in 1:length(row.names(G1_MEFs_in_SI12D_stage)) ){
  points(G1_MEFs_in_SI12D_stage[k,"data_dim_1"],G1_MEFs_in_SI12D_stage[k,"data_dim_2"],col="#0365C0",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(feeders_in_SI12D_stage)) ){
  points(feeders_in_SI12D_stage[k,"data_dim_1"],feeders_in_SI12D_stage[k,"data_dim_2"],col="#8B572A",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(RPG_I_all_in_SI12D_stage)) ){
  points(RPG_I_all_in_SI12D_stage[k,"data_dim_1"],RPG_I_all_in_SI12D_stage[k,"data_dim_2"],col="#01882B",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(feeders_II_III_in_SI12D_stage)) ){
  points(feeders_II_III_in_SI12D_stage[k,"data_dim_1"],feeders_II_III_in_SI12D_stage[k,"data_dim_2"],col="#8B572A",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(G3_XENs_in_SI12D_stage)) ){
  points(G3_XENs_in_SI12D_stage[k,"data_dim_1"],G3_XENs_in_SI12D_stage[k,"data_dim_2"],col="#DCBD23",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(G4_CiPSCs_in_SI12D_stage)) ){
  points(G4_CiPSCs_in_SI12D_stage[k,"data_dim_1"],G4_CiPSCs_in_SI12D_stage[k,"data_dim_2"],col="#C82506",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(G6_RPG_II_III_in_SI12D_stage)) ){
  points(G6_RPG_II_III_in_SI12D_stage[k,"data_dim_1"],G6_RPG_II_III_in_SI12D_stage[k,"data_dim_2"],col="#DF6A0F",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(G5_Ci2Cs_in_SI12D_stage)) ){
  points(G5_Ci2Cs_in_SI12D_stage[k,"data_dim_1"],G5_Ci2Cs_in_SI12D_stage[k,"data_dim_2"],col="#773F9B",pch=20,cex=0.3)
}
dev.off()

# XEN_stage
G1_MEFs_in_XEN_stage <- G1_MEFs[row.names(G1_MEFs) %in% row.names(is_XEN_stage),]
G2_RPG_I_in_XEN_stage <- G2_RPG_I[row.names(G2_RPG_I) %in% row.names(is_XEN_stage),]
G3_XENs_in_XEN_stage <- G3_XENs[row.names(G3_XENs) %in% row.names(is_XEN_stage),]
G4_CiPSCs_in_XEN_stage <- G4_CiPSCs[row.names(G4_CiPSCs) %in% row.names(is_XEN_stage),]
G5_Ci2Cs_in_XEN_stage <- G5_Ci2Cs[row.names(G5_Ci2Cs) %in% row.names(is_XEN_stage),]
G6_RPG_II_III_in_XEN_stage <- G6_RPG_II_III[row.names(G6_RPG_II_III) %in% row.names(is_XEN_stage),]

length(row.names(G1_MEFs_in_XEN_stage))
length(row.names(G2_RPG_I_in_XEN_stage))
length(row.names(G3_XENs_in_XEN_stage))
length(row.names(G4_CiPSCs_in_XEN_stage))
length(row.names(G5_Ci2Cs_in_XEN_stage))
length(row.names(G6_RPG_II_III_in_XEN_stage))

feeders_in_XEN_stage <- feeders[row.names(feeders) %in% row.names(is_XEN_stage),]
RPG_I_all_in_XEN_stage <- RPG_I_all[row.names(RPG_I_all) %in% row.names(is_XEN_stage),]
feeders_II_III_in_XEN_stage <- feeders_II_III[row.names(feeders_II_III) %in% row.names(is_XEN_stage),]

length(row.names(feeders_in_XEN_stage))
length(row.names(RPG_I_all_in_XEN_stage))

tiff(file = paste0("groups_in_","XEN_stage.tiff"), width = 7, height = 7, units = 'in', res = 300, compression = 'none')
par(mar = c(3, 3, 3, 3))
plot(m_tSNE[,c("data_dim_1","data_dim_2")], pch = 20, cex = 0.3, col = "#F2F2F2", frame.plot = F, xaxt = 'n', yaxt = 'n', ann = FALSE)
for ( k in 1:length(row.names(is_XEN_stage)) ){
  points(is_XEN_stage[k,1], is_XEN_stage[k,2], col = "#C82406", pch = 20, cex = 0.3)
}
for ( k in 1:length(row.names(G1_MEFs_in_XEN_stage)) ){
  points(G1_MEFs_in_XEN_stage[k,"data_dim_1"],G1_MEFs_in_XEN_stage[k,"data_dim_2"],col="#0365C0",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(feeders_in_XEN_stage)) ){
  points(feeders_in_XEN_stage[k,"data_dim_1"],feeders_in_XEN_stage[k,"data_dim_2"],col="#8B572A",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(RPG_I_all_in_XEN_stage)) ){
  points(RPG_I_all_in_XEN_stage[k,"data_dim_1"],RPG_I_all_in_XEN_stage[k,"data_dim_2"],col="#01882B",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(feeders_II_III_in_XEN_stage)) ){
  points(feeders_II_III_in_XEN_stage[k,"data_dim_1"],feeders_II_III_in_XEN_stage[k,"data_dim_2"],col="#8B572A",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(G3_XENs_in_XEN_stage)) ){
  points(G3_XENs_in_XEN_stage[k,"data_dim_1"],G3_XENs_in_XEN_stage[k,"data_dim_2"],col="#DCBD23",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(G4_CiPSCs_in_XEN_stage)) ){
  points(G4_CiPSCs_in_XEN_stage[k,"data_dim_1"],G4_CiPSCs_in_XEN_stage[k,"data_dim_2"],col="#C82506",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(G6_RPG_II_III_in_XEN_stage)) ){
  points(G6_RPG_II_III_in_XEN_stage[k,"data_dim_1"],G6_RPG_II_III_in_XEN_stage[k,"data_dim_2"],col="#DF6A0F",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(G5_Ci2Cs_in_XEN_stage)) ){
  points(G5_Ci2Cs_in_XEN_stage[k,"data_dim_1"],G5_Ci2Cs_in_XEN_stage[k,"data_dim_2"],col="#773F9B",pch=20,cex=0.3)
}
dev.off()

# SII8D_stage
G1_MEFs_in_SII8D_stage <- G1_MEFs[row.names(G1_MEFs) %in% row.names(is_SII8D_stage),]
G2_RPG_I_in_SII8D_stage <- G2_RPG_I[row.names(G2_RPG_I) %in% row.names(is_SII8D_stage),]
G3_XENs_in_SII8D_stage <- G3_XENs[row.names(G3_XENs) %in% row.names(is_SII8D_stage),]
G4_CiPSCs_in_SII8D_stage <- G4_CiPSCs[row.names(G4_CiPSCs) %in% row.names(is_SII8D_stage),]
G5_Ci2Cs_in_SII8D_stage <- G5_Ci2Cs[row.names(G5_Ci2Cs) %in% row.names(is_SII8D_stage),]
G6_RPG_II_III_in_SII8D_stage <- G6_RPG_II_III[row.names(G6_RPG_II_III) %in% row.names(is_SII8D_stage),]

length(row.names(G1_MEFs_in_SII8D_stage))
length(row.names(G2_RPG_I_in_SII8D_stage))
length(row.names(G3_XENs_in_SII8D_stage))
length(row.names(G4_CiPSCs_in_SII8D_stage))
length(row.names(G5_Ci2Cs_in_SII8D_stage))
length(row.names(G6_RPG_II_III_in_SII8D_stage))

feeders_in_SII8D_stage <- feeders[row.names(feeders) %in% row.names(is_SII8D_stage),]
RPG_I_all_in_SII8D_stage <- RPG_I_all[row.names(RPG_I_all) %in% row.names(is_SII8D_stage),]
feeders_II_III_in_SII8D_stage <- feeders_II_III[row.names(feeders_II_III) %in% row.names(is_SII8D_stage),]

length(row.names(feeders_in_SII8D_stage))
length(row.names(RPG_I_all_in_SII8D_stage))

tiff(file = paste0("groups_in_","SII8D_stage.tiff"), width = 7, height = 7, units = 'in', res = 300, compression = 'none')
par(mar = c(3, 3, 3, 3))
plot(m_tSNE[,c("data_dim_1","data_dim_2")], pch = 20, cex = 0.3, col = "#F2F2F2", frame.plot = F, xaxt = 'n', yaxt = 'n', ann = FALSE)
for ( k in 1:length(row.names(is_SII8D_stage)) ){
  points(is_SII8D_stage[k,1], is_SII8D_stage[k,2], col = "#C82406", pch = 20, cex = 0.3)
}
for ( k in 1:length(row.names(G1_MEFs_in_SII8D_stage)) ){
  points(G1_MEFs_in_SII8D_stage[k,"data_dim_1"],G1_MEFs_in_SII8D_stage[k,"data_dim_2"],col="#0365C0",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(feeders_in_SII8D_stage)) ){
  points(feeders_in_SII8D_stage[k,"data_dim_1"],feeders_in_SII8D_stage[k,"data_dim_2"],col="#8B572A",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(RPG_I_all_in_SII8D_stage)) ){
  points(RPG_I_all_in_SII8D_stage[k,"data_dim_1"],RPG_I_all_in_SII8D_stage[k,"data_dim_2"],col="#01882B",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(feeders_II_III_in_SII8D_stage)) ){
  points(feeders_II_III_in_SII8D_stage[k,"data_dim_1"],feeders_II_III_in_SII8D_stage[k,"data_dim_2"],col="#8B572A",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(G3_XENs_in_SII8D_stage)) ){
  points(G3_XENs_in_SII8D_stage[k,"data_dim_1"],G3_XENs_in_SII8D_stage[k,"data_dim_2"],col="#DCBD23",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(G4_CiPSCs_in_SII8D_stage)) ){
  points(G4_CiPSCs_in_SII8D_stage[k,"data_dim_1"],G4_CiPSCs_in_SII8D_stage[k,"data_dim_2"],col="#C82506",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(G6_RPG_II_III_in_SII8D_stage)) ){
  points(G6_RPG_II_III_in_SII8D_stage[k,"data_dim_1"],G6_RPG_II_III_in_SII8D_stage[k,"data_dim_2"],col="#DF6A0F",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(G5_Ci2Cs_in_SII8D_stage)) ){
  points(G5_Ci2Cs_in_SII8D_stage[k,"data_dim_1"],G5_Ci2Cs_in_SII8D_stage[k,"data_dim_2"],col="#773F9B",pch=20,cex=0.3)
}
dev.off()

# SII12D_stage
G1_MEFs_in_SII12D_stage <- G1_MEFs[row.names(G1_MEFs) %in% row.names(is_SII12D_stage),]
G2_RPG_I_in_SII12D_stage <- G2_RPG_I[row.names(G2_RPG_I) %in% row.names(is_SII12D_stage),]
G3_XENs_in_SII12D_stage <- G3_XENs[row.names(G3_XENs) %in% row.names(is_SII12D_stage),]
G4_CiPSCs_in_SII12D_stage <- G4_CiPSCs[row.names(G4_CiPSCs) %in% row.names(is_SII12D_stage),]
G5_Ci2Cs_in_SII12D_stage <- G5_Ci2Cs[row.names(G5_Ci2Cs) %in% row.names(is_SII12D_stage),]
G6_RPG_II_III_in_SII12D_stage <- G6_RPG_II_III[row.names(G6_RPG_II_III) %in% row.names(is_SII12D_stage),]

length(row.names(G1_MEFs_in_SII12D_stage))
length(row.names(G2_RPG_I_in_SII12D_stage))
length(row.names(G3_XENs_in_SII12D_stage))
length(row.names(G4_CiPSCs_in_SII12D_stage))
length(row.names(G5_Ci2Cs_in_SII12D_stage))
length(row.names(G6_RPG_II_III_in_SII12D_stage))

feeders_in_SII12D_stage <- feeders[row.names(feeders) %in% row.names(is_SII12D_stage),]
RPG_I_all_in_SII12D_stage <- RPG_I_all[row.names(RPG_I_all) %in% row.names(is_SII12D_stage),]
feeders_II_III_in_SII12D_stage <- feeders_II_III[row.names(feeders_II_III) %in% row.names(is_SII12D_stage),]

length(row.names(feeders_in_SII12D_stage))
length(row.names(RPG_I_all_in_SII12D_stage))

tiff(file = paste0("groups_in_","SII12D_stage.tiff"), width = 7, height = 7, units = 'in', res = 300, compression = 'none')
par(mar = c(3, 3, 3, 3))
plot(m_tSNE[,c("data_dim_1","data_dim_2")], pch = 20, cex = 0.3, col = "#F2F2F2", frame.plot = F, xaxt = 'n', yaxt = 'n', ann = FALSE)
for ( k in 1:length(row.names(is_SII12D_stage)) ){
  points(is_SII12D_stage[k,1], is_SII12D_stage[k,2], col = "#C82406", pch = 20, cex = 0.3)
}
for ( k in 1:length(row.names(G1_MEFs_in_SII12D_stage)) ){
  points(G1_MEFs_in_SII12D_stage[k,"data_dim_1"],G1_MEFs_in_SII12D_stage[k,"data_dim_2"],col="#0365C0",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(feeders_in_SII12D_stage)) ){
  points(feeders_in_SII12D_stage[k,"data_dim_1"],feeders_in_SII12D_stage[k,"data_dim_2"],col="#8B572A",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(RPG_I_all_in_SII12D_stage)) ){
  points(RPG_I_all_in_SII12D_stage[k,"data_dim_1"],RPG_I_all_in_SII12D_stage[k,"data_dim_2"],col="#01882B",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(feeders_II_III_in_SII12D_stage)) ){
  points(feeders_II_III_in_SII12D_stage[k,"data_dim_1"],feeders_II_III_in_SII12D_stage[k,"data_dim_2"],col="#8B572A",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(G3_XENs_in_SII12D_stage)) ){
  points(G3_XENs_in_SII12D_stage[k,"data_dim_1"],G3_XENs_in_SII12D_stage[k,"data_dim_2"],col="#DCBD23",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(G4_CiPSCs_in_SII12D_stage)) ){
  points(G4_CiPSCs_in_SII12D_stage[k,"data_dim_1"],G4_CiPSCs_in_SII12D_stage[k,"data_dim_2"],col="#C82506",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(G6_RPG_II_III_in_SII12D_stage)) ){
  points(G6_RPG_II_III_in_SII12D_stage[k,"data_dim_1"],G6_RPG_II_III_in_SII12D_stage[k,"data_dim_2"],col="#DF6A0F",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(G5_Ci2Cs_in_SII12D_stage)) ){
  points(G5_Ci2Cs_in_SII12D_stage[k,"data_dim_1"],G5_Ci2Cs_in_SII12D_stage[k,"data_dim_2"],col="#773F9B",pch=20,cex=0.3)
}
dev.off()

# SIII3D_stage
G1_MEFs_in_SIII3D_stage <- G1_MEFs[row.names(G1_MEFs) %in% row.names(is_SIII3D_stage),]
G2_RPG_I_in_SIII3D_stage <- G2_RPG_I[row.names(G2_RPG_I) %in% row.names(is_SIII3D_stage),]
G3_XENs_in_SIII3D_stage <- G3_XENs[row.names(G3_XENs) %in% row.names(is_SIII3D_stage),]
G4_CiPSCs_in_SIII3D_stage <- G4_CiPSCs[row.names(G4_CiPSCs) %in% row.names(is_SIII3D_stage),]
G5_Ci2Cs_in_SIII3D_stage <- G5_Ci2Cs[row.names(G5_Ci2Cs) %in% row.names(is_SIII3D_stage),]
G6_RPG_II_III_in_SIII3D_stage <- G6_RPG_II_III[row.names(G6_RPG_II_III) %in% row.names(is_SIII3D_stage),]

length(row.names(G1_MEFs_in_SIII3D_stage))
length(row.names(G2_RPG_I_in_SIII3D_stage))
length(row.names(G3_XENs_in_SIII3D_stage))
length(row.names(G4_CiPSCs_in_SIII3D_stage))
length(row.names(G5_Ci2Cs_in_SIII3D_stage))
length(row.names(G6_RPG_II_III_in_SIII3D_stage))

feeders_in_SIII3D_stage <- feeders[row.names(feeders) %in% row.names(is_SIII3D_stage),]
RPG_I_all_in_SIII3D_stage <- RPG_I_all[row.names(RPG_I_all) %in% row.names(is_SIII3D_stage),]
feeders_II_III_in_SIII3D_stage <- feeders_II_III[row.names(feeders_II_III) %in% row.names(is_SIII3D_stage),]

length(row.names(feeders_in_SIII3D_stage))
length(row.names(RPG_I_all_in_SIII3D_stage))

tiff(file = paste0("groups_in_","SIII3D_stage.tiff"), width = 7, height = 7, units = 'in', res = 300, compression = 'none')
par(mar = c(3, 3, 3, 3))
plot(m_tSNE[,c("data_dim_1","data_dim_2")], pch = 20, cex = 0.3, col = "#F2F2F2", frame.plot = F, xaxt = 'n', yaxt = 'n', ann = FALSE)
for ( k in 1:length(row.names(is_SIII3D_stage)) ){
  points(is_SIII3D_stage[k,1], is_SIII3D_stage[k,2], col = "#C82406", pch = 20, cex = 0.3)
}
for ( k in 1:length(row.names(G1_MEFs_in_SIII3D_stage)) ){
  points(G1_MEFs_in_SIII3D_stage[k,"data_dim_1"],G1_MEFs_in_SIII3D_stage[k,"data_dim_2"],col="#0365C0",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(feeders_in_SIII3D_stage)) ){
  points(feeders_in_SIII3D_stage[k,"data_dim_1"],feeders_in_SIII3D_stage[k,"data_dim_2"],col="#8B572A",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(RPG_I_all_in_SIII3D_stage)) ){
  points(RPG_I_all_in_SIII3D_stage[k,"data_dim_1"],RPG_I_all_in_SIII3D_stage[k,"data_dim_2"],col="#01882B",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(feeders_II_III_in_SIII3D_stage)) ){
  points(feeders_II_III_in_SIII3D_stage[k,"data_dim_1"],feeders_II_III_in_SIII3D_stage[k,"data_dim_2"],col="#8B572A",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(G3_XENs_in_SIII3D_stage)) ){
  points(G3_XENs_in_SIII3D_stage[k,"data_dim_1"],G3_XENs_in_SIII3D_stage[k,"data_dim_2"],col="#DCBD23",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(G4_CiPSCs_in_SIII3D_stage)) ){
  points(G4_CiPSCs_in_SIII3D_stage[k,"data_dim_1"],G4_CiPSCs_in_SIII3D_stage[k,"data_dim_2"],col="#C82506",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(G6_RPG_II_III_in_SIII3D_stage)) ){
  points(G6_RPG_II_III_in_SIII3D_stage[k,"data_dim_1"],G6_RPG_II_III_in_SIII3D_stage[k,"data_dim_2"],col="#DF6A0F",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(G5_Ci2Cs_in_SIII3D_stage)) ){
  points(G5_Ci2Cs_in_SIII3D_stage[k,"data_dim_1"],G5_Ci2Cs_in_SIII3D_stage[k,"data_dim_2"],col="#773F9B",pch=20,cex=0.3)
}
dev.off()

# SIII6D_stage
G1_MEFs_in_SIII6D_stage <- G1_MEFs[row.names(G1_MEFs) %in% row.names(is_SIII6D_stage),]
G2_RPG_I_in_SIII6D_stage <- G2_RPG_I[row.names(G2_RPG_I) %in% row.names(is_SIII6D_stage),]
G3_XENs_in_SIII6D_stage <- G3_XENs[row.names(G3_XENs) %in% row.names(is_SIII6D_stage),]
G4_CiPSCs_in_SIII6D_stage <- G4_CiPSCs[row.names(G4_CiPSCs) %in% row.names(is_SIII6D_stage),]
G5_Ci2Cs_in_SIII6D_stage <- G5_Ci2Cs[row.names(G5_Ci2Cs) %in% row.names(is_SIII6D_stage),]
G6_RPG_II_III_in_SIII6D_stage <- G6_RPG_II_III[row.names(G6_RPG_II_III) %in% row.names(is_SIII6D_stage),]

length(row.names(G1_MEFs_in_SIII6D_stage))
length(row.names(G2_RPG_I_in_SIII6D_stage))
length(row.names(G3_XENs_in_SIII6D_stage))
length(row.names(G4_CiPSCs_in_SIII6D_stage))
length(row.names(G5_Ci2Cs_in_SIII6D_stage))
length(row.names(G6_RPG_II_III_in_SIII6D_stage))

feeders_in_SIII6D_stage <- feeders[row.names(feeders) %in% row.names(is_SIII6D_stage),]
RPG_I_all_in_SIII6D_stage <- RPG_I_all[row.names(RPG_I_all) %in% row.names(is_SIII6D_stage),]
feeders_II_III_in_SIII6D_stage <- feeders_II_III[row.names(feeders_II_III) %in% row.names(is_SIII6D_stage),]

length(row.names(feeders_in_SIII6D_stage))
length(row.names(RPG_I_all_in_SIII6D_stage))

tiff(file = paste0("groups_in_","SIII6D_stage.tiff"), width = 7, height = 7, units = 'in', res = 300, compression = 'none')
par(mar = c(3, 3, 3, 3))
plot(m_tSNE[,c("data_dim_1","data_dim_2")], pch = 20, cex = 0.3, col = "#F2F2F2", frame.plot = F, xaxt = 'n', yaxt = 'n', ann = FALSE)
for ( k in 1:length(row.names(is_SIII6D_stage)) ){
  points(is_SIII6D_stage[k,1], is_SIII6D_stage[k,2], col = "#C82406", pch = 20, cex = 0.3)
}
for ( k in 1:length(row.names(G1_MEFs_in_SIII6D_stage)) ){
  points(G1_MEFs_in_SIII6D_stage[k,"data_dim_1"],G1_MEFs_in_SIII6D_stage[k,"data_dim_2"],col="#0365C0",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(feeders_in_SIII6D_stage)) ){
  points(feeders_in_SIII6D_stage[k,"data_dim_1"],feeders_in_SIII6D_stage[k,"data_dim_2"],col="#8B572A",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(RPG_I_all_in_SIII6D_stage)) ){
  points(RPG_I_all_in_SIII6D_stage[k,"data_dim_1"],RPG_I_all_in_SIII6D_stage[k,"data_dim_2"],col="#01882B",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(feeders_II_III_in_SIII6D_stage)) ){
  points(feeders_II_III_in_SIII6D_stage[k,"data_dim_1"],feeders_II_III_in_SIII6D_stage[k,"data_dim_2"],col="#8B572A",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(G3_XENs_in_SIII6D_stage)) ){
  points(G3_XENs_in_SIII6D_stage[k,"data_dim_1"],G3_XENs_in_SIII6D_stage[k,"data_dim_2"],col="#DCBD23",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(G4_CiPSCs_in_SIII6D_stage)) ){
  points(G4_CiPSCs_in_SIII6D_stage[k,"data_dim_1"],G4_CiPSCs_in_SIII6D_stage[k,"data_dim_2"],col="#C82506",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(G6_RPG_II_III_in_SIII6D_stage)) ){
  points(G6_RPG_II_III_in_SIII6D_stage[k,"data_dim_1"],G6_RPG_II_III_in_SIII6D_stage[k,"data_dim_2"],col="#DF6A0F",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(G5_Ci2Cs_in_SIII6D_stage)) ){
  points(G5_Ci2Cs_in_SIII6D_stage[k,"data_dim_1"],G5_Ci2Cs_in_SIII6D_stage[k,"data_dim_2"],col="#773F9B",pch=20,cex=0.3)
}
dev.off()

# SIII8D_stage
G1_MEFs_in_SIII8D_stage <- G1_MEFs[row.names(G1_MEFs) %in% row.names(is_SIII8D_stage),]
G2_RPG_I_in_SIII8D_stage <- G2_RPG_I[row.names(G2_RPG_I) %in% row.names(is_SIII8D_stage),]
G3_XENs_in_SIII8D_stage <- G3_XENs[row.names(G3_XENs) %in% row.names(is_SIII8D_stage),]
G4_CiPSCs_in_SIII8D_stage <- G4_CiPSCs[row.names(G4_CiPSCs) %in% row.names(is_SIII8D_stage),]
G5_Ci2Cs_in_SIII8D_stage <- G5_Ci2Cs[row.names(G5_Ci2Cs) %in% row.names(is_SIII8D_stage),]
G6_RPG_II_III_in_SIII8D_stage <- G6_RPG_II_III[row.names(G6_RPG_II_III) %in% row.names(is_SIII8D_stage),]

length(row.names(G1_MEFs_in_SIII8D_stage))
length(row.names(G2_RPG_I_in_SIII8D_stage))
length(row.names(G3_XENs_in_SIII8D_stage))
length(row.names(G4_CiPSCs_in_SIII8D_stage))
length(row.names(G5_Ci2Cs_in_SIII8D_stage))
length(row.names(G6_RPG_II_III_in_SIII8D_stage))

feeders_in_SIII8D_stage <- feeders[row.names(feeders) %in% row.names(is_SIII8D_stage),]
RPG_I_all_in_SIII8D_stage <- RPG_I_all[row.names(RPG_I_all) %in% row.names(is_SIII8D_stage),]
feeders_II_III_in_SIII8D_stage <- feeders_II_III[row.names(feeders_II_III) %in% row.names(is_SIII8D_stage),]

length(row.names(feeders_in_SIII8D_stage))
length(row.names(RPG_I_all_in_SIII8D_stage))

tiff(file = paste0("groups_in_","SIII8D_stage.tiff"), width = 7, height = 7, units = 'in', res = 300, compression = 'none')
par(mar = c(3, 3, 3, 3))
plot(m_tSNE[,c("data_dim_1","data_dim_2")], pch = 20, cex = 0.3, col = "#F2F2F2", frame.plot = F, xaxt = 'n', yaxt = 'n', ann = FALSE)
for ( k in 1:length(row.names(is_SIII8D_stage)) ){
  points(is_SIII8D_stage[k,1], is_SIII8D_stage[k,2], col = "#C82406", pch = 20, cex = 0.3)
}
for ( k in 1:length(row.names(G1_MEFs_in_SIII8D_stage)) ){
  points(G1_MEFs_in_SIII8D_stage[k,"data_dim_1"],G1_MEFs_in_SIII8D_stage[k,"data_dim_2"],col="#0365C0",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(feeders_in_SIII8D_stage)) ){
  points(feeders_in_SIII8D_stage[k,"data_dim_1"],feeders_in_SIII8D_stage[k,"data_dim_2"],col="#8B572A",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(RPG_I_all_in_SIII8D_stage)) ){
  points(RPG_I_all_in_SIII8D_stage[k,"data_dim_1"],RPG_I_all_in_SIII8D_stage[k,"data_dim_2"],col="#01882B",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(feeders_II_III_in_SIII8D_stage)) ){
  points(feeders_II_III_in_SIII8D_stage[k,"data_dim_1"],feeders_II_III_in_SIII8D_stage[k,"data_dim_2"],col="#8B572A",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(G3_XENs_in_SIII8D_stage)) ){
  points(G3_XENs_in_SIII8D_stage[k,"data_dim_1"],G3_XENs_in_SIII8D_stage[k,"data_dim_2"],col="#DCBD23",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(G4_CiPSCs_in_SIII8D_stage)) ){
  points(G4_CiPSCs_in_SIII8D_stage[k,"data_dim_1"],G4_CiPSCs_in_SIII8D_stage[k,"data_dim_2"],col="#C82506",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(G6_RPG_II_III_in_SIII8D_stage)) ){
  points(G6_RPG_II_III_in_SIII8D_stage[k,"data_dim_1"],G6_RPG_II_III_in_SIII8D_stage[k,"data_dim_2"],col="#DF6A0F",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(G5_Ci2Cs_in_SIII8D_stage)) ){
  points(G5_Ci2Cs_in_SIII8D_stage[k,"data_dim_1"],G5_Ci2Cs_in_SIII8D_stage[k,"data_dim_2"],col="#773F9B",pch=20,cex=0.3)
}
dev.off()

# SIII10D_stage
G1_MEFs_in_SIII10D_stage <- G1_MEFs[row.names(G1_MEFs) %in% row.names(is_SIII10D_stage),]
G2_RPG_I_in_SIII10D_stage <- G2_RPG_I[row.names(G2_RPG_I) %in% row.names(is_SIII10D_stage),]
G3_XENs_in_SIII10D_stage <- G3_XENs[row.names(G3_XENs) %in% row.names(is_SIII10D_stage),]
G4_CiPSCs_in_SIII10D_stage <- G4_CiPSCs[row.names(G4_CiPSCs) %in% row.names(is_SIII10D_stage),]
G5_Ci2Cs_in_SIII10D_stage <- G5_Ci2Cs[row.names(G5_Ci2Cs) %in% row.names(is_SIII10D_stage),]
G6_RPG_II_III_in_SIII10D_stage <- G6_RPG_II_III[row.names(G6_RPG_II_III) %in% row.names(is_SIII10D_stage),]

length(row.names(G1_MEFs_in_SIII10D_stage))
length(row.names(G2_RPG_I_in_SIII10D_stage))
length(row.names(G3_XENs_in_SIII10D_stage))
length(row.names(G4_CiPSCs_in_SIII10D_stage))
length(row.names(G5_Ci2Cs_in_SIII10D_stage))
length(row.names(G6_RPG_II_III_in_SIII10D_stage))

feeders_in_SIII10D_stage <- feeders[row.names(feeders) %in% row.names(is_SIII10D_stage),]
RPG_I_all_in_SIII10D_stage <- RPG_I_all[row.names(RPG_I_all) %in% row.names(is_SIII10D_stage),]
feeders_II_III_in_SIII10D_stage <- feeders_II_III[row.names(feeders_II_III) %in% row.names(is_SIII10D_stage),]

length(row.names(feeders_in_SIII10D_stage))
length(row.names(RPG_I_all_in_SIII10D_stage))

tiff(file = paste0("groups_in_","SIII10D_stage.tiff"), width = 7, height = 7, units = 'in', res = 300, compression = 'none')
par(mar = c(3, 3, 3, 3))
plot(m_tSNE[,c("data_dim_1","data_dim_2")], pch = 20, cex = 0.3, col = "#F2F2F2", frame.plot = F, xaxt = 'n', yaxt = 'n', ann = FALSE)
for ( k in 1:length(row.names(is_SIII10D_stage)) ){
  points(is_SIII10D_stage[k,1], is_SIII10D_stage[k,2], col = "#C82406", pch = 20, cex = 0.3)
}
for ( k in 1:length(row.names(G1_MEFs_in_SIII10D_stage)) ){
  points(G1_MEFs_in_SIII10D_stage[k,"data_dim_1"],G1_MEFs_in_SIII10D_stage[k,"data_dim_2"],col="#0365C0",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(feeders_in_SIII10D_stage)) ){
  points(feeders_in_SIII10D_stage[k,"data_dim_1"],feeders_in_SIII10D_stage[k,"data_dim_2"],col="#8B572A",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(RPG_I_all_in_SIII10D_stage)) ){
  points(RPG_I_all_in_SIII10D_stage[k,"data_dim_1"],RPG_I_all_in_SIII10D_stage[k,"data_dim_2"],col="#01882B",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(feeders_II_III_in_SIII10D_stage)) ){
  points(feeders_II_III_in_SIII10D_stage[k,"data_dim_1"],feeders_II_III_in_SIII10D_stage[k,"data_dim_2"],col="#8B572A",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(G3_XENs_in_SIII10D_stage)) ){
  points(G3_XENs_in_SIII10D_stage[k,"data_dim_1"],G3_XENs_in_SIII10D_stage[k,"data_dim_2"],col="#DCBD23",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(G4_CiPSCs_in_SIII10D_stage)) ){
  points(G4_CiPSCs_in_SIII10D_stage[k,"data_dim_1"],G4_CiPSCs_in_SIII10D_stage[k,"data_dim_2"],col="#C82506",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(G6_RPG_II_III_in_SIII10D_stage)) ){
  points(G6_RPG_II_III_in_SIII10D_stage[k,"data_dim_1"],G6_RPG_II_III_in_SIII10D_stage[k,"data_dim_2"],col="#DF6A0F",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(G5_Ci2Cs_in_SIII10D_stage)) ){
  points(G5_Ci2Cs_in_SIII10D_stage[k,"data_dim_1"],G5_Ci2Cs_in_SIII10D_stage[k,"data_dim_2"],col="#773F9B",pch=20,cex=0.3)
}
dev.off()

# SIII15D_stage
G1_MEFs_in_SIII15D_stage <- G1_MEFs[row.names(G1_MEFs) %in% row.names(is_SIII15D_stage),]
G2_RPG_I_in_SIII15D_stage <- G2_RPG_I[row.names(G2_RPG_I) %in% row.names(is_SIII15D_stage),]
G3_XENs_in_SIII15D_stage <- G3_XENs[row.names(G3_XENs) %in% row.names(is_SIII15D_stage),]
G4_CiPSCs_in_SIII15D_stage <- G4_CiPSCs[row.names(G4_CiPSCs) %in% row.names(is_SIII15D_stage),]
G5_Ci2Cs_in_SIII15D_stage <- G5_Ci2Cs[row.names(G5_Ci2Cs) %in% row.names(is_SIII15D_stage),]
G6_RPG_II_III_in_SIII15D_stage <- G6_RPG_II_III[row.names(G6_RPG_II_III) %in% row.names(is_SIII15D_stage),]

length(row.names(G1_MEFs_in_SIII15D_stage))
length(row.names(G2_RPG_I_in_SIII15D_stage))
length(row.names(G3_XENs_in_SIII15D_stage))
length(row.names(G4_CiPSCs_in_SIII15D_stage))
length(row.names(G5_Ci2Cs_in_SIII15D_stage))
length(row.names(G6_RPG_II_III_in_SIII15D_stage))

feeders_in_SIII15D_stage <- feeders[row.names(feeders) %in% row.names(is_SIII15D_stage),]
RPG_I_all_in_SIII15D_stage <- RPG_I_all[row.names(RPG_I_all) %in% row.names(is_SIII15D_stage),]
feeders_II_III_in_SIII15D_stage <- feeders_II_III[row.names(feeders_II_III) %in% row.names(is_SIII15D_stage),]

length(row.names(feeders_in_SIII15D_stage))
length(row.names(RPG_I_all_in_SIII15D_stage))

tiff(file = paste0("groups_in_","SIII15D_stage.tiff"), width = 7, height = 7, units = 'in', res = 300, compression = 'none')
par(mar = c(3, 3, 3, 3))
plot(m_tSNE[,c("data_dim_1","data_dim_2")], pch = 20, cex = 0.3, col = "#F2F2F2", frame.plot = F, xaxt = 'n', yaxt = 'n', ann = FALSE)
for ( k in 1:length(row.names(is_SIII15D_stage)) ){
  points(is_SIII15D_stage[k,1], is_SIII15D_stage[k,2], col = "#C82406", pch = 20, cex = 0.3)
}
for ( k in 1:length(row.names(G1_MEFs_in_SIII15D_stage)) ){
  points(G1_MEFs_in_SIII15D_stage[k,"data_dim_1"],G1_MEFs_in_SIII15D_stage[k,"data_dim_2"],col="#0365C0",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(feeders_in_SIII15D_stage)) ){
  points(feeders_in_SIII15D_stage[k,"data_dim_1"],feeders_in_SIII15D_stage[k,"data_dim_2"],col="#8B572A",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(RPG_I_all_in_SIII15D_stage)) ){
  points(RPG_I_all_in_SIII15D_stage[k,"data_dim_1"],RPG_I_all_in_SIII15D_stage[k,"data_dim_2"],col="#01882B",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(feeders_II_III_in_SIII15D_stage)) ){
  points(feeders_II_III_in_SIII15D_stage[k,"data_dim_1"],feeders_II_III_in_SIII15D_stage[k,"data_dim_2"],col="#8B572A",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(G3_XENs_in_SIII15D_stage)) ){
  points(G3_XENs_in_SIII15D_stage[k,"data_dim_1"],G3_XENs_in_SIII15D_stage[k,"data_dim_2"],col="#DCBD23",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(G4_CiPSCs_in_SIII15D_stage)) ){
  points(G4_CiPSCs_in_SIII15D_stage[k,"data_dim_1"],G4_CiPSCs_in_SIII15D_stage[k,"data_dim_2"],col="#C82506",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(G6_RPG_II_III_in_SIII15D_stage)) ){
  points(G6_RPG_II_III_in_SIII15D_stage[k,"data_dim_1"],G6_RPG_II_III_in_SIII15D_stage[k,"data_dim_2"],col="#DF6A0F",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(G5_Ci2Cs_in_SIII15D_stage)) ){
  points(G5_Ci2Cs_in_SIII15D_stage[k,"data_dim_1"],G5_Ci2Cs_in_SIII15D_stage[k,"data_dim_2"],col="#773F9B",pch=20,cex=0.3)
}
dev.off()

# SIII21D_stage
G1_MEFs_in_SIII21D_stage <- G1_MEFs[row.names(G1_MEFs) %in% row.names(is_SIII21D_stage),]
G2_RPG_I_in_SIII21D_stage <- G2_RPG_I[row.names(G2_RPG_I) %in% row.names(is_SIII21D_stage),]
G3_XENs_in_SIII21D_stage <- G3_XENs[row.names(G3_XENs) %in% row.names(is_SIII21D_stage),]
G4_CiPSCs_in_SIII21D_stage <- G4_CiPSCs[row.names(G4_CiPSCs) %in% row.names(is_SIII21D_stage),]
G5_Ci2Cs_in_SIII21D_stage <- G5_Ci2Cs[row.names(G5_Ci2Cs) %in% row.names(is_SIII21D_stage),]
G6_RPG_II_III_in_SIII21D_stage <- G6_RPG_II_III[row.names(G6_RPG_II_III) %in% row.names(is_SIII21D_stage),]

length(row.names(G1_MEFs_in_SIII21D_stage))
length(row.names(G2_RPG_I_in_SIII21D_stage))
length(row.names(G3_XENs_in_SIII21D_stage))
length(row.names(G4_CiPSCs_in_SIII21D_stage))
length(row.names(G5_Ci2Cs_in_SIII21D_stage))
length(row.names(G6_RPG_II_III_in_SIII21D_stage))

feeders_in_SIII21D_stage <- feeders[row.names(feeders) %in% row.names(is_SIII21D_stage),]
RPG_I_all_in_SIII21D_stage <- RPG_I_all[row.names(RPG_I_all) %in% row.names(is_SIII21D_stage),]
feeders_II_III_in_SIII21D_stage <- feeders_II_III[row.names(feeders_II_III) %in% row.names(is_SIII21D_stage),]

length(row.names(feeders_in_SIII21D_stage))
length(row.names(RPG_I_all_in_SIII21D_stage))

tiff(file = paste0("groups_in_","SIII21D_stage.tiff"), width = 7, height = 7, units = 'in', res = 300, compression = 'none')
par(mar = c(3, 3, 3, 3))
plot(m_tSNE[,c("data_dim_1","data_dim_2")], pch = 20, cex = 0.3, col = "#F2F2F2", frame.plot = F, xaxt = 'n', yaxt = 'n', ann = FALSE)
for ( k in 1:length(row.names(is_SIII21D_stage)) ){
  points(is_SIII21D_stage[k,1], is_SIII21D_stage[k,2], col = "#C82406", pch = 20, cex = 0.3)
}
for ( k in 1:length(row.names(G1_MEFs_in_SIII21D_stage)) ){
  points(G1_MEFs_in_SIII21D_stage[k,"data_dim_1"],G1_MEFs_in_SIII21D_stage[k,"data_dim_2"],col="#0365C0",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(feeders_in_SIII21D_stage)) ){
  points(feeders_in_SIII21D_stage[k,"data_dim_1"],feeders_in_SIII21D_stage[k,"data_dim_2"],col="#8B572A",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(RPG_I_all_in_SIII21D_stage)) ){
  points(RPG_I_all_in_SIII21D_stage[k,"data_dim_1"],RPG_I_all_in_SIII21D_stage[k,"data_dim_2"],col="#01882B",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(feeders_II_III_in_SIII21D_stage)) ){
  points(feeders_II_III_in_SIII21D_stage[k,"data_dim_1"],feeders_II_III_in_SIII21D_stage[k,"data_dim_2"],col="#8B572A",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(G3_XENs_in_SIII21D_stage)) ){
  points(G3_XENs_in_SIII21D_stage[k,"data_dim_1"],G3_XENs_in_SIII21D_stage[k,"data_dim_2"],col="#DCBD23",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(G4_CiPSCs_in_SIII21D_stage)) ){
  points(G4_CiPSCs_in_SIII21D_stage[k,"data_dim_1"],G4_CiPSCs_in_SIII21D_stage[k,"data_dim_2"],col="#C82506",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(G6_RPG_II_III_in_SIII21D_stage)) ){
  points(G6_RPG_II_III_in_SIII21D_stage[k,"data_dim_1"],G6_RPG_II_III_in_SIII21D_stage[k,"data_dim_2"],col="#DF6A0F",pch=20,cex=0.3)
}
for ( k in 1:length(row.names(G5_Ci2Cs_in_SIII21D_stage)) ){
  points(G5_Ci2Cs_in_SIII21D_stage[k,"data_dim_1"],G5_Ci2Cs_in_SIII21D_stage[k,"data_dim_2"],col="#773F9B",pch=20,cex=0.3)
}
dev.off()
