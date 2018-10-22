library(RColorBrewer)

load("/Users/m/tmp/data/_AGG11_170901/AGG11_mapped/monocle_filtered_data_hpc/monocle_raw_norm.Robj")

monocle_trajectory <- read.csv("t-SNE.csv", row.names = 1)
rownames(monocle_trajectory) <- monocle_trajectory$Barcode
# re-order
Barcode <- colnames(monocle_norm)
monocle_trajectory <- monocle_trajectory[Barcode,]
m_trajectory <- monocle_trajectory[,c("data_dim_1","data_dim_2")]

plot(m_trajectory)

monocle_genes <- read.csv("/Users/m/tmp/data/_AGG11_170901/AGG11_mapped/monocle_filtered_data_hpc/monocle_fData.csv", row.names = 1)

gene <- "Zscan4b"
genes <- c("Prrx1","Sox17","Sall4","Pou5f1","Nanog","Zscan4d","Thy1","Gata6","Dppa3","Pdgfra","Pecam1","Sox2","Tcstv1","Gata3")
# PPT Fig2e tSNE List
prefix <- "Fig2E"
genes <- c("Prrx1","Sox17","Sall4","Pou5f1","Nanog","Zscan4d")
# FigS2E
prefix <- "FigS2E"
genes <- c("Col6a1","Thy1","Fbn1","Twist2","Gata4","Gata6","Foxa2","Epcam","Tet1","Dppa4","Sox2","Zfp42")
# FigS2H
prefix <- "FigS2H"
genes <- c("Zscan4b","Zscan4c","Zscan4f")
# MEF
genes <- c("Thy1","Fbn1","Twist2","Col6a1")
# XEN
prefix <- "XEN"
genes <- c("Gata4", "Gata6", "Foxa2", "Epcam")
# Stemness
prefix <- "Stemness"
genes <- c("Sox2", "Zfp42", "Tet1", "Dppa4")
# 2C
prefix <- "2C"
genes <- c("Zscan4b", "Zscan4c", "Tcstv1", "Tcstv3", "Gm4340", "Zscan4f", "Gm5039")
# Add
prefix <- "Add"
genes <- c("Gm21761", "Tcstv1", "Gm4340", "Tmem92")
# Add2
prefix <- "Add2"
genes <- c("Tcstv1", "Tcstv3")

for (gene in genes) {
  id <- monocle_genes[monocle_genes$gene_short_name %in% gene,]
  id <- as.vector(as.matrix(id))
  id <- id[1]
  l <- apply(monocle_norm[id,] - .1,2,sum) + .1
  f <- l == 0
  l <- log2(l)
  l[f] <- NA
  mi <- min(l,na.rm=TRUE)
  ma <- max(l,na.rm=TRUE)
  ColorRamp <- colorRampPalette(c("#F2F2F2", "#F1E9E9", "#EEDBDB", "#F96060", "#FF0000"))(256)
  # plot(rep(1,256),col=ColorRamp, pch=19,cex=2)
  ColorLevels <- seq(mi, ma, length=length(ColorRamp))
  v <- round((l - mi)/(ma - mi)*255 + 1,0)
  tiff(file=paste0(prefix,"_tSNE_",gene,"_ColorRamp",".tiff"), width = 7, height = 7, units = 'in', res = 300, compression = 'none')
  par(mfrow=c(1,1))
  layout(matrix(data=c(1,3,2,4), nrow=2, ncol=2), widths=c(5,1,5,1), heights=c(5,1,1,1))
  plot(m_trajectory,main=gene,pch=20,cex=0.3,col="#F2F2F2",frame.plot=F, axes=FALSE, xlab="", ylab="")
  set.seed(42)
  rand_v <- 1:length(v)
  rand_v <- sample(rand_v)
  # for ( k in 1:length(v) ){
  for ( k in rand_v ){
    points(m_trajectory[k,1],m_trajectory[k,2],col=ColorRamp[v[k]],pch=20,cex=0.3)
  }
  # image(1, ColorLevels,
  #       matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
  #       col=ColorRamp,
  #       xlab="",ylab="",
  #       xaxt="n")
  dev.off()
}

# pdf
id <- monocle_genes[monocle_genes$gene_short_name %in% gene,]
id <- as.vector(as.matrix(id))
id <- id[1]
l <- apply(monocle_norm[id,] - .1,2,sum) + .1
f <- l == 0
l <- log2(l)
l[f] <- NA
mi <- min(l,na.rm=TRUE)
ma <- max(l,na.rm=TRUE)
ColorLevels <- seq(mi, ma, length=length(ColorRamp))
v <- round((l - mi)/(ma - mi)*255 + 1,0)
pdf(file=paste0("tSNE_",gene,".pdf"), width = 7, height = 7, onefile = F)
par(mfrow=c(1,1))
layout(matrix(data=c(1,3,2,4), nrow=2, ncol=2), widths=c(5,1,5,1), heights=c(5,1,1,1))
plot(m_trajectory[1,1],main=gene,pch=20,cex=0.3,col="#F2F2F2",frame.plot=F, axes=FALSE, xlab="", ylab="")
set.seed(42)
rand_v <- 1:length(v)
rand_v <- sample(rand_v)
# for ( k in 1:length(v) ){
# for ( k in rand_v ){
#   points(m_trajectory[k,1],m_trajectory[k,2],col=ColorRamp[v[k]],pch=20,cex=0.3)
# }
image(1, ColorLevels,
      matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
      col=ColorRamp,
      xlab="",ylab="",
      xaxt="n")
dev.off()
