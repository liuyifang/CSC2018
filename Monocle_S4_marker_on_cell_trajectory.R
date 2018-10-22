library(RColorBrewer)

load("/Users/m/tmp/data/_AGG11_170901/AGG11_mapped/monocle_filtered_data_hpc/monocle_raw_norm.Robj")

monocle_genes <- read.csv("/Users/m/tmp/data/_AGG11_170901/AGG11_mapped/monocle_filtered_data_hpc/monocle_fData.csv", row.names = 1)
genes <- c("Dppa3","Nanog","Sox17","Thy1","Sall4")
genes <- c("Gata6", "Pdgfra", "Sox2", "Pecam1")
genes <- c("Sox2", "Pou5f1")

genes <- c("Tcstv1","Tcstv3","Gm5039","Zscan4c","Zscan4d","Gm4340","Zscan4f","Tmem92","Gm21761")

# max_5
monocle_trajectory <- read.csv("monocle.csv", row.names = 1)
rownames(monocle_trajectory) <- monocle_trajectory$Barcode
Barcode <- colnames(monocle_norm)
monocle_trajectory <- monocle_trajectory[Barcode,]
m_trajectory <- monocle_trajectory[,c("data_dim_1","data_dim_2")]
for (gene in genes) {
  id <- monocle_genes[monocle_genes$gene_short_name %in% gene,]
  id <- as.vector(as.matrix(id))
  id <- id[1]
  l <- apply(monocle_norm[id,] - .1,2,sum) + .1
  f <- l == 0
  l <- log2(l)
  l[f] <- NA
  mi <- min(l, na.rm = TRUE)
  ma <- max(l, na.rm = TRUE)
  # ColorRamp <- colorRampPalette(rev(brewer.pal(n = 11,name = "RdYlBu")))(256)
  ColorRamp <- colorRampPalette(brewer.pal(n = 9, name = "YlOrRd"))(256)
  ColorLevels <- seq(mi, ma, length = length(ColorRamp))
  v <- round((l - mi)/(ma - mi)*255 + 1, 0)
  tiff(file = paste0("",gene,"_monocle_",".tiff"), width = 7, height = 7, units = 'in', res = 300, compression = 'none')
  par(mfrow = c(1,1))
  layout(matrix(data=c(1,3,2,4), nrow=2, ncol=2), widths = c(5,1,5,1), heights = c(5,1,1,1))
  plot(m_trajectory, main = gene, pch = 20, cex = 1.5, col = "#F2F2F2", frame.plot = F, axes = FALSE, xlab = "", ylab = "")
  set.seed(42)
  rand_v <- 1:length(v)
  rand_v <- sample(rand_v)
  # for ( k in 1:length(v) ){
  for ( k in rand_v ){
    points(m_trajectory[k,1],m_trajectory[k,2],col=ColorRamp[v[k]],pch=20,cex=1.5)
  }
  # image(1, ColorLevels,
  #       matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
  #       col=ColorRamp,
  #       xlab="",ylab="",
  #       xaxt="n")
  dev.off()
}
