library(monocle)

load(file = paste0("t-SNE.Robj"))

cds <- AGG
x=1
y=2
color_by="Cluster"
markers=NULL
show_cell_names=FALSE
cell_size=1.5
cell_name_size=2
if (is.null(cds@reducedDimA) | length(pData(cds)$Cluster) == 0){
  stop("Error: Clustering is not performed yet. Please call clusterCells() before calling this function.")
}

gene_short_name <- NULL
sample_name <- NULL
data_dim_1 <- NULL
data_dim_2 <- NULL

lib_info <- pData(cds)
tSNE_dim_coords <- reducedDimA(cds)
data_df <- data.frame(t(tSNE_dim_coords[c(x,y),]))
colnames(data_df) <- c("data_dim_1", "data_dim_2")
data_df$sample_name <- colnames(cds)
data_df <- merge(data_df, lib_info, by.x="sample_name", by.y="row.names")

write.csv(data_df, file = paste0("t-SNE.csv"), quote = F)
