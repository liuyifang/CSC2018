.libPaths("/lustre1/lch3000_pkuhpc/liuyf/R/R-3.4.1/library")

library(viridis)
library(igraph)
library(monocle)

load("/lustre1/lch3000_pkuhpc/liuyf/denglab/zhaoting/monocle/AGG_estimate_size_and_dispersion.Robj")

m <- 5
DEG_genes <- read.csv("AGG4_DEG_v2_top1518.csv", row.names = 1)
DEG_genes <- DEG_genes$id
DEG_genes <- as.character(DEG_genes)
AGG <- setOrderingFilter(AGG, DEG_genes)
plot_pc_variance_explained(AGG, return_all = T)
AGG <- reduceDimension(AGG,
                       max_components = m,
                       norm_method = 'log',
                       verbose = T,
                       scaling = T)
AGG <- orderCells(AGG)

levels(pData(AGG)$Batch)
Batch_vec <- c("MEF", "SI5D", "SI12D", "XEN",
               "SII8D", "SII12D",
               "SIII3D", "SIII6D", "SIII8D", "SIII10D", "SIII15D", "SIII21D")
Batch_cols <- viridis(12, option = "viridis")
names(Batch_cols) <- Batch_vec
Batch_labels <- c("MEF", "SI D5", "SI D12", "XEN",
                  "SII D8", "SII D12",
                  "SIII D3", "SIII D6", "SIII D8", "SIII D10", "SIII D15", "SIII D21")

png(file = paste0("cell_trajectory_by_Batch_max_",m,".png"), width = 7, height = 7, units = "in", res = 300)
plot_cell_trajectory(AGG, color_by = "as.factor(Batch)", cell_size = 0.3, show_branch_points=TRUE) + scale_color_manual(labels = Batch_labels, values = Batch_cols, name = "Batch") + guides(col = guide_legend(nrow = 3, byrow = TRUE))
dev.off()

png(file = paste0("cell_trajectory_by_Batch_facet_max_",m,".png"), width = 7, height = 7, units = "in", res = 300)
plot_cell_trajectory(AGG, color_by = "as.factor(Batch)", cell_size = 0.3, show_branch_points=FALSE) + facet_wrap(~Batch, nrow = 3) + scale_color_manual(labels = Batch_labels, values = Batch_cols, name = "Batch") + guides(col = guide_legend(nrow = 3, byrow = TRUE))
dev.off()

cds <- AGG
x=1
y=2
color_by="State"
show_tree=TRUE
show_backbone=TRUE
backbone_color="black"
markers=NULL
show_cell_names=FALSE
cell_size=1.5
cell_link_size=0.75
cell_name_size=2
show_branch_points=TRUE
theta = 0
gene_short_name <- NA
sample_name <- NA
data_dim_1 <- NA
data_dim_2 <- NA

lib_info_with_pseudo <- pData(cds)

if (is.null(cds@dim_reduce_type)){
  stop("Error: dimensionality not yet reduced. Please call reduceDimension() before calling this function.")
}

if (cds@dim_reduce_type == "ICA"){
  reduced_dim_coords <- reducedDimS(cds)
}else if (cds@dim_reduce_type %in% c("simplePPT", "DDRTree") ){
  reduced_dim_coords <- reducedDimK(cds)
}else {
  stop("Error: unrecognized dimensionality reduction method.")
}

if (is.null(reduced_dim_coords)){
  stop("You must first call reduceDimension() before using this function")
}

ica_space_df <- data.frame(Matrix::t(reduced_dim_coords[c(x,y),]))
colnames(ica_space_df) <- c("prin_graph_dim_1", "prin_graph_dim_2")

ica_space_df$sample_name <- row.names(ica_space_df)
#ica_space_with_state_df <- merge(ica_space_df, lib_info_with_pseudo, by.x="sample_name", by.y="row.names")
#print(ica_space_with_state_df)
dp_mst <- minSpanningTree(cds)

if (is.null(dp_mst)){
  stop("You must first call orderCells() before using this function")
}

edge_list <- as.data.frame(get.edgelist(dp_mst))
colnames(edge_list) <- c("source", "target")

edge_df <- merge(ica_space_df, edge_list, by.x="sample_name", by.y="source", all=TRUE)
#edge_df <- ica_space_df
edge_df <- plyr::rename(edge_df, c("prin_graph_dim_1"="source_prin_graph_dim_1", "prin_graph_dim_2"="source_prin_graph_dim_2"))
edge_df <- merge(edge_df, ica_space_df[,c("sample_name", "prin_graph_dim_1", "prin_graph_dim_2")], by.x="target", by.y="sample_name", all=TRUE)
edge_df <- plyr::rename(edge_df, c("prin_graph_dim_1"="target_prin_graph_dim_1", "prin_graph_dim_2"="target_prin_graph_dim_2"))

S_matrix <- reducedDimS(cds)
data_df <- data.frame(t(S_matrix[c(x,y),]))
colnames(data_df) <- c("data_dim_1", "data_dim_2")
data_df$sample_name <- row.names(data_df)
data_df <- merge(data_df, lib_info_with_pseudo, by.x="sample_name", by.y="row.names")

return_rotation_mat <- function(theta) {
  theta <- theta / 180 * pi
  matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow = 2)
}

tmp <- return_rotation_mat(theta) %*% t(as.matrix(data_df[, c(2, 3)]))
data_df$data_dim_1 <- tmp[1, ]
data_df$data_dim_2 <- tmp[2, ]

write.csv(data_df, file = paste0("max_",m,".csv"), quote = F)
save(AGG, file = paste0("max_",m,".Robj"))

sessionInfo()
