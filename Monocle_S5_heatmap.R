# heatmap_global_ward.D2_n_4

library(viridis)
library(igraph)
library(monocle)
library(pheatmap)
library(RColorBrewer)

load(file = paste0("hm_complete.Robj"))

add_annotation_row <- NULL
num_clusters <- 4
scale_max <- 3
scale_min <- -3

annotation_colors <- hm$annotation_colors
annotation_col <- hm$annotation_col
annotation_col_B <- data.frame(annotation_col[101:200,])
colnames(annotation_col_B) <- "Cell Type"

BranchB_exprs <- hm$BranchB_exprs
heatmap_matrix <- BranchB_exprs

heatmap_matrix=heatmap_matrix[!apply(heatmap_matrix, 1, sd)==0,]
heatmap_matrix=Matrix::t(scale(Matrix::t(heatmap_matrix),center=TRUE))
heatmap_matrix=heatmap_matrix[is.na(row.names(heatmap_matrix)) == FALSE,]
heatmap_matrix[is.nan(heatmap_matrix)] = 0
heatmap_matrix[heatmap_matrix>scale_max] = scale_max
heatmap_matrix[heatmap_matrix<scale_min] = scale_min

heatmap_matrix_ori <- heatmap_matrix

row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix)))/2)
row_dist[is.na(row_dist)] <- 1

exp_rng <- range(heatmap_matrix) #bks is based on the expression range
bks <- seq(exp_rng[1] - 0.1, exp_rng[2] + 0.1, by=0.1)
hmcols <- colorRampPalette(rev(brewer.pal(n = 5, name = "RdBu")))(length(bks) - 1)

# "ward.D", "ward.D2", "single", "complete", 
# "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC)

hclust_method <- "ward.D2"
ph <- pheatmap(heatmap_matrix, 
               useRaster = T,
               cluster_cols=FALSE, 
               cluster_rows=TRUE, 
               show_rownames=F, 
               show_colnames=F, 
               #scale="row",
               clustering_distance_rows=row_dist,
               clustering_method = hclust_method,
               cutree_rows=num_clusters,
               silent=TRUE,
               filename=NA,
               breaks=bks,
               color=hmcols
               #color=hmcols#,
               # filename="expression_pseudotime_pheatmap.pdf",
)
colnames(heatmap_matrix) <- c(1:ncol(heatmap_matrix))

num_clusters <- 4
annotation_row <- data.frame(Cluster=factor(cutree(ph$tree_row, num_clusters)))
ph_res <- pheatmap(heatmap_matrix[, ], #ph$tree_row$order
                   useRaster = T,
                   cluster_cols=FALSE, 
                   cluster_rows=TRUE, 
                   # show_rownames=show_rownames, 
                   show_colnames=F, 
                   show_rownames=F, 
                   border_color=NA,
                   #scale="row",
                   clustering_distance_rows=row_dist, #row_dist
                   clustering_method = hclust_method, #ward.D2
                   cutree_rows=num_clusters,
                   # cutree_cols = 2,
                   annotation_row=annotation_row,
                   annotation_col=annotation_col_B,
                   annotation_colors=annotation_colors,
                   # gaps_col = col_gap_ind,
                   # treeheight_row = 20, 
                   breaks=bks,
                   fontsize = 6,
                   color=hmcols, 
                   silent=TRUE)
tiff(file = paste0("heatmap_global_",hclust_method,"_n_",num_clusters,".tiff"), width = 7, height = 7, units = 'in', res = 300, compression = 'none')
# pdf(file = paste0("heatmap_global_complete.pdf"), width = 7, height = 7, onefile = F)
grid::grid.rect(gp=grid::gpar("fill", col=NA))
grid::grid.draw(ph_res$gtable)
dev.off()

heatmap_matrix_clust <- cbind(heatmap_matrix, 
                              cluster = cutree(ph_res$tree_row, 
                                               k = num_clusters))

res2 <- heatmap_matrix_clust[c(ph_res$tree_row[["order"]]),]

genes <- read.table("/Users/m/tmp/data/_AGG11_170901/AGG11_mapped/outs/filtered_gene_bc_matrices_mex/mm10/genes.tsv", row.names = 1)
genes <- data.frame(genes[row.names(res2), ,drop=FALSE])

res3 <- cbind(res2, genes)

write.csv(res3, "heatmap_matrix.csv", quote = F)
write.csv(genes, "heatmap_matrix_genes.csv", quote = F)
