library(viridis)
library(monocle)

load('t-SNE.Robj')

Batch_vec <- c('MEF', 'SI5D', 'SI12D', 'XEN', 'SII8D', 'SII12D', 'SIII3D', 'SIII6D', 'SIII8D', 'SIII10D', 'SIII15D', 'SIII21D')
Batch_cols <- viridis(12, option = 'viridis')
names(Batch_cols) <- Batch_vec
Batch_labels <- c('MEF', 'SI D5', 'SI D12', 'SI D16', 'SII D8', 'SII D12', 'SIII D3', 'SIII D6', 'SIII D8', 'SIII D10', 'SIII D15', 'SIII D21')

tiff(file = paste0('t-SNE','_ordering','.tiff'), width = 7, height = 7, units = 'in', res = 300, compression = 'none')
plot_ordering_genes(AGG)
dev.off()

tiff(file = paste0('t-SNE','_cell_clusters_Batch','.tiff'), width = 7, height = 7, units = 'in', res = 300, compression = 'none')
plot_cell_clusters(AGG, color_by = 'as.factor(Batch)', cell_size = 0.3) + scale_color_manual(labels = Batch_labels, values = Batch_cols, name = 'Batch') + guides(col = guide_legend(nrow = 3, byrow = TRUE))
dev.off()

pdf(file = paste0('t-SNE','_cell_clusters_Batch','.pdf'), width = 7, height = 7)
plot_cell_clusters(AGG, color_by = 'as.factor(Batch)', cell_size = 0.3) + scale_color_manual(labels = Batch_labels, values = Batch_cols, name = 'Batch') + guides(col = guide_legend(nrow = 3, byrow = TRUE))
dev.off()


tiff(file = paste0('t-SNE','_cell_clusters_Batch_facet','.tiff'), width = 7, height = 7, units = 'in', res = 300, compression = 'none')
plot_cell_clusters(AGG, color_by = 'as.factor(Batch)', cell_size = 0.3) + facet_wrap(~Batch) + scale_color_manual(labels = Batch_labels, values = Batch_cols, name = 'Batch') + guides(col = guide_legend(nrow = 3, byrow = TRUE))
dev.off()

sessionInfo()
