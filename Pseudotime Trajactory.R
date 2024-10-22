
library("monocle3")                                           
library(Seurat)
library("ggplot2")
cds <- new_cell_data_set(as(as.matrix(pbmc@assays$RNA@counts),"sparseMatrix"),
                         cell_metadata = pbmc@meta.data,
                         gene_metadata = data.frame(gene_short_name = row.names(pbmc),row.names = row.names(as.matrix(pbmc@assays$RNA@counts))))

cds <- preprocess_cds(cds,num_dim = 100)
cds <- align_cds(cds,alignment_group = "orig.ident")
cds <- reduce_dimension(cds,cores = 5,reduction_method = "UMAP")
cds <- cluster_cells(cds,reduction_method = "UMAP")
cds <- learn_graph(cds)
cds <- order_cells(cds)
cds_DGT_pseudotimegenes <- differentialGeneTest(cds,fullModelFormulaStr = "~sm.ns(Pseudotime)")
plot_cells(cds)
p <- plot_cells(cds,color_cells_by = "pseudotime",
           label_cell_groups = F, 
           label_leaves = F,
           label_branch_points = F,
           label_roots = T,
           graph_label_size = 5)+
  scale_color_continuous(name = "Pseudotime",low = "blue3",high = "orange")+
  theme(legend.position=c(0.9,0.8))

table1 <- data.frame(p$data$sample_name,p$data$cells,p$data$data_dim_1,p$data$data_dim_2,p$data$cell_color)
colnames(table1) <- c("cell ID","cluster","dim1","dim2","pseudotime data")
table1$cluster <- factor(table1$cluster,levels = c("TELCs","VCT like cell","unknown"))
table1 <- table1[order(table1$cluster),]

table2 <- data.frame(p$data$sample_name,p$data$cells,p$data$data_dim_1,p$data$data_dim_2,exprData["MXRA5",])
colnames(table2) <- c("cell ID","cluster","dim1","dim2","MXRA5")
table2$cluster <- factor(table2$cluster,levels = c("TELCs","VCT like cell","unknown"))
table2 <- table2[order(table2$cluster),]

plot_cells(cds,color_cells_by = "cells",
           show_trajectory_graph = F,
           label_cell_groups = F, 
           label_leaves = F,
           label_branch_points = F,
           label_roots = F,
           graph_label_size = 5)+
  scale_colour_manual(values=c("blue3","orange","gray"))+
  theme(legend.position=c(0.9,0.8))

exprData <- cds@assays@data@listData$counts
exprData_normalize <- LogNormalize(exprData) 
exprData_translog <- log10(exprData)

cds$MXRA5 <- exprData_normalize["MXRA5",]
cds$VGLL1 <- exprData_normalize["VGLL1",]
cds$NRP2 <- exprData_normalize["NRP2",]
ciliated_genes <- c("MXRA5","VGLL1","NRP2")
plot_genes_in_pseudotime(cds[ciliated_genes],color_cells_by = "cells",min_expr=0.5,)+
  scale_colour_manual(name = "",values=c("blue3","orange","gray"))+
  scale_y_continuous(limits=c(0,3),breaks = c(0,1,2,3))+
  theme(legend.position="right")
  
plot_cells(cds,color_cells_by = "MXRA5",
           cell_size=1 ,
           label_cell_groups = F, 
           label_leaves = F,
           label_branch_points = F,
           label_roots = T,
           graph_label_size = 5)+
  scale_color_gradient(name="value",low = "grey",high = "#006400")+
  theme(legend.position=c(0.9,0.8))