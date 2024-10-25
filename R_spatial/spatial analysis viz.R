library(imcRtools)
spe <- buildSpatialGraph(spe, img_id = "sample_id", type = "knn", k = 20)
spe <- buildSpatialGraph(spe, img_id = "sample_id", type = "expansion", threshold = 20)
spe <- buildSpatialGraph(spe, img_id = "sample_id", type = "delaunay", max_dist = 50)

library(ggplot2)
library(viridis)

# steinbock interaction graph 
plotSpatial(spe[,spe$sample_id == "Patient3_001"], 
            node_color_by = "celltype", 
            img_id = "sample_id", 
            draw_edges = TRUE, 
            colPairName = "neighborhood", 
            nodes_first = FALSE, 
            edge_color_fix = "grey") + 
  scale_color_manual(values = metadata(spe)$color_vectors$celltype) +
  ggtitle("steinbock interaction graph")



# knn interaction graph 
plotSpatial(spe[,spe$sample_id == "Patient3_001"], 
            node_color_by = "celltype", 
            img_id = "sample_id", 
            draw_edges = TRUE, 
            colPairName = "knn_interaction_graph", 
            nodes_first = FALSE,
            edge_color_fix = "grey") + 
  scale_color_manual(values = metadata(spe)$color_vectors$celltype) +
  ggtitle("knn interaction graph")


plotSpatial(spe[,spe$sample_id == "Patient3_001"], 
            node_color_by = "Ecad", 
            assay_type = "exprs",
            img_id = "sample_id", 
            draw_edges = TRUE, 
            colPairName = "expansion_interaction_graph", 
            nodes_first = FALSE, 
            node_size_by = "area", 
            directed = FALSE,
            edge_color_fix = "grey") + 
  scale_size_continuous(range = c(0.1, 2)) +
  ggtitle("E-cadherin expression")



plotSpatial(spe, 
            node_color_by = "celltype", 
            img_id = "sample_id", 
            node_size_fix = 0.5) + 
  scale_color_manual(values = metadata(spe)$color_vectors$celltype)





#Community detection
library(igraph)

set.seed(220819)

# Spatial community detection - tumor
tumor_spe <- spe[,spe$celltype == "Tumor"]

gr <- graph_from_data_frame(as.data.frame(colPair(tumor_spe, "neighborhood")), 
                            directed = FALSE, 
                            vertices = data.frame(index = seq_len(ncol(tumor_spe))))

cl_comm <- cluster_louvain(gr)
comm_tumor <- paste0("Tumor_", membership(cl_comm))
comm_tumor[membership(cl_comm) %in% which(sizes(cl_comm) < 10)] <- NA
names(comm_tumor) <- colnames(tumor_spe)

# Spatial community detection - non-tumor
stroma_spe <- spe[,spe$celltype != "Tumor"]

gr <- graph_from_data_frame(as.data.frame(colPair(stroma_spe, "knn_interaction_graph")), 
                            directed = FALSE, 
                            vertices = data.frame(index = seq_len(ncol(stroma_spe))))

cl_comm <- cluster_louvain(gr)
comm_stroma <- paste0("Stroma_", membership(cl_comm))
comm_stroma[membership(cl_comm) %in% which(sizes(cl_comm) < 10)] <- NA
names(comm_stroma) <- colnames(stroma_spe)

comm <- c(comm_tumor, comm_stroma)

spe$spatial_community <- comm[colnames(spe)]



library(igraph)

set.seed(220819)

# Spatial community detection - tumor
tumor_spe <- spe[,spe$celltype == "Tumor"]

gr <- graph_from_data_frame(as.data.frame(colPair(tumor_spe, "neighborhood")), 
                            directed = FALSE, 
                            vertices = data.frame(index = seq_len(ncol(tumor_spe))))

cl_comm <- cluster_louvain(gr)
comm_tumor <- paste0("Tumor_", membership(cl_comm))
comm_tumor[membership(cl_comm) %in% which(sizes(cl_comm) < 10)] <- NA
names(comm_tumor) <- colnames(tumor_spe)

# Spatial community detection - non-tumor
stroma_spe <- spe[,spe$celltype != "Tumor"]

gr <- graph_from_data_frame(as.data.frame(colPair(stroma_spe, "neighborhood")), 
                            directed = FALSE, 
                            vertices = data.frame(index = seq_len(ncol(stroma_spe))))

cl_comm <- cluster_louvain(gr)
comm_stroma <- paste0("Stroma_", membership(cl_comm))
comm_stroma[membership(cl_comm) %in% which(sizes(cl_comm) < 10)] <- NA
names(comm_stroma) <- colnames(stroma_spe)

comm <- c(comm_tumor, comm_stroma)

spe$spatial_community <- comm[colnames(spe)]

# get number of communities
as.factor(spe[,spe$celltype == "Tumor"]$spatial_community)

#the fraction of cell types within each spatial stromal community is displayed.
library(pheatmap)
library(viridis)


for_plot <- prop.table(table(spe[,spe$celltype != "Tumor"]$spatial_community, spe[,spe$celltype != "Tumor"]$celltype), margin = 1)
pheatmap(for_plot, color = viridis(100), show_rownames = FALSE)


#neigbhorhood analysis
# By celltypes
spe <- aggregateNeighbors(spe, colPairName = "knn_interaction_graph", 
                          aggregate_by = "metadata", count_by = "celltype")

set.seed(220705)

cn_1 <- kmeans(spe$aggregatedNeighbors, centers = 6)
spe$cn_celltypes <- as.factor(cn_1$cluster)

plotSpatial(spe, 
            node_color_by = "cn_celltypes", 
            img_id = "sample_id", 
            node_size_fix = 0.5) +
  scale_color_brewer(palette = "Set3")




# By expression
spe <- aggregateNeighbors(spe, colPairName = "knn_interaction_graph", 
                          aggregate_by = "expression", assay_type = "exprs",
                          subset_row = rowData(spe)$use_channel)
cn_2 <- kmeans(spe$mean_aggregatedExpression, centers = 6)
spe$cn_expression <- as.factor(cn_2$cluster)

plotSpatial(spe, 
            node_color_by = "cn_expression", 
            img_id = "sample_id", 
            node_size_fix = 0.5) +
  scale_color_brewer(palette = "Set3")


## LISA clusters
library(lisaClust)
library(spicyR)

cells <- data.frame(row.names = colnames(spe))
cells$ObjectNumber <- spe$ObjectNumber
cells$ImageNumber <- spe$sample_id
cells$AreaShape_Center_X <- spatialCoords(spe)[,"Pos_X"]
cells$AreaShape_Center_Y <- spatialCoords(spe)[,"Pos_Y"]
cells$cellType <- spe$celltype

lisa_sc <- SegmentedCells(cells, cellProfiler = TRUE)

lisa_sc

lisaCurves <- lisa(lisa_sc, Rs = c(10, 20, 50))

# Set NA to 0
lisaCurves[is.na(lisaCurves)] <- 0

lisa_clusters <- kmeans(lisaCurves, centers = 6)$cluster

spe$lisa_clusters <- as.factor(lisa_clusters)

plotSpatial(spe, 
            node_color_by = "lisa_clusters", 
            img_id = "sample_id", 
            node_size_fix = 0.5) +
  scale_color_brewer(palette = "Set3")


