library(Rphenograph)
library(igraph)
library(dittoSeq)
library(viridis)
library(bluster)
library(BiocParallel)


set.seed(220619)
cur_cells <- sample(seq_len(ncol(spe)), 2000)


mat <- t(assay(spe, "exprs")[rowData(spe)$use_channel,])

out <- Rphenograph(mat, k = 45)

clusters <- factor(membership(out[[2]]))

spe$pg_clusters <- clusters


##plot UMAP
dittoDimPlot(spe, var = "pg_clusters", 
             reduction.use = "UMAP", size = 0.2,
             do.label = TRUE) +
  ggtitle("Phenograph clusters expression on UMAP")


##Heatmap
dittoHeatmap(spe[,cur_cells], 
             genes = rownames(spe)[rowData(spe)$use_channel],
             assay = "exprs", scale = "none",
             heatmap.colors = viridis(100), 
             annot.by = c("pg_clusters", "patient_id"),
             annot.colors = c(dittoColors(1)[1:length(unique(spe$pg_clusters))],
                              metadata(spe)$color_vectors$patient_id))


## Integrated low dim embedding data

mat <- reducedDim(spe, "fastMNN")

out <- Rphenograph(mat, k = 45)

clusters <- factor(membership(out[[2]]))

spe$pg_clusters_corrected <- clusters

dittoDimPlot(spe, var = "pg_clusters_corrected", 
             reduction.use = "UMAP_mnnCorrected", size = 0.2,
             do.label = TRUE) +
  ggtitle("Phenograph clusters expression on UMAP, integrated cells")


sam_cells <- sample(seq_len(ncol(spe)), 10000)
mat <- t(assay(spe, "exprs")[rowData(spe)$use_channel,])



# tuning clustering settings
x <- clusterSweep(mat[sam_cells,], BLUSPARAM=KNNGraphParam(),
             k=c( 100L, 200L, 300L, 400L, 700L, 1000L), 
             type = "jaccard", 
             cluster.fun="louvain",
             BPPARAM = MulticoreParam(RNGseed = 220427))


sil <- vapply(as.list(x$clusters), 
              function(x) mean(approxSilhouette(mat[sam_cells,], x)$width), 
              0)



ggplot(data.frame(method = names(sil),
                  sil = sil)) +
  geom_point(aes(method, sil)) +
  theme_classic(base_size = 15) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Cluster parameter combination") +
  ylab("Average silhouette width")
