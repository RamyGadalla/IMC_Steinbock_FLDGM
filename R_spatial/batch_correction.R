library(batchelor)
library(scater)

set.seed(220228)

out <- fastMNN(spe, batch = spe$patient_id,
               auto.merge = TRUE,
               subset.row = rowData(spe)$use_channel,
               assay.type = "exprs")
spe

# Transfer the correction results to the main spe object
reducedDim(spe, "fastMNN") <- reducedDim(out, "corrected")

spe

#running UMAP for the corrected values and the uncorrected values
library(scater)
set.seed(220228)
spe <- runUMAP(spe, dimred= "fastMNN", name = "UMAP_mnnCorrected") 
spe <- runUMAP(spe, exprs_values="exprs" , name = "UMAP")

spe

saveRDS(spe, "/Volumes/GoogleDrive/My Drive/spatial/RDS/batch_corrected_spe_stein.rds")

