library(stringr)
library(openxlsx)
library(imcRtools)


read_steinbock("/Volumes/GoogleDrive/My Drive/spatial/R_spatial/steinbock/") -> spe_stein

meta <- read.xlsx("/Volumes/GoogleDrive/My Drive/spatial/sample_metadata.xlsx")

spe_stein$patient_id <- as.vector(str_extract_all(spe_stein$sample_id, "Patient[1-4]", simplify = TRUE))
spe_stein$ROI <- as.vector(str_extract_all(spe_stein$sample_id, "00[1-8]", simplify = TRUE))
spe_stein$indication <- meta$Indication[match(spe_stein$patient_id, meta$Sample.ID)]


colnames(spe_stein) <- paste0(spe_stein$sample_id, "_", spe_stein$ObjectNumber)


#add Boolean vector to the rowData markers to use are TRUE
rowData(spe_stein)$use_channel <- !grepl("DNA|Histone", rownames(spe_stein))


saveRDS(spe_stein,'/Volumes/GoogleDrive/My Drive/spatial/RDS/spe_steinbock_use.rds')
