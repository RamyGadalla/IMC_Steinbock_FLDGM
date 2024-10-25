library(BiocParallel)
library(imcRtools)
library(CATALYST)

#read txt file from single metal spot agrose gel
readSCEfromTXT("./compensation/") -> sce
assay(sce, "exprs") <- asinh(counts(sce)/5)

#Filtering incorrectly assigned pixels
bc_key <- as.numeric(unique(sce$sample_mass))
bc_key <- bc_key[order(bc_key)]
sce <- assignPrelim(sce, bc_key = bc_key)
sce <- estCutoffs(sce)
sce <- applyCutoffs(sce)


#compute spillover matrix
sce <- computeSpillmat(sce)
isotope_list <- CATALYST::isotope_list
isotope_list$Ar <- 80
plotSpillmat(sce, isotope_list = isotope_list)
# Save spillover matrix in new object
sm <- metadata(sce)$spillover_matrix


#I am using here Stienbock pipeline output

#Single cell data compensation
#spe <- readRDS("data/spe.rds")
rowData(spe)$channel_name <- paste0(rowData(spe)$channel, "Di")
#this function add the compensation data to slot assay
spe <- compCytof(spe, sm, 
                 transform = TRUE, cofactor = 1,
                 isotope_list = isotope_list, 
                 overwrite = FALSE)
spe
#remove uncompensated data and replace it with the comp data keeping same name
assay(spe, "counts") <- assay(spe, "compcounts") 
assay(spe, "exprs") <- assay(spe, "compexprs") 
assay(spe, "compcounts") <- assay(spe, "compexprs") <- NULL


#Image compensation
#images <- readRDS("data/images.rds")
channelNames(images) <- rowData(spe)$channel_name
panel <- read.csv("")
adapted_sm <- adaptSpillmat(sm, paste0(panel$channel[panel$keep == 1], "Di"), 
                            isotope_list = isotope_list)
images_comp <- compImage(images, adapted_sm, BPPARAM = MulticoreParam())

# Before compensation
plotPixels(images[5], colour_by = "Yb173Di", 
           image_title = list(text = "Yb173 (Ecad) - before", position = "topleft"), 
           legend = NULL, bcg = list(Yb173Di = c(0, 4, 1)))
plotPixels(images[5], colour_by = "Yb174Di", 
           image_title = list(text = "Yb174 (CD303) - before", position = "topleft"), 
           legend = NULL, bcg = list(Yb174Di = c(0, 4, 1)))

# After compensation
plotPixels(images_comp[5], colour_by = "Yb173Di",
           image_title = list(text = "Yb173 (Ecad) - after", position = "topleft"), 
           legend = NULL, bcg = list(Yb173Di = c(0, 4, 1)))
plotPixels(images_comp[5], colour_by = "Yb174Di", 
           image_title = list(text = "Yb174 (CD303) - after", position = "topleft"),
           legend = NULL, bcg = list(Yb174Di = c(0, 4, 1)))


#for convenience, reset channelnames to their biological name.
channelNames(images_comp) <- rownames(spe)

#write out/save compensated image - takes time

dir.create("/Volumes/GoogleDrive/My Drive/spatial/R_spatial/comp_img")
lapply(names(images_comp), function(x){
  writeImage(as.array(images_comp[[x]])/(2^16 - 1), 
             paste0("", x, ".tiff"),
             bits.per.sample = 16)
})

#saving rds object of compensated single cell data
saveRDS(spe, "")
saveRDS(images_comp, "")






