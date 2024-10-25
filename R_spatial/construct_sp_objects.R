# packages
library(imcRtools)
library(cytomapper)
library(openxlsx)
library(dittoSeq)
library(patchwork)
library(stringr)
library(RColorBrewer)

#read data into spe object
spe <- read_cpout("")
rownames(spe) <- rowData(spe)$Clean_Target
spe

colnames(spe) <- paste0(spe$sample_id, "_", spe$ObjectNumber)

#read extra metadata and add it
meta <- read.xlsx("")

spe$indication <- meta$Indication[match(spe$Metadata_acname, meta$Sample.ID)]
colData(spe)$patient_id <- colData(spe)$Metadata_acname
colData(spe)$ROI <- str_extract(colData(spe)$Metadata_description, ".$")

colData(spe)


#arcsin trans cofactor 1
assay(spe, "exprs") <- asinh(counts(spe)/1)

#before tansformation and after
p1 <- dittoRidgePlot(spe, var = "CD3", group.by = "patient_id", assay = "counts") +
  ggtitle("CD3 - before transformation")
p2 <- dittoRidgePlot(spe, var = "CD3", group.by = "patient_id", assay = "exprs") +
  ggtitle("CD3 - after transformation")

p1+p2

#add Boolean vector to the rowData markers to use are TRUE
rowData(spe)$use_channel <- !grepl("DNA|Histone", rownames(spe))


#unique color for different metadata entries
color_vectors <- list()

ROI <- setNames(brewer.pal(length(unique(spe$ROI)), name = "BrBG"), 
                unique(spe$ROI))
patient_id <- setNames(brewer.pal(length(unique(spe$patient_id)), name = "Set1"), 
                       unique(spe$patient_id))
sample_id <- setNames(dittoColors(reps = 1)[seq_along(unique(spe$sample_id))], 
                      unique(spe$sample_id))
indication <- setNames(brewer.pal(length(unique(spe$indication)), name = "Set2"), 
                       unique(spe$indication))

color_vectors$ROI <- ROI
color_vectors$patient_id <- patient_id
color_vectors$sample_id <- sample_id
color_vectors$indication <- indication

## Images
#Add images to CytoImagelist (one for masks and one for multi channel tiff images )
images <- loadImages("")
masks <- loadImages("", as.is = TRUE)

channelNames(images) <- rownames(spe)
all.equal(channelNames(images), rownames(rowData(spe)))
images
masks


#make sure images in both images and masks are in the same order for the next step
names(images) -> names(masks)


#add metadata to CytoImageList
patient_id <- str_extract(names(images), "^.{0,8}")
indication <- meta$Indication[match(patient_id, meta$Sample.ID)]

mcols(images) <- mcols(masks) <- DataFrame(sample_id = names(images),
                                           patient_id = patient_id,
                                           indication = indication)



#alternatively you can generate sce directly from images through measureobject () in cytomapper
#cytomapper_sce <- measureObjects(masks, image = images, img_id = "sample_id")


#save RDS object
saveRDS(spe,"")
saveRDS(masks, "")
saveRDS(images, "")



