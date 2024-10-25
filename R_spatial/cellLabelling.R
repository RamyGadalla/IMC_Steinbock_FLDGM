library(cytomapper)
celltype <- setNames(c("#3F1B03", "#F4AD31", "#894F36", "#1C750C", "#EF8ECC", 
                       "#6471E2", "#4DB23B", "grey", "#F4800C", "#BF0A3D", "#066970"),
                     c("Tumor", "Stroma", "Myeloid", "CD8", "Plasma_cell", 
                       "Treg", "CD4", "undefined", "BnTcell", "Bcell", "Neutrophil"))

metadata(spe)$color_vectors$celltype <- celltype


library(SingleCellExperiment)
label_files <- list.files("/Volumes/GoogleDrive/My Drive/spatial/R_spatial/steinbock/gated_cells", 
                          full.names = TRUE, pattern = ".rds$")

# Read in SPE objects
spes <- lapply(label_files, readRDS)

# Merge SPE objects
concat_spe <- do.call("cbind", spes)


cur_tab <- unclass(table(colnames(concat_spe), 
                         concat_spe$cytomapper_CellLabel))
cur_labels <- rep("doublets", nrow(cur_tab))
names(cur_labels) <- rownames(cur_tab)

# Single assignments
single_index <- rowSums(cur_tab) == 1
cur_labels[single_index] <- colnames(cur_tab)[apply(cur_tab[single_index,], 1, 
                                                    which.max)]

# Double assignment within the tumor
double_index <- rowSums(cur_tab) == 2 & cur_tab[,"Tumor"] == 1
no_tumor <- cur_tab[,colnames(cur_tab) != "Tumor"]
cur_labels[double_index] <- colnames(no_tumor)[apply(no_tumor[double_index,], 1, 
                                                     which.max)]

# Remove doublets
cur_labels <- cur_labels[cur_labels != "doublets"]
table(cur_labels)




# Transfer labels to SPE object
spe_labels <- rep("unlabeled", ncol(spe))
names(spe_labels) <- colnames(spe)
spe_labels[names(cur_labels)] <- cur_labels
spe$cell_labels <- spe_labels

# Number of cells labeled per patient
table(spe$cell_labels, spe$patient_id)






library(caret)

# Split between labeled and unlabeled cells
lab_spe <- spe[,spe$cell_labels != "unlabeled"]
unlab_spe <- spe[,spe$cell_labels == "unlabeled"]

# Randomly split into train and test data
set.seed(220925)
trainIndex <- createDataPartition(factor(lab_spe$cell_labels), p = 0.75)
train_spe <- lab_spe[,trainIndex$Resample1]
test_spe <- lab_spe[,-trainIndex$Resample1]

# Specify train parameters for 5-fold cross validation
fitControl <- trainControl(method = "cv",
                           number = 5)

# Select the data for training
cur_mat <- t(assay(train_spe, "exprs")[rowData(train_spe)$use_channel,])

# Train a random forest model for predicting cell labels
# This call also performs parameter tuning
rffit <- train(x = cur_mat, 
               y = factor(train_spe$cell_labels),
               method = "rf", ntree = 1000,
               tuneLength = 5,
               trControl = fitControl)

rffit

# Select test data
cur_mat <- t(assay(test_spe, "exprs")[rowData(test_spe)$use_channel,])

# Predict cell phenotypes in test data
cur_pred <- predict(rffit, 
                    newdata = cur_mat)


cm <- confusionMatrix(data = cur_pred, 
                      reference = factor(test_spe$cell_labels), 
                      mode = "everything")

cm

# Select unlabeled data
cur_mat <- t(assay(unlab_spe, "exprs")[rowData(unlab_spe)$use_channel,])

# Predict cell phenotypes
cell_class <- as.character(predict.train(rffit, 
                                         newdata = cur_mat, 
                                         type = "raw"))
names(cell_class) <- rownames(cur_mat)

table(cell_class)


# Extract classification probabilities
cell_prob <- predict.train(rffit, 
                           newdata = cur_mat, 
                           type = "prob")

# Label undefined cells
cell_class[rowMax(as.matrix(cell_prob)) < 0.4] <- "undefined"

# Store labels in SpatialExperiment onject
cell_labels <- spe$cell_labels
cell_labels[colnames(unlab_spe)] <- cell_class
spe$celltype <- cell_labels 

table(spe$celltype, spe$patient_id)