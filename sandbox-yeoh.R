library(Biocomb)

library(reshape2)
library(rgl)
library(ggplot2)
library(cowplot)
library(reshape2)
library(RColorBrewer)
library(dplyr)
library(pheatmap)
# library(Rtsne)
# library(dendextend)
# library(cluster)

library(sva)  # ComBat
library(scran)  # MNN
library(Harman)

# library(xtable)
# library(gPCA)

source("../functions.R")
source("bin/bcm.R")

# theme_set(theme_dark())
# theme_set(theme_gray())
theme_set(theme_cowplot())

# FUNCTIONS ---------------------------------------------------------------
plot_mean <- function(df, batch_vec1) {
  # Melt dataframe
  melt_df <- melt(df, variable.name = "ID")
  print(head(melt_df))
  # Trimmed mean probe intensities for each chip
  mean_tibble <- melt_df %>% group_by(ID) %>%
    summarise(mean = mean(value))
  mean_batch_tibble <- cbind(mean_tibble,
                             batch_vec1 = batch_vec1[mean_tibble$ID])
  
  mean_scatter <- ggplot(mean_batch_tibble, aes(x = ID, y = mean)) +
    geom_point(aes(col = factor(batch_vec1)),
               show.legend = F, size = 3) +
    facet_wrap(factor(batch_vec1), scales = "free_x") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  return(mean_scatter)
}

# Selecting drug responsive genes between D0 and D8
# using paired t-test and logfc
selectFeatures <- function(df1, metadata_df,
                           alpha = 0.05, logfc_threshold = 1) {
  # Subset df according to D0 and D8
  class_info <- metadata_df[colnames(df1), "class_info"]
  df_d0 <- df1[,class_info == "D0"]
  df_d8 <- df1[,class_info == "D8"]
  print(head(colnames(df_d0)))
  print(head(colnames(df_d8)))
  stopifnot(ncol(df_d0) == ncol(df_d8))
  
  # Identify drug responsive probesets
  ttest_pvalue <- calc_ttest(cbind(df_d0, df_d8), ncol(df_d0), is_paired = T)
  log_fc <- rowMeans(df_d8) - rowMeans(df_d0)
  ttest_probesets <- names(ttest_pvalue)[ttest_pvalue <= alpha]
  fc_probesets <- names(log_fc)[log_fc > logfc_threshold]
  intersect_probesets <- intersect(ttest_probesets, fc_probesets)
  print(paste("T-test:", length(ttest_probesets)))
  print(paste("Log fold change:", length(fc_probesets)))
  print(paste("Intersection:", length(intersect_probesets)))
  return(intersect_probesets)
}

# All dataframes have samples in rows and features in columns
# D0 centroid used to define D0-Normal vector
calcERM <- function(response_df, normal_df, labels_df) {
  # Split response df into D0 and D8 df
  n <- nrow(response_df)/2
  d0_df <- response_df[1:n,]
  d8_df <- response_df[-(1:n),]
  stopifnot(substring(rownames(d8_df),1,4) == substring(rownames(d0_df),1,4))
  
  # Calculate centroids
  leuk_centroid <- apply(d0_df, 2, median)
  normal_centroid <- apply(normal_df, 2, median)
  
  # Calculate leuk-normal unit vector
  leuk_normal <- normal_centroid - leuk_centroid
  unit_leuk_normal <- leuk_normal/calcL2Norm(leuk_normal)
  
  # Assume that patients from top rows match correspondingly with bottom rows
  # Calculate vector by: D8-D0
  d0_d8_hstack <- d8_df - d0_df
  # Multiplication of erm_factor is propagated through every column
  ### ERM1 ###
  erm1 <- colSums(t(d0_d8_hstack) * unit_leuk_normal)
  # Vertical stack of individual D0-Normal vectors
  d0_normal_vstack <- normal_centroid - t(d0_df)
  ### D0-Normal projection ###
  d0_normal_proj <- colSums(d0_normal_vstack * unit_leuk_normal)
  ### ERM1 Ratio ###
  erm1_ratio <- erm1/d0_normal_proj
  
  d8_normal_vstack <- normal_centroid - t(d8_df)
  ### D8-Normal projection ###
  d8_normal_proj <- colSums(d8_normal_vstack * unit_leuk_normal)
  
  stopifnot(identical(names(erm1), names(erm1_ratio)))
  
  # Calculate vstack of unit D0-Normal vectors
  l2norm_d0_normal <- apply(d0_normal_vstack, 2, calcL2Norm)
  unit_d0_normal_vstack <- sweep(d0_normal_vstack, 2, l2norm_d0_normal, "/")
  
  ### ERM2 ###
  erm2 <- colSums(t(d0_d8_hstack) * unit_d0_normal_vstack)
  erm2_ratio <- erm2/l2norm_d0_normal
  
  stopifnot(identical(names(erm2), names(erm2_ratio)))
  
  ### ERM3 ###
  erm3 <- d8_df[,1] - d0_df[,1]
  # Divide by D0-Normal along PC1
  erm3_ratio <- erm3/(median(normal_df[,1]) - d0_df[,1])
  
  stopifnot(identical(names(erm3), names(erm3_ratio)))
  
  ### l2norm ###
  l2norm_d0_d8 <- apply(d0_d8_hstack, 1, calcL2Norm)
  l2norm_d0 <- apply(d0_df, 1, calcL2Norm)
  l2norm_d8 <- apply(d8_df, 1, calcL2Norm)
  diff_l2norm <- l2norm_d8 - l2norm_d0
  
  ### Angle between D0 and D8 ###
  angle_d0_d8 <- mapply(calcAngleVectors,
                        data.frame(t(d0_df)), data.frame(t(d8_df)))
  
  ### Angle between D0-D8 and Leuk-Normal ###
  angle_d0d8_normal <- apply(
    d0_d8_hstack, 1, function(row_vec) calcAngleVectors(row_vec, leuk_normal)
  )
  
  ### Angle between D0 and normal ###
  angle_d0_normal <- apply(
    d0_df, 1, function(row_vec) calcAngleVectors(row_vec, normal_centroid)
  )
  
  ### Angle between D8 and Normal ###
  angle_d8_normal <- apply(
    d8_df, 1, function(row_vec) calcAngleVectors(row_vec, normal_centroid)
  )
  
  ### Concatenate all features ###
  features_df <- cbind(erm1, erm1_ratio, erm2, erm2_ratio, erm3, erm3_ratio,
                       d0_normal_proj, d8_normal_proj,
                       l2norm_d0_d8, diff_l2norm, angle_d0_d8,
                       angle_d0d8_normal, angle_d0_normal, angle_d8_normal)
  row_idxstr <- substring(rownames(features_df),1,4)
  
  ### RESULTS DF ###
  # Concatenate features with labels
  results_df <- cbind(features_df, labels_df[row_idxstr, c(1,5)])
  rownames(results_df) <- row_idxstr
  
  return(results_df)
}

# 3D PCA plot
plotPCA3D <- function(df, colour, pch, pc_labels = NULL,
                      ratio_list = list(2,1,1)) {
  if (is.null(pc_labels)) {
    print("PCA performed!")
    pca_obj <- prcomp(t(df), center = T, scale. = F)
    pca_df <- as.data.frame(pca_obj$x[,1:3])
    eigenvalues <- (pca_obj$sdev)^2
    var_pc <- eigenvalues[1:3]/sum(eigenvalues)
    print(var_pc)
    pc_labels <- sprintf("PC%d (%.2f%%)", 1:3, var_pc*100)
  } else {
    print("No PCA performed!")
    pca_df <- as.data.frame(df)
  }
  
  # RGL plot parameters
  rgl.open()
  rgl.bg(color="white")
  rgl.viewpoint(zoom = 0.8)
  # rgl.viewpoint(theta = 110, phi = 5, zoom = 0.8)
  par3d(windowRect = c(50, 20, 500, 500))
  with(pca_df, pch3d(PC1, PC2, PC3, bg = colour,
                     pch = pch, cex = 0.5, lwd = 1.5))
  box3d(col = "black")
  title3d(xlab = pc_labels[1], ylab = pc_labels[2],
          zlab = pc_labels[3], col = "black")
  # Plot aspect ratios of axis according to variance
  do.call(aspect3d, ratio_list)
}

# Plot PCA before selecting features
# Batch information of all the timepoints
plotPCA3DYeoh <- function(df1, metadata_df) {
  batch_info <- metadata_df[colnames(df1), "batch_info"]
  generate_colour <- colorRampPalette(c("lightblue", "darkblue"))
  batch_palette <- generate_colour(10)
  # batch_palette <- brewer.pal(10, "Set3")
  batch_colour <- batch_palette[batch_info]
  # Shape of all timepoints
  class_info <- metadata_df[colnames(df1), "class_info"]
  levels(class_info) <- 21:23
  timepoint_shape <- as.numeric(as.character(class_info))
  plotPCA3D(df1, batch_colour, timepoint_shape)
}

# Plot PCA before selecting features
# Batch information of all the timepoints
plotPCA3DYeoh1 <- function(df1, metadata_df) {
  batch_info <- metadata_df[colnames(df1), "batch_info"]
  batch_factor <- droplevels(as.factor(batch_info))
  print(batch_factor)
  print(levels(batch_factor))
  levels(batch_factor) <- 21:22
  pch <- as.numeric(as.character(batch_factor))
  # generate_colour <- colorRampPalette(c("lightblue", "darkblue"))
  # batch_palette <- generate_colour(10)
  
  # Shape of all timepoints
  class_info <- metadata_df[colnames(df1), "subtype"]
  palette <- brewer.pal(10, "Set3")
  col <- palette[class_info]
  
  plotPCA3D(df1, col, pch)
}

plot.vectors_3d <- function(vectors_arr, colour_code, shape_vec,
                            pc_labels = NULL, ratio_list = list(2,1,1)) {
  if (is.null(pc_labels)) {
    pca_obj <- prcomp(t(df))
    pca_df <- as.data.frame(pca_obj$x[,1:3])
    eigenvalues <- (pca_obj$sdev)^2
    var_pc <- eigenvalues[1:3]/sum(eigenvalues)
    pc_labels <- sprintf("PC%d (%.2f%%)", 1:3, var_pc*100)
  } else pca_df <- as.data.frame(df)
  # RGL plot parameters
  rgl.open()
  rgl.bg(color="white")
  rgl.viewpoint(zoom = 0.8)
  # rgl.viewpoint(theta = 110, phi = 5, zoom = 0.8)
  # par3d(windowRect = c(50, 20, 500, 500))
  label_colour <- ifelse(vectors_arr[,7] == 0, "black", "red")
  arrows_df <- data.frame(t(vectors_arr[,1:6]))
  mapply(arrow3d, arrows_df[1:3,], arrows_df[4:6,],
         s = 0.05, type = "line", col = label_colour)
  box3d(col = "black", expand = 0.5)
  title3d(xlab = pc_labels[1], ylab = pc_labels[2], zlab = pc_labels[3],
          col = "black")
  # Plot aspect ratios of axis according to variance
  do.call(aspect3d, ratio_list)
}

# Plot vectors 2D
plot_vectors <- function(df, centroid_df, pc_labels,
                         batch_colour, subtype_colour) {
  pca_1 <- ggplot(data = df) +
    geom_point(aes(x = PC1_A, y = PC2_A), colour = subtype_colour, shape = 21,
               size = 5, stroke = 2, fill = batch_colour, show.legend = F) +
    geom_point(aes(x = PC1_B, y = PC2_B), colour = subtype_colour, shape = 22,
               size = 5, stroke = 2, fill = batch_colour, show.legend = F) +
    geom_point(data = centroid_df, aes(x = PC1, y = PC2),
               size = 5, shape = 17,
               colour = c("tomato3", "orange", "hotpink", "pink")) +
    geom_segment(aes(x = PC1_A, y = PC2_A, xend = PC1_B,
                     yend = PC2_B, colour = as.factor(labels)),
                 arrow = arrow(length = unit(0.3, "cm")),
                 alpha = 0.5, show.legend = F) +
    scale_color_manual(values = c("black",  "red")) +
    xlab(pc_labels[1]) + ylab(pc_labels[2])
  
  pca_2 <- ggplot(data = df) +
    geom_point(aes(x = PC3_A, y = PC2_A), colour = subtype_colour, shape = 21,
               size = 5, stroke = 2, fill = batch_colour, show.legend = F) +
    geom_point(aes(x = PC3_B, y = PC2_B), colour = subtype_colour, shape = 22, 
               size = 5, stroke = 2, fill = batch_colour, show.legend = F) +
    geom_point(data = centroid_df, aes(x = PC2, y = PC3),
               size = 5, shape = 17,
               colour = c("tomato3", "orange", "hotpink", "pink")) +
    geom_segment(aes(x = PC3_A, y = PC2_A, xend = PC3_B, yend = PC2_B,
                     colour = as.factor(labels)),
                 arrow = arrow(length = unit(0.3, "cm")), alpha = 0.5,
                 show.legend = F) +
    scale_color_manual(values = c("black",  "red")) +
    xlab(pc_labels[3]) + ylab(pc_labels[2])
  
  multiplot <- plot_grid(pca_1, pca_2, ncol = 2)
  return(multiplot)
}

# IMPORT DATA -------------------------------------------------------------
## Subset of original data
# Removed outliers, patients with timepoints from different batches and batch 5
SUBSET_RPATH <- "data/GSE67684/processed/subset_yeoh.tsv"
subset_yeoh <- read.table(SUBSET_RPATH, sep = "\t")

## Metadata
# Preprocessed metadata
METADATA_RPATH <- "data/GSE67684/processed/metadata/full_metadata.tsv"
metadata_df <- read.table(METADATA_RPATH, sep = "\t")

# BATCH_RPATH <- "data/GSE67684/processed/metadata/metadata-batch.tsv"
# LABEL_RPATH <- "data/GSE67684/processed/metadata/metadata-label_mrd_subtype.tsv"
# yeoh_batch <- read.table(BATCH_RPATH, sep = "\t", header = T, row.names = 1)
# yeoh_label <- read.table(LABEL_RPATH, sep = "\t", header = T, row.names = 1)

# SCALE->REMOVE->FILTER->LOG
scaled_yeoh <- removeProbesets(normaliseMeanScaling(subset_yeoh))
data_yeoh <- log2_transform(filterProbesets(scaled_yeoh, 0.7, metadata_df))

# # Filter out all rows with zero values
# logi_idx <- rowSums(data_yeoh == 0) == 0
# filtered_yeoh <- data_yeoh[logi_idx,]

# QUANTILE ---------------------------------------------------------
# ### QUANTILE (ALL)
# quantile_yeoh <- normaliseQuantile(subset_yeoh)
colnames(subset_yeoh[403:405])
### QUANTILE (TIMEPOINT)
quantile_d0 <- normaliseQuantile(subset_yeoh[, 1:201])
quantile_d8 <- normaliseQuantile(subset_yeoh[, 202:402])
quantile_normal <- normaliseQuantile(subset_yeoh[, 403:405])
quantile_yeoh <- cbind(quantile_d0, quantile_d8, quantile_normal)

write.table(quantile_yeoh,
            "dump/quantile_yeoh.tsv",
            quote = F, sep = "\t")

# # Evaluation
# calc.var_preservation(log_yeoh, quantile_yeoh)
# metrics <- eval.batch_effects(quantile_yeoh, batch_info, class_numeric)
# plotPCA3DYeoh(quantile_yeoh, metadata_df)
# rgl.postscript("dump/pca_3d-cs_quantile.pdf", "pdf")

intersect_probesets <- selectFeatures(quantile_yeoh, metadata_df)
selected_quantile <- quantile_yeoh[intersect_probesets,]
results_df <- calcERM_PCA(selected_quantile, metadata_df)

# Quantile (Class) - ComBat -----------------------------------------------
colnames(log_yeoh)

# Obtaining batch information of selected_yeoh df
# Rows of metadata are to be in same order as columns of edata
batch_info <- yeoh_batch[colnames(log_yeoh), "batch"]
timepoint_info <- as.factor(rep(c("D0","D8","N"), c(208,208,45)))

yeoh_metadata <- data.frame(batch = as.factor(batch_info),
                            timepoint = timepoint_info)
# Place adjustment/confounding variables in model.matrix (e.g. age)
# Do not put batch variables in model.matrix
# # # Put batch variables directly in combat function
# model_combat <- model.matrix(~1, data = yeoh_metadata)
# Include biological variable of interest as covariate
model_combat <- model.matrix(~timepoint, data = yeoh_metadata)
combat_yeoh <- ComBat(data.matrix(quantile_yeoh), batch_info, model_combat)

# Replacing negative values with 0
combat_yeoh[combat_yeoh < 0] <- 0

# Plot PCA before selecting features
plotPCA3DYeoh(combat_yeoh, metadata_df)
rgl.postscript("dump/pca_3d-cs_quantile_combat_cov.pdf", "pdf")

calc.var_preservation(log_yeoh, combat_yeoh)
metrics <- eval.batch_effects(combat_yeoh, batch_info, class_numeric)

# Selecting drug responsive genes between D0 and D8
ttest_pvalue <- calc_ttest(combat_yeoh[,1:416], 208, is_paired = T)
combat_d0 <- combat_yeoh[,1:208]
combat_d8 <- combat_yeoh[,209:416]
log_fc <- rowMeans(combat_d8) - rowMeans(combat_d0)
pvalue_probesets <- names(ttest_pvalue)[ttest_pvalue <= 0.05]
fc_probesets <- names(log_fc)[log_fc > 1]
intersect_probesets <- fc_probesets[fc_probesets %in% pvalue_probesets]
print(length(intersect_probesets))

corrected_df <- combat_yeoh[intersect_probesets,]

# MNN ---------------------------------------------------------------------
correctMNN <- function(X, metadata) {
  # Splitting up data into batches
  batch_info <- metadata[colnames(X), "batch_info"]
  batch_levels <- unique(batch_info)
  batch_list <- lapply(batch_levels, function(i) X[, batch_info == i])
  arr_list <- lapply(batch_list, data.matrix)
  
  colnum_list <- sapply(arr_list, ncol)
  batch_order <- order(colnum_list, decreasing = T)
  # Order batches according to number of samples
  ordered_arr <- lapply(batch_order, function(i) arr_list[[i]])
  
  # MNN correction
  mnn_yeoh_obj <- do.call(mnnCorrect, c(ordered_arr, k = 3, cos.norm.out = F))
  mnn_yeoh <- do.call(cbind, mnn_yeoh_obj$corrected)
  
  # Column names for matrix arranged in above order
  colnames_vec <- unlist(sapply(arr_list[batch_order], colnames))
  colnames(mnn_yeoh) <- colnames_vec
  # Order columns according to name and timepoint
  ordered_mnn_yeoh <- mnn_yeoh[, colnames(X)]
  
  # Replace negative values with zero
  ordered_mnn_yeoh[ordered_mnn_yeoh < 0] <- 0
  return(ordered_mnn_yeoh)
}

yeoh_mnn <- correctMNN(subset_yeoh, metadata_df)
write.table(yeoh_mnn, "dump/yeoh_mnn.tsv",
            quote = F, sep = "\t")

# Plot PCA before selecting features
plotPCA3DYeoh(yeoh_mnn, metadata_df)
rgl.postscript("dump/pca_3d-MNN_K3.pdf", "pdf")

calc.var_preservation(log_yeoh, ordered_mnn_yeoh)
metrics <- eval.batch_effects(ordered_mnn_yeoh, batch_info, class_numeric)

mnn_d0 <- ordered_mnn_yeoh[,1:208]
mnn_d8 <- ordered_mnn_yeoh[,209:416]

# Selecting drug responsive genes between D0 and D8
ttest_pvalue <- calc_ttest(cbind(mnn_d0, mnn_d8), 208, is_paired = T)
log_fc <- rowMeans(mnn_d8) - rowMeans(mnn_d0)
pvalue_probesets <- names(ttest_pvalue)[ttest_pvalue <= 0.05]
fc_probesets <- names(log_fc)[log_fc > 1]
intersect_probesets <- fc_probesets[fc_probesets %in% pvalue_probesets]
print(length(intersect_probesets))

# Filter and log transform
corrected_df <- ordered_mnn_yeoh[intersect_probesets,]

# BCM ------------------------------------------------------------------
timepoint_info <- rep(c("D0","D8","N"), c(208,208,45))
table(batch_info, timepoint_info)
order_batch <- c(10,9,8,7,6,4,3,1,2,5)
cbc_yeoh <- norm.CBC(log_yeoh, batch_info, timepoint_info, order_batch)

quantile_d0 <- norm_quantile(log_yeoh[, 1:208])
quantile_d8 <- norm_quantile(log_yeoh[, 209:416])
quantile_normal <- norm_quantile(log_yeoh[, 417:461])
quantile_yeoh <- cbind(quantile_d0, quantile_d8, quantile_normal)
cbc_yeoh <- norm.BCM(quantile_yeoh, batch_info, timepoint_info, 2)

plotPCA3DYeoh(cbc_yeoh, metadata_df)
rgl.postscript("dump/pca_3d-bcm_ref.pdf", "pdf")

calc.var_preservation(log_yeoh, cbc_yeoh)
metrics <- eval.batch_effects(cbc_yeoh, batch_info, class_numeric)

# Selecting drug responsive genes between D0 and D8
cbc_d0 <- cbc_yeoh[,1:208]
cbc_d8 <- cbc_yeoh[,209:416]
ttest_pvalue <- calc_ttest(cbc_yeoh[,1:416], 208, is_paired = T)
log_fc <- rowMeans(cbc_d8) - rowMeans(cbc_d0)
pvalue_probesets <- names(ttest_pvalue)[ttest_pvalue <= 0.05]
fc_probesets <- names(log_fc)[log_fc > 1]
intersect_probesets <- fc_probesets[fc_probesets %in% pvalue_probesets]
print(length(intersect_probesets))

corrected_df <- cbc_yeoh[intersect_probesets,]
# BCM - SUBSET -------------------------------------------------------
id_notin5 <- rownames(yeoh_batch)[yeoh_batch$batch != 5]

# Remove batch 5 from df and P107_D8
subset_yeoh <- log_yeoh[, colnames(log_yeoh) %in% id_notin5 & colnames(log_yeoh) != "P107_D8"]
subset_batch <- yeoh_batch[colnames(subset_yeoh), "batch"]
subset_class <- substring(colnames(subset_yeoh), 6, 8)
subset_class[415:417] <- "N"

table(batch_info, timepoint_info)
subset_order <- c(10,9,8,7,6,4,3,1,2)
subset_cbc_yeoh <- norm.CBC(subset_yeoh, subset_batch, subset_class, subset_order)

subset_batch_colour <- batch_palette[subset_batch]
# Shape of all timepoints
subset_class_shape <- rep(21:23, c(207,208,3))
plot.pca_3d(subset_cbc_yeoh, subset_batch_colour, subset_class_shape)
rgl.postscript("dump/pca_3d-cbc_yeoh_no5.pdf", "pdf")

# PCA
pca_obj <- prcomp(t(subset_cbc_yeoh))
# PCA: Eigenvalues
eigenvalues <- (pca_obj$sdev)^2
var_pc <- eigenvalues[1:5]/sum(eigenvalues)
pc_labels <- sprintf("PC%d (%.2f%%)", 1:5, var_pc*100)

# PCA: Coordinates
pca_coord <- pca_obj$x[,1:4]
# Response df and normal df
response_df <- pca_coord[1:414, 1:3]
normal_df <- pca_coord[415:417, 1:3]

# Calculating ERM distance
features_df1 <- calc_erm1(response_df, normal_df)
features_df2 <- calc_erm2(response_df, normal_df)
features_df3 <- calc_erm3(response_df, normal_df)

# # Using probesets values directly
# t_subset_cbc <- t(subset_cbc_yeoh)
# response_df <- t_subset_cbc[1:414,]
# normal_df <- t_subset_cbc[415:417,]

# ### TRIAL PLOTTING
# relapse_id <- rownames(yeoh_label)[yeoh_label$label == 1]
# relapse_logvec <- substring(colnames(subset_cbc_yeoh), 1, 4) %in% relapse_id
# pca_relapse <- pca_coord[relapse_logvec,]
# pca_norelapse <- pca_coord[!relapse_logvec,]
# 
# plot.pca_3d(pca_relapse, "skyblue", rep(21:22, each = 5), pc_labels)
# plot.pca_3d(pca_norelapse, "skyblue", rep(21:23, c(156,156,3)), pc_labels)

# BCM - DIFF --------------------------------------------------------------
id_notin5 <- rownames(yeoh_batch)[yeoh_batch$batch != 5]
# Remove batch 5 from df
subset_yeoh <- log_yeoh[, colnames(log_yeoh) %in% id_notin5]
subset_batch <- yeoh_batch[colnames(subset_yeoh), "batch"]
subset_class <- substring(colnames(subset_yeoh), 6, 8)
subset_class[416:418] <- "N"

table(batch_info, timepoint_info)
subset_order <- c(10,9,8,7,6,4,3,1,2)
# Perform BMC correction without batch 5
subset_cbc_yeoh <- norm.CBC(subset_yeoh, subset_batch, subset_class, subset_order)

# Batch 5 df
id_in5 <- rownames(yeoh_batch)[yeoh_batch$batch == 5]
batch5_yeoh <- log_yeoh[, colnames(log_yeoh) %in% id_in5]

kantan_df <- cbind(subset_cbc_yeoh, batch5_yeoh)
kantan_batch <- rep(1:2, c(ncol(subset_cbc_yeoh), ncol(batch5_yeoh)))
kantan_class <- substring(colnames(kantan_df), 6, 8)
kantan_class[416:418] <- "N"

cbc_kantan_df <- norm.CBC(kantan_df, kantan_batch, kantan_class, 1:2,
                          "dump/diff_kantan_vectors.tsv")
# Sort columns of corrected df
cbc_kantan_df <- cbc_kantan_df[,colnames(log_yeoh)]

kantan_real_batch <- yeoh_batch[colnames(cbc_kantan_df), "batch"]
kantan_batch_colour <- batch_palette[kantan_real_batch]
kantan_class_factor <- as.factor(kantan_class)
levels(kantan_class_factor) <- c(21,23,22,24)
kantan_class_shape <- as.numeric(levels(kantan_class_factor))[kantan_class_factor]

# Shape of all timepoints
plot.pca_3d(cbc_kantan_df, kantan_batch_colour, kantan_class_shape)
rgl.postscript("dump/pca_3d-cbc_yeoh.pdf", "pdf")

# Selecting drug responsive genes between D0 and D8
cbc_d0 <- cbc_kantan_df[,1:208]
cbc_d8 <- cbc_kantan_df[,209:416]
ttest_pvalue <- calc_ttest(cbc_kantan_df[,1:416], 208, is_paired = T)
log_fc <- rowMeans(cbc_d8) - rowMeans(cbc_d0)
pvalue_probesets <- names(ttest_pvalue)[ttest_pvalue <= 0.05]
fc_probesets <- names(log_fc)[log_fc > 1]
intersect_probesets <- fc_probesets[fc_probesets %in% pvalue_probesets]
print(length(intersect_probesets))

corrected_df <- cbc_kantan_df[intersect_probesets,]
# ComBat ------------------------------------------------------------------
# Obtaining batch information of selected_yeoh df
# Rows of metadata are to be in same order as columns of edata
batch_info <- yeoh_batch[colnames(log_yeoh), "batch"]
timepoint_info <- as.factor(rep(c("D0","D8","N"), c(208,208,45)))
yeoh_metadata <- data.frame(batch = as.factor(batch_info),
                            timepoint = timepoint_info)

# Place adjustment/confounding variables in model.matrix (e.g. age)
# Do not put batch variables in model.matrix
# ### Put batch variables directly in combat function
model_combat <- model.matrix(~1, data = yeoh_metadata)

# ### Include biological variable of interest as covariate
# model_combat <- model.matrix(~timepoint, data = yeoh_metadata)

# ComBat batch correction
combat_yeoh <- ComBat(data.matrix(log_yeoh), batch_info, model_combat)

# Replacing negative values with 0
combat_yeoh[combat_yeoh < 0] <- 0

# Plot PCA before selecting features
plotPCA3DYeoh(combat_yeoh, metadata_df)
rgl.postscript("dump/pca_3d-combat_cov.pdf", "pdf")

# Quantitative evaluation
calc.var_preservation(log_yeoh, combat_yeoh)
metrics <- eval.batch_effects(combat_yeoh, batch_info, class_numeric)

# Selecting drug responsive genes between D0 and D8
combat_d0 <- combat_yeoh[,1:208]
combat_d8 <- combat_yeoh[,209:416]
ttest_pvalue <- calc_ttest(cbind(combat_d0, combat_d8), 208, is_paired = T)
log_fc <- rowMeans(combat_d8) - rowMeans(combat_d0)
pvalue_probesets <- names(ttest_pvalue)[ttest_pvalue <= 0.05]
fc_probesets <- names(log_fc)[log_fc > 1]
intersect_probesets <- fc_probesets[fc_probesets %in% pvalue_probesets]
print(length(intersect_probesets))

corrected_yeoh <- combat_yeoh[intersect_probesets,]

# ComBat (Class-specific) ------------------------------------------------------------------
# Special filtering for class-specific ComBat
# Filters out probesets that have zero variance within subset df in each class
selected_probesets1 <- filter_probesets(scaled_yeoh, 0.3, class_info,
                                        logical_func = all)
log_yeoh1 <- log_yeoh[selected_probesets1,]
# Quantile by timepoint
quantile_d0 <- norm_quantile(log_yeoh1[, 1:208])
quantile_d8 <- norm_quantile(log_yeoh1[, 209:416])
quantile_normal <- norm_quantile(log_yeoh1[, 417:461])
quantile_yeoh <- cbind(quantile_d0, quantile_d8, quantile_normal)

list_yeoh_df <- split.default(quantile_yeoh, class_info)
list_yeoh_arr <- lapply(list_yeoh_df, data.matrix)
print(yeoh_batch)
list_batch_info <- lapply(list_yeoh_df, function(df1) yeoh_batch[colnames(df1), "batch"])

# Creation of list of model.matrix
list_class_info <- split(class_info, class_info)
list_metadata <- lapply(list_class_info, data.frame)
# All samples in batch are the same class, hence no need for covariate
list_model_matrix <- lapply(list_metadata, model.matrix, object = ~1)

list_combat_arr <- mapply(ComBat, list_yeoh_arr, list_batch_info, list_model_matrix)
combat_arr <- do.call(cbind, list_combat_arr)
# Reorder columns according to initial df
cs_combat_yeoh <- data.frame(combat_arr[,colnames(quantile_yeoh)])

# Replacing negative values with 0
cs_combat_yeoh[cs_combat_yeoh < 0] <- 0

# Plot PCA before selecting features
plotPCA3DYeoh(cs_combat_yeoh, metadata_df)
rgl.postscript("dump/pca_3d-quantile_cs_combat.pdf", "pdf")

# Quantitative evaluation
calc.var_preservation(log_yeoh, combat_yeoh)
metrics <- eval.batch_effects(combat_yeoh, batch_info, class_numeric)

# Selecting drug responsive genes between D0 and D8
combat_d0 <- cs_combat_yeoh[,1:208]
combat_d8 <- cs_combat_yeoh[,209:416]
ttest_pvalue <- calc_ttest(cbind(combat_d0, combat_d8), 208, is_paired = T)
log_fc <- rowMeans(combat_d8) - rowMeans(combat_d0)
pvalue_probesets <- names(ttest_pvalue)[ttest_pvalue <= 0.05]
fc_probesets <- names(log_fc)[log_fc > 1]
intersect_probesets <- fc_probesets[fc_probesets %in% pvalue_probesets]
print(length(intersect_probesets))

corrected_df <- cs_combat_yeoh[intersect_probesets,]

# Harman ------------------------------------------------------------------
harman_obj <- harman(log_yeoh, class_info, batch_info, limit = 0.95)
harman_yeoh <- data.frame(reconstructData(harman_obj))

# Batch information of all the timepoints
plotPCA3DYeoh(harman_yeoh, metadata_df)
rgl.postscript("dump/pca_3d-harman_yeoh.pdf", "pdf")

calc.var_preservation(log_yeoh, harman_yeoh)
metrics <- eval.batch_effects(harman_yeoh, batch_info, class_numeric)

# Selecting drug responsive genes between D0 and D8
harman_d0 <- harman_yeoh[,1:208]
harman_d8 <- harman_yeoh[,209:416]
ttest_pvalue <- calc_ttest(harman_yeoh[,1:416], 208, is_paired = T)
log_fc <- rowMeans(harman_d8) - rowMeans(harman_d0)
pvalue_probesets <- names(ttest_pvalue)[ttest_pvalue <= 0.05]
fc_probesets <- names(log_fc)[log_fc > 1]
intersect_probesets <- fc_probesets[fc_probesets %in% pvalue_probesets]
print(length(intersect_probesets))

corrected_df <- harman_yeoh[intersect_probesets,]

# Scanorama ---------------------------------------------------------------
# Order samples according to batches
batch_list <- lapply(1:10, function(i) log_yeoh[, batch_info == i])
batch_ordered_yeoh <- do.call(cbind, batch_list)

# # Write log_maqc for numpy array
# write.table(batch_ordered_yeoh, "data/scanorama/batch_ordered_yeoh.tsv",
#             quote = F, sep = "\t", row.names = F, col.names = F)
# write(rownames(log_yeoh), file = "data/scanorama/yeoh-probesets.txt")

scanorama_yeoh <- read.table("data/scanorama/yeoh-scanorama_k10.tsv",
                             sep = "\t", row.names = 1)
colnames(scanorama_yeoh) <- colnames(batch_ordered_yeoh)
# Reorder columns according to log_yeoh
ordered_scanorama <- scanorama_yeoh[,colnames(log_yeoh)]

# Plot PCA before selecting features
# Batch information of all the timepoints
plotPCA3DYeoh(ordered_scanorama*100, metadata_df)
rgl.postscript("dump/pca_3d-scanorama_K20.pdf", "pdf")

calc.var_preservation(log_yeoh, ordered_scanorama)
metrics <- eval.batch_effects(ordered_scanorama, batch_info, class_numeric)
colnames(ordered_scanorama)

scanorama_d0 <- ordered_scanorama[,1:208]
scanorama_d8 <- ordered_scanorama[,209:416]

# Selecting drug responsive genes between D0 and D8
ttest_pvalue <- calc_ttest(cbind(scanorama_d0, scanorama_d8), 208, is_paired = T)
# bh_qvalue <- p.adjust(ttest_pvalue, method = "BH")
log_fc <- calc_logfc(scanorama_d0, scanorama_d8)
fc_probesets <- names(log_fc)[log_fc > 1]
pvalue_probesets <- names(ttest_pvalue)[ttest_pvalue <= 0.05]
length(pvalue_probesets)
intersect_probesets <- fc_probesets[fc_probesets %in% pvalue_probesets]
print(length(intersect_probesets))
# Filter and log transform
corrected_df <- ordered_scanorama[intersect_probesets,]

# PLOT: ROC -----------------------------------------------------------------
# Plot ROC
# Visualise current results now
head(results_df)
labels_vec <- results_df[, 10]
labels_vec
par(mar = rep(5,4))
line_labels <- c("ERM1", "ERM1-Ratio",
                 "ERM2", "ERM2-Ratio",
                 "PC1", "PC1-Ratio")
# line_labels <- c("Quantile", "CS-Quantile",
#                  "ComBat", "Harman",
#                  "MNN [k=3]", "BCM")

plotROC(results_df[,1:6], labels_vec, name_vec = line_labels)
results_roc <- recordPlot()

save_fig(results_roc, "dump/roc-telaml1_bcm_qpsp.pdf",
         width = 9, height = 9)

# # Process results from all batch correction methods
# results_rpath_vec <- list.files("dump/fig/erm", full.names = T)
# list_results <- lapply(results_rpath_vec, read.table, sep = "\t")
# names(list_results) <- substring(results_rpath_vec, 22)
# # lapply(list_results, rownames)
# 
# all_erm <- lapply(list_results, function(df) df[,1])
# # Reorder list of erm vectors
# reordered_erm <- all_erm[c(6,3,2,4,5,1)]
# str(reordered_erm)
# labels_vec <- list_results[[1]][,7]
# names(list_results)
# # Plot ROC
# # Visualise current results now
# head(results_df)
# labels_vec <- results_df[, 7]
# labels_vec
# par(mar = rep(5,4))
# line_labels <- c("Quantile", "CS-Quantile",
#                  "ComBat", "Harman",
#                  "MNN [k=3]", "BCM")
# plot_roc(results_df[,1:6], labels_vec,
#          name_vec = line_labels)
# results_roc <- recordPlot()

# Subtype analysis
head(results_df)
head(yeoh_label)
analysis_df <- cbind(results_df[,c(1,10)],
                     yeoh_label[rownames(results_df),"subtype", drop = F])
head(analysis_df)
plot_subtype <- ggplot(analysis_df) +
  geom_jitter(aes(x = subtype, y = erm1, color = as.factor(label)),
              position = position_jitter(0.1),
              show.legend = F) +
  scale_color_manual(values = c("darkolivegreen3", "tomato3"))

plot_subtype

ggsave("dump/subtype_analysis-qpsp_nea.pdf", plot_subtype,
       width = 8, height = 5)

# Load results
# results_df <- read.table("dump/results-quantile_timepoint_d33_projected_d8.tsv",
#                          sep = "\t")
# labels_vec <- results_df[, 7]
# head(results_df[,1:6])

# # Visualise all except subtype: Others
# head(filtered_subtype_results)
# labels_vec <- results_df[, 7]
# labels_vec
# par(mar = rep(5,4))
# line_labels <- c("ERM1", "ERM1-Ratio",
#                  "ERM2", "ERM2-Ratio",
#                  "PC1", "PC1-Ratio", "MRD33")
# plot_roc(filtered_subtype_results[,c(1:6,9)], filtered_subtype_results[,7],
#          name_vec = line_labels, is_bigger_better_vec = c(rep(F,6),T))
# 
# results_roc <- recordPlot()
# save_fig(results_roc, "dump/roc-quantile_timepoint_no_others.pdf",
#          width = 9, height = 9)

# VISUALISATION SUBTYPES ---------------------------------------------------
# Plot PCA before selecting features
# Batch information of all the timepoints

plotPCA3DYeoh1 <- function(df1, metadata_df) {
  # 3D PCA plot
  plotPCA3D1 <- function(df, colour, pch, pc_labels = NULL,
                         ratio_list = list(2,1,1)) {
    if (is.null(pc_labels)) {
      print("PCA performed!")
      pca_obj <- prcomp(t(df), center = T, scale. = F)
      pca_df <- as.data.frame(pca_obj$x[,1:3])
      eigenvalues <- (pca_obj$sdev)^2
      var_pc <- eigenvalues[1:3]/sum(eigenvalues)
      print(var_pc)
      pc_labels <- sprintf("PC%d (%.2f%%)", 1:3, var_pc*100)
    } else {
      print("No PCA performed!")
      pca_df <- as.data.frame(df)
    }
    
    # RGL plot parameters
    rgl.open()
    rgl.bg(color="white")
    rgl.viewpoint(zoom = 0.8)
    # rgl.viewpoint(theta = 110, phi = 5, zoom = 0.8)
    par3d(windowRect = c(50, 20, 500, 500))
    with(pca_df, pch3d(PC1, PC2, PC3, col = colour,
                       pch = pch, cex = 0.5, lwd = 2))
    box3d(col = "black")
    title3d(xlab = pc_labels[1], ylab = pc_labels[2],
            zlab = pc_labels[3], col = "black")
    # Plot aspect ratios of axis according to variance
    do.call(aspect3d, ratio_list)
  }
  # COLOUR
  class_info <- metadata_df[colnames(df1), "class_info"]
  generate_colour <- colorRampPalette(c("lightblue", "darkblue"))
  palette <- generate_colour(3)
  colour <- palette[class_info]
  # SHAPE
  batch_info <- metadata_df[colnames(df1), "batch_info"]
  levels(batch_info) <- 0:(nlevels(batch_info)-1)
  pch <- as.numeric(as.character(batch_info))
  # PLOT
  plotPCA3D1(df1, colour, pch)
}

plotPCA3DYeoh1(subset_yeoh, metadata_df)
plotPCA3DYeoh(subset_yeoh, metadata_df)

# Plot PCA-3D
plot.pca_3d(pca_coord, batch_colour, timepoint_shape, pc_labels)
rgl.postscript("dump/pca_3d_selected-bcm_ref.pdf", "pdf")

# Plot PCA-2D
# plot_pca <- pca_all(as.data.frame(plot_arr), batch_colour,
#                     batch_shape, pc_labels)
# plot_pca
# ggsave("dump/pca-quantile_d0leuk_d8_normal.pdf", plot_pca,
#        width = 9, height = 9)

# Visualising vectors
# Subtype sizes
num_patient <- nrow(response_df)/2
row_index <- substring(rownames(response_df)[1:num_patient],1,4)
label_yeoh <- yeoh_label[row_index, "label"]
subtype_yeoh <- yeoh_label[row_index, "subtype"]
generate_greens <- colorRampPalette(c("greenyellow", "darkgreen"))
subtype_palette <- generate_greens(8)
subtype_colour <- subtype_palette[subtype_yeoh]

# # Legend [levels(subtype_yeoh)]
# plot(rep(0, 8), 1:8, pch = 19, cex = 2,
#      xlim = c(-0.3,0.7),
#      col = subtype_palette, axes = F, ann = F)
# text(rep(0.1, 8), 1:8, levels(subtype_yeoh), pos = 4)
# plot_legend <- recordPlot()
# save_fig(plot_legend, "dump/legend-subtype.pdf",
#          width = 4, height = 6)

# For vectors plot
# patient_subtype <- as.numeric(subtype_yeoh)
# subtype_size <- c(rep(as.numeric(subtype_yeoh)/10, 2), rep(0.5, 45))

# Coordinates for arrows
arrows_arr <- cbind(response_df[1:num_patient,],
                    response_df[-(1:num_patient),],
                    matrix(label_yeoh))
colnames(arrows_arr) <- c(paste(colnames(arrows_arr)[1:6],
                                rep(LETTERS[1:2], each = 3),
                                sep = "_"),
                          "labels")
# Coordinates for centroids
centroid_arr <- rbind(apply(response_df[1:num_patient,], 2, median),
                      apply(response_df[-(1:num_patient),], 2, median),
                      # apply(plot_arr[1:750, 1:3], 2, median),
                      apply(normal_df[1:42,], 2, median),
                      apply(normal_df[43:45,], 2, median))

patient_colour <- batch_colour[1:num_patient]
vectors_plot <- plot_vectors(as.data.frame(arrows_arr),
                             as.data.frame(centroid_arr),
                             pc_labels,
                             patient_colour,
                             subtype_colour)
vectors_plot
ggsave("dump/vectors-qpsp_ref.pdf", vectors_plot,
       width = 12, height = 6)

# # Subtype colours
# num_subtype <- table(substring(rownames(plot_arr)[1:750], 1, 1))
# red_palette <- brewer.pal(9, "Reds")[2:9]
# subtype_colour <- rep(red_palette, num_subtype)

# # Batch information: Associating colour
# num_subtype <- table(substring(colnames(mile_data), 1, 1))
# red_palette <- brewer.pal(9, "Reds")[2:9]
# mile_palette <- c(red_palette, "darkolivegreen3")
# colour_code <- rep(mile_palette, num_subtype)

# # Batch information of yeoh is encoded
# batch_info <- yeoh_metadata[rownames(response_df),6]
# blue_palette <- brewer.pal(9, "Blues")
# batch_colour <- c(blue_palette[batch_info],
#                   "tomato3", "darkolivegreen3")
# # Colour information
# all_colour <- c(rep(c("steelblue4", "turquoise3"), c(210, 210)),
#                 rep("tomato3", 750), rep("darkolivegreen3", 74))

# SUBTYPE -------------------------------------------------------------------
analysis_df <- cbind(results_df,
                     yeoh_label[rownames(results_df),"subtype", drop = F])
plot_2 <- ggplot(analysis_df) +
  geom_point(aes(x = d0d8_l2norm, y = erm1, col = factor(label)), show.legend = F) +
  scale_color_manual(values = c("darkolivegreen3", "tomato3")) +
  facet_wrap(~subtype, nrow = 2, ncol = 4)

ggsave("dump/quantile_pca_subtype.pdf", width = 12, height = 6)

# Quantile %>% Feature selection
quantile_feat_selection <- t(corrected_df)
response_df <- quantile_feat_selection[1:416,]
normal_df <- quantile_feat_selection[417:461,]

# Calculating ERM distance
features_df1 <- calc_erm1(response_df, normal_df)
features_df2 <- calc_erm2(response_df, normal_df)
features_df3 <- calc_erm3(response_df, normal_df)
features_df4 <- calc_erm4(response_df)

results_df <- collate_results(features_df1,
                              features_df2,
                              features_df3,
                              features_df4,
                              yeoh_label)
head(results_df)


multiplot <- plot.angle(results_df)
ggsave("dump/angle_selected_2.pdf", multiplot, width = 8, height = 4)

# Quantile
t_quantile <- t(quantile_yeoh)
response_df <- t_quantile[1:416,]
normal_df <- t_quantile[417:461,]

# Calculating ERM distance
features_df1 <- calc_erm1(response_df, normal_df)
features_df2 <- calc_erm2(response_df, normal_df)
features_df3 <- calc_erm3(response_df, normal_df)
features_df4 <- calc_erm4(response_df)

results_df <- collate_results(features_df1,
                              features_df2,
                              features_df3,
                              features_df4,
                              yeoh_label)
head(results_df)

multiplot <- plot.angle(results_df)
multiplot
ggsave("dump/angle_quantile_2.pdf", multiplot, width = 8, height = 4)

# Plot density curves separated by class
colour_info <- rep(c("darkolivegreen3", "tomato3","gold"), c(208,208,45))
class_colour <- rep(colour_info, each = nrow(quantile_yeoh))
quantile_density <- ggplot(cbind(stack(quantile_yeoh), class_colour)) +
  geom_density(aes(x = values, col = ind), show.legend = F) + 
  facet_wrap(~class_colour)

# Plot 3D scatter
plot.pch_3d <- function(df, colour_code = "lightblue", shape_vec = 21) {
  rgl.open()
  rgl.bg(color="white")
  rgl.viewpoint(zoom = 1)
  # rgl.viewpoint(theta = 110, phi = 5, zoom = 0.8)
  par3d(windowRect = c(50, 20, 500, 500))
  # Plot of MILE dataset
  pch3d(df[,1], df[,2], df[,3], bg = colour_code,
        pch = shape_vec, cex = 0.3, lwd = 1.5)
  # axes3d(c('x-', 'y+', 'z+'),
  #        col = "gray8", labels = F, tick = F)
  box3d(col = "black")
  title3d(xlab = colnames(df)[1], ylab = colnames(df)[2],
          zlab = colnames(df)[3], col = "black")
}

head(results_df[,c(1,8,9)])
colour_label <- ifelse(results_df[,10] == 1, "tomato3", "darkolivegreen3")
plot.pch_3d(results_df[,c(1,8,9)], colour_label)


# GLOBAL BCM (D0) ----------------------------------------------
# Subtype of subsetted patients
table(yeoh_label[unique(substring(colnames(subset_yeoh), 1, 4)), "subtype"])
xtable(table(yeoh_batch[colnames(subset_yeoh), "batch"],
      yeoh_label[substring(colnames(subset_yeoh), 1, 4), "subtype"]))
subset_subtype <- yeoh_label[substring(colnames(subset_yeoh), 1, 4), "subtype"]
# Adds NA as a level
subset_subtype <- addNA(subset_subtype)
levels(subset_subtype)[9] <- "Normal"
# Split df by subtypes into list
list_subtype_df <- split.default(subset_yeoh, subset_subtype)

# Check batch-class configuration of subsetted data
# Combine both TEL-AML1 and normal patients for batch correction
telaml1_normal <- cbind(list_subtype_df$`TEL-AML1`, list_subtype_df$Normal)
telaml1_metadata <- metadata_df[colnames(telaml1_normal),]
library(xtable)
xtable(table(telaml1_metadata$batch_info, telaml1_metadata$class_info))
bcm_telaml1 <- normGlobalBCM(telaml1_normal, yeoh_batch)

### TEL-AML1 PLOT ###
# Plot PCA before selecting features
# Batch information of all the timepoints
generate_colour <- colorRampPalette(c("lightblue", "darkblue"))
batch_palette <- generate_colour(10)
subset_colour <- batch_palette[telaml1_metadata$batch_info]
# Shape of all timepoints
subset_shape <- rep(21:23, table(telaml1_metadata$class_info))
plot.pca_3d(telaml1_normal, subset_colour, subset_shape)
plot.pca_3d(bcm_telaml1, subset_colour, subset_shape)
rgl.postscript("dump/pca_3d-tel_aml1_bcm.pdf", "pdf")

# Selecting drug responsive genes between D0 and D8
bcm_d0 <- bcm_telaml1[,1:39]
bcm_d8 <- bcm_telaml1[,40:78]
print(colnames(bcm_d0)); print(colnames(bcm_d8))
ttest_pvalue <- calc_ttest(cbind(bcm_d0, bcm_d8), 39, is_paired = T)
log_fc <- rowMeans(bcm_d8) - rowMeans(bcm_d0)
pvalue_probesets <- names(ttest_pvalue)[ttest_pvalue <= 0.05]
fc_probesets <- names(log_fc)[log_fc > 1]
intersect_probesets <- fc_probesets[fc_probesets %in% pvalue_probesets]
print(length(intersect_probesets))
corrected_df <- bcm_telaml1[intersect_probesets,]
# corrected_df <- selected_quantile_yeoh

# PCA
pca_obj <- prcomp(t(corrected_df))
# PCA: Eigenvalues
eigenvalues <- (pca_obj$sdev)^2
variance_proportion <- eigenvalues/sum(eigenvalues)
# Use till PC15 to account for at least 70% variance
cumsum(variance_proportion)[1:50]
# Identify PC that has just above 70% variance
pc_ind <- which.max(cumsum(variance_proportion)[1:50] > 0.70)
pc_ind # PC20
# pc_labels <- sprintf("PC%d (%.2f%%)", 1:3, (var_pc*100)[1:3])

# PCA: Coordinates
pca_coord <- pca_obj$x[,1:pc_ind]
# Response df and normal df
response_df <- pca_coord[1:78,]
normal_df <- pca_coord[79:81,]
print(rownames(response_df))
print(rownames(normal_df))

# Calculating ERM distance
features_df1 <- calc_erm1(response_df, normal_df)
features_df2 <- calc_erm2(response_df, normal_df)
features_df3 <- calc_erm3(response_df, normal_df)
features_df4 <- calc_erm4(response_df)

### USING ALL SUBTYPES ###

# Plot PCA before selecting features
# Batch information of all the timepoints
generate_colour <- colorRampPalette(c("lightblue", "darkblue"))
batch_palette <- generate_colour(10)
subset_colour <- batch_palette[subset_batch]
# Shape of all timepoints
subset_shape <- rep(21:23, table(subset_class))
plot.pca_3d(bcm_df, subset_colour, subset_shape)
rgl.postscript("dump/pca_3d-bcm_d0.pdf", "pdf")

# Selecting drug responsive genes between D0 and D8
bcm_d0 <- bcm_df[,1:204]
bcm_d8 <- bcm_df[,205:408]
print(colnames(bcm_d0)); print(colnames(bcm_d8))
ttest_pvalue <- calc_ttest(cbind(bcm_d0, bcm_d8), 204, is_paired = T)
log_fc <- rowMeans(bcm_d8) - rowMeans(bcm_d0)
pvalue_probesets <- names(ttest_pvalue)[ttest_pvalue <= 0.05]
fc_probesets <- names(log_fc)[log_fc > 1]
intersect_probesets <- fc_probesets[fc_probesets %in% pvalue_probesets]
print(length(intersect_probesets))
corrected_df <- bcm_df[intersect_probesets,]

### All features
colnames(without5_df)
t_without5 <- t(without5_df)
response_df <- t_without5[1:414,]
normal_df <- t_without5[415:417,]
features_df1 <- calc_erm1(response_df, normal_df)
features_df2 <- calc_erm2(response_df, normal_df)
features_df3 <- calc_erm3(response_df, normal_df)
features_df4 <- calc_erm4(response_df)

### Quantile %>% Feature selection
# Selecting drug responsive genes between D0 and D8
bcm_d0 <- without5_df[,1:207]
bcm_d8 <- without5_df[,208:414]
ttest_pvalue <- calc_ttest(cbind(bcm_d0, bcm_d8), 207, is_paired = T)
log_fc <- rowMeans(bcm_d8) - rowMeans(bcm_d0)
pvalue_probesets <- names(ttest_pvalue)[ttest_pvalue <= 0.05]
fc_probesets <- names(log_fc)[log_fc > 1]
intersect_probesets <- fc_probesets[fc_probesets %in% pvalue_probesets]
print(length(intersect_probesets))
corrected_df <- without5_df[intersect_probesets,]

t_selected_df <- t(corrected_df)
response_df <- t_selected_df[1:414,]
normal_df <- t_selected_df[415:417,]
# Calculating ERM distance
features_df1 <- calc_erm1(response_df, normal_df)
features_df2 <- calc_erm2(response_df, normal_df)
features_df3 <- calc_erm3(response_df, normal_df)
features_df4 <- calc_erm4(response_df)

results_df <- collate_results(features_df1,
                              features_df2,
                              features_df3,
                              features_df4,
                              yeoh_label)

# PLOT ANGLE --------------------------------------------------------------
plot_angle <- function(df1, flag) {
  if (flag == 1) {
    plot_1 <- ggplot(df1) +
      geom_point(aes(x = angle_d0d8, y = d0d8_l2norm, col = factor(label)),
                 show.legend = F) +
      scale_color_manual(values = c("darkolivegreen3", "tomato3"))
    plot_2 <- ggplot(df1) +
      geom_point(aes(x = angle_d0d8, y = erm1, col = factor(label)),
                 show.legend = F) +
      scale_color_manual(values = c("darkolivegreen3", "tomato3"))
    multiplot1 <- plot_grid(plot_1, plot_2, ncol = 2)
    return(multiplot1)
  } else {
    plot_3 <- ggplot(df1) +
      geom_point(aes(x = angle_d0d8, y = diff_l2norm, col = factor(label)),
                 show.legend = F) +
      scale_color_manual(values = c("darkolivegreen3", "tomato3"))
    plot_4 <- ggplot(df1) +
      geom_point(aes(x = d0d8_l2norm, y = erm1, col = factor(label)),
                 show.legend = F) +
      scale_color_manual(values = c("darkolivegreen3", "tomato3"))
    multiplot2 <- plot_grid(plot_3, plot_4, ncol = 2)
  }
}

scatter_angle <- plot_angle(results_df1)
scatter_angle
ggsave("dump/angle-telaml1_bcm_qpsp_2.pdf",
       scatter_angle, width = 8, height = 4)

plotMRD <- function(df1) {
  # Log10 transform MRD33 values
  df1$d33_mrd <- -log10(df1$d33_mrd)
  df1$d33_mrd[is.infinite(df1$d33_mrd)] <- -1
  print(head(df1$d33_mrd))
  
  plot_mrd1 <- ggplot(df1) +
    geom_point(aes(x = d33_mrd, y = erm1, col = factor(label)),
               show.legend = F) +
    scale_color_manual(values = c("darkolivegreen3", "tomato3")) +
    labs(x = bquote(-log[10](mrd_d33)))
  plot_mrd2 <- ggplot(df1) +
    geom_point(aes(x = d33_mrd, y = angle_d0d8, col = factor(label)),
               show.legend = F) +
    scale_color_manual(values = c("darkolivegreen3", "tomato3")) +
    labs(x = bquote(-log[10](mrd_d33)))
  plot_mrd3 <- ggplot(df1) +
    geom_point(aes(x = d33_mrd, y = d0d8_l2norm, col = factor(label)),
               show.legend = F) +
    scale_color_manual(values = c("darkolivegreen3", "tomato3")) +
    labs(x = bquote(-log[10](mrd_d33)))
  
  plot_mrd <- plot_grid(plot_mrd1, plot_mrd2, plot_mrd3, ncol = 3)
  
  return(plot_mrd)
}

mrd_plot <- plotMRD(results_df1)
mrd_plot

ggsave(sprintf("dump/mrd-%s_qpsp_bcmD0.pdf", SUBTYPE),
       mrd_plot, width = 12, height = 4)

# EXPLORATION --------------------------------------------------
# Split df by batches into list
list_batch_df <- split.default(log_yeoh, batch_info)

generate_colour <- colorRampPalette(c("lightblue", "darkblue"))
batch_palette <- generate_colour(10)
colour_code <- c("darkolivegreen3", "tomato3", "darkorchid4")

### Investigate whether batch imbalance changes the biological hyperplane
pair_batch <- cbind(list_batch_df[[2]][,1:45], list_batch_df[[3]][,15:42])

pair_batch <- cbind(list_batch_df[[2]], list_batch_df[[3]])
# For specific dataframe
class_batch <- metadata_df[colnames(pair_batch), "class"]
class_colour <- colour_code[as.factor(class_batch)]

batch_info2 <- metadata_df[colnames(pair_batch), "batch"]
batch_factor2 <- as.factor(batch_info2)
levels(batch_factor2) <- 21:22
batch_pch <- as.numeric(as.character(batch_factor2))

plot_title <- sprintf("B2 vs. B%d", 3)
plot_pca_2d(pair_batch, class_colour, shape_vec = batch_pch)
# Batch imbalance does not affect parallel alignment between batches
# More variation due to batch (Batches have different compositions?)

table(batch_info2, class_batch)

### All possible pairwise combinations
list_plot <- list()
j <- 1

for (i in c(1,3:10)) {
  pair_batch <- cbind(list_batch_df[[2]], list_batch_df[[i]])
  # For specific dataframe
  class_batch <- metadata_df[colnames(pair_batch), "class"]
  class_colour <- colour_code[as.factor(class_batch)]
  
  batch_info2 <- metadata_df[colnames(pair_batch), "batch"]
  batch_factor2 <- droplevels(as.factor(batch_info2))
  levels(batch_factor2) <- 21:22
  batch_pch <- as.numeric(as.character(batch_factor2))
  
  plot_title <- sprintf("B2 vs. B%d", 3)
  list_plot[[j]] <- plot_pca_2d(pair_batch, class_colour, shape_vec = batch_pch)
  j <- j + 1
}

pairwise_pca <- do.call(plot_grid, c(list_plot, nrow = 3, ncol = 3))
pairwise_pca
ggsave("dump/pairwise_pca1.pdf", pairwise_pca, width = 11, height = 11)



# Investigate by batch
pairwise_dist <- dist(t(pair_batch))
pairwise_mat <- as.matrix(pairwise_dist)

hcluster <- hclust(pairwise_dist)
dendo_obj <- as.dendrogram(hcluster)
# sample_id <- labels(dendo_obj)
# nodePar <- list(lab.cex = 0.3, pch = c(NA, NA),
#                 cex = 0.5, col = "blue")
# Settings of dendogram
dendo_obj <- set(dendo_obj, "labels_cex", 0.4)
plot(dendo_obj, horiz = F)
colored_bars(batch_pch, dendo_obj,
             rowLabels = "Batch ", y_shift = -150, y_scale = 30,
             sort_by_labels_order = T)
colored_bars(class_colour, dendo_obj,
             rowLabels = "Class ", y_shift = -100, y_scale = 30,
             sort_by_labels_order = T)
dendogram <- recordPlot()
# save_fig(dendogram, "dump/hclust-mnn_yeoh.pdf",
#          width = 12, height = 6)

### Measure cosine distance by cosine normalising data
cosine_yeoh <- norm_cosine(log_yeoh)
pairwise_dist1 <- dist(t(cosine_yeoh))
pairwise_mat1 <- as.matrix(pairwise_dist)

hcluster <- hclust(pairwise_dist1)
dendo_obj1 <- as.dendrogram(hcluster)
# sample_id <- labels(dendo_obj)
# nodePar <- list(lab.cex = 0.3, pch = c(NA, NA),
#                 cex = 0.5, col = "blue")

generate_colour <- colorRampPalette(c("lightblue", "darkblue"))
batch_palette <- generate_colour(10)
batch_colour <- batch_palette[batch_info]
colour_code <- c("darkolivegreen3", "tomato3", "darkorchid4")
class_colour <- colour_code[as.factor(class_info)]

# Settings of dendogram
dendo_obj1 <- set(dendo_obj1, "labels_cex", 0.4)
plot(dendo_obj1, horiz = F)
colored_bars(batch_colour, dendo_obj1,
             rowLabels = "Batch ", y_shift = -0.1, y_scale = 0.05,
             sort_by_labels_order = T)
colored_bars(class_colour, dendo_obj1,
             rowLabels = "Class ", y_shift = -0.15, y_scale = 0.05,
             sort_by_labels_order = T)

dendogram <- recordPlot()
save_fig(dendogram, "dump/hclust-mnn_yeoh.pdf",
         width = 12, height = 6)

# LOADINGS ----------------------------------------------------------------
### PCA: Entire dataset ###
generate_colour <- colorRampPalette(c("lightblue", "darkblue"))
batch_palette <- generate_colour(10)
batch_colour <- batch_palette[metadata_df$batch_info]
# Shape of all timepoints
timepoint_shape <- rep(21:23, c(208,208,45))

entire_prcomp <- prcomp(t(log_yeoh))
plot.pca_3d(entire_prcomp$x, batch_colour, timepoint_shape, 1:3)
entire_pca <- entire_prcomp$x

# All eigenvectors are of unit norm even when original variance..
# is not standardised
apply(entire_prcomp$rotation, 2, calcL2Norm)

# Histogram of PC1 and PC2 loading coefficients
sort(abs(entire_prcomp$rotation[,2]), decreasing = F)[10000:10100]
hist(abs(entire_prcomp$rotation[,1]), breaks = 20, probability = F)
hist(abs(entire_prcomp$rotation[,2]), breaks = 20, probability = F, add = T)

### PAIRWISE PCA ###
# Split df by batches into list
list_batch_df <- split.default(log_yeoh, batch_info)

### All possible pairwise combinations
list_plot <- list()
list_bio_vec <- list()
list_batch_vec <- list()
j <- 1
for (i in c(1,3:10)) {
  ref_df <- list_batch_df[[2]]
  other_df <- list_batch_df[[i]]
  pair_batch <- cbind(ref_df, other_df)
  pair_prcomp <- prcomp(t(pair_batch), scale. = F)
  # Subset samples from ref_df
  ref_pca <- pair_prcomp$x[1:ncol(ref_df),1:2]
  other_pca <- pair_prcomp$x[-(1:ncol(ref_df)),1:2]
  
  # Linear regression to get gradient
  ref_coef <- coef(lm(ref_pca[,2] ~ ref_pca[,1]))
  other_coef <- coef(lm(other_pca[,2] ~ other_pca[,1]))
  mean_gradient <- mean(c(ref_coef[2], other_coef[2]))
  
  # Biological vector in PCA space (unit norm)
  bio_vec_pca <- c(1, mean_gradient) / calcL2Norm(c(1, mean_gradient))
  rotation_90 <- calcRotationMatrix(pi/2)
  # Batch effect vector in PCA space (unit norm)
  batch_vec_pca <- rotation_90 %*% bio_vec_pca
  # Biological vector in gene expression space (unit norm)
  bio_mat <- pair_prcomp$rotation[,1:2] %*% bio_vec_pca
  list_bio_vec[[j]] <- setNames(as.vector(bio_mat), rownames(bio_mat)) 
  # Batch effect vector in PCA space (unit norm)
  batch_mat <- pair_prcomp$rotation[,1:2] %*% batch_vec_pca
  list_batch_vec[[j]] <- setNames(as.vector(batch_mat), rownames(batch_mat)) 
  print(calcAngleVectors(list_bio_vec[[j]],
                         list_batch_vec[[j]]))
  
  ## Plotting parameters
  CLASS_COLOUR_PALETTE <- c("darkolivegreen3", "tomato3", "darkorchid4")
  plot_title <- sprintf("B2 vs. B%d", i)
  eigenvalues <- (pair_prcomp$sdev)^2
  var_pct <- eigenvalues[1:5]/sum(eigenvalues)
  pc_labels <- sprintf("PC%d (%.2f%%)", 1:5, var_pct*100)
  
  class_batch <- metadata_df[colnames(pair_batch), "class"]
  class_colour <- CLASS_COLOUR_PALETTE[as.factor(class_batch)]
  batch_factor2 <- metadata_df[colnames(pair_batch), "batch"]
  batch_factor2 <- droplevels(batch_factor2)
  levels(batch_factor2) <- 21:22
  batch_pch <- as.numeric(as.character(batch_factor2))
  
  list_plot[[j]] <- ggplot(data.frame(pair_prcomp$x), aes(x = PC1, y = PC2)) +
    geom_point(size = 3, fill = class_colour, shape = batch_pch, show.legend = F) + 
    geom_vline(xintercept = 0, color = "black", alpha = 0.5) +
    geom_hline(yintercept = 0, color = "black", alpha = 0.5) +
    geom_abline(slope = ref_coef[2], intercept = ref_coef[1],
                color = "blue", alpha = 0.5) +
    geom_abline(slope = other_coef[2], intercept = other_coef[1],
                color = "blue", alpha = 0.5) +
    geom_abline(slope = mean_gradient,
                color = "orange", alpha = 0.5) +
    geom_abline(slope = batch_vec_pca[2]/batch_vec_pca[1],
                color = "orange", alpha = 0.5) +
    labs(x = pc_labels[1], y = pc_labels[2], title = plot_title) +
    theme(plot.title = element_text(hjust = 0.5)) +
    coord_fixed(ratio = 1)
  
  j <- j + 1
}

pairwise_pca <- do.call(plot_grid, c(list_plot, nrow = 3, ncol = 3))
pairwise_pca
ggsave("dump/pairwise_pca-gradient.pdf", pairwise_pca, width = 11, height = 11)

### PAIRWISE COMPARISONS OF BIO AND BATCH VECTORS
pairwise_comb <- combn(1:9,2)
pairwise_comb

angles_bio_vec <- apply(
  pairwise_comb, 2,
  function(cvec) calcAngleVectors(list_bio_vec[[cvec[1]]],
                                  list_bio_vec[[cvec[2]]])
)
distance_bio_vec <- apply(
  pairwise_comb, 2,
  function(cvec) calcL2Norm(list_bio_vec[[cvec[1]]] - list_bio_vec[[cvec[2]]])
)


angles_batch_vec <- apply(
  pairwise_comb, 2,
  function(cvec) calcAngleVectors(list_batch_vec[[cvec[1]]],
                                  list_batch_vec[[cvec[2]]])
)
distance_batch_vec <- apply(
  pairwise_comb, 2,
  function(cvec) calcL2Norm(list_batch_vec[[cvec[1]]] -
                              list_batch_vec[[cvec[2]]])
)

dev.new()
par(mfrow=c(1,2))
hist(angles_bio_vec, breaks = 30, main = NA,
     xlab = "Pairwise angles (biological vectors)")
hist(angles_batch_vec, breaks = 30, main = NA,
     xlab = "Pairwise angles (batch vectors)")
compare_hist <- recordPlot()
save_fig(compare_hist, "dump/hist-pairwise_angles.pdf",
         width = 8, height = 4)

sort(angles_bio_vec)

hist(distance_bio_vec, breaks = 30)
hist(distance_batch_vec, breaks = 30)

# Identifying genes that are most responsible for biological variation
list_bio_genes <- lapply(
  list_bio_vec,
  function(vec) names(vec)[order(abs(vec), decreasing = T)[1:100]]
)
bio_genes <- do.call(c, list_bio_genes)
# Out of 900 genes identified
length(unique(bio_genes))
barplot(table(table(bio_genes)),
        xlab = "Frequency of genes in 9 biological vectors")

# Identifying genes that are most responsible
list_batch_genes <- lapply(
  list_batch_vec,
  function(vec) names(vec)[order(abs(vec), decreasing = T)[1:100]]
)
batch_genes <- do.call(c, list_batch_genes)
# Out of 900 genes identified
length(unique(batch_genes))
barplot(table(table(batch_genes)),
        xlab = "Frequency of genes in 9 batch vectors")

par(mfrow=c(1,2))
compare_bar <- recordPlot()
save_fig(compare_bar, "dump/count-frequency_top100_loadings.pdf",
         width = 8, height = 4)

## Plot: Parallel coordinates
# Choose most significant bio genes
sig_bio_genes <- names(table(bio_genes))[table(bio_genes) == 9]

plot(log_yeoh[sig_bio_genes,1], ylim = c(0,15), type = "l")
for (i in 2:ncol(log_yeoh)) {
  lines(log_yeoh[sig_bio_genes,i], type = "l",
        col = metadata_df$batch_info[i])
}

subset_sigbio_mat <- data.matrix(log_yeoh[sig_bio_genes,])

# Annotation parameters
pheatmap(subset_sigbio_mat,
         treeheight_row = 0,
         annotation_col = full_metadata_df,
         show_colnames = F,
         fontsize_row = 8)
heatmap_sigbio <- recordPlot()
save_fig(heatmap_sigbio, "dump/heatmap_sigbio.pdf",
         width = 12, height = 6)
# Different set of genes perturbed for different batches
# Same set of genes are responsible for biological variation

# Histogram of loading coefficients
hist(abs(list_batch_vec[[1]]), breaks = 50)
hist_loadings <- recordPlot()
save_fig(hist_loadings, "dump/hist-loadings_batch.pdf")

# Try to perform BCM only along these subset of genes
batch_genes1 <- names(head(sort(abs(list_batch_vec[[1]]), decreasing = T), 100))
correction_vec <- list_batch_df[[1]] - list_batch_df[[2]]
correction_vec[, -batch_genes1]

# Plot log_yeoh filtered with most significant bio genes
label_factor <- full_metadata_df[colnames(log_yeoh), "label_info"]
label_factor <- addNA(label_factor)
levels(label_factor)[3] <- 0
levels(label_factor) <- c(4,6)
label_size <- as.numeric(levels(label_factor))[label_factor]
plot_pca1(log_yeoh[sig_bio_genes,], batch_colour, timepoint_shape, label_size)

# Calculate pairwise distance between batches -> To determine reference batch
tail(full_metadata_df, 60)

b2n8_d0 <- log_yeoh[, batch_info %in% c(2,8) & class_info == "D0"]
b2n8_batch <- full_metadata_df[colnames(b2n8_d0), "batch"]
b2n8_shape <- ifelse(as.character(b2n8_batch) == "2", 21, 22)
b2n8_label <- full_metadata_df[colnames(b2n8_d0), "subtype_info"]
b2n8_label
plot_pca1(b2n8_d0, 21)

# PAIRWISE PCA ------------------------------------------------------------
# Dataset wo batch 5 and ...
bcm_obj <- correctSVDBCM(data_yeoh, metadata_df, 2)
pairwise_pca <- do.call(plot_grid, c(bcm_obj$plot, nrow = 3, ncol = 3))
corrected_pca <- do.call(plot_grid, c(bcm_obj$corrected_plot,
                                      nrow = 3, ncol = 3))
pairwise_pca

ggsave("dump/pairwise_pca-foo.pdf",
       pairwise_pca, width = 11, height = 11)


bcm_df <- bcm_obj$data

x <- data.frame(PC1 = 1:5, PC2 = 1:5)
y <- factor(LETTERS[1:5])
baby <- y[c(2,2,2,1,1)]
ggplot(x, aes(x=PC1, y=PC2)) +
  geom_point(aes(fill = baby), pch = 21, size = 3)

plotPCA3DYeoh(bcm_df, metadata_df)
rgl.postscript("dump/pca_3d-bcm_pairwise_wo_batch5.pdf", "pdf")

intersect_probesets <- selectFeatures(bcm_df, metadata_df)
results_df <- calcERM_PCA(bcm_df[intersect_probesets,], metadata_df)

# QPSP --------------------------------------------------------------------
# Import NEA
NEA_RPATH <- paste0("../diff_expr/data/subnetwork/nea-hsa/",
                    "ovarian_cancer/geneset-nea_kegg_ovarian.tsv")
nea_df <- read.table(NEA_RPATH, sep = "\t", header = T, stringsAsFactors = F)
subnetwork_nea <- split(as.character(nea_df$gene_id), nea_df$subnetwork_id)

# 1. Removes affymetrix ambiguous and control probesets
# 2. Map probesets to IDs
# Removes one-to-many probesets and probesets with no ID
# Selects maximum if two probesets match to same gene
# CHECK: What microarray platform is the data from?
ENTREZ_GPL570 <- "../info/microarray/HG-U133_Plus_2/annot_entrez-GPL570.tsv"
SYMBOL_GPL570 <- "../info/microarray/HG-U133_Plus_2/annot_genesymbol-GPL570.tsv"
entrez_yeoh <- affy2id(data_yeoh, ENTREZ_GPL570)
symbol_yeoh <- affy2id(data_yeoh, SYMBOL_GPL570)
head(symbol_yeoh)

# Calculate QPSP profiles
gfs_yeoh <- normaliseGFS(entrez_yeoh, num_intervals = 4)
qpsp_yeoh <- calcQPSP(gfs_yeoh, subnetwork_nea)
dim(qpsp_yeoh)
dim(subset_yeoh)

# No need for feature selection as features have been reduced
SUBTYPE <- levels(full_metadata_df$subtype)[6]
print(SUBTYPE)

normal_pid <- paste0("N0", c(1,2,4))
subtype_pid <- rownames(subset(full_metadata_df,
                               subtype == SUBTYPE & class_info != "N" & 
                                 rownames(full_metadata_df) %in% 
                                 colnames(qpsp_yeoh)))
subset_pid <- c(subtype_pid, normal_pid)
# subtype_metadata <- full_metadata_df[subset_pid,]

# Subtype and normal samples
subtype_qpsp <- qpsp_yeoh[,subset_pid]

### DOWNSAMPLE HETEROGENEOUS SUBTYPES ###
# Sample size = 40 (33 remission, 7 relapse)
subtype_metadata <- full_metadata_df[colnames(subtype_qpsp),]
subtype_metadata1 <- subtype_metadata[endsWith(rownames(subtype_metadata), "_D0"),]
pid_0 <- rownames(subtype_metadata1)[subtype_metadata1$label == 0]
pid_1 <- rownames(subtype_metadata1)[subtype_metadata1$label == 1]

sample_pid_0 <- sample(pid_0, 33)
sample_pid_1 <- sample(pid_1, 7)
sample_pid <- substring(c(sample_pid_0, sample_pid_1), 1, 4)

logi_idx <- substring(colnames(subtype_qpsp), 1, 4) %in% sample_pid
logi_idx[171:173] <- TRUE
sampled_qpsp <- subtype_qpsp[,logi_idx]
writeLines(colnames(sampled_qpsp), con = "dump/sampled_pid.txt")

# Check that there are D0 samples in batch 2
table(metadata_df[colnames(sampled_qpsp),])
###

# BCM on QPSP transformed data
# bcm_qpsp <- correctGlobalBCM(subtype_qpsp, metadata_df)
### OPTION 2
bcm_qpsp <- correctGlobalBCM(sampled_qpsp, metadata_df)

# PCA
pca_obj <- prcomp(t(bcm_qpsp))
# PCA: Eigenvalues
eigenvalues <- (pca_obj$sdev)^2
var_pc <- eigenvalues/sum(eigenvalues)
# Use till PC15 to account for at least 70% variance
cum_var <- cumsum(var_pc)
cum_var
# Identify PC that has just above 70% variance
pc_ind <- which.max(cum_var > 0.70)
print(pc_ind)

# PCA: Coordinates
pca_coord <- pca_obj$x[,1:pc_ind]
# Response df and normal df
rownames(pca_coord)
n_idx <- nrow(pca_coord) - 3
response_df <- pca_coord[1:n_idx,]
normal_df <- pca_coord[-(1:n_idx),]
print(rownames(response_df))
# CHECK
print(rownames(normal_df))

# # PLOT
# pc_labels <- sprintf("PC%d (%.2f%%)", 1:3, (var_pc*100)[1:3])
# color1 <- rep("tomato3", 81)
# color1[c(10,49)] <- "darkolivegreen3"
# pch1 <- rep(21:23, c(39,39,3))
# plotPCA3D(pca_coord[,1:3], color1, pch1, pc_labels)

# Collate MRD results as well
results_df <- calcERM(response_df, normal_df, yeoh_label)
results_df
SUBTYPE
write.table(results_df, sprintf("dump/features2-%s.tsv", SUBTYPE),
            sep = "\t", quote = F)

# par(mar=c(8,3,2,3))
# # Standardise each feature
# calibrated_features <- cbind(
#   results_df[,1:10],
#   results_df[,11:14]/60, # angle
#   -log10(results_df[,15])/2 # mrd
# )

# # Plotting parameters
# line_color <- ifelse(results_df[,16] == 1, "tomato3", "darkolivegreen3")
# feature_names <- c(colnames(results_df)[1:11], "ang_d0d8_normal",
#                    colnames(results_df)[13:14], "-log10(mrd)")
# 
# # Plot parallel coordinates
# plot(0, xlim = c(1,15), ylim = c(-2,3),
#      type = "n", xaxt = "n", ann = F)
# axis(1, at = 1:15, feature_names, cex = 0.5)
# for (r in 1:nrow(calibrated_features)) {
#   lines(1:15, calibrated_features[r,], col = line_color[r])
# }
# 
# subtype_parallel <- recordPlot()
# save_fig(subtype_parallel, sprintf("dump/parallel-%s.pdf", SUBTYPE),
#          width = 10, height = 7)

# LIKELIHOOD RATIO --------------------------------------------------------
curve(dnorm(x,3), xlim=c(-3,3))

restrictor <- function(x) x/(x+1)
calcLR <- function(x) restrictor(dnorm(x,-3)/dnorm(x,3))
curve(calcLR, xlim = c(-3,3))

x <- rnorm(10000)
hist(x, breaks = 100)

full_metadata_df[full_metadata_df$subtype == "Hypodiploid",]

##### SANDBOX: Same or different #####
# # Cosine normalise data
# normalised_yeoh <- normaliseCosine(data_yeoh)

### Subset all D0 patients
data_d0 <- data_yeoh[,endsWith(colnames(data_yeoh), "0")]
d0_metadata <- metadata_df[colnames(data_d0),]
# # Order metadata
# d0_metadata_ord <- d0_metadata[order(d0_metadata$batch_info),]
# data_d0_ordered <- data_d0[,rownames(d0_metadata_ord)]

### FEATURE SELECTION (CHI2) ###
# Add class labels to each individual sample (row)
## Run the chi2 algorithm n times (one vs rest)
t_d0 <- t(data_d0)
subtype <- d0_metadata[rownames(t_d0), "subtype"] -> original_subtype
subtype <- original_subtype

chi2_features <- character()
for (i in 1:9) {
  subtype <- original_subtype
  levels(subtype)[-i] <- "Rest"
  weka_d0 <- data.frame(t_d0, subtype)
  selected_d0 <- chi2.algorithm(weka_d0, numeric(), threshold = 0)
  # Selected features
  chi2_features <- c(chi2_features, colnames(selected_d0[[1]]))
}

writeLines(substring(chi2_features, 2), con = "dump/chi2_features.txt")


X <- data_d0[substring(chi2_features, 2),]
y <- metadata_df[colnames(X),]
X <- X[,order(y$subtype)]

subtype_df <- metadata_df[colnames(data_d0), "subtype", drop = F]
pheatmap(X, col = brewer.pal(9, "Blues"),
         legend = T, border_color = NA, scale = "none",
         cluster_rows = T, cluster_cols = F,
         show_colnames = F, show_rownames = F,
         annotation_col = subtype_df)
heatmap <- recordPlot()
save_fig(heatmap, "dump/heatmap-chi2.pdf")]

# ## Save only D0 patients for WEKA
# # Add class labels to each individual sample (row)
# t_d0 <- t(data_d0)
# subtype <- d0_metadata[rownames(t_d0), "subtype"]
# weka_d0 <- cbind(t_d0, subtype)
# 
# write.table(weka_d0, "temp/weka_d0.csv", sep = ",", quote = F,
#             row.names = F, col.names = T)

## Batch-subtype configuration
table(d0_metadata_ord$batch_info, d0_metadata_ord$subtype)

## Entire batch 1 D0
pw_dist <- dist(t(data_d0_ordered), method = "euclidean")
dist_mat <- as.matrix(pw_dist)

i <- 1
x_i <- dist_mat[i,]
x_i1 <- x_i[x_i != 0]
hist(x_i1, breaks = 30)
x_i1

# Batch colours
batch_ordered <- d0_metadata[colnames(data_d0_ordered), "batch_info"]
batch_factor <- as.factor(batch_ordered)
levels(batch_factor) <- brewer.pal(10, "Set3")
batch_col <- as.character(batch_factor)

heatmap(dist_mat, Rowv = NA, Colv = "Rowv", scale = "none",
        labRow = NA, labCol = NA, col = brewer.pal(9, "Blues"),
        ColSideColors = batch_col, RowSideColors = batch_col)

# heatmap_plot <- recordPlot()
# save_fig(heatmap_plot, "dump/heatmap-dist_all_batches.pdf",
#          width = 8, height = 8)

# Choosing the batches
pid_1 <- rownames(d0_metadata_ord)[d0_metadata_ord$batch_info == 8]
pid_2 <- rownames(d0_metadata_ord)[d0_metadata_ord$batch_info == 9]
all_pid <- c(pid_1, pid_2)

# Batch 1 vs Batch 2
pair_dist <- dist_mat[pid_1, pid_2]
pheatmap(pair_dist, col = brewer.pal(9, "Blues"),
         legend = F, border_color = NA, scale = "none",
         cluster_rows = F, cluster_cols = F)

heatmap_b8b9 <- recordPlot()
save_fig(heatmap_b8b9, "dump/heatmap-dist_b8b9.pdf",
         width = 12, height = 8)

# Subtype breakdown for pair of batches
pair_metadata <- d0_metadata_ord[c(pid_1, pid_2),]
table(pair_metadata$subtype, pair_metadata$batch_info)

getNN <- function(pairwise_mat, flag) {
  # Distance: Choose smallest distance
  # Correlation: Choose largest correlation
  if(flag == "dist") FUNC <- which.min
  else FUNC <- which.max
  
  # 1-NN for batch 1 samples
  for(p1 in rownames(pairwise_mat)) {
    p1_subtype <- as.character(d0_metadata[p1,"subtype"])
    
    p2 <- names(FUNC(pairwise_mat[p1,]))
    p2_subtype <- as.character(d0_metadata[p2,"subtype"])
    
    msg <- sprintf("%s (%s): %s (%s)\n", p1, p1_subtype, p2, p2_subtype)
    cat(msg)
    
    # dist_incr <- sort(pairwise_mat[p1,])
    # names(dist_incr) <- metadata_df[names(dist_incr), "subtype"]
    # print(dist_incr)
    
    # hist(pairwise_mat[p1,], breaks = 10, main = p1)
    # abline(v = min(pairwise_mat[p1,]), col = "red")
  }
  
  cat("\n")
  
  # 1-NN for batch 1 samples
  for(p2 in colnames(pairwise_mat)) {
    p2_subtype <- as.character(d0_metadata[p2,"subtype"])
    
    p1 <- names(FUNC(pairwise_mat[,p2]))
    p1_subtype <- as.character(d0_metadata[p1,"subtype"])
    
    msg <- sprintf("%s (%s): %s (%s)\n", p2, p2_subtype, p1, p1_subtype)
    cat(msg)
    
    # dist_incr <- sort(pairwise_mat[,p2])
    # names(dist_incr) <- metadata_df[names(dist_incr), "subtype"]
    # print(dist_incr)
    
    # hist(pairwise_mat[,p2], breaks = 10, main = p2)
    # abline(v = min(pairwise_mat[,p2]), col = "red")
  }  
}

getNN(pair_dist, "dist")
getNN(pair_cor, "cor")
getNN(pair_cor1, "cor")

pair_cor1[1,]



plotPCA3DYeoh1(data_yeoh[,all_pid], metadata_df)
# rgl.postscript("dump/pca_3d-b7b8_d0.pdf", "pdf")

### Chi-square selected genes from yeoh_2002 ###
# Neither probesets nor gene symbols correspond with current technology
CHI_GENES <- "data/yeoh_2002/README/chi_square_probesets/combined.tsv"
chi_genes <- read.delim(CHI_GENES)
chi_symbols <- levels(chi_genes$gene_symbol)[-1]

SYMBOL_GPL570 <- "../info/microarray/HG-U133_Plus_2/annot_genesymbol-GPL570.tsv"
affy_yeoh <- affy2id(data_yeoh, SYMBOL_GPL570)
intersect(chi_symbols, rownames(affy_yeoh))

plotPCA3DYeoh1(data_yeoh[chi_genes,all_pid], metadata_df)

### SELECTING GENE SIGNATURES OF SUBTYPES ###
# Using only D0 samples
# T-test based on the ranks
getTTestFeatures <- function(data, subtype) {
  col_idx <- which(metadata_df[colnames(data), "subtype"] == subtype)
  ordered_data <- cbind(data[,col_idx], data[,-col_idx])
  pvalue <- calc_ttest(ordered_data, length(col_idx), is_paired = F)
  names(head(sort(pvalue), 20))
}

filtered_d0 <- filterProbesets(data_d0, 0.95, metadata_df)

list_probesets <- sapply(levels(metadata_df$subtype),
                    getTTestFeatures, data = filtered_d0)
# Hypodiploid and Normal are not included
probesets <- unique(unlist(list_probesets, use.names = F))

# Ordering data_d0 according to subtype
subtype_info <- metadata_df[colnames(data_d0), "subtype"]
col_idx <- order(subtype_info)

subtype_df <- metadata_df[colnames(data_d0), "subtype", drop = F]
ordered_d0 <- data_d0[probesets, col_idx]
ps_median <- apply(ordered_d0, 1, median)
ps_mad <- apply(ordered_d0, 1, mad)
znorm_d0 <- (ordered_d0 - ps_median)/ ps_mad
znorm_d0[1:5,1:5]
cust_pal <- c(rep("#67001F", 45),
              brewer.pal(11, "RdBu"),
              rep("#053061", 5))

# FOR PROGNOSTIC MARKERS. USUALLY UNIQUELY HIGH EXPRESSION..
# ..FOR THAT SUBTYPE, REST ARE ZERO VALUES

pheatmap(znorm_d0, color = cust_pal,
         legend = T, border_color = NA, scale = "none",
         cluster_rows = F, cluster_cols = F,
         annotation_col = subtype_df)

# Correlation ignoring all zeros
# Correlated data
a <- rnorm(100)
b <- rnorm(100) + a
plot(a,a*3)

# All zeroes
# Both zeroes does not affect correlation
c <- rnorm(100)
d <- c(rep(0,100), c)
e <- c(rep(0,100), rnorm(100))
# cbind(d,e)
plot(d,e)
cor(d,e)

# Zero vs RV
# More zero-signal values decrease correlation
f <- rnorm(100) + 5
g <- rep(f, 21)
h <- c(rep(0, 2000), f)
plot(g, h, xlim = c(0,10))
cor(g,h)

# Pearson correlation: Zero-inflated data
cor(a, b, method = "pearson")
cor(a, b, method = "spearman")
cor(a, b, method = "kendall")

## TEST
# Batch 8 and 9 metadata
d0_metadata_ord[all_pid,]
p1 <- data_d0_ordered[,"P103_D0"]
p2 <- data_d0_ordered[,"P151_D0"] # CORRECT
p3 <- data_d0_ordered[,"P137_D0"] # WRONG

calcCor(p1, p2)
calcCor(p1, p3)

# corr_p1 <- sapply(pid_1, function(pid) cor(p1, data_d0_ordered[,pid]))
# which.max(corr_p1)

fig1 <- ggplot(data.frame(p1, p2), aes(x = p1, y = p2)) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon")

fig2 <- ggplot(data.frame(p1, p3), aes(x = p1, y = p3)) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon")

density_plot <- plot_grid(fig1, fig2, ncol = 2)
density_plot
 
par(mfrow=c(1,2))
plot(p1, p2) # Correct
plot(p1, p3) # Wrong
cor_plot <- recordPlot()
save_fig(cor_plot, "dump/scatter-cor.pdf", 12, 6)

# Investigating patterns in zero values
# No pattern in dispersion of values that accompany zeros
logi1 <- p1 == 0
logi2 <- p2 == 0
logi3 <- p3 == 0
xor_logi1 <- xor(logi1, logi2)
xor_logi2 <- xor(logi1, logi3)

sum(xor_logi1)
sum(logi1 & logi2)

sum(xor_logi2)
sum(logi1 & logi3)
# P3 has disproportionately little zeros? Try other patientss

mat1 <- cbind(p1, p2)[xor_logi1,]
mat2 <- cbind(p1, p3)[xor_logi2,]

hist(mat1[mat1[,1] == 0, 2], breaks = 15,
     xlab = NULL, main = "Zero-signal values [P1 (0) vs P2]")
hist(mat2[mat2[,1] == 0, 2], breaks = 16,
     xlab = NULL, main = "Zero-signal values [P1 (0) vs P3]")
hist_plot <- recordPlot()
save_fig(hist_plot, "dump/hist-xzero_ysignal.pdf", 12, 6)

cor(p1, p2)
cor(p2, p1)

#' @param X df with samples as columns and features as rows
calcPairwise <- function(X, pid_grp1, pid_grp2, FUNC, ...) {
  sapply(pid_grp2, function(p2) {
    sapply(pid_grp1, function(p1) FUNC(X[,p2], X[,p1], ...))
  })
}

calcCor <- function(v1, v2) {
  logi_idx <- v1 != 0 & v2 != 0
  cor(v1[logi_idx], v2[logi_idx])
}

pair_cor <- calcPairwise(data_d0_ordered, pid_1, pid_2, cor)
pair_cor1 <- calcPairwise(data_d0_ordered, pid_1, pid_2, calcCor)

pheatmap(pair_cor, col = brewer.pal(9, "Blues"),
         legend = F, border_color = NA, scale = "none",
         cluster_rows = F, cluster_cols = F)

heatmap_b8b9_cor <- recordPlot()
save_fig(heatmap_b8b9_cor, "dump/heatmap-cor_b8b9.pdf",
         width = 12, height = 8)

epsilon <- rnorm(3000, 10, 1000)

sampleError <- function(n, k = 0.2, beta = 0.001) {
  c(rgamma(n/2, k, beta), -1 * rgamma(n/2, k, beta))  
}

x <- 2^rgamma(5000, shape = 24.2, rate = 2.7)
epsilon <- sampleError(5000)
plot(log2(x), log2(x + epsilon),
     xlim = c(0, 20), ylim = c(0,20))

orig_p1 <- 2^p1
p1_p <- p1[p1 != 0]
hist(orig_p1[orig_p1 != 1], breaks = 20)

hist(p1_p, breaks = 20)
hist(rgamma(5000, 24.2, 2.7), breaks = 20)

hist(p1, breaks = 30, prob = T)
curve(dgamma(x, 24.2, 2.7), xlim = c(0,15), add = T)

curve(dgamma(x, 4, 1), xlim = c(0,15))
# All timepoints D0, D8, Normal
subset_metadata <- metadata_df[colnames(data_yeoh),]
logi_idx <- subset_metadata$batch_info == 8 | subset_metadata$batch_info == 9
pair_data <- filtered_yeoh[,logi_idx]
plotPCA3DYeoh1(pair_data, metadata_df)
rgl.postscript("dump/pca_3d-b7b8.pdf", "pdf")

# ## Cosine normalisation
# normalised_yeoh <- normaliseCosine(pair_data)
# plotPCA3DYeoh(normalised_yeoh, metadata_df)

hcluster <- hclust(pw_dist)
dendo_obj1 <- as.dendrogram(hcluster)

##### PLOTS #####
# Choosing the batches
pid_1 <- rownames(d0_metadata)[d0_metadata$batch_info == 1]
pid_2 <- rownames(d0_metadata)[d0_metadata$batch_info == 2]
all_pid <- c(pid_1, pid_2)
batch_1_2 <- data_d0[,all_pid]

## B1 vs B2
plotPCA3DYeoh1(batch_1_2, metadata_df)
pca_obj <- prcomp(t(batch_1_2))
rot_mat <- pca_obj$rotation

# PC1
desc_ps <- names(sort(abs(rot_mat[,1]), decreasing = TRUE))
b_1_2_pc1_1 <- batch_1_2[desc_ps[1:21],]
b_1_2_pc1_2 <- batch_1_2[desc_ps[101:121],]

# PC2
desc_ps <- names(sort(abs(rot_mat[,2]), decreasing = TRUE))
b_1_2_pc2_1 <- batch_1_2[desc_ps[1:21],]
b_1_2_pc2_2 <- batch_1_2[desc_ps[5001:5021],]

x <- b_1_2_pc2_2
figname <- "b_1_2_pc2_2"
plotMulti(x, metadata_df)
ps_plot <- recordPlot()
save_fig(ps_plot, sprintf("dump/%s-scatter.pdf", figname), 11, 15)

plotMultiHist(x, metadata_df)
ps_hist <- recordPlot()
save_fig(ps_hist, sprintf("dump/%s-hist.pdf", figname), 11, 15)



## B8 vs B9
plotPCA3DYeoh1(batch_8_9, metadata_df)
pca_obj <- prcomp(t(batch_8_9))
rot_mat <- pca_obj$rotation

# PC1
desc_ps <- names(sort(abs(rot_mat[,1]), decreasing = TRUE))
b_8_9_pc1_1 <- batch_8_9[desc_ps[1:21],]
b_8_9_pc1_2 <- batch_8_9[desc_ps[5001:5021],]

# PC2
desc_ps <- names(sort(abs(rot_mat[,2]), decreasing = TRUE))
b_8_9_pc2_1 <- batch_8_9[desc_ps[1:21],]
b_8_9_pc2_2 <- batch_8_9[desc_ps[5001:5021],]

x <- b_8_9_pc2_2
figname <- "b_8_9_pc2_2"
plotMulti(x, metadata_df)
ps_plot <- recordPlot()
save_fig(ps_plot, sprintf("dump/%s-scatter.pdf", figname), 11, 15)

plotMultiHist(x, metadata_df)
ps_hist <- recordPlot()
save_fig(ps_hist, sprintf("dump/%s-hist.pdf", figname), 11, 15)

# Batch effect
ps4 <- "214250_at" #b12
ps1 <- "210116_at" #b12
ps5 <- "218895_at" #b12

# Class effect
ps8 <- "211005_at" #b12

# No effect
ps2 <- "203209_at" #b12
ps3 <- "219935_at" #b12

# Both
ps6 <- "221566_s_at" #maqc
ps7 <- "201608_s_at" #maqc

rownames(subset_maqc)
x <- batch_1_2[c(ps4, ps1, ps5),]

plotSample(x, metadata_df)
ps_sample <- recordPlot()
save_fig(ps_sample, sprintf("dump/ps-batch.pdf"), 10, 9)

plotMulti <- function(X, metadata_df, subplot=c(7,3)) {
  cat(dim(X))
  n_feat <- nrow(X)
  n <- ncol(X)
  ### Plot parameters
  # Shape: Batch
  batch_info <- metadata_df[colnames(X), "batch_info"]
  batch_factor <- droplevels(as.factor(batch_info))
  print(batch_factor)
  print(levels(batch_factor))
  levels(batch_factor) <- 16:17
  pch <- as.numeric(as.character(batch_factor))
  # Colour: Subtype
  class_info <- metadata_df[colnames(X), "subtype"]
  palette <- brewer.pal(9, "Set1")
  col <- palette[class_info]
  
  par(mfrow=subplot, mar=rep(2.3,4))
  for (i in 1:n_feat) {
    plot(1:n, X[i,], cex = 1.5, col = col, pch = pch, main = rownames(X)[i])
  }
  par(mfrow=c(1,1))
}

plotMultiHist <- function(X, metadata_df, subplot=c(7,3)) {
  cat(dim(X))
  n_feat <- nrow(X)
  n <- ncol(X)
  ### Plot parameters
  # Shape: Batch
  batch_info <- metadata_df[colnames(X), "batch_info"]
  batch_factor <- droplevels(as.factor(batch_info))
  print(batch_factor)
  print(levels(batch_factor))
  levels(batch_factor) <- 16:17
  pch <- as.numeric(as.character(batch_factor))
  # Colour: Subtype
  class_info <- metadata_df[colnames(X), "subtype"]
  palette <- brewer.pal(9, "Set1")
  col <- palette[class_info]
  
  par(mfrow=subplot, mar=rep(2.3,4))
  for (i in 1:n_feat) {
    hist(data.matrix(X[i, batch_factor == 16]),
         main = rownames(X)[i],
         seq(0,16,0.5), col = rgb(1,0,0,0.2))
    hist(data.matrix(X[i, batch_factor == 17]),
         seq(0,16,0.5), col = rgb(0,1,0,0.2), add = T)
  }
  par(mfrow=c(1,1))
}

plotSample <- function(X, metadata_df, subplot=c(3,2)) {
  cat(dim(X))
  n_feat <- nrow(X)
  n <- ncol(X)
  ### Plot parameters
  # Shape: Batch
  batch_info <- metadata_df[colnames(X), "batch_info"]
  batch_factor <- droplevels(as.factor(batch_info))
  print(batch_factor)
  print(levels(batch_factor))
  levels(batch_factor) <- 16:17
  pch <- as.numeric(as.character(batch_factor))
  # Colour: Subtype
  class_info <- metadata_df[colnames(X), "subtype"]
  palette <- brewer.pal(9, "Set1")
  col <- palette[class_info]
  
  par(mfrow=subplot, mar=rep(2.5,4))
  for (i in 1:n_feat) {
    plot(1:n, X[i,], cex = 1.5, col = col, pch = pch, main = rownames(X)[i])
    hist(data.matrix(X[i, batch_factor == 16]),
         main = rownames(X)[i],
         seq(0,16,0.5), col = rgb(1,0,0,0.2))
    hist(data.matrix(X[i, batch_factor == 17]),
         seq(0,16,0.5), col = rgb(0,1,0,0.2), add = T)
  }
  par(mfrow=c(1,1))
}


