library(xtable)
# library(Rtsne)
# library(dendextend)
# library(gPCA)
# library(cluster)
# library(pheatmap)
library(sva)
# library(scran)
library(Harman)
library(reshape2)
library(rgl)
library(ggplot2)
library(cowplot)
library(reshape2)
library(RColorBrewer)
library(dplyr)

source("../functions.R")
source("bcm.R")

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

# Plot PCA before selecting features
# Batch information of all the timepoints
plotPCA3DYeoh <- function(df1, metadata_df) {
  batch_info <- metadata_df[colnames(df1), "batch"]
  generate_colour <- colorRampPalette(c("lightblue", "darkblue"))
  batch_palette <- generate_colour(10)
  batch_colour <- batch_palette[batch_info]
  # Shape of all timepoints
  class_info <- metadata_df[colnames(df1), "class"]
  levels(class_info) <- 21:23
  timepoint_shape <- as.numeric(as.character(class_info))
  plotPCA3D(df1, batch_colour, timepoint_shape)
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
    geom_point(aes(x = PC1_A, y = PC2_A), size = 5, stroke = 2,
               colour = subtype_colour, shape = 21, fill = batch_colour, show.legend = F) +
    geom_point(aes(x = PC1_B, y = PC2_B), size = 5, stroke = 2,
               colour = subtype_colour, shape = 22, fill = batch_colour, show.legend = F) +
    geom_point(data = centroid_df, aes(x = PC1, y = PC2),
               size = 5, shape = 17,
               colour = c("tomato3", "orange", "hotpink", "pink")) +
    geom_segment(aes(x = PC1_A, y = PC2_A, xend = PC1_B, yend = PC2_B, colour = as.factor(labels)),
                 arrow = arrow(length = unit(0.3, "cm")), alpha = 0.5, show.legend = F) +
    scale_color_manual(values = c("black",  "red")) +
    xlab(pc_labels[1]) + ylab(pc_labels[2])
  pca_2 <- ggplot(data = df) +
    geom_point(aes(x = PC3_A, y = PC2_A), size = 5, stroke = 2,
               colour = subtype_colour, shape = 21, fill = batch_colour, show.legend = F) +
    geom_point(aes(x = PC3_B, y = PC2_B), size = 5, stroke = 2,
               colour = subtype_colour, shape = 22, fill = batch_colour, show.legend = F) +
    geom_point(data = centroid_df, aes(x = PC2, y = PC3),
               size = 5, shape = 17,
               colour = c("tomato3", "orange", "hotpink", "pink")) +
    geom_segment(aes(x = PC3_A, y = PC2_A, xend = PC3_B, yend = PC2_B, colour = as.factor(labels)),
                 arrow = arrow(length = unit(0.3, "cm")), alpha = 0.5, show.legend = F) +
    scale_color_manual(values = c("black",  "red")) +
    xlab(pc_labels[3]) + ylab(pc_labels[2])
  multiplot <- plot_grid(pca_1, pca_2, ncol = 2)
  return(multiplot)
}

# IMPORT DATA -------------------------------------------------------------
yeoh_d0d8 <- read.table("data/GSE67684/processed/mas5_ordered.tsv",
                        sep = "\t", header = T, row.names = 1)
yeoh_d33 <- read.table("data/leuk_D33/processed/mas5_filtered.tsv",
                       sep = "\t", header = T, row.names = 1)
yeoh_normal <- read.table("data/leuk_normal/processed/mas5_filtered.tsv",
                          sep = "\t", header = T, row.names = 1)
# Metadata
yeoh_batch <- read.table("data/GSE67684/processed/metadata_combined-batch.tsv",
                         sep = "\t", header = T, row.names = 1)
yeoh_label <- read.table("data/GSE67684/processed/metadata_combined-label_subtype_edited.tsv",
                         sep = "\t", header = T, row.names = 1)

# Edit yeoh_label: MRD33
# Change all equivalent strings to the same
levels(yeoh_label$d33_mrd)[1:4] <- levels(yeoh_label$d33_mrd)[1]
# Assign <10E-04
levels(yeoh_label$d33_mrd)[1] <- 0.00001
yeoh_label$d33_mrd <- addNA(yeoh_label$d33_mrd)
# Assign NA <- 0
levels(yeoh_label$d33_mrd)[35] <- 0
yeoh_label$d33_mrd <- as.numeric(as.character(yeoh_label$d33_mrd))

# Removal of outlier samples and their associated pairs
# Patients "P198_D8", "P186_D0" are associated pairs
outlier_samples <- c("P198_D0", "P186_D8", "N03", "P198_D8", "P186_D0")
yeoh_d0d8 <- yeoh_d0d8[, !(colnames(yeoh_d0d8) %in% outlier_samples)]
yeoh_normal <- yeoh_normal[, !(colnames(yeoh_normal) %in% outlier_samples)]

# Selection of d33 remission samples
d33_label <- yeoh_label[substring(colnames(yeoh_d33), 1, 4), "label", drop = F]
# D33 patients that experience remission
d33_remission <- rownames(d33_label)[d33_label == 0 & !is.na(d33_label)]
yeoh_remission <- yeoh_d33[, paste0(d33_remission, "_D33")]

# SCALE & FILTER & LOG -------------------------------------------------------
yeoh_combined <- cbind(yeoh_d0d8, yeoh_remission, yeoh_normal)

# Create metadata df
batch_info <- as.factor(yeoh_batch[colnames(yeoh_combined), "batch"])
class_info <- rep(c("D0","D8","N"), c(208,208,45))
class_numeric <- rep(1:3, c(208,208,45))
metadata_df <- data.frame(batch_info, class_info)
rownames(metadata_df) <- colnames(yeoh_combined)

# Add subtype info to metadata
subtype <- yeoh_label[substr(colnames(yeoh_combined), 1, 4), "subtype"]
label <- as.factor(yeoh_label[substr(colnames(yeoh_combined), 1, 4),
                                   "label"])
full_metadata_df <- cbind(metadata_df, subtype, label)

### SFL ###
processed_yeoh <- log2_transform(
  filterProbesets(normaliseMeanScaling(yeoh_combined), 0.3, metadata_df)
)

# SUBSET DATA -------------------------------------------------------------
# # Identify patients with D0 and D8 profiles from different batches
# for (i in seq(1,420,2)) {
#   if (yeoh_batch[i,"batch"] != yeoh_batch[i+1,"batch"]) {
#     print(yeoh_batch[i:(i+1),])
#   }
# }

patients_diffbatch <- c("P107", "P110", "P112", "P113", "P114", "P118", "P168")

# REMOVE BATCH 5 AND SAMPLES FIRST
# Taking out batch 5 samples and P107_D0 (its pair present in batch5)
not_diffbatch <-
  !(substring(colnames(processed_yeoh), 1, 4) %in% patients_diffbatch)
not_batch5 <- metadata_df[colnames(processed_yeoh), "batch_info"] != 5
subset_yeoh <- processed_yeoh[, not_batch5 & not_diffbatch]

# QUANTILE ---------------------------------------------------------
# ### QUANTILE (ALL)
# quantile_yeoh <- normaliseQuantile(subset_yeoh)
colnames(subset_yeoh[403:405])
### QUANTILE (TIMEPOINT)
quantile_d0 <- normaliseQuantile(subset_yeoh[, 1:201])
quantile_d8 <- normaliseQuantile(subset_yeoh[, 202:402])
quantile_normal <- normaliseQuantile(subset_yeoh[, 403:405])
quantile_yeoh <- cbind(quantile_d0, quantile_d8, quantile_normal)

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
# Quantile by timepoint
quantile_d0 <- norm_quantile(log_yeoh[, 1:208])
quantile_d8 <- norm_quantile(log_yeoh[, 209:416])
quantile_normal <- norm_quantile(log_yeoh[, 417:461])
quantile_yeoh <- cbind(quantile_d0, quantile_d8, quantile_normal)

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
# Splitting up yeoh data into batches
batch_info <- yeoh_batch[colnames(log_yeoh), "batch"]
batch_list <- lapply(1:10, function(i) log_yeoh[, batch_info == i])
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
ordered_mnn_yeoh <- mnn_yeoh[, colnames(log_yeoh)]

# Replace negative values with zero
ordered_mnn_yeoh[ordered_mnn_yeoh < 0] <- 0

# Plot PCA before selecting features
plotPCA3DYeoh(ordered_mnn_yeoh, metadata_df)
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

# PCA PLOT ----------------------------------------------------------------
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
bcm_obj <- correctSVDBCM(subset_yeoh, metadata_df, 2)
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
nea_df <- read.table(paste0("../diff_expr/data/subnetwork/nea-hsa/",
                            "ovarian_cancer/geneset-nea_kegg_ovarian.tsv"),
                     sep = "\t", header = T, stringsAsFactors = F)
subnetwork_nea <- split(as.character(nea_df$gene_id), nea_df$subnetwork_id)

# 1. Removes affymetrix ambiguous and control probesets
# 2. Map probesets to IDs
# Removes one-to-many probesets and probesets with no ID
# Selects maximum if two probesets match to same gene
# CHECK: What microarray platform is the data from?
ANNOT_PROBESET_RPATH <- "../info/microarray/HG-U133A/annot_entrez-GPL96.tsv"
entrez_yeoh <- affy2id(removeProbesets(subset_yeoh), ANNOT_PROBESET_RPATH)

# # Perform BCM before QPSP
# bcm_yeoh <- correctGlobalBCM(entrez_yeoh, metadata_df)

# Calculate QPSP profiles
gfs_yeoh <- normaliseGFS(entrez_yeoh, num_intervals = 4)
qpsp_yeoh <- calcQPSP(gfs_yeoh, subnetwork_nea)

dim(qpsp_yeoh)
dim(subset_yeoh)
# # Plot PCA before selecting features
# plotPCA3DYeoh(qpsp_yeoh, metadata_df)
# rgl.postscript("dump/pca_3d-qpsp_nea.pdf", "pdf")

# No need for feature selection as features have been reduced
SUBTYPE <- levels(full_metadata_df$subtype)[9]
SUBTYPE
normal_pid <- paste0("N0", c(1,2,4))

subtype_pid <- rownames(subset(full_metadata_df,
                               subtype == SUBTYPE & class_info != "N" & 
                                 rownames(full_metadata_df) %in% 
                                 colnames(qpsp_yeoh)))
subset_pid <- c(subtype_pid, normal_pid)

subtype_metadata <- full_metadata_df[subset_pid,]

# Subtype and normal samples
subtype_qpsp <- qpsp_yeoh[,subset_pid]

# BCM on QPSP transformed data
bcm_qpsp <- correctGlobalBCM(subtype_qpsp, metadata_df)

# PCA
pca_obj <- prcomp(t(selected_quantile))
# PCA: Eigenvalues
eigenvalues <- (pca_obj$sdev)^2
var_pc <- eigenvalues/sum(eigenvalues)
# Use till PC15 to account for at least 70% variance
cum_var <- cumsum(var_pc)[1:100]
cum_var
# Identify PC that has just above 70% variance
pc_ind <- which.max(cum_var > 0.70)
pc_ind # PC20

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
SUBTYPE <- "cs_quantile"
write.table(results_df, sprintf("dump/features-%s.tsv", SUBTYPE),
            sep = "\t", quote = F)

par(mar=c(8,3,2,3))
# Standardise each feature
calibrated_features <- cbind(
  results_df[,1:10],
  results_df[,11:14]/60, # angle
  -log10(results_df[,15])/2 # mrd
)

# Plotting parameters
line_color <- ifelse(results_df[,16] == 1, "tomato3", "darkolivegreen3")
feature_names <- c(colnames(results_df)[1:11], "ang_d0d8_normal",
                   colnames(results_df)[13:14], "-log10(mrd)")

# Plot parallel coordinates
plot(0, xlim = c(1,15), ylim = c(-2,3),
     type = "n", xaxt = "n", ann = F)
axis(1, at = 1:15, feature_names, cex = 0.5)
for (r in 1:nrow(calibrated_features)) {
  lines(1:15, calibrated_features[r,], col = line_color[r])
}

subtype_parallel <- recordPlot()
save_fig(subtype_parallel, sprintf("dump/parallel-%s.pdf", SUBTYPE),
         width = 10, height = 7)
