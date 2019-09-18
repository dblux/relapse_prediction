# library(Rtsne)
# library(dendextend)
library(sva)
library(scran)
library(rgl)
library(ggplot2)
library(cowplot)
library(reshape2)
library(RColorBrewer)
source("../functions.R")
source("class_batch_correction.R")

# theme_set(theme_dark())
theme_set(theme_cowplot())
theme_set(theme_gray())
# FUNCTIONS ---------------------------------------------------------------
# Filter probes with too many zeros
filter_probesets <- function(df, percent_threshold) {
  logical_df <- df != 0
  selected_rows <- rowSums(logical_df) > percent_threshold * ncol(df)
  return(selected_rows)
}

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

# Dataframe has probesets as rows and patients as columns
# Dataframe has patients with D0 and D8 data
select_probes <- function(df) {
  # Row vector is halved into two vectors
  # Elements in the two vectors are assumed to form a pair according to their index!!!
  row_pvalue <- function(row_vec) {
    half_index <- length(row_vec)/2
    # Wilcoxon signed rank takes x-y
    wilcox_obj <- wilcox.test(row_vec[1:half_index],
                              row_vec[-(1:half_index)],
                              paired = T)
    return(wilcox_obj$p.value)
  }
  # Calculating significance of probeset using Wilcoxon signed rank test
  signedrank_pvalue <- apply(df, 1, row_pvalue)
  
  # Mean of D8 - Mean of D0 for every probeset
  effect_size <- rowMeans(df[,-(1:210)]) - rowMeans(df[,1:210])
  print(paste("No.of effect size == 0:", sum(effect_size == 0)))
  effect_direction <- ifelse(effect_size >= 0, "up", "down")

  # Selection of top 500 down-regulated probesets
  probes_df <- data.frame(signedrank_pvalue, effect_direction)
  sorted_probes <- probes_df[order(probes_df$signedrank_pvalue),]
  downreg_probes <- rownames(sorted_probes)[sorted_probes$effect_direction == "down"][1:500]
  return(downreg_probes)
}

# All dataframes have patients as rows and probesets/features as columns
# D0 centroid used to define D0-Normal vector
calc_erm1 <- function(response_df, normal_df, leukemia_df = NA, flag = "replace") {
  num_patients <- nrow(response_df)/2
  # Calculation of centroids and ERM factor
  if (flag == "replace") {
    leukemia_centroid <- apply(response_df[1:num_patients,], 2, median)  
  } else if (flag == "original") {
    # Replace leukemia centroid with D0 centroid
    leukemia_centroid <- apply(leukemia_df, 2, median)
  } else {
    stop("Incorrect flag passed")
  }
  normal_centroid <- apply(normal_df, 2, median)
  
  # Calculation of constant vector factor
  leukemia_normal <- normal_centroid - leukemia_centroid
  unit_leukemia_normal <- leukemia_normal/l2_norm(leukemia_normal)
  
  # Assume that patients from top rows match correspondingly with bottom rows
  # Calculate vector by: D8-D0
  d0d8_hstack <- response_df[-(1:num_patients),] - response_df[1:num_patients,]
  # Multiplication of erm_factor is propagated through every column
  erm <- colSums(t(d0d8_hstack) * unit_leukemia_normal)
  # Vertical stack of individual D0-Normal vectors
  d0normal_vec_vstack <- normal_centroid - t(response_df[1:num_patients,])
  d0normal_scalar_proj <- colSums(d0normal_vec_vstack * unit_leukemia_normal)
  erm_ratio <- erm/d0normal_scalar_proj
  return(cbind(erm, erm_ratio))
}

# All dataframes have patients as rows and probesets/features as columns
calc_erm2 <- function(response_df, normal_df, leukemia_df = NA) {
  num_patients <- nrow(response_df)/2
  normal_centroid <- apply(normal_df, 2, median)
  l2_norm <- function(vec) sqrt(sum(vec^2))
  # D0-Normal vectors for each patient
  d0normal_vec_vstack <- normal_centroid - t(response_df[1:num_patients,])
  l2norm_d0normal <- apply(d0normal_vec_vstack, 2, l2_norm)
  l2norm_arr <- matrix(l2norm_d0normal, 3, length(l2norm_d0normal), byrow = T)
  # Divided by l2norm
  unit_d0normal_vstack <- d0normal_vec_vstack/l2norm_arr

  # Assume that patients from top rows match correspondingly with bottom rows
  # Calculate vector by: D8-D0
  d0d8_hstack <- response_df[-(1:num_patients),] - response_df[1:num_patients,]
  # Element-wise multiplication of two df
  erm <- colSums(t(d0d8_hstack) * unit_d0normal_vstack)
  erm_ratio <- erm/l2norm_d0normal
  return(cbind(erm, erm_ratio))
}

# Calculating ERM distance (PC1)
calc_erm3 <- function(response_df, normal_df) {
  num_patient <- nrow(response_df)/2
  response_pc1 <- response_df[,1]
  normal_pc1 <- normal_df[,1]
  d0d8_pc1 <- response_pc1[-(1:num_patient)] - response_pc1[1:num_patient]
  d0normal_pc1 <- median(normal_pc1) - response_pc1[1:num_patient]
  return(cbind(erm = d0d8_pc1, erm_ratio = d0d8_pc1/d0normal_pc1))
}

# Plots ROC and calculates AUC in a primitive fashion (i.e. ROC is step function)
# Does not resolve ties in the score
# Assumption: Lower score will be labelled preferentially as 1, ROC is step function
# Assumption that score vec and label vecs are corresponding
plot_roc <- function(score_list, label_vec,
                     name_vec = NULL, plot_title = NULL,
                     is_bigger_better_vec = rep(F, length(score_list))) {
  # Function to plot a single ROC curve and calculate AUC
  ROC_AUC <- function(score_vec, is_bigger_better, color) {
    # Sort label vector according to score vector in ascending order
    sort_label <- label_vec[order(score_vec, decreasing = is_bigger_better)]
    # Dataframe of TPR and FPR
    df <- data.frame(TPR = cumsum(sort_label)/sum(sort_label),
                     FPR = cumsum(!sort_label)/sum(!sort_label))
    # Insert 0, 0 at the first row
    roc_df <- rbind(c(0,0), df)
    # Calculates change in FPR at each step
    dFPR <- c(0, diff(roc_df$FPR))
    # Sum of individual rectangle steps
    AUC <- sum(roc_df$TPR * dFPR)
    lines(roc_df$FPR, roc_df$TPR,
          col = color, lwd = 3)
    return(AUC)
  }
  # Initialise plot
  colour_palette <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A",
                      "#66A61E", "#E6AB02", "#A6761D", "#666666")
  color_index <- colour_palette[1:length(score_list)]
  plot(NULL, main = plot_title,
       xlab = "1 - Specificity", ylab = "Sensitivity",
       xlim = c(0,1), ylim = c(0,1),
       xaxs = "i", yaxs = "i",
       cex.main = 1.8, cex.lab = 1.7, cex.axis = 1.6)
  abline(0, 1, lty = 5, lwd = 2)
  auc_vec <- mapply(ROC_AUC, score_list, is_bigger_better_vec,
                    color_index, SIMPLIFY = T)
  # If name_vec is not NULL display legend
  if (!is.null(name_vec)) {
    format_text <- function(name, auc) sprintf("%s (%.3f)", name, auc)
    legend_text <- mapply(format_text, name_vec, auc_vec)
    legend("bottomright", inset = 0.03, lty = 1, lwd = 2,
           cex = 1.4, bty = "o", text.width = 0.45,
           legend = legend_text, col = color_index)
  }
  return(auc_vec)
}

# Environment: yeoh_metadata
plot_pca <- function(df, metadata) {
  pca_obj <- prcomp(t(df))
  pca_df <- as.data.frame(pca_obj$x[,1:6])
  eig_value <- (pca_obj$sdev)^2
  var_pc <- eig_value[1:5]/sum(eig_value)
  pc_labels <- sprintf("PC%d (%.2f%%)", 1:5, var_pc*100)
  # Create plot_df by concatenating metadata labels to data
  plot_df <- cbind(metadata[rownames(pca_df), c(3,6)], pca_df)
  print(plot_df[1:3,])
  pc1_pc2 <- ggplot(plot_df, aes(x = PC1, y = PC2)) +
    geom_point(aes(col = factor(batch), shape = factor(time_point)),
               size = 3, show.legend = F) +
    xlab(pc_labels[1]) + ylab(pc_labels[2])
  pc2_pc3 <- ggplot(plot_df, aes(x = PC2, y = PC3)) +
    geom_point(aes(col = factor(batch), shape = factor(time_point)),
               size = 3, show.legend = F) +
    xlab(pc_labels[2]) + ylab(pc_labels[3])
  pc1_pc3 <- ggplot(plot_df, aes(x = PC1, y = PC3)) +
    geom_point(aes(col = factor(batch), shape = factor(time_point)),
               size = 3, show.legend = F) +
    xlab(pc_labels[1]) + ylab(pc_labels[3])
  multiplot <- plot_grid(pc1_pc2, pc2_pc3, pc1_pc3,
                         ncol = 3, nrow = 1)
  return(multiplot)
}

# 2D PCA plot
# Arguments: PCA-transformed df
pca_all <- function(df, colour_code, shape_vec, pc_labels) {
  pc1_pc2 <- ggplot(df, aes(x = PC1, y = PC2)) +
    geom_point(size = 3, fill = colour_code, colour = "black", shape = shape_vec, show.legend = F) +
    xlab(pc_labels[1]) + ylab(pc_labels[2])
  pc2_pc3 <- ggplot(df, aes(x = PC2, y = PC3)) +
    geom_point(size = 3, fill = colour_code, colour = "black", shape = shape_vec, show.legend = F) +
    xlab(pc_labels[2]) + ylab(pc_labels[3])
  pc1_pc3 <- ggplot(df, aes(x = PC1, y = PC3)) +
    geom_point(size = 3, fill = colour_code, colour = "black", shape = shape_vec, show.legend = F) +
    xlab(pc_labels[1]) + ylab(pc_labels[3])
  pc3_pc4 <- ggplot(df, aes(x = PC3, y = PC4)) +
    geom_point(size = 3, fill = colour_code, colour = "black", shape = shape_vec, show.legend = F) +
    xlab(pc_labels[3]) + ylab(pc_labels[4])
  multiplot <- plot_grid(pc1_pc2, pc2_pc3, pc1_pc3, pc3_pc4,
                         ncol = 2, nrow = 2)
  return(multiplot)
}

# 3D PCA plot
plot.pca_3d <- function(df, colour_code, shape_vec, pc_labels = NULL, ratio_list = list(2,1,1)) {
  if (is.null(pc_labels)) {
    pca_obj <- prcomp(t(df))
    pca_df <- as.data.frame(pca_obj$x[,1:3])
    eig_value <- (pca_obj$sdev)^2
    var_pc <- eig_value[1:3]/sum(eig_value)
    pc_labels <- sprintf("PC%d (%.2f%%)", 1:3, var_pc*100)
  } else {
    pca_df <- as.data.frame(df)
  }
  # RGL plot parameters
  rgl.open()
  rgl.bg(color="white")
  rgl.viewpoint(zoom = 0.8)
  # rgl.viewpoint(theta = 110, phi = 5, zoom = 0.8)
  par3d(windowRect = c(50, 20, 500, 500))
  # Plot of MILE dataset
  with(pca_df, pch3d(PC1, PC2, PC3, bg = colour_code,
                     pch = shape_vec, cex = 0.5, lwd = 1.5))
  box3d(col = "black")
  title3d(xlab = pc_labels[1], ylab = pc_labels[2],
          zlab = pc_labels[3], col = "black")
  # Plot aspect ratios of axis according to variance
  do.call(aspect3d, ratio_list)
}

plot.vectors_3d <- function(vectors_arr, colour_code, shape_vec, pc_labels = NULL, ratio_list = list(2,1,1)) {
  if (is.null(pc_labels)) {
    pca_obj <- prcomp(t(df))
    pca_df <- as.data.frame(pca_obj$x[,1:3])
    eig_value <- (pca_obj$sdev)^2
    var_pc <- eig_value[1:3]/sum(eig_value)
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
plot_vectors <- function(df, centroid_df, pc_labels, batch_colour, subtype_colour) {
  pca_1 <- ggplot(data = df) +
    geom_point(aes(x = PC1_A, y = PC2_A), size = 5, stroke = 2,
               colour = subtype_colour, shape = 21, fill = batch_colour, show.legend = F) +
    geom_point(aes(x = PC1_B, y = PC2_B), size = 5, stroke = 2,
               colour = subtype_colour, shape = 22, fill = batch_colour, show.legend = F) +
    geom_point(data = centroid_df, aes(x = PC1, y = PC2),
               size = 5, shape = 17,
               colour = c("tomato3", "orange", "hotpink", "pink")) +
    geom_segment(aes(x = PC1_A, y = PC2_A, xend = PC1_B, yend = PC2_B, colour = as.factor(labels)),
                 arrow = arrow(length = unit(0.3, "cm")), alpha = 0.8, show.legend = F) +
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
                 arrow = arrow(length = unit(0.3, "cm")), alpha = 0.8, show.legend = F) +
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
yeoh_batch <- read.table("data/GSE67684/processed/metadata_combined-batch.tsv",
                         sep = "\t", header = T, row.names = 1)
yeoh_label <- read.table("data/GSE67684/processed/metadata_combined-label_subtype_edited.tsv",
                         sep = "\t", header = T, row.names = 1)

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

# mile_data <- read.table("data/GSE13204/processed/mas5_ordered.tsv",
#                         sep = "\t", header = T, row.names = 1)
# mile_metadata <- read.table("data/GSE13204/processed/metadata.tsv",
#                             sep = "\t", header = T, row.names = 1)

# SCALE & FILTER ----------------------------------------------------------
yeoh_combined <- cbind(yeoh_d0d8, yeoh_remission, yeoh_normal)
scaled_yeoh <- norm.mean_scaling(yeoh_combined)
# Filtering of probesets
selected_probesets <- filter_probesets(scaled_yeoh, 0.1)
selected_yeoh <- scaled_yeoh[selected_probesets,]
# Log2_transform
log_yeoh <- log2_transform(selected_yeoh)

# VISUALISATION -----------------------------------------------------------
yeoh_combined <- cbind(yeoh_d0d8, yeoh_d33, yeoh_normal)
batch_info <- yeoh_batch[colnames(yeoh_combined), "batch"]

generate_colour <- colorRampPalette(c("lightblue", "darkblue"))
batch_palette <- generate_colour(10)
batch_colour <- batch_palette[batch_info]
shape_2d <- rep(21:24, c(210,210,59,4))

# Plot un-normalised
plot_pca2d <- plot.pca_batch(yeoh_combined, batch_colour, shape_2d)
ggsave("dump/pca_2d-yeoh_unnorm.pdf", plot_pca2d,
       width = 9, height = 9)

plot.pca_3d(yeoh_combined, batch_colour, shape_2d)
rgl.viewpoint(zoom = 0.8)

# Save viewpoint parameters
saved_matrix <- par3d()$userMatrix
saved_zoom <- par3d()$zoom
# Reload viewpoint parameters
view3d(userMatrix = saved_matrix, zoom = saved_zoom)
rgl.postscript("dump/pca_3d-yeoh_unnorm.pdf", "pdf")

scaled_combined_yeoh <- norm.mean_scaling(yeoh_combined)
plot_pca2d <- plot.pca_batch(scaled_combined_yeoh, batch_colour, shape_2d)
plot_pca2d
plot.pca_3d(scaled_combined_yeoh, batch_colour, shape_2d)
rgl.postscript("dump/pca_3d-yeoh_scaled.pdf", "pdf")
mean_plot <- plot_mean(scaled_combined_yeoh, batch_info)
ggsave("dump/mean_yeoh.pdf", mean_plot,
       width = 12, height = 8)
col_mean <- colMeans(scaled_combined_yeoh)
which(col_mean > 1100)

outlier_samples <- c("P198_D0", "P186_D8", "N03")
outlier_index <- which(colnames(scaled_combined_yeoh) %in% outlier_samples)
selected_scaled_yeoh <- scaled_combined_yeoh[,-outlier_index]
shape_2d1 <- rep(21:24, c(209,209,59,3))
batch_info1 <- yeoh_batch[colnames(selected_scaled_yeoh), "batch"]
batch_colour1 <- batch_palette[batch_info1]
plot.pca_3d(selected_scaled_yeoh, batch_colour1, shape_2d1, list(2,1,1))
rgl.postscript("dump/pca_3d-yeoh_scaled_selected.pdf", "pdf")

# Quantile-normalised (old)
quantile_yeoh <- norm.quantile(selected_scaled_yeoh)
plot.pca_3d(quantile_yeoh, batch_colour1, shape_2d1, list(2,1,1))
rgl.postscript("dump/pca_3d-yeoh_quantile.pdf", "pdf")

# GFS
gfs_yeoh <- norm.gfs(selected_scaled_yeoh)
plot.pca_3d(gfs_yeoh, batch_colour1, shape_2d1, list(2,1,1))
rgl.postscript("dump/pca_3d-yeoh_gfs.pdf", "pdf")

# Quantile-normalised (old)
quantile_d0 <- norm.quantile(selected_scaled_yeoh[,1:209])
quantile_d8 <- norm.quantile(selected_scaled_yeoh[,210:418])
quantile_d33 <- norm.quantile(selected_scaled_yeoh[,419:480])
quantile1_yeoh <- cbind(quantile_d0, quantile_d8, quantile_d33)
plot.pca_3d(quantile1_yeoh, batch_colour1, shape_2d1, list(2,1,1))
rgl.postscript("dump/pca_3d-quantile_new.pdf", "pdf")


# QUANTILE (TIMEPOINT) ---------------------------------------------------------
# # Quantile (all time points together)
# quantile_yeoh <- norm.quantile(scaled_yeoh[selected_probesets,])
# quantile_d0 <- quantile_yeoh[,1:208]
# quantile_d8 <- quantile_yeoh[,209:416]
# quantile_d33 <- quantile_yeoh[,417:458]
# quantile_normal <- quantile_yeoh[,459:461]

# # Quantile by timepoint
# quantile_d0 <- norm.quantile(log_yeoh[, 1:208])
# quantile_d8 <- norm.quantile(log_yeoh[, 209:416])
# quantile_d33 <- norm.quantile(log_yeoh[, 417:458])
# quantile_normal <- norm.quantile(log_yeoh[, 459:461])
# quantile_yeoh <- cbind(quantile_d0, quantile_d8, quantile_d33, quantile_normal)
# colnames(log_yeoh)[459:461]
# dim(log_yeoh)

# Quantile by timepoint
quantile_d0 <- norm.quantile(selected_yeoh[, 1:208])
quantile_d8 <- norm.quantile(selected_yeoh[, 209:416])
quantile_d33 <- norm.quantile(selected_yeoh[, 417:458])
quantile_normal <- norm.quantile(selected_yeoh[, 459:461])
quantile_yeoh <- cbind(quantile_d0, quantile_d8, quantile_d33, quantile_normal)
colnames(selected_yeoh)[459:461]
dim(selected_yeoh)

# Plot PCA before selecting features
# Batch information of all the timepoints
batch_info <- yeoh_batch[colnames(quantile_yeoh), "batch"]
generate_colour <- colorRampPalette(c("lightblue", "darkblue"))
batch_palette <- generate_colour(10)
batch_colour <- batch_palette[batch_info]
# Shape of all timepoints
timepoint_shape <- rep(21:24, c(208,208,42,3))
plot.pca_3d(quantile_yeoh, batch_colour, timepoint_shape)
rgl.postscript("dump/pca_3d-quantile_timepoint_log_all_feature.pdf", "pdf")

# Selecting drug responsive genes between D0 and D8
ttest_pvalue <- calc_ttest(cbind(quantile_d0, quantile_d8), 208, is_paired = T)
log_fc <- calc_logfc(quantile_d0, quantile_d8)
pvalue_probesets <- names(ttest_pvalue)[ttest_pvalue <= 0.05]
fc_probesets <- names(log_fc)[log_fc > 1]
intersect_probesets <- fc_probesets[fc_probesets %in% pvalue_probesets]
print(length(intersect_probesets))

log_d0 <- log2_transform(quantile_d0[intersect_probesets,])
log_d8 <- log2_transform(quantile_d8[intersect_probesets,])
log_normal <- log2_transform(quantile_normal[intersect_probesets,])
log_d33 <- log2_transform(quantile_d33[intersect_probesets,])
log_combined <- cbind(log_d0, log_d8, log_d33, log_normal)

# log_d0 <- quantile_d0
# log_d8 <- quantile_d8
# log_normal <- quantile_normal
# log_d33 <- quantile_d33

# # MILE data
# scaled_mile <- norm.mean_scaling(mile_data[selected_probesets, 751:824])
# quantile_mile <- norm.quantile(scaled_mile)
# log_mile <- log2_transform(quantile_mile[intersect_probesets,])
# transposed_df <- t(cbind(log_d0, log_mile))

# PCA
pca_obj <- prcomp(t(log_combined))
# PCA: Eigenvalues
eig_value <- (pca_obj$sdev)^2
var_pc <- eig_value[1:5]/sum(eig_value)
pc_labels <- sprintf("PC%d (%.2f%%)", 1:5, var_pc*100)

# PCA: Coordinates
pca_coord <- pca_obj$x[,1:4]
# Response df and normal df
response_df <- pca_coord[1:416, 1:3]
normal_df <- pca_coord[417:461, 1:3]

# Calculating ERM distance
features_df1 <- calc_erm1(response_df, normal_df)
features_df2 <- calc_erm2(response_df, normal_df)
features_df3 <- calc_erm3(response_df, normal_df)

# # PCA basis: D0 and normal
# transposed_df <- t(cbind(log_d0, log_d33, log_normal))
# pca_obj <- prcomp(transposed_df, center = T, scale. = T)
# # PCA coordinates
# pca_basis <- pca_obj$x[,1:3]
# # PCA eigenvalues
# eig_value <- (pca_obj$sdev)^2
# var_pc <- eig_value[1:5]/sum(eig_value)
# pc_labels <- sprintf("PC%d (%.2f%%)", 1:5, var_pc*100)
# # Projection of D8 data
# add_data <- t(log_d8)
# pca_add <- predict(pca_obj, add_data)[,1:3]
# # Final PCA coordinates arr
# plot_arr <- rbind(pca_basis[1:208,],
#                   pca_add,
#                   pca_basis[-(1:208),])
# 
# # Response df for D0 and D8 samples and D33 and normal samples for centroid
# response_df <- plot_arr[1:416, 1:3]
# normal_df <- plot_arr[417:461, 1:3]
# rownames(plot_arr)[417:490]
# rownames(response_df)
# rownames(normal_df)

# # Calculating ERM distance
# features_df1 <- calc_erm1(response_df, normal_df)
# features_df2 <- calc_erm2(response_df, normal_df)
# features_df3 <- calc_erm3(response_df, normal_df)

# QUANTILE (ALL) ---------------------------------------------------------
# Quantile (all time points together)
quantile_yeoh <- norm.quantile(scaled_yeoh[selected_probesets,])
quantile_d0 <- quantile_yeoh[,1:208]
quantile_d8 <- quantile_yeoh[,209:416]
quantile_d33 <- quantile_yeoh[,417:458]
quantile_normal <- quantile_yeoh[,459:461]

# Plot PCA before selecting features
# Batch information of all the timepoints
batch_info <- yeoh_batch[colnames(quantile_yeoh), "batch"]
generate_colour <- colorRampPalette(c("lightblue", "darkblue"))
batch_palette <- generate_colour(10)
batch_colour <- batch_palette[batch_info]
# Shape of all timepoints
timepoint_shape <- rep(21:24, c(208,208,42,3))
plot.pca_3d(quantile_yeoh, batch_colour, timepoint_shape)
rgl.postscript("dump/pca_3d-quantile_timepoint_all_feature.pdf", "pdf")

# Selecting drug responsive genes between D0 and D8
ttest_pvalue <- calc_ttest(cbind(quantile_d0, quantile_d8), 208, is_paired = T)
log_fc <- calc_logfc(quantile_d0, quantile_d8)
pvalue_probesets <- names(ttest_pvalue)[ttest_pvalue <= 0.05]
fc_probesets <- names(log_fc)[log_fc > 1]
intersect_probesets <- fc_probesets[fc_probesets %in% pvalue_probesets]
print(length(intersect_probesets))

log_d0 <- log2_transform(quantile_d0[intersect_probesets,])
log_d8 <- log2_transform(quantile_d8[intersect_probesets,])
log_normal <- log2_transform(quantile_normal[intersect_probesets,])
log_d33 <- log2_transform(quantile_d33[intersect_probesets,])
log_combined <- cbind(log_d0, log_d8, log_d33, log_normal)

# PCA
pca_obj <- prcomp(t(log_combined))
# PCA: Eigenvalues
eig_value <- (pca_obj$sdev)^2
var_pc <- eig_value[1:5]/sum(eig_value)
pc_labels <- sprintf("PC%d (%.2f%%)", 1:5, var_pc*100)

# PCA: Coordinates
pca_coord <- pca_obj$x[,1:4]
# Response df and normal df
response_df <- pca_coord[1:416, 1:3]
normal_df <- pca_coord[417:461, 1:3]

# Calculating ERM distance
features_df1 <- calc_erm1(response_df, normal_df)
features_df2 <- calc_erm2(response_df, normal_df)
features_df3 <- calc_erm3(response_df, normal_df)

# GFS ---------------------------------------------------------------------
yeoh_combined <- cbind(yeoh_d0d8, yeoh_remission, yeoh_normal)
scaled_yeoh <- norm.mean_scaling(yeoh_combined)
# Filtering of probesets
selected_probesets <- filter_probesets(scaled_yeoh, 0.1)

filtered_yeoh <- yeoh_data[selected_probes,]


# MILE data
filtered_mile <- mile_data[selected_probes,]
# log_mile <- log2_transform(filtered_mile)

### GFS normalisation and selection of probesets
gfs_yeoh <- norm_gfs(filtered_yeoh)

# plot_tsne <- Rtsne(t(gfs_yeoh), dims = 2, perplexity = 50)
# names_index <- colnames(gfs_yeoh)
# batch_info <- yeoh_metadata[names_index,6]
# batch_palette <- brewer.pal(9, "Set1")
# batch_colour <- c(batch_palette[batch_info],
#                   batch_palette[batch_info])
# plot(plot_tsne$Y, col = batch_colour, pch = rep(c(17,19), each = 210))
# pca_gfs_yeoh <- plot_pca(gfs_yeoh)
# ggsave("dump/pca-yeoh_gfs.pdf", pca_gfs_yeoh,
#        width = 12, height = 8)
gfs_mile <- norm_gfs(filtered_mile)

# # Selection of probesets according to paired t-test between D0 and D8
# paired_pvalues <- calc_ttest(gfs_yeoh, 210, is_paired = T)
# top_probesets <- names(sort(paired_pvalues, na.last = NA))[1:500]
# select_gfs_yeoh <- gfs_yeoh[top_probesets,]

# Selection of probesets according to unpaired t-test between normal and D0
d0_normal <- cbind(gfs_yeoh[,1:210], gfs_mile[,751:824])
# Filtering of new combined df
selected_probes1 <- filter_probesets(d0_normal, 0.2)
unpaired_pvalues <- calc_ttest(d0_normal[selected_probes1,], 210, is_paired = F)
top_probesets <- names(sort(unpaired_pvalues, na.last = NA))[1:500]
select_gfs_yeoh <- gfs_yeoh[top_probesets,]
select_gfs_mile <- gfs_mile[top_probesets,]

# # Selection of probesets with most variance
# probeset_var <- apply(gfs_yeoh, 1, var)
# top_probesets <- names(sort(probeset_var, decreasing = T)[1:500])
# select_gfs_yeoh <- gfs_yeoh[top_probesets,]

# plot_select_gfs <- plot_pca(select_gfs_yeoh, yeoh_metadata)
# plot_select_gfs
# ggsave("dump/pca-yeoh_gfs_ttest.pdf", plot_select_gfs,
#        width = 12, height = 4)

# PCA of all data (centroid, training, test)
pca_obj <- prcomp(t(cbind(select_gfs_yeoh, select_gfs_mile)))
# Top 3 PCs without batch effects
pca_arr <- pca_obj$x[,1:3]
response_df <- pca_arr[1:420,]
leukemia_df <- pca_arr[421:1170,]
normal_df <- pca_arr[1171:1244,]

# PCA of D0 and normal data
basis_data <- t(cbind(select_gfs_yeoh[,1:210], select_gfs_mile[,751:824]))
rownames(basis_data)
pca_obj <- prcomp(basis_data)
pca_basis <- pca_obj$x[,1:4]
# Projection of D8 and leukemia data
add_data <- t(cbind(select_gfs_yeoh[,-(1:210)], select_gfs_mile[,-(751:824)]))
pca_add <- predict(pca_obj, add_data)[,1:4]
plot_arr <- rbind(pca_basis[1:210,],
                  pca_add[1:210,],
                  pca_basis[-(1:210),])

pca_arr <- rbind(pca_basis[1:210, c(1,2,4)],
                 pca_add[1:210, c(1,2,4)],
                 pca_basis[-(1:210), c(1,2,4)])

pca_arr <- rbind(pca_basis[1:210, 1:3],
                 pca_add[1:210, 1:3],
                 pca_basis[-(1:210), 1:3])

colnames(pca_arr)
response_df <- pca_arr[1:420,]
normal_df <- pca_arr[421:494,]

# PCA: Eigenvalues
eig_value <- (pca_obj$sdev)^2
var_pc <- eig_value[1:5]/sum(eig_value)
pc_labels <- sprintf("PC%d (%.2f%%)", 1:5, var_pc*100)

erm <- calc_erm(response_df, normal_df)

# MNN ---------------------------------------------------------------------
yeoh_combined <- cbind(yeoh_d0d8, yeoh_remission, yeoh_normal)
scaled_yeoh <- norm.mean_scaling(yeoh_combined)
# Filtering of probesets
selected_probesets <- filter_probesets(scaled_yeoh, 0.1)
selected_yeoh <- scaled_yeoh[selected_probesets,]
# Splitting up yeoh data into batches
batch_info <- yeoh_batch[colnames(selected_yeoh), "batch"]
batch_list <- lapply(1:10, function(i) selected_yeoh[, batch_info == i])
arr_list <- lapply(batch_list, data.matrix)

colnum_list <- sapply(arr_list, ncol)
batch_order <- order(colnum_list, decreasing = T)
print(batch_order)

# Ordered according to number of samples in batch
mnn_yeoh_obj <- mnnCorrect(arr_list[[2]],
                           arr_list[[9]],
                           arr_list[[10]],
                           arr_list[[8]],
                           arr_list[[5]],
                           arr_list[[3]],
                           arr_list[[1]],
                           arr_list[[6]],
                           arr_list[[4]],
                           arr_list[[7]])
mnn_yeoh <- do.call(cbind, mnn_yeoh_obj$corrected)

# Column names for matrix arranged in above order
colnames_vec <- unlist(sapply(arr_list[batch_order], colnames))
colnames(mnn_yeoh) <- colnames_vec
# Order columns according to name and timepoint
ordered_mnn_yeoh <- mnn_yeoh[,order(colnames(mnn_yeoh))]
ordered_mnn_yeoh1 <- cbind(ordered_mnn_yeoh[,endsWith(colnames(ordered_mnn_yeoh), "0")],
                           ordered_mnn_yeoh[,endsWith(colnames(ordered_mnn_yeoh), "8")],
                           ordered_mnn_yeoh[,endsWith(colnames(ordered_mnn_yeoh), "33")],
                           ordered_mnn_yeoh[,startsWith(colnames(ordered_mnn_yeoh), "N")])
# Replace negative values with zero
ordered_mnn_yeoh1[ordered_mnn_yeoh1 < 0] <- 0
# ordered_mnn_yeoh1[1:10,1:10]

mnn_d0 <- ordered_mnn_yeoh1[,1:208]
mnn_d8 <- ordered_mnn_yeoh1[,209:416]
mnn_d33 <- ordered_mnn_yeoh1[,417:458]
mnn_normal <- ordered_mnn_yeoh1[,459:461]

# Plot PCA before selecting features
# Batch information of all the timepoints
batch_info <- yeoh_batch[colnames(ordered_mnn_yeoh1), "batch"]
generate_colour <- colorRampPalette(c("lightblue", "darkblue"))
batch_palette <- generate_colour(10)
batch_colour <- batch_palette[batch_info]
# Shape of all timepoints
timepoint_shape <- rep(21:24, c(208,208,42,3))
plot.pca_3d(ordered_mnn_yeoh1*20, batch_colour, timepoint_shape)
rgl.postscript("dump/pca_3d-mnn_all_feature.pdf", "pdf")

# Selecting drug responsive genes between D0 and D8
ttest_pvalue <- calc_ttest(cbind(mnn_d0, mnn_d8), 208, is_paired = T)
log_fc <- calc_logfc(mnn_d0, mnn_d8)
pvalue_probesets <- names(ttest_pvalue)[ttest_pvalue <= 0.05]
fc_probesets <- names(log_fc)[log_fc > 1]
intersect_probesets <- fc_probesets[fc_probesets %in% pvalue_probesets]
print(length(intersect_probesets))

plot.pca_3d(selected_yeoh, batch_colour, timepoint_shape)
rgl.postscript("dump/pca_3d-before_mnn.pdf", "pdf")

# Filter and log transform
log_mnn_yeoh <- log2_transform(ordered_mnn_yeoh1[intersect_probesets,])

# PCA of D0 and normal data
basis_data <- t(log_mnn_yeoh[, c(1:208, 417:461)])
rownames(basis_data)
pca_obj <- prcomp(basis_data)
# PCA: Coordinates
pca_basis <- pca_obj$x[,1:4]
# PCA: Eigenvalues
eig_value <- (pca_obj$sdev)^2
var_pc <- eig_value[1:5]/sum(eig_value)
pc_labels <- sprintf("PC%d (%.2f%%)", 1:5, var_pc*100)
# Projection of D8 data
add_data <- t(log_mnn_yeoh[, 209:416])
rownames(add_data)
pca_add <- predict(pca_obj, add_data)[,1:4]
# Final PCA coordinates
plot_arr <- rbind(pca_basis[1:208,],
                  pca_add,
                  pca_basis[-(1:208),])
# Response df and normal df
response_df <- plot_arr[1:416, 1:3]
normal_df <- plot_arr[417:461, 1:3]
rownames(normal_df)

# Calculating ERM distance
features_df1 <- calc_erm1(response_df, normal_df)
features_df2 <- calc_erm2(response_df, normal_df)
features_df3 <- calc_erm3(response_df, normal_df)

# MANUAL ------------------------------------------------------------------
# Batch information of all the timepoints
batch_info <- yeoh_batch[colnames(scaled_yeoh), "batch"]
generate_colour <- colorRampPalette(c("lightblue", "darkblue"))
batch_palette <- generate_colour(10)
batch_colour <- batch_palette[batch_info]
# Shape of all timepoints
timepoint_shape <- rep(21:24, c(208,208,42,3))
plot.pca_3d(selected_yeoh, batch_colour, timepoint_shape)
plot.pca_3d(log2_transform(selected_yeoh), batch_colour, timepoint_shape)

timepoint_info <- rep(c("D0","D8","Normal"), c(208,208,45))

table(batch_info, timepoint_info)
order_batch <- c(10,9,8,7,6,4,3,1,2,5)

log_yeoh <- log2_transform(selected_yeoh)
cbc_yeoh <- norm.CBC(log_yeoh, batch_info, timepoint_info, order_batch)
plot.pca_3d(cbc_yeoh, batch_colour, timepoint_shape)
# plot.pca_3d(log2_transform(cbc_yeoh), batch_colour, timepoint_shape)
rgl.postscript("dump/pca_3d-cbc_log_all", "pdf")

colnames(cbc_yeoh)
# Selecting drug responsive genes between D0 and D8
ttest_pvalue <- calc_ttest(cbc_yeoh[,1:416], 208, is_paired = T)
log_fc <- calc_logfc(2^cbc_yeoh[,1:208], 2^cbc_yeoh[,209:416])
pvalue_probesets <- names(ttest_pvalue)[ttest_pvalue <= 0.05]
fc_probesets <- names(log_fc)[log_fc > 1]
intersect_probesets <- fc_probesets[fc_probesets %in% pvalue_probesets]
print(length(intersect_probesets))

# Selecting drug responsive genes between D0 and D33
ttest_pvalue <- calc_ttest(cbind(cbc_yeoh[,1:208], cbc_yeoh[,417:461]),
                           size_a = 208, is_paired = F)
log_fc <- calc_logfc(2^cbc_yeoh[,1:208], 2^cbc_yeoh[,417:461])
pvalue_probesets <- names(ttest_pvalue)[ttest_pvalue <= 0.05]
fc_probesets <- names(log_fc)[log_fc > 1]
intersect_probesets <- fc_probesets[fc_probesets %in% pvalue_probesets]
print(length(intersect_probesets))

selected_cbc_yeoh <- cbc_yeoh[intersect_probesets,]

# PCA
pca_obj <- prcomp(t(selected_cbc_yeoh))
# PCA: Eigenvalues
eig_value <- (pca_obj$sdev)^2
var_pc <- eig_value[1:5]/sum(eig_value)
pc_labels <- sprintf("PC%d (%.2f%%)", 1:5, var_pc*100)

# PCA: Coordinates
pca_coord <- pca_obj$x[,1:4]
# Response df and normal df
response_df <- pca_coord[1:416, 1:3]
normal_df <- pca_coord[417:461, 1:3]

# Calculating ERM distance
features_df1 <- calc_erm1(response_df, normal_df)
features_df2 <- calc_erm2(response_df, normal_df)
features_df3 <- calc_erm3(response_df, normal_df)

# ComBat ------------------------------------------------------------------
yeoh_combined <- cbind(yeoh_d0d8, yeoh_remission, yeoh_normal)
# Mean scaling of data
scaled_yeoh <- norm.mean_scaling(yeoh_combined)
# Filtering of probesets
selected_probesets <- filter_probesets(scaled_yeoh, 0.1)
selected_yeoh <- scaled_yeoh[selected_probesets,]

# Obtaining batch information of selected_yeoh df
# Rows of metadata are to be in same order as columns of edata
batch_info <- yeoh_batch[colnames(selected_yeoh), "batch"]
yeoh_metadata <- data.frame(batch_info)
# Place adjustment/confounding variables in model.matrix (e.g. age)
# Do not put batch variables in model.matrix
# Put batch variables directly in combat function
model_combat <- model.matrix(~1, data = yeoh_metadata)
combat_yeoh <- ComBat(data.matrix(selected_yeoh), batch_info, model_combat)
# Replacing negative values with 0
combat_yeoh[combat_yeoh < 0] <- 0

# Plot PCA before selecting features
# Batch information of all the timepoints
generate_colour <- colorRampPalette(c("lightblue", "darkblue"))
batch_palette <- generate_colour(10)
batch_colour <- batch_palette[batch_info]
# Shape of all timepoints
timepoint_shape <- rep(21:24, c(208,208,42,3))

plot.pca_3d(combat_yeoh, batch_colour, timepoint_shape)
rgl.postscript("dump/pca_3d-combat_all_feature.pdf", "pdf")
plot.pca_3d(log2_transform(combat_yeoh), batch_colour, timepoint_shape)

# Selecting drug responsive genes between D0 and D8
combat_d0 <- combat_yeoh[,1:208]
combat_d8 <- combat_yeoh[,209:416]
ttest_pvalue <- calc_ttest(cbind(combat_d0, combat_d8), 208, is_paired = T)
log_fc <- calc_logfc(combat_d0, combat_d8)
pvalue_probesets <- names(ttest_pvalue)[ttest_pvalue <= 0.05]
fc_probesets <- names(log_fc)[log_fc > 1]
intersect_probesets <- fc_probesets[fc_probesets %in% pvalue_probesets]
print(length(intersect_probesets))

# Filter and log transform
log_combat_yeoh <- log2_transform(combat_yeoh[intersect_probesets,])

# PCA of D0 and normal data
basis_data <- t(log_combat_yeoh[, c(1:208, 417:461)])
rownames(basis_data)
pca_obj <- prcomp(basis_data)
# PCA: Coordinates
pca_basis <- pca_obj$x[,1:4]
# PCA: Eigenvalues
eig_value <- (pca_obj$sdev)^2
var_pc <- eig_value[1:5]/sum(eig_value)
pc_labels <- sprintf("PC%d (%.2f%%)", 1:5, var_pc*100)
# Projection of D8 data
# CHANGE
add_data <- t(log_combat_yeoh[, 209:416])
rownames(add_data)
pca_add <- predict(pca_obj, add_data)[,1:4]
# Final PCA coordinates
plot_arr <- rbind(pca_basis[1:208,],
                  pca_add,
                  pca_basis[-(1:208),])

# Response df and normal df
response_df <- plot_arr[1:416, 1:3]
normal_df <- plot_arr[417:461, 1:3]
rownames(normal_df)

# Calculating ERM distance
features_df1 <- calc_erm1(response_df, normal_df)
features_df2 <- calc_erm2(response_df, normal_df)
features_df3 <- calc_erm3(response_df, normal_df)

# Hierachical clustering --------------------------------------------------
pairwise_dist <- dist(t(ordered_mnn_yeoh1[,1:210]))
hcluster <- hclust(pairwise_dist)
dendo_obj <- as.dendrogram(hcluster)
# sample_id <- labels(dendo_obj)
# nodePar <- list(lab.cex = 0.3, pch = c(NA, NA),
#                 cex = 0.5, col = "blue")

colnames_data <- colnames(ordered_mnn_yeoh1)[1:210]
batch_info <- yeoh_metadata[colnames_data, 6]
batch_palette <- brewer.pal(9, "Set1")
batch_colour <- batch_palette[batch_info]

# Settings of dendogram
dendo_obj <- set(dendo_obj, "labels_cex", 0.4)
plot(dendo_obj, horiz = F)
colored_bars(batch_colour, dendo_obj,
             rowLabels = "Batch ", y_shift = -0.15,
             sort_by_labels_order = F)
dendogram <- recordPlot()
save_fig(dendogram, "dump/hclust-mnn_yeoh.pdf",
         width = 12, height = 6)

# PCA PLOT ----------------------------------------------------------------
# Plot PCA-3D
plot.pca_3d(pca_coord, batch_colour, timepoint_shape, pc_labels)
rgl.postscript("dump/pca_3d-cdc_log_select.pdf", "pdf")

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
ggsave("dump/vectors-cbc_log_select.pdf", vectors_plot,
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

# # Batch information of yeoh is encoded
# names_index <- rownames(plot_arr)[751:(751+209)]
# batch_info <- yeoh_metadata[names_index,6]
# blue_palette <- brewer.pal(9, "Blues")
# batch_colour <- c(subtype_colour,
#                   rep(blue_palette[batch_info], 2),
#                   rep("darkolivegreen3", 74))

# SAVE RESULTS -----------------------------------------------------------------
### Extracting truth labels and training/test
features_df <- cbind(features_df1, features_df2, features_df3)
colnames(features_df) <- c("erm1", "erm1_ratio", "erm2", "erm2_ratio", "erm3", "erm3_ratio")
row_index <- substring(rownames(features_df),1,4)
row_index
labels_yeoh <- yeoh_label[row_index, 5:6]
results_df <- cbind(features_df, labels_yeoh)
rownames(results_df) <- row_index
# results_df[order(results_df$erm),]
write.table(results_df, "dump/results-quantile_timepoint_log_allsamples.tsv",
            quote = F, sep = "\t", row.names = T, col.names = T)

# Plot ROC
# Load results
# results_df <- read.table("dump/results-quantile_timepoint_d33_projected_d8.tsv",
#                          sep = "\t")
# labels_vec <- results_df[, 7]
# head(results_df[,1:6])

# Visualise current results now
head(results_df)
labels_vec <- results_df[, 7]
labels_vec
par(mar = rep(5,4))
line_labels <- c("ERM1", "ERM1-Ratio",
                 "ERM2", "ERM2-Ratio",
                 "PC1", "PC1-Ratio")
plot_roc(results_df[,1:6], labels_vec,
         name_vec = line_labels)

results_roc <- recordPlot()
save_fig(results_roc, "dump/roc-quantile_timepoint_log_allsamples.pdf",
         width = 9, height = 9)

# Visualise all except subtype: Others
head(filtered_subtype_results)
labels_vec <- results_df[, 7]
labels_vec
par(mar = rep(5,4))
line_labels <- c("ERM1", "ERM1-Ratio",
                 "ERM2", "ERM2-Ratio",
                 "PC1", "PC1-Ratio", "MRD33")
plot_roc(filtered_subtype_results[,c(1:6,9)], filtered_subtype_results[,7],
         name_vec = line_labels, is_bigger_better_vec = c(rep(F,6),T))

results_roc <- recordPlot()
save_fig(results_roc, "dump/roc-quantile_timepoint_no_others.pdf",
         width = 9, height = 9)

View(filtered_subtype_results)

# Subtype analysis
head(results_df)
head(yeoh_label)
analysis_df <- cbind(results_df[,c(1,7)],
                     yeoh_label[rownames(results_df),"subtype", drop = F])
head(analysis_df[order(analysis_df$erm1),])

plot_subtype <- ggplot(analysis_df) +
  geom_jitter(aes(x = subtype, y = erm1, color = as.factor(label)),
              position = position_jitter(0.1),
              show.legend = F) +
  scale_color_manual(values = c("darkolivegreen3", "tomato3"))
plot_subtype
ggsave("dump/subtype_analysis.pdf", plot_subtype,
       width = 8, height = 5)

# Join results_df and yeoh_label
merged_results <- merge(results_df, yeoh_label[,c(1,7)], by = "row.names")
merged_results <- merged_results[,-c(1)]
# Filter out MRD33 missing values
filtered_results<- merged_results[!is.na(merged_results$d33_mrd),]
# Convert to numeric
filtered_results$d33_mrd <- as.numeric(as.character(filtered_results$d33_mrd))
# Change NA values (<0.0001) to 5e-5
filtered_results$d33_mrd[is.na(filtered_results$d33_mrd)] <- 5e-5
# Filter out subtype: Others
filtered_subtype_results <- filtered_results[filtered_results$subtype != "Others",]

write.table(filtered_subtype_results, "dump/results_quantile_timepoint_mrd33_no_others.tsv",
            quote = F, sep = "\t")


# MILE: Normal samples ----------------------------------------------------
# PCA
plot.pca_batch <- function(untransformed_df, colour_vec, shape_vec = 19) {
  colour_vec <- as.factor(colour_vec)
  pca_obj <- prcomp(t(untransformed_df))
  df <- data.frame(pca_obj$x[,1:5])
  eig_value <- (pca_obj$sdev)^2
  var_pc <- eig_value[1:5]/sum(eig_value)
  pc_labels <- sprintf("PC%d (%.2f%%)", 1:5, var_pc*100)
  # Scatter plots
  # pc1_pc2 <- ggplot(df, aes(x = PC1, y = PC2)) +
  #   geom_point(size = 2, col = colour_vec, shape = shape_vec, show.legend = F) +
  #   xlab(pc_labels[1]) + ylab(pc_labels[2])
  # pc2_pc3 <- ggplot(df, aes(x = PC2, y = PC3)) +
  #   geom_point(size = 2, col = colour_vec, shape = shape_vec, show.legend = F) +
  #   xlab(pc_labels[2]) + ylab(pc_labels[3])
  # pc1_pc3 <- ggplot(df, aes(x = PC1, y = PC3)) +
  #   geom_point(size = 2, col = colour_vec, shape = shape_vec, show.legend = F) +
  #   xlab(pc_labels[1]) + ylab(pc_labels[3])
  # pc3_pc4 <- ggplot(df, aes(x = PC3, y = PC4)) +
  #   geom_point(size = 2, col = colour_vec, shape = shape_vec, show.legend = F) +
  #   xlab(pc_labels[3]) + ylab(pc_labels[4])
  
  pc1_pc2 <- ggplot(df, aes(x = PC1, y = PC2)) +
    geom_point(size = 3, fill = colour_vec, shape = shape_vec, show.legend = F) +
    xlab(pc_labels[1]) + ylab(pc_labels[2])
  pc2_pc3 <- ggplot(df, aes(x = PC2, y = PC3)) +
    geom_point(size = 3, fill = colour_vec, shape = shape_vec, show.legend = F) +
    xlab(pc_labels[2]) + ylab(pc_labels[3])
  pc1_pc3 <- ggplot(df, aes(x = PC1, y = PC3)) +
    geom_point(size = 3, fill = colour_vec, shape = shape_vec, show.legend = F) +
    xlab(pc_labels[1]) + ylab(pc_labels[3])
  pc3_pc4 <- ggplot(df, aes(x = PC3, y = PC4)) +
    geom_point(size = 3, fill = colour_vec, shape = shape_vec, show.legend = F) +
    xlab(pc_labels[3]) + ylab(pc_labels[4])
  
  multiplot <- plot_grid(pc1_pc2, pc2_pc3, pc1_pc3, pc3_pc4,
                         ncol = 2, nrow = 2)
  return(multiplot)
}

normal_mile <- mile_data[,751:824]
scaled_mile <- norm.mean_scaling(normal_mile)
leuk_yeoh <- yeoh_data[,1:210]

colour_info <- rep(c("black", "red"), c(74, 210))
plot.pca_batch(cbind(normal_mile, leuk_yeoh), colour_info)
plot.pca_3d(cbind(normal_mile, leuk_yeoh), colour_info)

# QUANTILE (MISC) ---------------------------------------------------------------
# Selecting genes based on subtypes in MILE data
leukemia_subtype <- substring(colnames(qnorm_leukemia),1,1)
union_probesets <- character()
for (subtype in LETTERS[1:8]){
  print(subtype)
  logfc <- calc_logfc(qnorm_normal, qnorm_leukemia[,leukemia_subtype == subtype])
  probesets <- rownames(qnorm_normal)[logfc > 2]
  union_probesets <- union(union_probesets, probesets)
}

# # Selecting drug responsive genes between D0 and D8 (Quantile: ALL)
# ttest_pvalue <- calc_ttest(qnorm_yeoh, 210, is_paired = T)
# log_fc <- calc_logfc(qnorm_yeoh[,1:210], qnorm_yeoh[,-(1:210)])
# pvalue_probesets <- names(ttest_pvalue)[ttest_pvalue <= 0.05]
# fc_probesets <- names(log_fc)[log_fc > 1]
# intersect_probesets <- fc_probesets[fc_probesets %in% pvalue_probesets]
# print(length(intersect_probesets))

log_d0 <- log2_transform(qnorm_d0[intersect_probesets,])
log_d8 <- log2_transform(qnorm_d8[intersect_probesets,])
log_normal <- log2_transform(qnorm_normal[intersect_probesets,])
log_leukemia <- log2_transform(qnorm_leukemia[intersect_probesets,])
dim(log_leukemia)

