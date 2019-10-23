# library(Rtsne)
# library(dendextend)
library(gPCA)
library(xtable)
library(cluster)
library(sva)
library(scran)
library(Harman)
library(rgl)
library(ggplot2)
library(cowplot)
library(reshape2)
library(RColorBrewer)
source("../functions.R")
source("class_batch_correction.R")

# theme_set(theme_dark())
theme_set(theme_cowplot())
# theme_set(theme_gray())

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
  print(tail(names(erm)))
  print(tail(names(erm_ratio)))
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
  print(tail(names(erm)))
  print(tail(names(erm_ratio)))
  return(cbind(erm, erm_ratio))
}

# Calculating ERM distance (PC1)
calc_erm3 <- function(response_df, normal_df) {
  num_patient <- nrow(response_df)/2
  response_pc1 <- response_df[,1]
  normal_pc1 <- normal_df[,1]
  d0d8_pc1 <- response_pc1[-(1:num_patient)] - response_pc1[1:num_patient]
  d0normal_pc1 <- median(normal_pc1) - response_pc1[1:num_patient]
  print(tail(names(d0d8_pc1)))
  print(tail(names(d0d8_pc1/d0normal_pc1)))
  return(cbind(erm = d0d8_pc1, erm_ratio = d0d8_pc1/d0normal_pc1))
}


# Calculates l2norm of D0D8 vector
# Calculates l2norm(D8) - l2norm(D0)
calc_erm4 <- function(response_df) {
  num_patient <- nrow(response_df)/2
  d0_df <- response_df[1:num_patient,]
  d8_df <- response_df[-(1:num_patient),]
  
  stopifnot(nrow(d8_df) == nrow(d0_df))
  d0d8_df <- d8_df - d0_df
  d0d8_l2norm <- apply(d0d8_df, 1, calc.l2_norm)
  
  d0_l2norm <- apply(d0_df, 1, calc.l2_norm)
  d8_l2norm <- apply(d8_df, 1, calc.l2_norm)
  diff_l2norm <- d8_l2norm - d0_l2norm
  
  # Angle between D0 and D8
  cos_similarity <- mapply(calc.cos_sim, data.frame(t(d0_df)), data.frame(t(d8_df)))
  angle_d0d8 <- calc.rad2degree(acos(cos_similarity))
  print(tail(names(d0d8_l2norm)))
  print(tail(names(diff_l2norm)))
  print(tail(names(angle_d0d8)))
  return(data.frame(d0d8_l2norm, diff_l2norm, angle_d0d8))
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
           cex = 1.2, bty = "o", text.width = 0.35,
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
    print("PCA performed!")
    pca_obj <- prcomp(t(df), center = T, scale. = F)
    pca_df <- as.data.frame(pca_obj$x[,1:3])
    eig_value <- (pca_obj$sdev)^2
    var_pc <- eig_value[1:3]/sum(eig_value)
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

# Metadata
batch_info <- yeoh_batch[colnames(yeoh_combined), "batch"]
class_info <- rep(c("D0","D8","N"), c(208,208,45))
class_numeric <- rep(1:3, c(208,208,45))

scaled_yeoh <- norm.mean_scaling(yeoh_combined)
# Filtering of probesets
selected_probesets <- filter_probesets(scaled_yeoh, 0.3, class_info)

filtered_yeoh <- scaled_yeoh[selected_probesets,]
# Log2_transform
log_yeoh <- log2_transform(filtered_yeoh)

# # Plot PCA before selecting features
# # Batch information of all the timepoints
# batch_info <- yeoh_batch[colnames(log_yeoh), "batch"]
# generate_colour <- colorRampPalette(c("lightblue", "darkblue"))
# batch_palette <- generate_colour(10)
# batch_colour <- batch_palette[batch_info]
# # Shape of all timepoints
# timepoint_shape <- rep(21:23, c(208,208,45))
# 
# plot.pca_3d(log_yeoh, batch_colour, timepoint_shape)
# rgl.postscript("dump/fig/pca_3d-uncorrected_data.pdf", "pdf")
# 
# # Original batch effects metric
# eval.batch_effects(log_yeoh, batch_info, class_numeric)

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

# QUANTILE ---------------------------------------------------------
# ### QUANTILE (ALL)
# quantile_yeoh <- norm.quantile(log_yeoh)
# quantile_d0 <- quantile_yeoh[,1:208]
# quantile_d8 <- quantile_yeoh[,209:416]
# quantile_normal <- quantile_yeoh[,417:461]

### QUANTILE (TIMEPOINT)
quantile_d0 <- norm.quantile(log_yeoh[, 1:208])
quantile_d8 <- norm.quantile(log_yeoh[, 209:416])
quantile_normal <- norm.quantile(log_yeoh[, 417:461])
quantile_yeoh <- cbind(quantile_d0, quantile_d8, quantile_normal)

calc.var_preservation(log_yeoh, quantile_yeoh)
metrics <- eval.batch_effects(quantile_yeoh, batch_info, class_numeric)

# Plot PCA before selecting features
# Batch information of all the timepoints
batch_info <- yeoh_batch[colnames(quantile_yeoh), "batch"]
generate_colour <- colorRampPalette(c("lightblue", "darkblue"))
batch_palette <- generate_colour(10)
batch_colour <- batch_palette[batch_info]


# Shape of all timepoints
timepoint_shape <- rep(21:23, c(208,208,45))
plot.pca_3d(quantile_yeoh, batch_colour, timepoint_shape)
rgl.postscript("dump/pca_3d-cs_quantile.pdf", "pdf")

# Selecting drug responsive genes between D0 and D8
ttest_pvalue <- calc_ttest(cbind(quantile_d0, quantile_d8), 208, is_paired = T)
log_fc <- rowMeans(quantile_d8) - rowMeans(quantile_d0)
pvalue_probesets <- names(ttest_pvalue)[ttest_pvalue <= 0.05]
fc_probesets <- names(log_fc)[log_fc > 1]
intersect_probesets <- fc_probesets[fc_probesets %in% pvalue_probesets]
print(length(intersect_probesets))

corrected_df <- quantile_yeoh[intersect_probesets,]

# # MILE data
# scaled_mile <- norm.mean_scaling(mile_data[selected_probesets, 751:824])
# quantile_mile <- norm.quantile(scaled_mile)
# log_mile <- log2_transform(quantile_mile[intersect_probesets,])
# transposed_df <- t(cbind(log_d0, log_mile))
# QNORM - PROTOTYPE -------------------------------------------------------
### QUANTILE (TIMEPOINT)
quantile_d0 <- norm.quantile(subset_yeoh[, 1:207])
quantile_d8 <- norm.quantile(subset_yeoh[, 208:414])
quantile_normal <- norm.quantile(subset_yeoh[, 415:417])
quantile_yeoh <- cbind(quantile_d0, quantile_d8, quantile_normal)
colnames(subset_yeoh)[415:417]

calc.var_composition(quantile_yeoh, batch_info, class_info)
calc.var_preservation(log_yeoh, quantile_yeoh)

# Plot PCA before selecting features
# Batch information of all the timepoints
batch_info <- yeoh_batch[colnames(quantile_yeoh), "batch"]
generate_colour <- colorRampPalette(c("lightblue", "darkblue"))
batch_palette <- generate_colour(10)
batch_colour <- batch_palette[batch_info]
# Shape of all timepoints
timepoint_shape <- rep(21:23, c(207,207,3))
plot.pca_3d(quantile_yeoh, batch_colour, timepoint_shape)
rgl.postscript("dump/pca_3d-quantile_timepoint_log_all_feature.pdf", "pdf")

# Selecting drug responsive genes between D0 and D8
ttest_pvalue <- calc_ttest(cbind(quantile_d0, quantile_d8), 207, is_paired = T)
log_fc <- rowMeans(quantile_d8) - rowMeans(quantile_d0)
pvalue_probesets <- names(ttest_pvalue)[ttest_pvalue <= 0.05]
fc_probesets <- names(log_fc)[log_fc > 1]
intersect_probesets <- fc_probesets[fc_probesets %in% pvalue_probesets]
print(length(intersect_probesets))

selected_quantile_yeoh <- quantile_yeoh[intersect_probesets,]

corrected_df <- selected_quantile_yeoh
# PCA
pca_obj <- prcomp(t(corrected_df))
# PCA: Eigenvalues
eig_value <- (pca_obj$sdev)^2
var_pc <- eig_value[1:5]/sum(eig_value)
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

# Quantile (Class) - ComBat -----------------------------------------------
colnames(log_yeoh)
# Quantile by timepoint
quantile_d0 <- norm.quantile(log_yeoh[, 1:208])
quantile_d8 <- norm.quantile(log_yeoh[, 209:416])
quantile_normal <- norm.quantile(log_yeoh[, 417:461])
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
# Batch information of all the timepoints
generate_colour <- colorRampPalette(c("lightblue", "darkblue"))
batch_palette <- generate_colour(10)
batch_colour <- batch_palette[batch_info]
# Shape of all timepoints
timepoint_shape <- rep(21:23, c(208,208,45))

plot.pca_3d(combat_yeoh, batch_colour, timepoint_shape)
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
# Batch information of all the timepoints
batch_info <- yeoh_batch[colnames(ordered_mnn_yeoh), "batch"]
generate_colour <- colorRampPalette(c("lightblue", "darkblue"))
batch_palette <- generate_colour(10)
batch_colour <- batch_palette[batch_info]
# Shape of all timepoints
timepoint_shape <- rep(21:23, c(208,208,45))
plot.pca_3d(ordered_mnn_yeoh, batch_colour, timepoint_shape)
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


quantile_d0 <- norm.quantile(log_yeoh[, 1:208])
quantile_d8 <- norm.quantile(log_yeoh[, 209:416])
quantile_normal <- norm.quantile(log_yeoh[, 417:461])
quantile_yeoh <- cbind(quantile_d0, quantile_d8, quantile_normal)
cbc_yeoh <- norm.BCM(quantile_yeoh, batch_info, timepoint_info, 2)

# Batch information of all the timepoints
batch_info <- yeoh_batch[colnames(cbc_yeoh), "batch"]
generate_colour <- colorRampPalette(c("lightblue", "darkblue"))
batch_palette <- generate_colour(10)
batch_colour <- batch_palette[batch_info]
# Shape of all timepoints
timepoint_shape <- rep(21:23, c(208,208,45))

plot.pca_3d(cbc_yeoh, batch_colour, timepoint_shape)
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
eig_value <- (pca_obj$sdev)^2
var_pc <- eig_value[1:5]/sum(eig_value)
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

cbc_kantan_df <- norm.CBC(kantan_df, kantan_batch, kantan_class, 1:2, "dump/diff_kantan_vectors.tsv")
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
# Batch information of all the timepoints
generate_colour <- colorRampPalette(c("lightblue", "darkblue"))
batch_palette <- generate_colour(10)
batch_colour <- batch_palette[batch_info]
# Shape of all timepoints
timepoint_shape <- rep(21:23, c(208,208,45))

plot.pca_3d(combat_yeoh, batch_colour, timepoint_shape)
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
selected_probesets1 <- filter_probesets(scaled_yeoh, 0.3, class_info, logical_func = all)
log_yeoh1 <- log_yeoh[selected_probesets1,]
# Quantile by timepoint
quantile_d0 <- norm.quantile(log_yeoh1[, 1:208])
quantile_d8 <- norm.quantile(log_yeoh1[, 209:416])
quantile_normal <- norm.quantile(log_yeoh1[, 417:461])
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
# Batch information of all the timepoints
generate_colour <- colorRampPalette(c("lightblue", "darkblue"))
batch_palette <- generate_colour(10)
batch_colour <- batch_palette[batch_info]
# Shape of all timepoints
timepoint_shape <- rep(21:23, c(208,208,45))

plot.pca_3d(cs_combat_yeoh, batch_colour, timepoint_shape)
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
batch_info <- yeoh_batch[colnames(harman_yeoh), "batch"]
generate_colour <- colorRampPalette(c("lightblue", "darkblue"))
batch_palette <- generate_colour(10)
batch_colour <- batch_palette[batch_info]
# Shape of all timepoints
timepoint_shape <- rep(21:23, c(208,208,45))

plot.pca_3d(harman_yeoh, batch_colour, timepoint_shape)
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
batch_info <- yeoh_batch[colnames(ordered_scanorama), "batch"]
generate_colour <- colorRampPalette(c("lightblue", "darkblue"))
batch_palette <- generate_colour(10)
batch_colour <- batch_palette[batch_info]
# Shape of all timepoints
timepoint_shape <- rep(21:23, c(208,208,45))
plot.pca_3d(ordered_scanorama*100, batch_colour, timepoint_shape)
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

# ERM CALCULATION ---------------------------------------------------------
# corrected_df <- selected_quantile_yeoh
# PCA
pca_obj <- prcomp(t(corrected_df))
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
features_df4 <- calc_erm4(response_df)

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

# SAVE RESULTS -----------------------------------------------------------------
### Extracting truth labels and training/test
collate_results <- function(features_df1, features_df2, features_df3, features_df4, yeoh_label) {
  features_df <- cbind(features_df1, features_df2, features_df3, features_df4)
  colnames(features_df) <- c("erm1", "erm1_ratio",
                             "erm2", "erm2_ratio",
                             "erm3", "erm3_ratio",
                             "d0d8_l2norm", "diff_l2norm","angle_d0d8")
  row_index <- substring(rownames(features_df),1,4)
  labels_yeoh <- yeoh_label[row_index, 5:6]
  results_df <- cbind(features_df, labels_yeoh)
  rownames(results_df) <- row_index
  return(results_df)
}

results_df <- collate_results(features_df1, features_df2, features_df3, features_df4, yeoh_label)
head(results_df)

write.table(results_df, "dump/results-mnn_k3.tsv",
            quote = F, sep = "\t", row.names = T, col.names = T)



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
plot_roc(results_df[,1:6], labels_vec,
         name_vec = line_labels)
results_roc <- recordPlot()

save_fig(results_roc, "dump/roc-quantile_cs_combat.pdf",
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

# # Join results_df and yeoh_label
# merged_results <- merge(results_df, yeoh_label[,c(1,7)], by = "row.names")
# merged_results <- merged_results[,-c(1)]
# # Filter out MRD33 missing values
# filtered_results<- merged_results[!is.na(merged_results$d33_mrd),]
# # Convert to numeric
# filtered_results$d33_mrd <- as.numeric(as.character(filtered_results$d33_mrd))
# # Change NA values (<0.0001) to 5e-5
# filtered_results$d33_mrd[is.na(filtered_results$d33_mrd)] <- 5e-5
# # Filter out subtype: Others
# filtered_subtype_results <- filtered_results[filtered_results$subtype != "Others",]
# 
# write.table(filtered_subtype_results, "dump/results_quantile_timepoint_mrd33_no_others.tsv",
#             quote = F, sep = "\t")

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

# # Batch information of yeoh is encoded
# names_index <- rownames(plot_arr)[751:(751+209)]
# batch_info <- yeoh_metadata[names_index,6]
# blue_palette <- brewer.pal(9, "Blues")
# batch_colour <- c(subtype_colour,
#                   rep(blue_palette[batch_info], 2),
#                   rep("darkolivegreen3", 74))

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

# QPSP --------------------------------------------------------------------
library(NetProt)
library(genefilter)

# Import CORUM df
raw_corum <- read.table("../info/CORUM/entrezId.txt",
                        sep = "\t", header = T, row.names = 1, stringsAsFactors = F)
# Only human complexes
human_corum <- raw_corum[raw_corum$Organism == "Human",]
list_corum <- strsplit(human_corum[,2], ';')
names(list_corum) <- rownames(human_corum)
head(list_corum)

# Import NEA
nea_df <- read.table("../diff_expr/data/subnetwork/nea-hsa/ovarian_cancer/geneset-nea_kegg_ovarian.tsv",
                     sep = "\t", header = T, stringsAsFactors = F)
subnetwork_nea <- split(as.character(nea_df$gene_id), nea_df$subnetwork_id)

#' Removes ambiguous and AFFY probesets from dataframe
#' Rowname of affymetrix probesets
remove_probesets <- function(df) {
  logical_vec <- grepl("[0-9]_at", rownames(df)) & !startsWith(rownames(df), "AFFX")
  print(paste0("No. of ambiguous and AFFY probesets removed: ",
               nrow(df) - sum(logical_vec)))
  return(df[logical_vec, , drop=F])
}

processed_yeoh <- remove_probesets(log_yeoh)

# Map probesets to IDs
# Removes ambiguous probesets and probesets with no ID
# Selects maximum if two probesets match to same gene
# CHECK: What microarray platform is the data from?
ANNOT_PROBESET_RPATH <- "../info/microarray/HG-U133A/annot_entrez-GPL96.tsv"
entrez_yeoh <- affy2id(processed_yeoh, ANNOT_PROBESET_RPATH)

gfs_yeoh <- norm.gfs(entrez_yeoh, upper = 0.1, lower = 0.2, num_intervals = 4)

# QPSP
calc.qpsp <- function (rank_weight_matrix, complex_list) {
  qpsp_matrix <- c()
  for (j in 1:length(complex_list)) {
    if (length(rownames(rank_weight_matrix)[which(rownames(rank_weight_matrix) %in% 
                                                  complex_list[[j]])]) > 1) {
      qpsp_matrix <- rbind(qpsp_matrix, colSums(rank_weight_matrix[rownames(rank_weight_matrix)[which(rownames(rank_weight_matrix) %in% 
                                                                                                        complex_list[[j]])], ]))
    }
    else if (length(rownames(rank_weight_matrix)[which(rownames(rank_weight_matrix) %in% 
                                                       complex_list[[j]])]) == 1) {
      qpsp_matrix <- rbind(qpsp_matrix, rank_weight_matrix[rownames(rank_weight_matrix)[which(rownames(rank_weight_matrix) %in% 
                                                                                                complex_list[[j]])], ])
    }
    else {
      qpsp_matrix <- rbind(qpsp_matrix, c(rep(0, ncol(rank_weight_matrix))))
    }
  }
  rownames(qpsp_matrix) <- names(complex_list)
  return(qpsp_matrix)
}

#' Rownames of df has to be the same annotation type as list of protein complexes
calc.qpsp1 <- function(df, list_complex) {
  # Filter out protein complexes with proteins that are not measured in df
  df_proteins <- rownames(df)
  # Logvec of whether entire complex is present
  complex_logvec <- sapply(list_complex, function(proteins_vec) all(proteins_vec %in% df_proteins))
  print(sum(complex_logvec))
  # Check for any NA in complex_logvec
  if (anyNA(complex_logvec)) stop("NA present in logvec")
  sublist_complex <- list_complex[complex_logvec]
  
  # Driver function that takes calculates QPSP profile for single sample
  # @param col_vec Column vector representing GFS-transformed single sample
  calc.qpsp_profile <- function(col_vec) {
    sapply(sublist_complex, function(proteins_vec) mean(col_vec[proteins_vec]))
  }
  return(apply(df, 2, calc.qpsp_profile))
}

qpsp_yeoh <- calc.qpsp1(gfs_yeoh, subnetwork_nea)
tail(qpsp_yeoh)
length(list_corum)

nea_size <- sapply(subnetwork_nea, length)

# Plot PCA before selecting features
# Batch information of all the timepoints
batch_info <- yeoh_batch[colnames(qpsp_yeoh), "batch"]
generate_colour <- colorRampPalette(c("lightblue", "darkblue"))
batch_palette <- generate_colour(10)
batch_colour <- batch_palette[batch_info]
# Shape of all timepoints
timepoint_shape <- rep(21:23, c(208,208,45))
plot.pca_3d(qpsp_yeoh, batch_colour, timepoint_shape)
rgl.postscript("dump/pca_3d-qpsp_nea.pdf", "pdf")

# Selecting drug responsive genes between D0 and D8
ttest_pvalue <- calc_ttest(qpsp_yeoh[,1:416], 208, is_paired = T)
# log_fc <- rowMeans(quantile_d8) - rowMeans(quantile_d0)
pvalue_probesets <- names(ttest_pvalue)[ttest_pvalue <= 0.05 & !is.nan(ttest_pvalue)]
# intersect_probesets <- fc_probesets[fc_probesets %in% pvalue_probesets]
print(length(pvalue_probesets))
selected_qpsp_yeoh <- qpsp_yeoh[pvalue_probesets,]

corrected_df <- qpsp_yeoh




# ANGLE -------------------------------------------------------------------
plot.angle <- function(df) {
  plot_1 <- ggplot(df) +
    geom_point(aes(x = angle_d0d8, y = diff_l2norm, col = factor(label)), show.legend = F) +
    scale_color_manual(values = c("darkolivegreen3", "tomato3"))

  plot_2 <- ggplot(df) +
    geom_point(aes(x = angle_d0d8, y = erm, col = factor(label)), show.legend = F) +
    scale_color_manual(values = c("darkolivegreen3", "tomato3"))

  multiplot <- plot_grid(plot_1, plot_2, ncol = 2)
  return(multiplot)
}

# Quantile %>% Feature selection %>% PCA
results_df <- collate_results(features_df1,
                              features_df2,
                              features_df3,
                              features_df4,
                              yeoh_label)

multiplot <- plot.angle(results_df)
multiplot
ggsave("dump/angle_pca_2.pdf", multiplot, width = 8, height = 4)

analysis_df <- cbind(results_df,
                     yeoh_label[rownames(results_df),"subtype", drop = F])
plot_2 <- ggplot(analysis_df) +
  geom_point(aes(x = d0d8_l2norm, y = erm1, col = factor(label)), show.legend = F) +
  scale_color_manual(values = c("darkolivegreen3", "tomato3")) +
  facet_wrap(~subtype, nrow = 2, ncol = 4)

ggsave("dump/quantile_pca_subtype.pdf", width = 12, height = 6)

# Subtype analysis


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

plot.angle <- function(df) {
  plot_1 <- ggplot(df) +
    geom_point(aes(x = angle_d0d8, y = d0d8_l2norm, col = factor(label)), show.legend = F) +
    scale_color_manual(values = c("darkolivegreen3", "tomato3"))
  
  plot_2 <- ggplot(df) +
    geom_point(aes(x = d0d8_l2norm, y = erm1, col = factor(label)), show.legend = F) +
    scale_color_manual(values = c("darkolivegreen3", "tomato3"))
  
  multiplot <- plot_grid(plot_1, plot_2, ncol = 2)
  return(multiplot)
}

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

