# library(Rtsne)
library(scran)
library(ggplot2)
library(cowplot)
library(dendextend)
library(rgl)
library(reshape2)
library(RColorBrewer)
source("../functions.R")
# FUNCTIONS ---------------------------------------------------------------
# Filter probes with too many zeros
filter_probesets <- function(df, percent_threshold) {
  logical_df <- df != 0
  selected_rows <- rowSums(logical_df) > percent_threshold * ncol(df)
  return(selected_rows)
}

plot_mean <- function(df, batch_vec1, batch_vec2) {
  # Melt dataframe
  melt_df <- melt(df, variable.name = "ID")
  print(head(melt_df))
  # Trimmed mean probe intensities for each chip
  mean_tibble <- melt_df %>% group_by(ID) %>%
    summarise(mean = mean(value, trim_percentage = 0.02))
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
calc_erm <- function(response_df, normal_df, leukemia_df = NA, flag = "replace") {
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
  leukemia_normal <- normal_centroid - leukemia_centroid
  l2_norm <- function(vec) sqrt(sum(vec^2))
  erm_factor <- leukemia_normal/l2_norm(leukemia_normal)
  
  # Assume that patients from top rows match correspondingly with bottom rows
  patient_arr <- cbind(response_df[1:num_patients,],
                       response_df[-(1:num_patients),])
  # Calculate vector by: D8-D0
  calc_response_vec <- function(row) row[-(1:num_dim)] - row[1:num_dim]
  num_dim <- ncol(patient_arr)/2
  response_vec_vstack <- apply(patient_arr, 1, calc_response_vec)
  # Multiplication of erm_factor is propagated through every column
  erm <- colSums(response_vec_vstack * erm_factor)
  return(erm)
}

# Plots ROC and calculates AUC in a primitive fashion (i.e. ROC is step function)
# Does not resolve ties in the score
# Assumption: Lower score will be labelled preferentially as 1, ROC is step function
# Assumption that score vec and label vecs are corresponding
plot_roc <- function(score_list, label_vec,
                     is_bigger_better_vec = rep(F, length(score_list)),
                     name_vec = NULL) {
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
  color_index <- (1:length(score_list)) + 1
  plot(NULL,
       xlab = "1 - Specificity", ylab = "Sensitivity",
       xlim = c(0,1), ylim = c(0,1),
       xaxs = "i", yaxs = "i",
       cex.main = 1.8, cex.lab = 1.7, cex.axis = 1.6)
  abline(0, 1, lty = 5, lwd = 2)
  auc_vec <- mapply(ROC_AUC, score_list, is_bigger_better_vec,
                    color_index, SIMPLIFY = T)
  # If name_vec is not NULL display legend
  if (!is.null(name_vec)) {
    format_text <- function(name, auc) {
      return(sprintf("%s (%.3f)", name, auc))
    }
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
    geom_point(size = 2, col = colour_code, shape = shape_vec, show.legend = F) +
    xlab(pc_labels[1]) + ylab(pc_labels[2])
  pc2_pc3 <- ggplot(df, aes(x = PC2, y = PC3)) +
    geom_point(size = 2, col = colour_code, shape = shape_vec, show.legend = F) +
    xlab(pc_labels[2]) + ylab(pc_labels[3])
  pc1_pc3 <- ggplot(df, aes(x = PC1, y = PC3)) +
    geom_point(size = 2, col = colour_code, shape = shape_vec, show.legend = F) +
    xlab(pc_labels[1]) + ylab(pc_labels[3])
  pc3_pc4 <- ggplot(df, aes(x = PC3, y = PC4)) +
    geom_point(size = 2, col = colour_code, shape = shape_vec, show.legend = F) +
    xlab(pc_labels[3]) + ylab(pc_labels[4])
  multiplot <- plot_grid(pc1_pc2, pc2_pc3, pc1_pc3, pc3_pc4,
                         ncol = 2, nrow = 2)
  return(multiplot)
}

# IMPORT DATA -------------------------------------------------------------
yeoh_data <- read.table("data/GSE67684/processed/mas5_ordered.tsv",
                        sep = "\t", header = T, row.names = 1)
yeoh_metadata <- read.table("data/GSE67684/processed/metadata_batch.tsv",
                            sep = "\t", header = T, row.names = 1)

mile_data <- read.table("data/GSE13204/processed/mas5_ordered.tsv",
                        sep = "\t", header = T, row.names = 1)
mile_metadata <- read.table("data/GSE13204/processed/metadata.tsv",
                            sep = "\t", header = T, row.names = 1)
colnames(mile_data)

yeoh_mean <- plot_mean(yeoh_data,
                       yeoh_metadata$batch,
                       yeoh_metadata$time_point)
ggsave("dump/trimmed_mean.pdf", yeoh_mean_plot,
       width = 9, height = 9)

mile_mean <- plot_mean(mile_data,
                       substring(mile_metadata$subtype, 17))
ggsave("dump/trimmed_mean-mile.pdf", mile_mean,
       width = 9, height = 9)

# DIFENG ------------------------------------------------------------------
# Trimmed-mean scaling
scaled_yeoh <- norm_mean_scaling(yeoh_data)
scaled_mile <- norm_mean_scaling(mile_data)
# Pre-processing of yeoh_data
selected_probes <- filter_probesets(scaled_yeoh, 0.2)
filtered_yeoh <- scaled_yeoh[selected_probes,]
# MILE data
filtered_mile <- scaled_mile[selected_probes,]
# # Quantile normalise before log transform to calculate logfc
# qnorm_data <- norm_quantile(cbind(filtered_yeoh, filtered_mile))
# qnorm_yeoh <- qnorm_data[1:ncol(filtered_yeoh)]
# qnorm_mile <- qnorm_data[-(1:ncol(filtered_yeoh))]

# Quantile normalise by time point
colnames(filtered_mile[,751:824])
qnorm_d0 <- norm_quantile(filtered_yeoh[,1:210])
qnorm_d8 <- norm_quantile(filtered_yeoh[,-(1:210)])
qnorm_normal <- norm_quantile(filtered_mile[,751:824])

# Selecting drug responsive genes between D0 and D8
ttest_pvalue <- calc_ttest(cbind(qnorm_d0, qnorm_d8), 210, is_paired = T)
log_fc <- calc_logfc(qnorm_d0, qnorm_d8)
pvalue_probesets <- names(ttest_pvalue)[ttest_pvalue <= 0.05]
fc_probesets <- names(log_fc)[log_fc > 1]
intersect_probesets <- fc_probesets[fc_probesets %in% pvalue_probesets]
print(length(intersect_probesets))

# Selecting drug responsive genes between D0 and normal
ttest_pvalue <- calc_ttest(cbind(qnorm_d0, qnorm_normal), 210, is_paired = F)
log_fc <- calc_logfc(qnorm_d0, qnorm_normal)
fc_probesets <- names(log_fc)[log_fc > 2]
intersect_probesets <- names(sort(ttest_pvalue[fc_probesets])[1:500])
print(length(intersect_probesets))

log_d0 <- log2_transform(qnorm_d0[intersect_probesets,])
log_d8 <- log2_transform(qnorm_d8[intersect_probesets,])
log_normal <- log2_transform(qnorm_normal[intersect_probesets,])

# PCA basis: D0 and normal
transposed_df <- t(cbind(log_d0, log_normal))
pca_obj <- prcomp(transposed_df, center = T, scale. = T)
pca_basis <- pca_obj$x[,1:4]

# Projection of D8 data
add_data <- t(log_d8)
pca_add <- predict(pca_obj, add_data)[,1:4]
plot_arr <- rbind(pca_basis[1:210,],
                  pca_add,
                  pca_basis[-(1:210),])

# PCA: Eigenvalues
eig_value <- (pca_obj$sdev)^2
var_pc <- eig_value[1:5]/sum(eig_value)
pc_labels <- sprintf("PC%d (%.2f%%)", 1:5, var_pc*100)

pca_arr <- plot_arr[,1:3]
colnames(pca_arr)
dim(pca_arr)
response_df <- pca_arr[1:420,]
normal_df <- pca_arr[421:494,]

erm <- calc_erm(response_df, normal_df)

# YEOH --------------------------------------------------------------------
# Trimmed-mean scaling
scaled_yeoh <- norm_mean_scaling(yeoh_data)
scaled_mile <- norm_mean_scaling(mile_data)

# Pre-processing of yeoh_data
selected_probes <- filter_probesets(scaled_yeoh, 0.2)
filtered_yeoh <- scaled_yeoh[selected_probes,]
log_yeoh <- log2_transform(filtered_yeoh)

# MILE data
filtered_mile <- scaled_mile[selected_probes,]
log_mile <- log2_transform(filtered_mile)

qnorm_data <- norm_quantile(cbind(log_yeoh, log_mile))
qnorm_yeoh <- qnorm_data[1:ncol(log_yeoh)]
qnorm_mile <- qnorm_data[-(1:ncol(log_yeoh))]

# Selection of down-regulated probes
downreg_probes <- select_probes(qnorm_yeoh)
# Check if all down-regulated probesets are present in MILE data
sum(downreg_probes %in% rownames(qnorm_mile))
leukemia_mile <- qnorm_mile[downreg_probes, 1:750]
normal_mile <- qnorm_mile[downreg_probes, -(1:750)]
yeoh_500_df <- qnorm_yeoh[downreg_probes,]
erm <- calc_erm(t(leukemia_mile), t(normal_mile), t(yeoh_500_df))

# GFS ---------------------------------------------------------------------
# Pre-processing of yeoh_data
selected_probes <- filter_probesets(yeoh_data, 0.2)
filtered_yeoh <- yeoh_data[selected_probes,]
# log_yeoh <- log2_transform(filtered_yeoh)

# MILE data
filtered_mile <- mile_data[selected_probes,]
# log_mile <- log2_transform(filtered_mile)

### GFS normalisation and selection of probesets
gfs_yeoh <- norm_gfs(filtered_yeoh)

plot_tsne <- Rtsne(t(gfs_yeoh), dims = 2, perplexity = 50)
names_index <- colnames(gfs_yeoh)
batch_info <- yeoh_metadata[names_index,6]
batch_palette <- brewer.pal(9, "Set1")
batch_colour <- c(batch_palette[batch_info],
                  batch_palette[batch_info])
plot(plot_tsne$Y, col = batch_colour, pch = rep(c(17,19), each = 210))
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
basis_data <- t(cbind(select_gfs_yeoh[1:210], select_gfs_mile[751:824]))
rownames(basis_data)
pca_obj <- prcomp(basis_data)
pca_basis <- pca_obj$x[,1:4]
# Projection of D8 and leukemia data
add_data <- t(cbind(select_gfs_yeoh[-(1:210)], select_gfs_mile[-(751:824)]))
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
# Trimmed-mean scaling
scaled_yeoh <- norm_mean_scaling(yeoh_data)
scaled_mile <- norm_mean_scaling(mile_data)
# Filter probesets
selected_probes <- filter_probesets(scaled_yeoh, 0.2)
filtered_yeoh <- yeoh_data[selected_probes,]
filtered_mile <- mile_data[selected_probes,]
# Feature selection using D0 and D8 data
ttest_pvalue <- calc_ttest(filtered_yeoh, 210, is_paired = T)
pvalue_probesets <- names(ttest_pvalue)[ttest_pvalue <= 0.05]
log_fc <- calc_logfc(filtered_yeoh[,1:210], filtered_yeoh[,-(1:210)])
fc_probesets <- names(log_fc)[log_fc > 1]
intersect_probesets <- fc_probesets[fc_probesets %in% pvalue_probesets]
print(length(intersect_probesets))
# Log-transform
log_yeoh <- log2_transform(filtered_yeoh[intersect_probesets,])
log_mile <- log2_transform(filtered_mile[intersect_probesets,])  

batch_info <- yeoh_metadata[colnames(log_yeoh),6]
batch_list <- lapply(1:9, function(i) log_yeoh[, batch_info == i])
arr_list <- lapply(batch_list, data.matrix)

colnum_list <- sapply(arr_list, ncol)
batch_order <- order(colnum_list, decreasing = T)
print(batch_order)

# Ordered according to number of samples in batch
mnn_yeoh_obj <- mnnCorrect(arr_list[[3]],
                           arr_list[[2]],
                           arr_list[[8]],
                           arr_list[[9]],
                           arr_list[[7]],
                           arr_list[[1]],
                           arr_list[[6]],
                           arr_list[[4]],
                           arr_list[[5]])
mnn_yeoh <- do.call(cbind, mnn_yeoh_obj$corrected)

# Column names for matrix arranged in above order
colnames_vec <- unlist(sapply(arr_list[batch_order], colnames))
colnames(mnn_yeoh) <- colnames_vec
ordered_mnn_yeoh <- mnn_yeoh[,order(colnames(mnn_yeoh))]
ordered_mnn_yeoh1 <- cbind(ordered_mnn_yeoh[,endsWith(colnames(ordered_mnn_yeoh), "0")],
                           ordered_mnn_yeoh[,endsWith(colnames(ordered_mnn_yeoh), "8")])

plot_mnn <- plot_pca(ordered_mnn_yeoh1, yeoh_metadata)
ggsave("dump/pca-yeoh_mnn.pdf", plot_mnn,
       width = 12, height = 4)

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
# # Batch information: Associating colour
# num_subtype <- table(substring(colnames(mile_data), 1, 1))
# red_palette <- brewer.pal(9, "Reds")[2:9]
# mile_palette <- c(red_palette, "darkolivegreen3")
# colour_code <- rep(mile_palette, num_subtype)

# Batch information of yeoh is encoded
batch_info <- yeoh_metadata[rownames(response_df),6]
blue_palette <- brewer.pal(9, "Blues")
batch_colour <- c(blue_palette[batch_info],
                  "tomato3", "darkolivegreen3")
# Colour information
all_colour <- c(rep(c("steelblue4", "turquoise3"), c(210, 210)),
                rep("tomato3", 750), rep("darkolivegreen3", 74))
# Colour information
centroid_colour <- c(rep(c("steelblue4", "turquoise3"), c(210, 210)),
                     rep("tomato3", 1), rep("darkolivegreen3", 1))

# Batch information of yeoh is encoded
names_index <- rownames(pca_arr)[1:210]
head(names_index)
batch_info <- yeoh_metadata[names_index,6]
batch_palette <- brewer.pal(9, "Blues")
batch_colour <- c(batch_palette[batch_info],
                  batch_palette[batch_info],
                  rep("darkolivegreen3", 74))
batch_shape <- c(rep(17, 210),
                 rep(19, 210+74))

# Dataframe to be visualised
visualise_arr <- rbind(response_df,
                       apply(leukemia_df, 2, median),
                       apply(normal_df, 2, median))

timepoint_vec <- c(rep(c(19, 17), each = 210), 19, 19)

rownames(plot_arr)
plot_pca <- pca_all(as.data.frame(plot_arr), batch_colour,
                    batch_shape, pc_labels)
plot_pca
dev.new()
ggsave("dump/pca_basis-quantile_d0normal.pdf", plot_pca,
       width = 9, height = 9)

# Visualising vectors
arrows_df <- cbind(response_df[1:210,], response_df[-(1:210),])
colnames(arrows_df)
colnames(arrows_df) <- paste(colnames(arrows_df),
                             rep(LETTERS[1:2], each = 3),
                             sep = "_")
centroid_df <- rbind(apply(response_df[1:210,], 2, median),
                     apply(response_df[-(1:210),], 2, median),
                     # apply(leukemia_df, 2, median),
                     apply(normal_df, 2, median))
plot_vectors <- function(df, centroid_df, pc_labels) {
  pc1_pc2 <- ggplot(data = df) +
    geom_point(aes(x = PC1_A, y = PC2_A), size = 3,
               col = "steelblue4", show.legend = F) +
    geom_point(aes(x = PC1_B, y = PC2_B), size = 3,
               col = "turquoise3", show.legend = F) +
    geom_point(data = centroid_df, aes(x = PC1, y = PC2),
               size = 5, shape = 17, colour = c("purple4", "violet", "darkolivegreen3")) +
    geom_segment(aes(x = PC1_A, y = PC2_A,
                     xend = PC1_B, yend = PC2_B),
                 arrow = arrow(length = unit(0.3, "cm")),
                 alpha = 0.5) +
    xlab(pc_labels[1]) + ylab(pc_labels[2])
  pc2_pc4 <- ggplot(data = df) +
    geom_point(aes(x = PC2_A, y = PC4_A), size = 3,
               col = "steelblue4", show.legend = F) +
    geom_point(aes(x = PC2_B, y = PC4_B), size = 3,
               col = "turquoise3", show.legend = F) +
    geom_point(data = centroid_df, aes(x = PC2, y = PC4),
               size = 5, shape = 17, colour = c("purple4", "violet", "darkolivegreen3")) +
    geom_segment(aes(x = PC2_A, y = PC4_A,
                     xend = PC2_B, yend = PC4_B),
                 arrow = arrow(length = unit(0.3, "cm")),
                 alpha = 0.5) +
    xlab(pc_labels[2]) + ylab(pc_labels[4])
  multiplot <- plot_grid(pc1_pc2, pc2_pc4, ncol = 2)
  return(multiplot)
}

vectors_plot <- plot_vectors(as.data.frame(arrows_df),
                             as.data.frame(centroid_df),
                             pc_labels)
vectors_plot

ggsave("dump/vectors-qnorm_d0normal_pca_basis_choose3.pdf", vectors_plot,
       width = 12, height = 6)

# RESULTS -----------------------------------------------------------------
### Extracting truth labels and training/test
labels_yeoh <- yeoh_metadata[names(erm), c(6, 8, 12)]
results_df <- cbind(erm, labels_yeoh)
head(results_df)
# Change truth labels with value 2 to 1
results_df[results_df$event_code == 2, 4] <- 1
# results_df <- results_df[order(results_df$erm),]

write.table(results_df, "dump/results-qnorm_d0d8_pca_basis_top3.tsv",
            quote = F, sep = "\t", row.names = T, col.names = T)

write.table(results_df, "dump/remove_centroid/results-gfs.tsv",
            quote = F, sep = "\t", row.names = T, col.names = T)

# Investigate
results_list <- split(results_df, results_df[,2])
results_list[[1]]

ROC_TITLE <- "GFS - PCA"
ROC_WPATH <- "dump/remove_centroid/roc-gfs.pdf"

plot_roc(list(results_df$erm),
         label_vec = results_df$event_code,
         name_vec = ROC_TITLE)
results_roc <- recordPlot()
save_fig(results_roc, ROC_WPATH,
         width = 9, height = 9)

# RETROSPECTIVE -----------------------------------------------------------
results_1 <- read.table("dump/results-gfs_ttest_pca_all.tsv",
                        sep = "\t", header = T, row.names = 1)
results_2 <- read.table("dump/results-gfs_ttest_d0d8_pca_basis_choose3.tsv",
                        sep = "\t", header = T, row.names = 1)
results_3 <- read.table("dump/results-gfs_ttest_d0d8_pca_basis_top3.tsv",
                        sep = "\t", header = T, row.names = 1)
results_4 <- read.table("dump/results-gfs_ttest_d0normal_pca_basis_choose3.tsv",
                        sep = "\t", header = T, row.names = 1)
results_5 <- read.table("dump/results-gfs_ttest_d0normal_pca_basis_top3.tsv",
                        sep = "\t", header = T, row.names = 1)

labels_vec <- results_1[order(rownames(results_1)), 4]
results1_vec <- results_1[order(rownames(results_1)), 1]
results2_vec <- results_2[order(rownames(results_2)), 1]
results3_vec <- results_3[order(rownames(results_3)), 1]
results4_vec <- results_4[order(rownames(results_4)), 1]
results5_vec <- results_5[order(rownames(results_5)), 1]
length(results5_vec)

par(mar = rep(5,4))
plot_roc(list(results1_vec, results2_vec, results3_vec, results4_vec, results5_vec), labels_vec,
         name_vec = c("PCA (All)",
                      "D0-D8, 3 PCs", "D0-D8, Top 3 PCs",
                      "D0-Normal, 3 PCs", "D0-Normal, Top 3 PCs"))
results_roc <- recordPlot()
save_fig(results_roc, "dump/results-gfs.pdf",
         width = 9, height = 9)

# yeoh_metadata <- read.table("data/GSE67684/processed/metadata_batch.tsv",
#                             sep = "\t", header = T, row.names = 1)
# results_list2 <- split(results_df2, results_df2[,2])
# for (i in 1:length(results_list)) {
#   print(results_list[[i]])
#   print(results_list2[[i]])
# }