library(reshape2)
library(ggplot2)
library(cowplot)
library(scran)
library(sva)
library(Harman)
source("../functions.R")
theme_set(theme_cowplot())
# Assumes that dataframe has been log-transformed
plot_batch <- function(df, batch_info, shape_info) {
  batch_factor <- factor(batch_info)
  shape_factor <- factor(shape_info)
  # Melt dataframe
  melt_df <- melt(df, variable.name = "ID")
  
  # Mean probe intensities for each chip
  mean_tibble <- melt_df %>% group_by(ID) %>%
    summarise(mean = mean(value))
  mean_jitter <- ggplot(mean_tibble, aes(x = ID,
                                         y = mean,
                                         col = batch_factor,
                                         shape = shape_factor)) +
    geom_point(show.legend = F, size = 3) + 
    theme(axis.text.x = element_text(size = 5, angle = 90, hjust = 1))
  
  # Principal component analysis
  col_logical <- apply(t(df), 2, sum) != 0 & apply(t(df), 2, var) != 0
  pca_df <- t(df)[, col_logical]
  pca_obj <- prcomp(pca_df, center = T, scale. = T)
  # Eigenvalues
  eig_value <- (pca_obj$sdev)^2
  var_pc <- eig_value[1:5]/sum(eig_value)
  pc_labels <- sprintf("PC%d (%.2f%%)", 1:5, var_pc*100)
  
  top_pc <- as.data.frame(pca_obj$x[,1:4])
  pc1_pc2 <- ggplot(top_pc, aes(x = PC1, y = PC2, col = batch_factor, shape = shape_factor)) +
    geom_point(size = 3, show.legend = F) +
    xlab(pc_labels[1]) + ylab(pc_labels[2])
  pc3_pc4 <- ggplot(top_pc, aes(x = PC3, y = PC4, col = batch_factor, shape = shape_factor)) +
    geom_point(size = 3, show.legend = F) +
    xlab(pc_labels[3]) + ylab(pc_labels[4])
  # Plot all graphs
  pca <- plot_grid(pc1_pc2, pc3_pc4)
  multiplot <- plot_grid(mean_jitter, pca, nrow = 2)
  return(multiplot)
}

# Returns: ggplot2 of PCA plot
plot_pca <- function(df, batch_info) {
  # Principal component analysis
  # Removes columns with all zeroes
  col_logical <- apply(t(df), 2, var) != 0
  pca_df <- t(df)[, col_logical]
  pca_obj <- prcomp(pca_df, center = T, scale. = T)
  top_pc <- as.data.frame(pca_obj$x[,1:2])
  pc1_pc2 <- ggplot(top_pc, aes(x = PC1, y = PC2, col = factor(batch_info))) +
    geom_point(size = 2, show.legend = F) +
    geom_label(label = rownames(pca_obj$x),
               nudge_x = 1, nudge_y = 2, size = 4,
               show.legend = F)
  return(pc1_pc2)
}
# IMPORT DATA -------------------------------------------------------------
raw_maqc <- read.table("../diff_expr/data/MAQC-I/processed/mas5_original-ref.tsv",
                       sep = "\t", header = T, row.names = 1)
colnames(raw_maqc)
batch_info <- rep(1:6, each = 10)
shape_info <- rep(rep(1:2, each = 5), 6)

# EXPLORE -----------------------------------------------------------------
# No normalisation
plot_before <- plot_batch(raw_maqc, batch_info, shape_info)
plot_before
ggsave("dump/plot_before.pdf", plot_before,
       width = 6, height = 6)

# Select only rows with all zeros
no_zero_maqc <- raw_maqc[rowSums(raw_maqc == 0) == 0,]
filtered_plot <- plot_batch(no_zero_maqc, batch_info, shape_info)
# Selecting only probesets without zeros accentuates the batch effects
filtered_plot

batch_1_A <- raw_maqc[,1:5]
batch_2_A <- raw_maqc[,11:15]
batch_1_2_A <- cbind(batch_1_A, batch_2_A)
# Select only rows with all zeros
filtered_1_2_A <- batch_1_2_A[rowSums(batch_1_2_A == 0) == 0,]
rank_1_2_A <- apply(filtered_1_2_A, 2, rank)
rank_1_A <- rank_1_2_A[,1:5]
rank_2_A <- rank_1_2_A[,6:10]

batch_1_A_sd <- apply(rank_1_A, 1, sd)
batch_2_A_sd <- apply(rank_2_A, 1, sd)
hist(batch_1_A_sd, breaks = 60)
hist(batch_2_A_sd, breaks = 60)
low_ps <- names(batch_1_A_sd[batch_1_A_sd < 100])
high_ps <- names(batch_1_A_sd[batch_1_A_sd > 2000])
# Most of the probesets with low sd in ranks are highly expressed
rank_1_A[high_ps,]
# 22155 probesets
dim(rank_1_A)
# Ratio of technical variation to signal
# If ratio passed threshold of one (MAYBE)

batch_1_B <- raw_maqc[,6:10]
batch_2_B <- raw_maqc[,16:20]
batch_1_2_B <- cbind(batch_1_B, batch_2_B)
# Select only rows with all zeros
filtered_1_2_B <- batch_1_2_B[rowSums(batch_1_2_B == 0) == 0,]
rank_1_2_B <- apply(filtered_1_2_B, 2, rank)
rank_1_B <- rank_1_2_B[,1:5]
rank_2_B <- rank_1_2_B[,6:10]

batch_1_B_sd <- apply(rank_1_B, 1, sd)
batch_2_B_sd <- apply(rank_2_B, 1, sd)
hist(batch_1_B_sd, breaks = 60)
hist(batch_2_B_sd, breaks = 60)
rank_diff_B <- rowMeans(rank_1_B) - rowMeans(rank_2_B)
hist(rank_diff_B, breaks = 60)
little_change_B <- names(rank_diff_B[rank_diff_B > -1 & rank_diff_B < 1])
rank_1_B[little_change_B,]
rank_2_B[little_change_B,]

# Selection of probesets
# Do not select probesets with high variance (Technical and biological variance)
rank_diff_A <- rowMeans(rank_1_A) - rowMeans(rank_2_A)
batch_1_A_sd[1:10] 
batch_2_A_sd[1:10]
mean(rank_diff)
sd(rank_diff)
hist(rnorm(20000,0,1688), breaks = 80, col = "lightblue")
hist(rank_diff, breaks = 80, col = "tomato")


# PCA loadings ------------------------------------------------------------
col_logical <- apply(t(df), 2, var) != 0
pca_df <- t(df)[, col_logical]
pca_obj <- prcomp(pca_df, center = T, scale. = T)
top_pc <- as.data.frame(pca_obj$x[,1:2])

df <- scaled_maqc
col_logical <- apply(t(df), 2, var) != 0
pca_df <- t(df)[, col_logical]
pca_obj <- prcomp(pca_df, center = T, scale. = T)
top_pc <- as.data.frame(pca_obj$x[,1:2])
hist(pca_obj$rotation[,3], breaks = 50)
loading_1 <- pca_obj$rotation[,1]
head(loading_1)
sorted_loading_1 <- sort(abs(loading_1), decreasing = T)
head(sorted_loading_1)

loading_2 <- pca_obj$rotation[,2]
sorted_loading_2 <- sort(abs(loading_2), decreasing = T)
ps_2 <- names(head(sorted_loading_2, 2000))

loading_3 <- pca_obj$rotation[,3]
sorted_loading_3 <- sort(abs(loading_3), decreasing = T)
ps_3 <- names(head(sorted_loading_3, 2000))
union_ps <- union(ps_2, ps_3)
loading_3[intersect_ps]
loading_2[intersect_ps]

subset_maqc <- scaled_maqc[union_ps,]
sample_maqc <- scaled_maqc[,1]
dim_maqc <- nrow(scaled_maqc)
quantile_maqc <- 1:dim_maqc/dim_maqc

View(subset_maqc)
# Filter out probesets that are heavily loaded in PC2 and PC3
filtered_scaled <- scaled_maqc[!(rownames(scaled_maqc) %in% union_ps),]
plot_filtered <- plot_batch(filtered_scaled, batch_info, shape_info)
plot_filtered

scaled_maqc <- norm.mean_scaling(raw_maqc)
plot_scaled <- plot_batch(scaled_maqc, batch_info, shape_info)
plot_scaled


# Local Batch Effects? ----------------------------------------------------
# Correct batch 5 and 6
batch_5_6 <- scaled_maqc[40:60]
batch_info1 <- rep(1:2, each = 10)
class_info1 <- rep(rep(LETTERS[1:2], each = 5), 2)

# factor_batch <- as.factor(batch_info1)
# factor_class <- as.factor(class_info1)
# # Metadata dataframe
# metadata_df <- data.frame(factor_batch, factor_class)
# head(metadata_df)
# 
# levels(factor_class)
# 
# # Returns vector of batches for each class
# get_batch <- function(class) metadata_df[metadata_df$factor_class == class, "factor_batch"]
# class_batch_info <- lapply(levels(factor_class), get_batch)
# names(class_batch_info) <- levels(factor_class)

# Correct batch 5 and 6
batch_5_6 <- scaled_maqc[40:60]
batch_info1 <- rep(1:2, each = 10)
class_info1 <- rep(rep(LETTERS[1:2], each = 5), 2)

batch_5_A <- scaled_maqc[41:45]
batch_5_B <- scaled_maqc[46:50]
batch_6_A <- scaled_maqc[51:55]
batch_6_B <- scaled_maqc[56:60]
# Two correction vectors: Class A & B
# Anchor: Batch 5
correction_A <- rowMeans(batch_6_A) - rowMeans(batch_5_A) 
correction_B <- rowMeans(batch_6_B) - rowMeans(batch_5_B)
# Correction performed
corrected_6_A <- batch_6_A - correction_A
corrected_6_B <- batch_6_B - correction_B
# Change negative values to 0
corrected_6_A[corrected_6_A < 0] <- 0
corrected_6_B[corrected_6_B < 0] <- 0

corrected_maqc <- cbind(scaled_maqc[,1:50], corrected_6_A, corrected_6_B)

plot_batch(corrected_maqc, batch_info, shape_info)

# COMPARISON (EQUAL) ---------------------------------------------------------------
# No normalisation
plot_before <- plot_batch(raw_maqc, batch_info, shape_info)
plot_before
ggsave("dump/plot_before.pdf", plot_before,
       width = 6, height = 6)

selected_probesets <- filter_probesets(raw_maqc, 0.98)
filtered_maqc <- raw_maqc[selected_probesets,]
nrow(filtered_maqc)
filtered_plot <- plot_batch(filtered_maqc, batch_info, shape_info)
filtered_plot


# Mean-scaling
scaled_maqc <- norm_mean_scaling(raw_maqc)
plot_scaled <- plot_batch(scaled_maqc, batch_info, shape_info)
plot_scaled
ggsave("dump/plot_scaled.pdf", plot_scaled,
       width = 6, height = 6)

# CDF
rank_maqc <- norm.cdf(raw_maqc)
plot_batch(rank_maqc, batch_info, shape_info)

# GFS
gfs_maqc <- norm_gfs(raw_maqc)
plot_gfs <- plot_batch(gfs_maqc, batch_info, shape_info)
ggsave("dump/plot_gfs.pdf", plot_gfs,
       width = 6, height = 6)

# MNN
raw_arr <- data.matrix(raw_maqc)
# Ordered according to number of samples in batch
mnn_maqc_obj <- mnnCorrect(raw_arr[,1:10],
                           raw_arr[,11:20],
                           raw_arr[,21:30],
                           raw_arr[,31:40],
                           raw_arr[,41:50],
                           raw_arr[,51:60])
                           
mnn_maqc <- do.call(cbind, mnn_maqc_obj$corrected)
# Column names for matrix arranged in above order
colnames(mnn_maqc) <- colnames(raw_arr)
plot_mnn <- plot_batch(data.frame(mnn_maqc), batch_info, shape_info)
ggsave("dump/plot_mnn.pdf", plot_mnn,
       width = 6, height = 6)

# MNN - Scaled
scaled_arr <- data.matrix(norm.mean_scaling(raw_maqc))
# Ordered according to number of samples in batch
mnn_maqc_obj <- mnnCorrect(scaled_arr[,1:10],
                           scaled_arr[,11:20],
                           scaled_arr[,21:30],
                           scaled_arr[,31:40],
                           scaled_arr[,41:50],
                           scaled_arr[,51:60])

mnn_maqc <- do.call(cbind, mnn_maqc_obj$corrected)
# Column names for matrix arranged in above order
colnames(mnn_maqc) <- colnames(raw_arr)
plot_mnn_scaled <- plot_batch(data.frame(mnn_maqc), batch_info, shape_info)
ggsave("dump/plot_mnn_scaled.pdf", plot_mnn_scaled,
       width = 6, height = 6)

# Quantile normalisation
col_indx <- substring(colnames(raw_maqc), 7, 7) == "A"
qnorm_a <- norm_quantile(raw_maqc[,col_indx])
qnorm_b <- norm_quantile(raw_maqc[,!col_indx])
qnorm_maqc <- cbind(qnorm_a, qnorm_b)

colnames(qnorm_maqc)
batch_qnorm <- rep(rep(1:6, each = 5), 2)
shape_qnorm <- rep(1:2, each = 30)
plot_qnorm <- plot_batch(qnorm_maqc, batch_qnorm, shape_qnorm)
ggsave("dump/plot_qnorm.pdf", plot_qnorm,
       width = 6, height = 6)

# Harman
scaled_maqc <- norm.mean_scaling(raw_maqc)
log_maqc <- log2_transform(scaled_maqc)

harman_obj <- harman(log_maqc, shape_info, batch_info)
harman_maqc <- data.frame(reconstructData(harman_obj))

plot_harman <- plot_batch(harman_maqc, batch_info, shape_info)
ggsave("dump/plot_harman.pdf", plot_harman,
       width = 6, height = 6)

# Scanorama
library(reticulate)
scanorama <- import('scanorama')
# ComBat ------------------------------------------------------------------
# Creation of pData dataframe (metadata)
class <- as.factor(substring(colnames(raw_maqc), 7, 7))
batch <- as.factor(substring(colnames(raw_maqc), 5, 5))
maqc_metadata <- data.frame(class, batch)
# Rownames of metadata are same as colnames of data df
rownames(maqc_metadata) <- colnames(raw_maqc)
maqc_metadata

# ComBat assumes that data has been normalised and probesets have been filtered
# Error if probesets are not filtered as rows have 0 variance
scaled_maqc <- norm.mean_scaling(raw_maqc)
selected_probetsets <- filter_probesets(scaled_maqc, 0.1)
filtered_maqc <- scaled_maqc[selected_probetsets,]

# Place adjustment/confounding variables in model.matrix (e.g. age)
# Do not put batch variables in model.matrix
# Put batch variables directly in combat function
model_combat <- model.matrix(~1, data = maqc_metadata)
combat_maqc <- ComBat(data.matrix(filtered_maqc),
                      batch = maqc_metadata$batch,
                      mod = model_combat)
combat_maqc_df <- data.frame(combat_maqc)
# Replace negative values with 0
combat_maqc_df[combat_maqc_df < 0] <- 0

plot_combat_scaled <- plot_batch(combat_maqc_df, batch_info, shape_info)
plot_combat_scaled

plot_log_combat_scaled <- plot_batch(log2_transform(combat_maqc_df), batch_info, shape_info)
plot_log_combat_scaled
ggsave("dump/plot_log_combat_scaled.pdf", plot_log_combat_scaled,
       width = 6, height = 6)

# Log-transform data after mean-scaling
log_maqc <- log2_transform(filtered_maqc)
combat_log_maqc <- ComBat(data.matrix(log_maqc),
                          batch = maqc_metadata$batch,
                          mod = model_combat)
combat_log_maqc_df <- data.frame(combat_log_maqc)
# Replace negative values with 0
combat_log_maqc_df[combat_log_maqc_df < 0] <- 0

plot_combat_log_scaled <- plot_batch(combat_log_maqc_df, batch_info, shape_info)
ggsave("dump/plot_combat_log_scaled.pdf", plot_combat_log_scaled,
       width = 6, height = 6)

# CBC ---------------------------------------------------------------------
class_info <- ifelse(shape_info == 1, "A", "B")

# Perform trimmed mean scaling first
scaled_maqc <- norm.mean_scaling(raw_maqc)
cbc_maqc <- norm.CBC(scaled_maqc, batch_info, class_info, 1:6)

plot_scaled <- plot_batch(scaled_maqc, batch_info, class_info)
plot_log_scaled <- plot_batch(log2_transform(scaled_maqc), batch_info, class_info)

ggsave("dump/plot_scaled.pdf", plot_scaled,
       width = 6, height = 6)

ggsave("dump/plot_log_scaled.pdf", plot_log_scaled,
       width = 6, height = 6)

# Apply log first before scaling or normalisation
log_maqc <- log2_transform(raw_maqc)
scaled_log_maqc <- norm.mean_scaling(log_maqc)
cbc_scaled_log_maqc <- norm.CBC(scaled_log_maqc, batch_info, class_info, 1:6)
plot_log_raw <- plot_batch(log_maqc, batch_info, class_info)
plot_scaled_log <- plot_batch(scaled_log_maqc, batch_info, class_info)
plot_scaled_log
plot_cbc_scaled_log <- plot_batch(cbc_scaled_log_maqc, batch_info, class_info)

# Log-transform after scaling but before batch correction
log_scaled_maqc <- log2_transform(scaled_maqc)
cbc_log_scaled_maqc <- norm.CBC(log_scaled_maqc, batch_info, class_info, 1:6)
plot_cbc_log_scaled <- plot_batch(cbc_log_scaled_maqc, batch_info, class_info)
plot_log_cbc_log_scaled <- plot_batch(log2_transform(cbc_log_scaled_maqc), batch_info, class_info)

ggsave("dump/plot_cbc_log_scaled.pdf", plot_cbc_log_scaled,
       width = 6, height = 6)
ggsave("dump/plot_log_cbc_log_scaled.pdf", plot_log_cbc_log_scaled,
       width = 6, height = 6)


# Proportion Study --------------------------------------------------------
# Selects samples such that each batch has a diff proportion of classes
# Proportions rep(rep(LETTERS[1:2], 6), c(5,5,4,5,5,4,3,5,5,3,3,4))
selection_index <- c(1:10,11:14,16:20,21:25,26:29,31:33,36:40,41:45,46:48,51:53,56:59)
odd_maqc <- raw_maqc[,selection_index]
odd_batch_info <- rep(1:6, c(10,9,9,8,8,7))
odd_class_info <- substring(colnames(odd_maqc), 7, 7)
# Scale and log-transform odd_maqc
scaled_maqc <- norm.mean_scaling(odd_maqc)
log_maqc <- log2_transform(scaled_maqc)

# ComBat
# Creation of pData dataframe (metadata)
class <- as.factor(odd_class_info)
batch <- as.factor(odd_batch_info)
maqc_metadata <- data.frame(class, batch)
# Rownames of metadata are same as colnames of data df
rownames(maqc_metadata) <- colnames(odd_maqc)
head(maqc_metadata, 20)

# ComBat assumes that data has been normalised and probesets have been filtered
# Error if probesets are not filtered as rows have 0 variance
scaled_maqc <- norm.mean_scaling(odd_maqc)
selected_probetsets <- filter_probesets(scaled_maqc, 0.1)
filtered_maqc <- scaled_maqc[selected_probetsets,]
log_maqc <- log2_transform(filtered_maqc)
# Place adjustment/confounding variables in model.matrix (e.g. age)
# Do not put batch variables in model.matrix
# Put batch variables directly in combat function
model_combat <- model.matrix(~1, data = maqc_metadata)
combat_maqc <- ComBat(data.matrix(log_maqc),
                      batch = maqc_metadata$batch,
                      mod = model_combat)
combat_maqc_df <- data.frame(combat_maqc)
# Replace negative values with 0
combat_maqc_df[combat_maqc_df < 0] <- 0

plot_combat_log_scaled <- plot_batch(combat_maqc_df, odd_batch_info, odd_class_info)
plot_combat_log_scaled

ggsave("dump/plot_combat_log_scaled_odd.pdf", plot_combat_log_scaled,
       width = 6, height = 6)

# CBC
scaled_maqc <- norm.mean_scaling(odd_maqc)
log_maqc <- log2_transform(scaled_maqc)
cbc_maqc <- norm.CBC(log_maqc, odd_batch_info, odd_class_info, 1:6)
plot_cbc <- plot_batch(cbc_maqc, odd_batch_info, odd_class_info)
ggsave("dump/plot_cbc_log_scaled_odd.pdf", plot_cbc,
       width = 6, height = 6)

# Harman
harman_obj <- harman(log_maqc, odd_class_info, odd_batch_info, limit = 0.95)
harman_maqc <- data.frame(reconstructData(harman_obj))

plot_harman_odd <- plot_batch(harman_maqc, odd_batch_info, odd_class_info)
plot_harman_odd
ggsave("dump/plot_harman_odd.pdf", plot_harman,
       width = 6, height = 6)
