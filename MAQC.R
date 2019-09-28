library(reshape2)
library(ggplot2)
library(cowplot)
# library(scran)
library(sva)
# library(Harman)
source("../functions.R")
source("class_batch_correction.R")
theme_set(theme_cowplot())
# Assumes that dataframe has been log-transformed
plot_batch <- function(df, batch_info, shape_info) {
  batch_factor <- factor(batch_info)
  shape_factor <- factor(shape_info)
  # Principal component analysis
  col_logical <- apply(t(df), 2, sum) != 0 & apply(t(df), 2, var) != 0
  pca_df <- t(df)[, col_logical]
  pca_obj <- prcomp(pca_df, center = T, scale. = T)
  pca_arr <- pca_obj$x
  # Eigenvalues
  pca_df <- data.frame(pca_arr)
  eig_value <- (pca_obj$sdev)^2
  # Percentage variance of each PC
  var_pct <- eig_value/sum(eig_value)
  pc_labels <- sprintf("PC%d (%.2f%%)", 1:5, var_pct[1:5]*100)

  # # Argument: PC coordinates for single PC
  # calc.pc_var <- function(vec) {
  #   # Argument: ANOVA attributes; Calculates percentage of variance
  #   # SS_between/SS_between + SS_within
  #   calc.var_percentage <- function(vec) unname(vec[3]/(vec[3] + vec[4]))
  #   pc_metadata <- data.frame(pc = vec,
  #                             batch = as.factor(batch_info),
  #                             class = as.factor(class_info))
  #   batch_anova_attr <- unlist(summary(aov(pc~batch, data = pc_metadata)))
  #   class_anova_attr <- unlist(summary(aov(pc~class, data = pc_metadata)))
  #   return(c(calc.var_percentage(batch_anova_attr),
  #            calc.var_percentage(class_anova_attr)))
  # }
  # var_composition <- sapply(pca_df, calc.pc_var)
  # pca_var_df <- data.frame(var_pct,
  #                          batch_pct = var_composition[1,],
  #                          class_pct = var_composition[2,])
  # total_batch_pct <- sum(pca_var_df$var_pct * pca_var_df$batch_pct)*100
  # total_class_pct <- sum(pca_var_df$var_pct * pca_var_df$class_pct)*100
  
  # PCA scatter plot
  top_pc <- as.data.frame(pca_obj$x[,1:4])
  pc1_pc2 <- ggplot(top_pc, aes(x = PC1, y = PC2, col = batch_factor, shape = shape_factor)) +
    geom_point(size = 3, show.legend = F) +
    labs(x = pc_labels[1], y = pc_labels[2])
  
  pc3_pc4 <- ggplot(top_pc, aes(x = PC3, y = PC4, col = batch_factor, shape = shape_factor)) +
    geom_point(size = 3, show.legend = F) +
    labs(x = pc_labels[3], y = pc_labels[4])
  
  # PLOT MEANS
  # Melt dataframe
  melt_df <- melt(df, variable.name = "ID")
  # Mean probe intensities for each chip
  mean_tibble <- melt_df %>% group_by(ID) %>%
    summarise(mean = mean(value))
  mean_scatter <- ggplot(mean_tibble, aes(x = ID,
                                          y = mean,
                                          col = batch_factor,
                                          shape = shape_factor)) +
    geom_point(show.legend = F, size = 3) +
    # labs(title = sprintf("BE: %.2f%% \n CV: %.2f%%",
    #                      total_batch_pct, total_class_pct)) +
    theme(axis.text.x = element_text(size = 5, angle = 90, hjust = 1))
  
  # Plot all graphs
  pca <- plot_grid(pc1_pc2, pc3_pc4)
  multiplot <- plot_grid(mean_scatter, pca, nrow = 2)
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
filtered_maqc <- raw_maqc[filter_probesets(raw_maqc, 0.1),]
scaled_maqc <- norm.mean_scaling(filtered_maqc)
log_maqc <- log2_transform(scaled_maqc)

# MAQC metadata
batch_info <- rep(1:6, each = 10)
class_info <- rep(rep(1:2, each = 5), 6)

# COMPARISON (EQUAL) ---------------------------------------------------------------
# # No normalisation
# plot_before <- plot_batch(raw_maqc, batch_info, shape_info)
# plot_before
# ggsave("dump/plot_before.pdf", plot_before,
#        width = 6, height = 6)
# 
# selected_probesets <- filter_probesets(raw_maqc, 0.98)
# filtered_maqc <- raw_maqc[selected_probesets,]
# nrow(filtered_maqc)
# filtered_plot <- plot_batch(filtered_maqc, batch_info, shape_info)
# filtered_plot
# 
# # Mean-scaling
# scaled_maqc <- norm_mean_scaling(raw_maqc)
# plot_scaled <- plot_batch(scaled_maqc, batch_info, shape_info)
# plot_scaled
# ggsave("dump/plot_scaled.pdf", plot_scaled,
#        width = 6, height = 6)
# 
# # CDF
# rank_maqc <- norm.cdf(raw_maqc)
# plot_batch(rank_maqc, batch_info, shape_info)
# 
# # GFS
# gfs_maqc <- norm_gfs(raw_maqc)
# plot_gfs <- plot_batch(gfs_maqc, batch_info, shape_info)
# ggsave("dump/plot_gfs.pdf", plot_gfs,
#        width = 6, height = 6)

# MNN
# Split dataframe according to batch info
list_df <- split.default(log_maqc, batch_info)
# Convert df to arr
list_arr <- lapply(list_df, data.matrix)
list_args <- c(list_arr, k = 5, cos.norm.out = F)
mnn_maqc_obj <- do.call(mnnCorrect, list_args)
mnn_maqc_arr <- do.call(cbind, mnn_maqc_obj$corrected)
# Column names for matrix arranged in above order
colnames(mnn_maqc_arr) <- colnames(log_maqc)
plot_mnn <- plot_batch(data.frame(mnn_maqc_arr), batch_info, class_info)
plot_mnn

ggsave("dump/plot_mnn_log_scaled_k5.pdf", plot_mnn_log_scaled,
       width = 6, height = 6)

# Investigate MNNs ability to correct for multiplicative error

# Quantile normalisation
col_indx <- substring(colnames(log_maqc), 7, 7) == "A"
qnorm_a <- norm.quantile(log_maqc[,col_indx])
qnorm_b <- norm.quantile(log_maqc[,!col_indx])
qnorm_maqc <- cbind(qnorm_a, qnorm_b)
colnames(qnorm_maqc)
plot_qnorm <- plot_batch(qnorm_maqc, batch_info, class_info)
plot_qnorm
ggsave("dump/plot_class_qnorm.pdf", plot_qnorm,
       width = 6, height = 6)

# Harman
scaled_maqc <- norm.mean_scaling(raw_maqc)
log_maqc <- log2_transform(scaled_maqc)

harman_obj <- harman(scaled_maqc, class_info, batch_info)
harman_maqc <- data.frame(reconstructData(harman_obj))

plot_harman <- plot_batch(harman_maqc, batch_info, class_info)
plot_harman
ggsave("dump/plot_harman.pdf", plot_harman,
       width = 6, height = 6)

# Scanorama
# Write log_maqc for numpy array
write.table(log_maqc, "data/scanorama/log_maqc.tsv",
            quote = F, sep = "\t", row.names = F, col.names = F)
write(rownames(log_maqc), file = "data/scanorama/probeset_names.txt")

scanorama_maqc <- read.table("data/scanorama/scanorama_data_k10.tsv",
                             sep = "\t", row.names = 1)
colnames(scanorama_maqc) <- colnames(log_maqc)
plot_scanorama <- plot_batch(scanorama_maqc, batch_info, class_info)
plot_scanorama
ggsave("dump/plot_scanorama_k10.pdf", plot_scanorama,
       width = 6, height = 6)

# ComBat ------------------------------------------------------------------
# Creation of pData dataframe (metadata)
class <- as.factor(class_info)
batch <- as.factor(batch_info)
maqc_metadata <- data.frame(class, batch)
# Rownames of metadata are same as colnames of data df
rownames(maqc_metadata) <- colnames(raw_maqc)
head(maqc_metadata)
# Place adjustment/confounding variables in model.matrix (e.g. age)
# Do not put batch variables in model.matrix
# Put batch variables directly in combat function
# model_combat <- model.matrix(~1, data = maqc_metadata)
# Include biological variable of interest as covariate
model_combat <- model.matrix(~class, data = maqc_metadata)

combat_maqc <- ComBat(data.matrix(log_maqc),
                      batch = maqc_metadata$batch,
                      mod = model_combat)
combat_maqc_df <- data.frame(combat_maqc)
# Replace negative values with 0
combat_maqc_df[combat_maqc_df < 0] <- 0

plot_combat <- plot_batch(combat_maqc_df, batch_info, class_info)
plot_combat

ggsave("dump/plot_FSL_combatCovariate.pdf", plot_combat,
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
class_info <- ifelse(class_info == 1, "A", "B")
cbc_maqc <- norm.CBC(log_maqc, batch_info, class_info, 1:6, "dump/correction_log_maqc.tsv")

plot_cbc <- plot_batch(cbc_maqc, batch_info, class_info)
plot_cbc

ggsave("dump/plot_cbc.pdf", plot_cbc,
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

# PROPORTION DATASET --------------------------------------------------------
# Selects samples such that each batch has a diff proportion of classes
# Proportions rep(rep(LETTERS[1:2], 6), c(5,5,4,5,5,4,3,5,5,3,3,4))
selection_index <- c(1:10,11:14,16:20,21:25,26:29,31:33,36:40,41:45,46:48,51:53,56:59)
odd_maqc <- log_maqc[,selection_index]
odd_batch_info <- rep(1:6, c(10,9,9,8,8,7))
odd_class_info <- substring(colnames(odd_maqc), 7, 7)

# ComBat
# Creation of pData dataframe (metadata)
class <- as.factor(odd_class_info)
batch <- as.factor(odd_batch_info)
odd_maqc_metadata <- data.frame(class, batch)
# Rownames of metadata are same as colnames of data df
rownames(odd_maqc_metadata) <- colnames(odd_maqc)
head(odd_maqc_metadata, 20)

# ComBat assumes that data has been normalised and probesets have been filtered
# Error if probesets are not filtered as rows have 0 variance
selected_probesets <- filter_probesets(odd_log_maqc, 0.1)
odd_filtered_maqc <- odd_maqc[selected_probesets,]
# Place adjustment/confounding variables in model.matrix (e.g. age)
# Do not put batch variables in model.matrix
# Put batch variables directly in combat function
# # Only include intercept term in design matrix
# model_combat <- model.matrix(~1, data = maqc_metadata)
# Include biological variable of interest as covariate
model_combat <- model.matrix(~class, data = odd_maqc_metadata)

combat_maqc <- ComBat(data.matrix(odd_filtered_maqc),
                      batch = odd_maqc_metadata$batch,
                      mod = model_combat)
combat_maqc_df <- data.frame(combat_maqc)
# Replace negative values with 0
combat_maqc_df[combat_maqc_df < 0] <- 0

plot_combat_log_scaled <- plot_batch(combat_maqc_df, odd_batch_info, odd_class_info)
plot_combat_log_scaled

ggsave("dump/plot_FSL_combatCov-odd.pdf", plot_combat_log_scaled,
       width = 6, height = 6)

# CBC
scaled_maqc <- norm.mean_scaling(odd_maqc)
log_maqc <- log2_transform(scaled_maqc)
cbc_maqc <- norm.CBC(log_maqc, odd_batch_info, odd_class_info, 1:6)
plot_cbc <- plot_batch(cbc_maqc, odd_batch_info, odd_class_info)
ggsave("dump/plot_cbc_log_scaled_odd.pdf", plot_cbc,
       width = 6, height = 6)

# Harman
harman_obj <- harman(odd_log_maqc, odd_class_info, odd_batch_info, limit = 0.95)
harman_maqc <- data.frame(reconstructData(harman_obj))

plot_harman_odd <- plot_batch(harman_maqc, odd_batch_info, odd_class_info)
plot_harman_odd
ggsave("dump/plot_harman_odd.pdf", plot_harman_odd,
       width = 6, height = 6)

# Scanorama
# Write log_maqc for numpy array
write.table(odd_log_maqc, "data/scanorama/odd_log_maqc.tsv",
            quote = F, sep = "\t", row.names = F, col.names = F)
# Read scanorama corrected data
scanorama_maqc_odd <- read.table("data/scanorama/scanorama_data_odd_k5.tsv",
                                 sep = "\t", row.names = 1)
colnames(scanorama_maqc_odd) <- colnames(odd_maqc)
plot_scanorama_odd <- plot_batch(scanorama_maqc_odd,
                                 odd_batch_info, odd_class_info)
plot_scanorama_odd
ggsave("dump/plot_scanorama_odd_k5.pdf", plot_scanorama_odd,
       width = 6, height = 6)

# MNN
list_small_df <- split.default(odd_log_maqc, odd_batch_info)
list_small_arr <- lapply(list_small_df, data.matrix)

list_args <- c(list_small_arr, k = 3)
mnn_maqc_obj <- do.call(mnnCorrect, list_args)

mnn_maqc_odd <- do.call(cbind, mnn_maqc_obj$corrected)
# Column names for matrix arranged in above order
colnames(mnn_maqc_odd) <- colnames(odd_log_maqc)
plot_mnn_odd_k5 <- plot_batch(data.frame(mnn_maqc_odd),
                              odd_batch_info, odd_class_info)
plot_mnn_odd_k5

ggsave("dump/plot_mnn_odd_k5.pdf", plot_mnn_odd_k5,
       width = 6, height = 6)

# SMALL DATASET -----------------------------------------------------------
small_index <- 1:3 + rep(seq(0,55,5), each = 3)
small_maqc <- raw_maqc[,small_index]
small_batch_info <- rep(1:6, each = 6)
small_class_info <- substring(colnames(small_maqc), 7, 7)
# Scale and log-transform odd_maqc
small_scaled_maqc <- norm.mean_scaling(small_maqc)
small_log_maqc <- log2_transform(small_scaled_maqc)

# Scanorama
# Write small_log_maqc for numpy array
write.table(small_log_maqc, "data/scanorama/small_log_maqc.tsv",
            quote = F, sep = "\t", row.names = F, col.names = F)
# Read scanorama corrected data
scanorama_maqc_small <- read.table("data/scanorama/scanorama_data_small_k15.tsv",
                                   sep = "\t", row.names = 1)
colnames(scanorama_maqc_small) <- colnames(small_log_maqc)
plot_scanorama_small <- plot_batch(scanorama_maqc_small,
                                   small_batch_info, small_class_info)
plot_scanorama_small
ggsave("dump/plot_scanorama_small_k15.pdf", plot_scanorama_small,
       width = 6, height = 6)

# ComBat
# Creation of pData dataframe (metadata)
class <- as.factor(small_class_info)
batch <- as.factor(small_batch_info)
maqc_metadata <- data.frame(class, batch)
# Rownames of metadata are same as colnames of data df
rownames(maqc_metadata) <- colnames(small_maqc)
head(maqc_metadata, 20)

# ComBat assumes that data has been normalised and probesets have been filtered
# Error if probesets are not filtered as rows have 0 variance
selected_probetsets <- filter_probesets(small_log_maqc, 0.1)
filtered_maqc <- small_log_maqc[selected_probetsets,]
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

plot_combat_log_scaled_small <- plot_batch(combat_maqc_df, small_batch_info, small_class_info)
plot_combat_log_scaled_small

ggsave("dump/plot_combat_log_scaled_small.pdf", plot_combat_log_scaled_small,
       width = 6, height = 6)

# CBC
cbc_maqc_small <- norm.CBC(small_log_maqc, small_batch_info, small_class_info, 1:6)
plot_cbc_log_scaled_small <- plot_batch(cbc_maqc_small, small_batch_info, small_class_info)
plot_cbc_log_scaled_small
ggsave("dump/plot_cbc_log_scaled_small.pdf", plot_cbc_log_scaled_small,
       width = 6, height = 6)

# MNN - Scaled
small_log_arr <- data.matrix(small_log_maqc)
list_small_df <- split.default(small_log_maqc, rep(1:6, each = 6))
list_small_arr <- lapply(list_small_df, data.matrix)
list_args <- c(list_small_arr, k = 3, cos.norm.out = T)
mnn_maqc_obj <- do.call(mnnCorrect, list_args)

# # Ordered according to number of samples in batch
# mnn_maqc_obj <- mnnCorrect(scaled_arr[,1:10],
#                            scaled_arr[,11:20],
#                            scaled_arr[,21:30],
#                            scaled_arr[,31:40],
#                            scaled_arr[,41:50],
#                            scaled_arr[,51:60],
#                            k = 5)

mnn_maqc <- do.call(cbind, mnn_maqc_obj$corrected)
# Column names for matrix arranged in above order
colnames(mnn_maqc) <- colnames(small_log_maqc)

# Replace negative values with zeros
plot_mnn_log_scaled_small <- plot_batch(data.frame(mnn_maqc), small_batch_info, small_class_info)
plot_mnn_log_scaled_small

ggsave("dump/plot_mnn_log_scaled_small_k3.pdf", plot_mnn_log_scaled_small,
       width = 6, height = 6)

# Harman
harman_obj <- harman(small_log_maqc, small_class_info,
                     small_batch_info, limit = 0.95)
harman_maqc <- data.frame(reconstructData(harman_obj))

plot_harman_small <- plot_batch(harman_maqc, small_batch_info, small_class_info)
plot_harman_small
ggsave("dump/plot_harman_log_scaled_small.pdf", plot_harman_small,
       width = 6, height = 6)


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



# Scanorama ---------------------------------------------------------------
library(reticulate)
scanorama <- import('scanorama')
# # Scanorama requires matrices of samples x genes
# # Transposing converts df into matrix
# # Call split.data.frame explicity to work with matrices
# # List of batch matrices
# list_batch_mat <- split.data.frame(t(log_maqc), as.factor(batch_info))
# num_batches <- 6
# genes_list <- rep(list(rownames(log_maqc)), num_batches)
# 
# integrated_data <- scanorama$integrate(list_batch_mat, genes_list)
# corrected.data <- scanorama$correct(datasets, genes_list, return_dense=TRUE)
# integrated.corrected.data <- scanorama$correct(datasets, genes_list,
#                                                return_dimred=TRUE, return_dense=TRUE)


# MODELLING ---------------------------------------------------------------

plot_evaluation(log2_transform(raw_maqc), batch_info)

x <- log2_transform(norm.mean_scaling(raw_maqc))
colnames(x)
plot(x[,6], x[,7], main = "mean(raw_maqc) (B1B1 vs B1B2)")
plot(rowMeans(x[,1:5]), rowMeans(x[,11:15]), main = "log2(raw_maqc) (B1-MeanA vs B2-MeanA)")
abline(a=0, b = 1)
dev.copy(png, "dump/scatter_fig10.png")
dev.off()

plot_maqc <- plot_batch(x, batch_info, class_info)
plot_maqc

batch_factor <- as.factor(rep(1:2, each = 10))
class_factor <- as.factor(rep(rep(LETTERS[1:2], each = 5), 2))

corrected_df <- norm.CBC(x[,1:20], batch_factor, class_factor, c(1,2), "dump/log_mean_correction.tsv")

correction_vectors <- read.table("dump/log_mean_correction.tsv", sep = "\t")

plot_batch(corrected_df, batch_factor, class_factor)
plot(correction_vectors$A, correction_vectors$B)

A_avg <- rowMeans(x[,11:15])
B_avg <- rowMeans(x[,16:15])
observed_avg <- data.frame(A_avg, B_avg)
big_diff_ps <- names(A_avg)[abs(A_avg-B_avg) > 5 & A_avg != 0 & B_avg != 0]

par(mfrow=c(1,1))

for (i in 1:12) {
  probeset <- big_diff_ps[i]
  plot(unlist(observed_avg[probeset,]),
       unlist(correction_vectors[probeset,]))
  abline(a = 0, b = 1)
}

k_arr <- matrix(0, 54675, 2)
for (i in 1:54675) {
  k_arr[i,1] <- correction_vectors[i,1]/observed_avg[i,1]
  k_arr[i,2] <- correction_vectors[i,2]/observed_avg[i,2]
}

plot(k_arr[,1], k_arr[,2])

df <- data.frame(A_avg[1:10], A_avg[1:10])

head(correction_vectors)
head(observed_avg)

max(correction_vectors[,1])

x <- log2_transform(norm.mean_scaling(raw_maqc))
# B1A1 - B2A1
plot(x[,1]-x[,11], x[,4]-x[,15]) 
plot(x[,1]-x[,11], x[,1]-x[,11])

par(mfrow=c(1,1))

plot(x[,6], x[,16],
     main = "log2(mean(raw_maqc)) (B1B1 vs B2B1)")
abline(a=0, b = 1)
dev.copy(png, "dump/scatter_fig11.png")
dev.off()

plot(rowMeans(x[,1:5]), rowMeans(x[,11:15]),
     main = "log2(raw_maqc) (B1-MeanA vs B2-MeanA)")


plot_arr <- rowMeans(x[,1:5])
ggplot(data.frame(plot_arr)) +
  geom_jitter(aes(x = plot_arr, y = 0))

head(cbind(x[,1:5], x[,11:15]), 100)

# Select rows with all positive numbers
x_nozero <- x[rowSums(x == 0) == 0,]

nozero <- plot_batch(x_nozero, batch_info, class_info)
gotzero <- plot_batch(x, batch_info, class_info)
nozero
gotzero

# EVALUATION --------------------------------------------------------------
filtered_maqc <- raw_maqc[filter_probesets(raw_maqc, 0.1),]
scaled_maqc <- norm.mean_scaling(filtered_maqc)
log_maqc <- log2_transform(scaled_maqc)

plot_batch(log_maqc, batch_info, class_info)
# ANOVA
batch_info <- rep(LETTERS[1:6], each = 10)
df <- log_maqc
pca_obj <- prcomp(t(df), center = T, scale. = T)
pca_arr <- pca_obj$x
pca_df <- data.frame(pca_arr)
eig_value <- (pca_obj$sdev)^2
var_pct <- eig_value/sum(eig_value)

# Argument: PC coordinates for single PC
calc.pc_var <- function(vec) {
  # Argument: ANOVA attributes; Calculates percentage of variance
  # SS_between/SS_between + SS_within
  calc.var_percentage <- function(vec) unname(vec[3]/(vec[3] + vec[4]))
  pc_metadata <- data.frame(pc = vec,
                            batch = as.factor(batch_info),
                            class = as.factor(class_info))
  batch_anova_attr <- unlist(summary(aov(pc~batch, data = pc_metadata)))
  class_anova_attr <- unlist(summary(aov(pc~class, data = pc_metadata)))
  return(c(calc.var_percentage(batch_anova_attr),
           calc.var_percentage(class_anova_attr)))
}
var_composition <- sapply(pca_df, calc.pc_var)
pca_var_df <- data.frame(var_pct,
                         batch_pct = var_composition[1,],
                         class_pct = var_composition[2,])
total_batch_pct <- sum(pca_var_df$var_pct * pca_var_df$batch_pct)
total_class_pct <- sum(pca_var_df$var_pct * pca_var_df$class_pct)

# ANOVA (Odd)
odd_batch_info
filtered_odd_maqc <- odd_maqc[filter_probesets(odd_maqc, 0.1),]
df <- filtered_odd_maqc
pca_obj <- prcomp(t(df), center = T, scale. = T)
pca_arr <- pca_obj$x
odd_single_pc_df <- data.frame(pc = pca_arr[,1],
                           batch = as.factor(odd_batch_info))
odd_single_pc_df
odd_aov_obj <- aov(pc~batch, data = odd_single_pc_df)
summary(odd_aov_obj)
dim(odd_single_pc_df)

# Explore total variance of data matrix
hist(data.matrix(log_maqc[10,]), breaks = 20)
# Checking that SS_total = SS_between + SS_within
vector <- unlist(log_maqc[10,])
ss_total <- sum((vector-mean(vector))^2)
list_batches <- split(vector, rep(1:6, each = 10))
list_groupmeans <- lapply(list_batches, mean)
vec_groupmeans <- do.call(c, list_groupmeans)
ss_between <- sum((vec_groupmeans - mean(vector))^2)*10
within_var <- mapply(function(vec, mean) sum((vec-mean)^2), list_batches, list_groupmeans)
ss_within <- sum(within_var)
