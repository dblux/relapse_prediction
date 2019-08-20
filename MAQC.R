library(reshape2)
library(ggplot2)
library(cowplot)
library(scran)
library(sva)
source("../functions.R")

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

# Visualise ---------------------------------------------------------------
# No normalisation
plot_before <- plot_batch(raw_maqc, batch_info, shape_info)
plot_before
ggsave("dump/plot_before.pdf", plot_before,
       width = 6, height = 6)

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

# ComBat ------------------------------------------------------------------
# Creation of pData dataframe (metadata)
class <- as.factor(substring(colnames(raw_maqc), 7, 7))
batch <- as.factor(substring(colnames(raw_maqc), 5, 5))
maqc_metadata <- data.frame(shape_info, batch_info)
# Rownames of metadata are same as colnames of data df
rownames(maqc_metadata) <- colnames(raw_maqc)

# ComBat assumes that data has been normalised and probesets have been filtered
# Error if probesets are not filtered as rows have 0 variance
scaled_maqc <- norm.mean_scaling(raw_maqc)
selected_probetsets <- filter_probesets(scaled_maqc, 0.2)
filtered_maqc <- scaled_maqc[selected_probetsets,]

model_combat <- model.matrix(~1, data = maqc_metadata)
combat_maqc <- ComBat(data.matrix(filtered_maqc), maqc_metadata$batch_info, model_combat)

plot_combat <- plot_batch(data.frame(combat_maqc), batch_info, shape_info)
ggsave("dump/plot_combat.pdf", plot_combat,
       width = 6, height = 6)
