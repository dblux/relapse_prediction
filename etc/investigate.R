library(ggplot2)
library(cowplot)
library(rgl)
library(reshape2)
library(RColorBrewer)
source("../functions.R")

mile_data <- read.table("data/GSE13204/processed/mas5_ordered.tsv",
                        sep = "\t", header = T, row.names = 1)
yeoh_data <- read.table("data/GSE67684/processed/mas5_ordered.tsv",
                        sep = "\t", header = T, row.names = 1)

# Batch information: Associating colour
num_subtype <- table(substring(colnames(mile_data), 1, 1))
red_palette <- brewer.pal(9, "Reds")[2:9]
mile_palette <- c(red_palette, "darkolivegreen3")
colour_code <- rep(mile_palette, num_subtype)

# PCA: MILE data
selected_probes <- apply(mile_data, 1, var) != 0
sum(!col_logical)
transposed_df <- t(mile_data[selected_probes,])
pca_obj <- prcomp(transposed_df, center = T, scale. = T)
pca_arr <- pca_obj$x[,1:4]

pca_3d <- function(df, batch_info) {
  pca_obj <- prcomp(t(df))
  pca_df <- as.data.frame(pca_obj$x[,1:5])
  # RGL plot parameters
  rgl.open()
  rgl.bg(color="white")
  rgl.viewpoint(theta = 110, phi = 5, zoom = 0.8)
  par3d(windowRect = c(50, 20, 500, 500))
  aspect3d(1,1,1)
  # Plot of MILE dataset
  with(df, plot3d(PC1, PC2, PC3, col = colour_code,
                  type = "p", size = 5))
}

# PCA: yeoh_data
selected_data <- yeoh_data[selected_probes,]
nrow(selected_data)
allzero_logical <- !(apply(selected_data, 1, var) != 0)
sum(allzero_logical)
x <- selected_data[allzero_logical,]

# S3 method of class: predict.prcomp
# Projection into PCA space of MILE data
projected_arr <- predict(pca_obj, t(selected_data))
getS3method("predict", "prcomp")

incremental_df <- as.data.frame(rbind(pca_arr, projected_arr[,1:4]))

# Colour code
colour_code1 <- c(colour_code, rep(c("steelblue4", "turquoise3"), c(210, 210)))

with(incremental_df, plot3d(PC1, PC2, PC3,
                            col = colour_code1,
                            type = "p", size = 1))

rgl.postscript("dump/scatter-yeoh_projected_mile_s.pdf", "pdf")

# PCA of both MILE and yeoh_data
combined_data <- cbind(mile_data, yeoh_data)
selected_probes1 <- apply(combined_data, 1, var) != 0
pca_obj1 <- prcomp(t(combined_data[selected_probes1,]))
pca_df1 <- as.data.frame(pca_obj1$x[,1:4])

with(pca_df1, plot3d(PC1, PC2, PC3,
                     col = colour_code1,
                     type = "p", size = 5))
rgl.postscript("dump/scatter-yeoh_mile_all.pdf", "pdf")

# 2D PCA plot
pca_2d <- function(df, colour_code) {
  pc1_pc2 <- ggplot(df, aes(x = PC1, y = PC2)) +
    geom_point(size = 2, col = colour_code, show.legend = F)
  pc2_pc3 <- ggplot(df, aes(x = PC2, y = PC3)) +
    geom_point(size = 2, col = colour_code, show.legend = F)
  pc3_pc4 <- ggplot(df, aes(x = PC3, y = PC4)) +
    geom_point(size = 2, col = colour_code, show.legend = F)
  multiplot <- plot_grid(pc1_pc2, pc2_pc3, pc3_pc4, ncol = 3)
  return(multiplot)
}

multiplot_all <- pca_2d(pca_df1, colour_code1)
multiplot_inc <- pca_2d(incremental_df, colour_code1)
ggsave("dump/pca_all.pdf", multiplot_all,
       height = 4, width = 12)
ggsave("dump/pca_incremental.pdf", multiplot_inc,
       height = 4, width = 12)


# BATCH EFFECTS -----------------------------------------------------------
yeoh_data <- read.table("data/GSE67684/processed/mas5_ordered.tsv",
                        sep = "\t", header = T, row.names = 1)
yeoh_metadata <- read.table("data/GSE67684/processed/metadata_batch.tsv",
                            sep = "\t", header = T, row.names = 1)

mile_data <- read.table("data/GSE13204/processed/mas5_ordered.tsv",
                        sep = "\t", header = T, row.names = 1)
mile_metadata <- read.table("data/GSE13204/processed/metadata.tsv",
                            sep = "\t", header = T, row.names = 1)

# Environment: yeoh_metadata
plot_pca <- function(df) {
  pca_obj <- prcomp(t(df))
  pca_df <- as.data.frame(pca_obj$x[,1:6])
  eig_value <- (pca_obj$sdev)^2
  var_pc <- eig_value[1:5]/sum(eig_value)
  pc_labels <- sprintf("PC%d (%.2f%%)", 1:5, var_pc*100)
  # Create plot_df by concatenating metadata labels to data
  plot_df <- cbind(yeoh_metadata[rownames(pca_df), c(3,6)], pca_df)
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
  pc3_pc4 <- ggplot(plot_df, aes(x = PC3, y = PC4)) +
    geom_point(aes(col = factor(batch), shape = factor(time_point)),
               size = 3, show.legend = F) +
    xlab(pc_labels[3]) + ylab(pc_labels[4])
  pc1_pc4 <- ggplot(plot_df, aes(x = PC1, y = PC4)) +
    geom_point(aes(col = factor(batch), shape = factor(time_point)),
               size = 3, show.legend = F) +
    xlab(pc_labels[1]) + ylab(pc_labels[4])
  pc4_pc5 <- ggplot(plot_df, aes(x = PC4, y = PC5)) +
    geom_point(aes(col = factor(batch), shape = factor(time_point)),
               size = 3, show.legend = F) +
    xlab(pc_labels[4]) + ylab(pc_labels[5])
  multiplot <- plot_grid(pc1_pc2, pc2_pc3, pc1_pc3, pc3_pc4, pc4_pc5, pc1_pc4,
                         ncol = 3, nrow = 2)
  return(multiplot)
}

# Environment: yeoh_metadata
plot_pca1 <- function(df) {
  pca_obj <- prcomp(t(df))
  pca_df <- as.data.frame(pca_obj$x[,1:6])
  # Create plot_df by concatenating metadata labels to data
  plot_df <- cbind(mile_metadata[rownames(pca_df), 2, drop = F], pca_df)
  pc1_pc2 <- ggplot(plot_df, aes(x = PC1, y = PC2)) +
    geom_point(aes(col = factor(subtype)),
               size = 3, show.legend = F)
  pc2_pc3 <- ggplot(plot_df, aes(x = PC2, y = PC3)) +
    geom_point(aes(col = factor(subtype)),
               size = 3, show.legend = F)
  pc1_pc3 <- ggplot(plot_df, aes(x = PC1, y = PC3)) +
    geom_point(aes(col = factor(subtype)),
               size = 3, show.legend = F)
  pc3_pc4 <- ggplot(plot_df, aes(x = PC3, y = PC4)) +
    geom_point(aes(col = factor(subtype)),
               size = 3, show.legend = F)
  pc1_pc4 <- ggplot(plot_df, aes(x = PC1, y = PC4)) +
    geom_point(aes(col = factor(subtype)),
               size = 3, show.legend = F)
  pc4_pc5 <- ggplot(plot_df, aes(x = PC4, y = PC5)) +
    geom_point(aes(col = factor(subtype)),
               size = 3, show.legend = F)
  multiplot <- plot_grid(pc1_pc2, pc2_pc3, pc1_pc3, pc3_pc4, pc4_pc5, pc1_pc4,
                         ncol = 3, nrow = 2)
  return(multiplot)
}

filtered_mile <- filter_probesets(mile_data)
log_mile <- log2_transform(filtered_mile)
qnorm_mile <- norm_quantile(log_mile)
gfs_mile <- norm_gfs(log_mile)
pca_mile_gfs <- plot_pca1(gfs_mile)

ggsave("dump/pca-mile_gfs.pdf", pca_mile_gfs,
       width = 12, height = 8)

# Log transform data
log_data <- log2_transform(filtered_yeoh)
qnorm_data <- norm_quantile(log_data)
gfs_data <- norm_gfs(log_data)

pca_qnorm <- plot_pca(gfs_data)
ggsave("dump/pca-yeoh_gfs.pdf", pca_qnorm,
       width = 12, height = 8)

# Assumes that dataframe has been log-transformed
plot_eval <- function(df, batch_info1, batch_info2) {
  df <- log2_transform(df)
  # Melt dataframe
  melt_df <- melt(df, variable.name = "ID")
  
  # Plot density curve
  pdf <- ggplot(melt_df, aes(x = value, col = ID)) +
    geom_density(show.legend = F) +
    facet_wrap(rep(factor(batch_info2), nrow(df))) +
    xlim(0, 16)
  
  # Mean probe intensities for each chip
  mean_tibble <- melt_df %>% group_by(ID) %>%
    summarise(mean = mean(value))
  mean_scatter <- ggplot(mean_tibble, aes(x = ID, y = mean)) +
    geom_point(aes(col = factor(batch_info1),
                   shape = factor(batch_info2)),
               show.legend = F, size = 3) + 
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  return(list(mean_scatter, pdf))
}

# Retrieve batch info from metadata
head(yeoh_metadata)
batch_timepoint <- yeoh_metadata[colnames(yeoh_data), 3]
batch_assigned <- yeoh_metadata[colnames(yeoh_data), 6]
yeoh_before <- plot_eval(filtered_yeoh, batch_assigned, batch_timepoint)
yeoh_before[[2]]
ggsave("dump/eval.pdf", yeoh_before)

x <- data.frame(a = rep(5, 100), b = rep(10,100))
melt_x <- melt(x)

