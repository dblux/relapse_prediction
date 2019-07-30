library(ggplot2)
library(cowplot)
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
  row_signedrank <- function(row_vec) {
    half_index <- length(row_vec)/2
    # Wilcoxon signed rank takes x-y
    wilcox_obj <- wilcox.test(row_vec[1:half_index],
                              row_vec[-(1:half_index)],
                              paired = T)
    return(wilcox_obj$p.value)
  }
  # Calculating significance of probeset using Wilcoxon signed rank test
  signedrank_pvalue <- apply(df, 1, row_signedrank)
  
  # # Mean of D8 - Mean of D0 for every probeset
  # effect_size <- rowMeans(df[,-(1:210)]) - rowMeans(df[,1:210])
  # sum(effect_size == 0)
  # effect_direction <- ifelse(effect_size >= 0, "up", "down")
  # 
  # # Selection of top 500 down-regulated probesets
  # probes_df <- data.frame(signedrank_pvalue, effect_direction)
  # sorted_probes <- probes_df[order(probes_df$signedrank_pvalue),]
  # downreg_probes <- rownames(sorted_probes)[sorted_probes$effect_direction == "down"][1:500]
  return(signedrank_pvalue)
}

# All dataframes have patients as rows and probesets/features as columns
calc_erm <- function(leukemia_df, normal_df, d0_d8_df) {
  # Calculation of centroids
  leukemia_centroid <- apply(leukemia_df, 2, median)
  normal_centroid <- apply(normal_df, 2, median)
  
  leukemia_normal <- normal_centroid - leukemia_centroid
  l2_norm <- function(vec) sqrt(sum(vec^2))
  erm_factor <- leukemia_normal/l2_norm(leukemia_normal)
  
  # subset_yeoh <- qnorm_yeoh[downreg_probes,]
  # # Tranpose in order to have patients as rows and probesets as columns
  # transposed_yeoh <- t(subset_yeoh)
  # Assume that patients from rows 1:210 match together with rows 211:420
  print(dim(d0_d8_df))
  patient_arr <- cbind(d0_d8_df[1:210,], d0_d8_df[-(1:210),])
  # dim(patient_arr) = c(210, 1000)
  # Calculate vector by: D8-D0
  calc_d0_d8 <- function(row_vec) {
    dim_x <- length(row_vec)/2
    return(row_vec[-(1:dim_x)] - row_vec[1:dim_x])
  }
  d0_d8_vec_vstack <- apply(patient_arr, 1, calc_d0_d8)
  rownames(d0_d8_vec_vstack)
  # Multiplication of erm_factor is propagated through every column
  erm <- colSums(d0_d8_vec_vstack * erm_factor)
  return(erm)
}

# Plots ROC and calculates AUC in a primitive fashion (i.e. ROC is step function)
# Does not resolve ties in the score
# Assumption: Lower score the better, ROC is step function
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
          col = color, lwd = 2)
    return(AUC)
  }
  # Initialise plot
  color_index <- (1:length(score_list)) + 1
  plot(NULL,
       xlab = "1 - Specificity", ylab = "Sensitivity",
       xlim = c(0,1), ylim = c(0,1),
       xaxs = "i", yaxs = "i")
  abline(0, 1, lty = 5)
  auc_vec <- mapply(ROC_AUC, score_list, is_bigger_better_vec,
                    color_index, SIMPLIFY = T)
  # If name_vec is not NULL display legend
  if (!is.null(name_vec)) {
    format_text <- function(name, auc) {
      return(sprintf("%s (%.3f)", name, auc))
    }
    legend_text <- mapply(format_text, name_vec, auc_vec)
    legend("bottomright", inset = 0.03, lty = 1, lwd = 2,
           legend = legend_text, col = color_index)
  }
  return(auc_vec)
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
log_yeoh <- log2_transform(filtered_yeoh)

# MILE data
filtered_mile <- scaled_mile[selected_probes,]
log_mile <- log2_transform(filtered_mile)

# ### Quantile normalisation of both datasets
# qnorm_data <- norm_quantile(cbind(log_yeoh, log_mile))
# qnorm_yeoh <- qnorm_data[,1:420]
# qnorm_mile <- qnorm_data[,-(1:420)]

# Quantile normalise before log transform to calculate logfc
qnorm_data <- norm_quantile(cbind(filtered_yeoh, filtered_mile))
qnorm_yeoh <- qnorm_data[,1:420]
qnorm_mile <- qnorm_data[,-(1:420)]

# Selecting drug responsive genes
signedrank_pvalue <- select_probes(qnorm_yeoh)
log_fc <- calc_logfc(qnorm_yeoh[,1:210], qnorm_yeoh[,-(1:210)])
pvalue_probesets <- names(signedrank_pvalue)[signedrank_pvalue <= 0.0001]
fc_probesets <- names(log_fc)[log_fc > 1]
intersect_probesets <- fc_probesets[fc_probesets %in% pvalue_probesets]
print(length(intersect_probesets))

log_yeoh <- log2_transform(qnorm_yeoh[intersect_probesets,])
log_mile <- log2_transform(qnorm_mile[intersect_probesets,])

transposed_df <- t(cbind(log_yeoh, log_mile))
rownames(transposed_df)
pca_obj <- prcomp(transposed_df, center = T, scale. = T)
pca_arr <- pca_obj$x[,1:4]
response_df <- pca_arr[1:420,]
leukemia_df <- pca_arr[421:1170,]
normal_df <- pca_arr[1171:1244,]
erm <- calc_erm(leukemia_df, normal_df, response_df)
rownames(pca_arr)

# YEOH --------------------------------------------------------------------
# ORDERING OF YEOH DATA IS VERY IMPORTANT
# Check that columns are paired
# Note that all D0 data is before D8 data
colnames(qnorm_yeoh)

# Selection of down-regulated probes
downreg_probes <- select_probes(qnorm_yeoh)
# Check if all down-regulated probesets are present in MILE data
sum(downreg_probes %in% rownames(qnorm_mile))
leukemia_mile <- qnorm_mile[downreg_probes, 1:750]
normal_mile <- qnorm_mile[downreg_probes, -(1:750)]
yeoh_500_df <- qnorm_yeoh[downreg_probes,]
erm <- calc_erm(t(leukemia_mile), t(normal_mile), t(yeoh_500_df))

# GFS ---------------------------------------------------------------------
### GFS normalisation and selection of probesets
gfs_yeoh <- norm_gfs(filtered_yeoh)
pca_gfs_yeoh <- plot_pca(gfs_yeoh)
ggsave("dump/pca-yeoh_gfs.pdf", pca_gfs_yeoh,
       width = 12, height = 8)

# Selection of probesets with most variance
probeset_var <- apply(gfs_yeoh, 1, var)
top_probesets <- names(sort(probeset_var, decreasing = T)[1:500])
subset_gfs_yeoh <- gfs_yeoh[top_probesets,]
pca_subset_gfs <- plot_pca(subset_gfs_yeoh)
ggsave("dump/pca-yeoh_gfs_subset.pdf", pca_subset_gfs,
       width = 12, height = 8)

gfs_mile <- norm_gfs(log_mile)
subset_gfs_mile <- gfs_mile[top_probesets,]

# PCA of all data (centroid, training, test)
all_pca_obj <- prcomp(t(cbind(subset_gfs_yeoh, subset_gfs_mile)))
all_pca_df <- all_pca_obj$x

# Top 3 PCs
pca_mapping <- all_pca_df[,1:3]
rownames(pca_mapping[421:1170,])
leukemia_df <- pca_mapping[421:1170,]
normal_df <- pca_mapping[1171:1244,]
d0_d8_df <- pca_mapping[1:420,]
erm <- calc_erm(leukemia_df, normal_df, d0_d8_df)

# PCA PLOT ----------------------------------------------------------------
# # Batch information: Associating colour
# num_subtype <- table(substring(colnames(mile_data), 1, 1))
# red_palette <- brewer.pal(9, "Reds")[2:9]
# mile_palette <- c(red_palette, "darkolivegreen3")
# colour_code <- rep(mile_palette, num_subtype)
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
# Dataframe to be visualised
visualise_df <- rbind(response_df,
                      apply(leukemia_df, 2, median),
                      apply(normal_df, 2, median))

timepoint_vec <- c(rep(c(19, 17), each = 210), 19, 19)

# 2D PCA plot
pca_2d <- function(df, colour_code) {
  pc1_pc2 <- ggplot(df, aes(x = PC1, y = PC2)) +
    geom_point(size = 2, col = colour_code, show.legend = F)
  pc2_pc3 <- ggplot(df, aes(x = PC2, y = PC3)) +
    geom_point(size = 2, col = colour_code, show.legend = F)
  pc1_pc3 <- ggplot(df, aes(x = PC1, y = PC3)) +
    geom_point(size = 2, col = colour_code, show.legend = F)
  pc3_pc4 <- ggplot(df, aes(x = PC3, y = PC4)) +
    geom_point(size = 2, col = colour_code, show.legend = F)
  multiplot <- plot_grid(pc1_pc2, pc2_pc3, pc1_pc3, pc3_pc4,
                         ncol = 2, nrow = 2)
  return(multiplot)
}

plot_pca <- pca_2d(as.data.frame(pca_arr), all_colour)
plot_pca
ggsave("dump/pca-qnorm_select_pca3_centroids_batch.pdf", plot_pca,
       width = 9, height = 9)

# Visualising vectors
arrows_df <- cbind(response_df[1:210,], response_df[-(1:210),])
colnames(arrows_df) <- paste(colnames(arrows_df),
                             rep(LETTERS[1:2], each = 4),
                             sep = "_")
centroid_df <- rbind(apply(leukemia_df, 2, median),
                     apply(normal_df, 2, median))

plot_vectors <- function(df, centroid_df) {
  pc1_pc2 <- ggplot(data = df) +
    geom_point(aes(x = PC1_A, y = PC2_A), size = 3,
               col = "steelblue4", show.legend = F) +
    geom_point(aes(x = PC1_B, y = PC2_B), size = 3,
               col = "turquoise3", show.legend = F) +
    geom_point(data = centroid_df, aes(x = PC1, y = PC2),
               size = 4, colour = c("tomato3", "darkolivegreen3")) +
    geom_segment(aes(x = PC1_A, y = PC2_A,
                     xend = PC1_B, yend = PC2_B),
                 arrow = arrow(length = unit(0.3, "cm")),
                 alpha = 0.8)
  pc2_pc3 <- ggplot(data = df) +
    geom_point(aes(x = PC2_A, y = PC3_A), size = 3,
               col = "steelblue4", show.legend = F) +
    geom_point(aes(x = PC2_B, y = PC3_B), size = 3,
               col = "turquoise3", show.legend = F) +
    geom_point(data = centroid_df, aes(x = PC2, y = PC3),
               size = 4, colour = c("tomato3", "darkolivegreen3")) +
    geom_segment(aes(x = PC2_A, y = PC3_A,
                     xend = PC2_B, yend = PC3_B),
                 arrow = arrow(length = unit(0.3, "cm")),
                 alpha = 0.8)
  multiplot <- plot_grid(pc1_pc2, pc2_pc3, ncol = 2)
}

vectors_plot <- plot_vectors(as.data.frame(arrows_df), as.data.frame(centroid_df))
vectors_plot
ggsave("dump/vectors-qnorm_select_pc3.pdf", vectors_plot,
       width = 12, height = 6)

# RESULTS -----------------------------------------------------------------
### Extracting truth labels and training/test
labels_yeoh <- yeoh_metadata[names(erm), c(6, 8, 12)]
results_df <- cbind(erm, labels_yeoh)
# Change truth labels with value 2 to 1
results_df[results_df$event_code == 2, 4] <- 1
results_df <- results_df[order(results_df$erm),]

write.table(results_df, "dump/qnorm_select_pca3.tsv",
            quote = F, sep = "\t", row.names = T, col.names = T)
# Investigate
results_list <- split(results_df, results_df[,2])
results_list[[1]]

ROC_TITLE <- "Qnorm: 500 down_reg"
ROC_WPATH <- "dump/qnorm_500down_roc.pdf"

plot_roc(list(results_df$erm),
         label_vec = results_df$event_code,
         name_vec = ROC_TITLE)
results_roc <- recordPlot()
save_fig(results_roc, ROC_WPATH,
         width = 10, height = 10)

# RETROSPECTIVE -----------------------------------------------------------
results_df2 <- read.table("dump/qnorm_downreg.tsv",
                          sep = "\t", header = T, row.names = 1)
yeoh_metadata <- read.table("data/GSE67684/processed/metadata_batch.tsv",
                            sep = "\t", header = T, row.names = 1)

results_list2 <- split(results_df2, results_df2[,2])

for (i in 1:length(results_list)) {
  print(results_list[[i]])
  print(results_list2[[i]])
}
