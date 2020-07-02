# library(Biocomb)

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

# library(sva)  # ComBat
# library(scran)  # MNN
# library(Harman)

# library(xtable)
# library(gPCA)

source("../functions.R")
source("bin/bcm.R")

theme_set(theme_dark())
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

# Plot PCA 3D: Batch effects
# Plot batches in different colours and classes in different shapes
plotPCA3DBatchEffects <- function(df1, metadata_df) {
  # Batch and class annotations
  batch_factor <- metadata_df[colnames(df1), "batch"]
  batch_palette <- generateGgplotColours(length(unique(batch_factor)))
  batch_colour <- batch_palette[batch_factor]
  
  class_factor <- metadata_df[colnames(df1), "celltype"]
  all_pch <- 21:24 # TODO: See number of batches and celltypes
  # Error if there are more classes than pch symbols (> 5)
  stopifnot(length(unique(class_factor)) <= 5)
  class_pch <- all_pch[class_factor]
  plotPCA3D(df1, batch_colour, class_pch)
}

plotHeatmapSubtype <- function(X, metadata_df) {
  subtype_factor <-  [colnames(X), "subtype"]
  set3_pal <- brewer.pal(9, "Set3")
  subtype_col <- set3_pal[subtype_factor]
  
  par(mar = c(1,1,1,1))
  heatmap(data.matrix(X),
          col = brewer.pal(9, "Blues"),
          ColSideColors = subtype_col,
          scale = "none",
          labRow = NA, labCol = NA)
  
  legend(x = -.04, y = 1350, legend = levels(subtype_factor),
         col = set3_pal[factor(levels(subtype_factor))],
         pch = 15, cex = .7)
  heatmap_subtype <- recordPlot()
  par(mar = c(5.1, 4.1, 4.1, 2.1)) # Reset to defaults
  return(heatmap_subtype)
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

# STRATIFICATION ----------------------------------------------------------
# Provide all available labels to algorithm
head(metadata_df)

# B1 vs B2
meta_b1b2<- metadata_df[
  (metadata_df$batch_info == 1 | metadata_df$batch_info == 2) &
  metadata_df$class_info == "D0",]

list_batch <- split(meta_b1b2, meta_b1b2$batch_info)
lapply(list_batch, unique)

# Dataset 1 ---------------------------------------------------------------
library(Rtsne)

DATA1_WPATH <- "../scrna_seq/data/tran_2020/dataset1/dataset1_sm_uc3.txt"
META1_WPATH <- "../scrna_seq/data/tran_2020/dataset1/sample_sm_uc3.txt"
raw1 <- read.table(DATA1_WPATH, header = T, sep = "\t")
meta1 <- read.table(META1_WPATH, header = T, sep = "\t")

# Filter probes with too many zeros
#' @param df dataframe
#' @param percent_threshold percentage threshold of non-zero values
#' @param metadata_df df containing class labels of samples
#' @param logical_func a function that is either "all" or "any". (Either all or
#' just one class have to pass the threshold)
#' @return dataframe containing rows that meet threshold of non-zero
#' values
filterSparseRows <- function(df1, percent_threshold, metadata_df = NULL,
                            logical_func = any) {
  if (is.null(metadata_df)) {
    logical_df <- df1 != 0
    selected_logvec <- rowSums(logical_df) > percent_threshold * ncol(df1)
    print(paste("No. of probesets removed =",
                nrow(df1) - sum(selected_logvec)))
    return(df1[selected_logvec,])
  } else {
    class_factor <- metadata_df[colnames(df1), "celltype"]
    logical_df <- data.frame(df1 != 0)
    list_logical_df <- split.default(logical_df, class_factor)
    list_logvec <- lapply(
      list_logical_df,
      function(df1) rowSums(df1) > (percent_threshold * ncol(df1))
    )
    combined_log_df <- do.call(cbind,list_logvec)
    print(head(combined_log_df))
    selected_logvec <- apply(combined_log_df, 1, logical_func)
    # selected_logvec <- do.call(mapply, c(logical_func, list_logvec))
    print(paste("No. of probesets removed =",
                nrow(df1) - sum(selected_logvec)))
    return(df1[selected_logvec,])
  }
}
filtered1 <- filterSparseRows(raw1, 0.5, meta1)
log1 <- log2_transform(filtered1)
log1[log1 < 0] <- 0

# Plot - Individual feature
i <- 83
D <- data.frame(value = t(log1[i,]),
                celltype = meta1[colnames(log1), "celltype"],
                batch = meta1[colnames(log1), "batch"])
colnames(D)[1] <- "value"
single_gene <- ggplot(D, aes(x=celltype, y=value, colour=batch)) +
  geom_point(position = position_jitterdodge(), cex=3, show.legend = F)
ggsave("~/Dropbox/temp/tran1-single_gene.pdf", single_gene,
       width = 8, height = 5)

# Plot uncorrected data
set.seed(0)
log1_tsne_obj <- Rtsne(t(log1), perplexity = 30, theta = 0.0)
log1_tsne <- data.frame(log1_tsne_obj$Y,
                        celltype = meta1[colnames(log1), "celltype"],
                        batch = meta1[colnames(log1), "batch"])
colnames(log1_tsne)[1:2] <- c("TSNE1", "TSNE2")
tsne_filtered <- ggplot(log1_tsne, aes(x=TSNE1, y=TSNE2, colour=celltype, pch=batch)) +
  geom_point(cex=3, alpha = .7, show.legend = F) +
  scale_shape_manual(values=15:19)
ggsave("~/Dropbox/temp/tran1-tsne_log.pdf", tsne_filtered,
       width = 6, height = 6)

log1_pca_obj <- prcomp(t(log1))
log1_pca <- data.frame(log1_pca_obj$x[,1:2],
                       celltype = meta1[colnames(log1), "celltype"],
                       batch = meta1[colnames(log1), "batch"])

ggplot(log1_pca, aes(x=PC1, y=PC2, colour=celltype, pch=batch)) +
  geom_point(cex=3, alpha = .7, show.legend = F) +
  scale_shape_manual(values=15:19)

# Batch correction
batch1 <- meta1[colnames(log1), "batch"]
class1 <- meta1[colnames(log1), "celltype", drop=F]

bcm1 <- correctMultiBCM(log1, batch1, class1, ref_batch = "Batch1")

# Plots
set.seed(0)
bcm1_tsne_obj <- Rtsne(t(bcm1), perplexity = 30, theta = 0.0)
bcm1_tsne <- data.frame(bcm1_tsne_obj$Y,
                        celltype = meta1[colnames(bcm1), "celltype"],
                        batch = meta1[colnames(bcm1), "batch"])
colnames(bcm1_tsne)[1:2] <- c("TSNE1", "TSNE2")
tsne_bcm1 <- ggplot(bcm1_tsne, aes(x=TSNE1, y=TSNE2, colour=celltype, pch=batch)) +
  geom_point(cex=3, alpha=.7, show.legend = F) +
  scale_shape_manual(values=15:19)
tsne_bcm1
ggsave("~/Dropbox/temp/tran1-tsne_bcm.pdf", tsne_bcm1,
       width = 6, height = 6)
