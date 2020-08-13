library(reshape2)
library(dplyr)

## Plotting
library(ggplot2)
library(cowplot)
library(rgl)
library(RColorBrewer)
library(pheatmap)
library(UpSetR)
library(VennDiagram)
library(xtable)
library(Rtsne)
# library(dendextend)

library(sva)

## Custom
source("../functions.R")

theme_set(theme_gray())

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
calcERM <- function(response_df, normal_df) {
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
  ## ERM1 / projection of D0-N on L-N
  erm1_ratio1 <- erm1/d0_normal_proj
  
  d8_normal_vstack <- normal_centroid - t(d8_df)
  ### D8-Normal projection ###
  d8_normal_proj <- colSums(d8_normal_vstack * unit_leuk_normal)
  
  stopifnot(identical(names(erm1), names(erm1_ratio1)))
  
  # Calculate vstack of unit D0-Normal vectors
  l2norm_d0_normal <- apply(d0_normal_vstack, 2, calcL2Norm)
  unit_d0_normal_vstack <- sweep(d0_normal_vstack, 2, l2norm_d0_normal, "/")
  
  ### ERM2 ###
  ## Projection of D0-D8 on D0-N
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
  
  ### L2-norm between D8 and Normal ###
  l2norm_d8_normal <- apply(d8_normal_vstack, 2, calcL2Norm)
  
  ### L2-norm ratios
  l2norm_ratio1 <- l2norm_d0_d8/l2norm_d0_normal
  l2norm_ratio2 <- l2norm_d0_d8/l2norm_d8_normal
  l2norm_diff <- l2norm_d0_normal - l2norm_d8_normal
  l2norm_diff_ratio <- l2norm_diff/l2norm_d0_d8
  
  ### Ratios
  erm1_ratio2 <- erm1/abs(d8_normal_proj)
  erm1_ratio3 <- erm1/l2norm_d0_d8
  
  ### Concatenate all features ###
  features_df <- data.frame(
    erm1, erm1_ratio1, erm2, erm2_ratio, erm3, erm3_ratio,
    d0_normal_proj, d8_normal_proj, l2norm_d0_d8,
    diff_l2norm, angle_d0_d8, angle_d0d8_normal,
    angle_d0_normal, angle_d8_normal,
    l2norm_d0_normal, l2norm_d8_normal,
    l2norm_ratio1, l2norm_ratio2,
    l2norm_diff, l2norm_diff_ratio,
    erm1_ratio2, erm1_ratio3
  )
  return(features_df)
}

# 3D PCA plot
plotPCA3D <- function(df, colour, pch, pc_labels = NULL,
                      ratio_list = list(2,1,1)) {
  if (is.null(pc_labels)) {
    print("PCA performed!")
    pca_obj <- prcomp(t(df), center = T, scale. = F)
    pca_df <- as.data.frame(pca_obj$x[,1:3])
    eigenvalues <- (pca_obj$sdev)^2
    var_pc <- eigenvalues[1:3]/sum(eigenvalues)
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
  with(pca_df, pch3d(PC1, PC2, PC3, bg = colour,
                     pch = pch, cex = 0.5, lwd = 1.5))
  box3d(col = "black")
  title3d(xlab = pc_labels[1], ylab = pc_labels[2],
          zlab = pc_labels[3], col = "black")
  # Plot aspect ratios of axis according to variance
  do.call(aspect3d, ratio_list)
}

# Plot PCA before selecting features
# Batch information of all the timepoints
plotPCA3DYeoh <- function(df1, metadata_df) {
  batch_info <- metadata_df[colnames(df1), "batch_info"]
  generate_colour <- colorRampPalette(c("lightblue", "darkblue"))
  batch_palette <- generate_colour(10)
  # batch_palette <- brewer.pal(10, "Set3")
  batch_colour <- batch_palette[batch_info]
  # Shape of all timepoints
  class_info <- metadata_df[colnames(df1), "class_info"]
  print(levels(class_info))
  levels(class_info) <- 21:23
  timepoint_shape <- as.numeric(as.character(class_info))
  plotPCA3D(df1, batch_colour, timepoint_shape)
}

# Plot PCA before selecting features
# Batch information of all the timepoints
plotPCA3DYeoh1 <- function(df1, metadata_df) {
  batch_info <- metadata_df[colnames(df1), "batch_info"]
  batch_factor <- droplevels(as.factor(batch_info))
  print(batch_factor)
  print(levels(batch_factor))
  levels(batch_factor) <- 21:22
  pch <- as.numeric(as.character(batch_factor))
  # generate_colour <- colorRampPalette(c("lightblue", "darkblue"))
  # batch_palette <- generate_colour(10)
  
  # Shape of all timepoints
  class_info <- metadata_df[colnames(df1), "subtype"]
  palette <- brewer.pal(10, "Set3")
  col <- palette[class_info]
  
  plotPCA3D(df1, col, pch)
}

plotJitterYeoh <- function(X, metadata_df, n_pc = 10) {
  pca_obj <- prcomp(t(X))
  X_pca <- data.frame(pca_obj$x)
  batch <- as.factor(metadata_df[rownames(X_pca),"batch_info"])
  class <- as.factor(metadata_df[rownames(X_pca),"class_info"])
  X_meta <- cbind(batch, class, X_pca[,1:n_pc])
  X_long <- melt(X_meta, id = c("batch", "class"), variable.name="PC")
  
  ax_batch <- ggplot(X_long, aes(x=PC, y=value)) +
    # geom_boxplot(aes(fill=batch), alpha=0.3, outlier.shape=NA) +
    geom_point(aes(colour=batch), position=position_jitterdodge(),
               size = 1, alpha = 1.0)
  
  ax_class <- ggplot(X_long, aes(x=PC, y=value)) +
    # geom_boxplot(aes(fill=class), alpha=0.3, outlier.shape=NA) +
    geom_point(aes(colour=class), position=position_jitterdodge(),
               size = 1, alpha = 1.0)
  
  fig <- plot_grid(ax_batch, ax_class, nrow = 2)
  return(fig)  
}

plotPrediction <- function(results, metadata_df, yeoh_label) {
  y <- as.factor(metadata_df[rownames(results),"label"])
  features1 <- results[,c("erm1_ratio2", "l2norm_ratio2"), drop=F]
  features2 <- results[, "angle_d0d8_normal", drop=F]
  
  # D33 MRD
  pid_idx <- substr(rownames(results), 1, 4)
  d33_mrd <- yeoh_label[pid_idx, "d33_mrd"]
  mrd_rank <- rank(d33_mrd, na.last = T, ties.method="min")
  mrd_percent <- (mrd_rank-1)/sum(!is.na(d33_mrd))
  mrd_percent[is.na(d33_mrd)] <- NA
  
  # Two different ways of rankings
  features_rankdesc <- apply(-features1, 2, rank, ties.method="min")
  features_percentdesc <- (features_rankdesc-1)/nrow(features1)
  features_rankasc <- apply(features2, 2, rank, ties.method="min")
  features_percentasc <- (features_rankasc-1)/nrow(features2)
  features_percent <- cbind(features_percentdesc, features_percentasc,
                            mrd_percent)
  
  if (sum(is.na(features_percentasc), is.na(features_percentdesc)) > 0)
    warning("Features contain NA values")
  
  # Calculate p for features and avg_p
  avg_percent <- rowMeans(features_percent, na.rm = T)
  avgpercent_y <- data.frame(p = avg_percent, label = y)
  percent_y <- cbind(pid = rownames(features_percent),
                     avgpercent_y, features_percent)
  long_percent_y <- melt(percent_y, id = c("pid", "label"),
                         variable.name = "feature")
  
  pid_idx <- substr(rownames(avgpercent_y), 1, 4)
  avgpercent_mrd <- cbind(avgpercent_y,
                          d33_mrd = -log10(yeoh_label[pid_idx, "d33_mrd"]))
  
  # Features and avg_p
  features_y <- data.frame(features1, features2, p = avg_percent, label = y)
  long_features_y <- melt(features_y, id="label", variable.name = "feature")
  
  # Calculating relative risk
  sort_avgp <- avgpercent_y[order(avgpercent_y$p),]
  sort_avgp$label <- as.numeric(as.character(sort_avgp$label))
  sort_avgp$total_le <- rank(sort_avgp$p, ties.method = "max")
  sort_avgp$total_g <- nrow(sort_avgp) - sort_avgp$total_le
  sort_avgp$relapse_le <- sapply(sort_avgp$total_le,
                                 function(i) sum(sort_avgp$label[1:i]))
  sort_avgp$relapse_g <- sum(sort_avgp$label) - sort_avgp$relapse_le
  sort_avgp <- within(sort_avgp,
                      relative_risk <-
                        (relapse_le/total_le)/(relapse_g/total_g))
  sort_avgp <- within(sort_avgp,
                      odds_ratio <-
                        (relapse_le/(total_le-relapse_le))/
                        (relapse_g/(total_g-relapse_g)))
  
  # PLOT: FEATURES
  jitter_features <- ggplot(long_features_y,
                             aes(feature, value, colour = label)) +
    geom_point(position = position_jitterdodge(), cex = 3, show.legend = F) +
    scale_color_manual(values = c("darkolivegreen3", "tomato3")) +
    facet_wrap(~feature, nrow = 1, scales = "free") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_text(angle = 10, vjust = 0.5))

  emp_cdf <- ggplot(avgpercent_y, aes(x = p, colour = label)) +
    stat_ecdf(show.legend = F) +
    scale_color_manual(values = c("darkolivegreen3", "tomato3"))

  rel_risk <- ggplot(sort_avgp) +
    geom_step(aes(p, relative_risk, colour = "RR"), direction = "hv") + 
    geom_step(aes(p, odds_ratio, colour = "OR"), direction = "hv") +
    scale_color_manual("",
                       breaks = c("RR", "OR"),
                       values = c("RR" = "orange", "OR" = "steelblue3")) +
    theme(axis.title.y = element_blank())
  
  ax1 <- plot_grid(jitter_features, emp_cdf, rel_risk,
                   ncol = 3, rel_widths = c(2.5,1.2,1.3))
  
  parallel <- ggplot(long_percent_y) +
    geom_line(aes(feature, value, colour = label, group = pid),
              show.legend = F) +
    scale_color_manual(values = c("darkolivegreen3", "tomato3"))
  
  mrd_p <- ggplot(avgpercent_mrd) +
    geom_point(aes(p, d33_mrd, colour = label), cex = 3, show.legend = F) +
    scale_color_manual(values = c("darkolivegreen3", "tomato3"))
  
  ax2 <- plot_grid(parallel, mrd_p, ncol = 2, rel_widths = c(2.5,1))
  
  fig <- plot_grid(ax1, ax2, nrow = 2)
  return(fig)
}

plotFeatures <- function(results, metadata_df) {
  y <- as.factor(metadata_df[rownames(results),"label"])
  subset_features1 <- c("erm1", "angle_d0d8_normal", "l2norm_d0_d8",
                        "l2norm_d0_normal", "l2norm_d8_normal", "l2norm_diff",
                        "erm1_ratio1", "erm1_ratio2", "erm1_ratio3",
                        "l2norm_ratio1", "l2norm_ratio2", "l2norm_diff_ratio")
  
  features1 <- results[, subset_features1, drop=F]
  features1_y <- data.frame(features1, label = y)
  long_features1_y <- melt(features1_y, id="label", variable.name = "feature")
  
  # PLOT: FEATURES
  jitter_features1 <- ggplot(long_features1_y) +
    geom_point(aes(feature, value, colour = label),
               position = position_jitterdodge(), cex = 3,
               show.legend = F) +
    scale_color_manual(values = c("darkolivegreen3", "tomato3")) +
    facet_wrap(~feature, nrow = 2, ncol = 6,  scales = "free") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_text(angle = 10, vjust = 0.5))
  
  return(jitter_features1)
}

# Factor to split data
splitSubtype <- function(X, metadata_df) {
  subtype_factor <- as.factor(metadata_df[colnames(X), "subtype"])
  split.default(X, subtype_factor, drop = F) # Split by subtype
}

# IMPORT DATA -------------------------------------------------------------
## Subset of original data
# Removed outliers, patients with timepoints from different batches and batch 5
SUBSET_RPATH <- "data/GSE67684/processed/subset_yeoh.tsv"
raw_yeoh <- read.table(SUBSET_RPATH, sep = "\t")

## Metadata
# Preprocessed metadata
METADATA_RPATH <- "data/GSE67684/processed/metadata/full_metadata.tsv"
metadata_df <- read.table(METADATA_RPATH, sep = "\t")

BATCH_RPATH <- "data/GSE67684/processed/metadata/metadata-batch.tsv"
LABEL_RPATH <- "data/GSE67684/processed/metadata/metadata-label_mrd_subtype.tsv"
yeoh_batch <- read.table(BATCH_RPATH, sep = "\t", header = T, row.names = 1)
yeoh_label <- read.table(LABEL_RPATH, sep = "\t", header = T, row.names = 1)

# SCALE->REMOVE->FILTER->LOG
scaled_yeoh <- normaliseMeanScaling(raw_yeoh)
selected_yeoh <- removeProbesets(scaled_yeoh)
data_yeoh <- log2_transform(filterProbesets(selected_yeoh, 0.7, metadata_df))

# # Filter out all rows with zero values
# logi_idx <- rowSums(data_yeoh == 0) == 0
# filtered_yeoh <- data_yeoh[logi_idx,]

# Prediction (Batch genes) -------------------------------------------------------
## Batch genes
# Only D0 samples
pid_d0 <- rownames(metadata_df)[metadata_df$class_info == "D0"]
pid_telaml1 <- rownames(metadata_df)[metadata_df$subtype == "TEL-AML1"]
pid_remission <- rownames(metadata_df)[metadata_df$label == 0]

# Recursive intersect
pid_idx <- intersect(
  intersect(pid_remission, intersect(pid_d0, pid_telaml1)),
  colnames(data_yeoh)
)
d0_telaml1 <- data_yeoh[,pid_idx]
d0_batch <- metadata_df[colnames(d0_telaml1), "batch_info"]


d0_telaml1_t <- t(d0_telaml1)
#' @param X matrix with samples as rows and features as columns
calcBatchANOVA <- function(X, batch, method = "welch") {
  .featureANOVA <- function(vec, d0_batch, method) {
    X <- data.frame(gene = vec,
                    batch = as.factor(d0_batch))
    
    if (method == "welch") return(oneway.test(gene~batch, X)$p.value)
    else if (method == "aov") return(unname(unlist(summary(aov(gene~batch, data = X)))[9]))
    else if (method == "kruskal") return(kruskal.test(gene~batch, X)$p.value)
    else stop("option not available for argument: method")
  }
  
  pvalue <- sapply(data.frame(X), .featureANOVA, batch, method)
  names(pvalue) <- substring(names(pvalue), 2)
  n_nan <- sum(sapply(pvalue, is.na))
  print(c("No. of NaNs =", n_nan))
  return(pvalue)
}

aov_pvalue <- calcBatchANOVA(d0_telaml1_t, d0_batch, method = "aov")
# welch_pvalue <- calcBatchANOVA(d0_telaml1_t, d0_batch, method = "welch")
# kruskal_pvalue <- calcBatchANOVA(d0_telaml1_t, d0_batch, method = "kruskal")

# Selecting by pvalue threshold
batch_genes <- names(aov_pvalue)[aov_pvalue < 0.05 & !is.na(aov_pvalue)]
# welch_genes <- names(welch_pvalue)[welch_pvalue < 0.05 & !is.na(welch_pvalue)]
# kruskal_genes <- names(kruskal_pvalue)[kruskal_pvalue < 0.05 & !is.na(kruskal_pvalue)]
length(batch_genes)

### PLOTS ###
X_batch <- data_yeoh[batch_genes,]
pheatmap(X_batch, col = brewer.pal(9, "Blues"),
         legend = T, border_color = "black", scale = "none",
         cluster_method = "ward.D2", cluster_rows = T, cluster_cols = T,
         show_colnames = F, show_rownames = F,
         annotation_col = metadata_df)
heatmap_batch <- recordPlot()
save_fig(heatmap_batch, "dump/heatmap-batch_2565.pdf",
         width = 10, height = 10)

table(metadata_df$batch_info, metadata_df$subtype)

list_batch_genes <- list(anova = batch_genes, welch = welch_genes,
                         kruskal = kruskal_genes)

upset(fromList(list_batch_genes),
      nsets = length(list_selected),
      nintersects = NA,
      order.by = "freq")
upset_plot <- recordPlot()
upset_plot
save_fig(upset_plot, "dump/upset-batch_genes.pdf",
         width = 8, height = 8)

## Label genes
pid_d8 <- rownames(metadata_df)[metadata_df$class_info == "D8"]
pid_idx <- intersect(pid_d8, colnames(data_yeoh))
d8_yeoh <- data_yeoh[,pid_idx]

subtype_factor1 <- as.factor(metadata_df[colnames(d8_yeoh), "subtype"])
subtypes_d8 <- split.default(d8_yeoh, subtype_factor1, drop = F) # Split by subtype

# Prediction (Drug genes) --------------------------------------------
## Drug responsive genes
# Nonlocals: pid_remission
getLocalGenes <- function(X_subtype, pid_remission,
                          alpha = 0.05, EXPR = 6, N = 50, LOGFC = 1) {
  pid_idx <- intersect(pid_remission, colnames(X_subtype))
  print(pid_idx)
  X_subtype_remission <- X_subtype[,pid_idx, drop = F]
  print(c("Dimension:", dim(X_subtype_remission)))
  n_pairs <- ncol(X_subtype_remission)/2
  # print(colnames(X_subtype_remission)[1:n_pairs])
  # print(colnames(X_subtype_remission)[-(1:n_pairs)])
  
  # P-value
  pvalue <- calc_ttest(X_subtype_remission, n_pairs, is_paired = T) # nan values!
  
  # # Plot
  # hist(pvalue, breaks = 20, main = subtype)
  # hist_class <- recordPlot()
  # # HIST_WPATH <- sprintf("dump/hist_class-%s.pdf", subtype)
  # save_fig(hist_class, HIST_WPATH,
  #          width = 6, height = 6)
  
  # # Q-value
  # calc_qvalue <- function(p) length(p)*p/rank(p)
  # qvalue <- calc_qvalue(pvalue) # FDR threshold
  # hist(qvalue, breaks =20)
  
  # Median paired log-FC
  d0_mu <- rowMeans(X_subtype_remission[,1:n_pairs])
  d8_mu <- rowMeans(X_subtype_remission[,-(1:n_pairs)])
  paired_logfc <- X_subtype_remission[,-(1:n_pairs)] -
    X_subtype_remission[,1:n_pairs] # D8 - D0
  median_logfc <- apply(paired_logfc, 1, median)
  print(sprintf("No. of NaN values in log-fc = %d",
                 sum(is.na(median_logfc))))
  median_logfc1 <- median_logfc[!is.na(median_logfc)]
  selected_median_logfc <- median_logfc1[d0_mu > EXPR | d8_mu > EXPR]
  print(sprintf("No. of probesets excluded by expr threshold = %d",
                length(median_logfc1) - length(selected_median_logfc)))
  # feat_top_median_logfc <- names(head(sort(selected_median_logfc), N))
  
  # # Custom t-statistic
  # deviation_median <- sweep(paired_logfc, 1, median_logfc, "-")
  # median_abs_dev <- apply(abs(deviation_median), 1, median)
  # test_stat <- median_logfc/(median_abs_dev/n_pairs^0.5)
  # pvalue <- pt(abs(test_stat)*-1, n_pairs-1)
  # hist(pvalue, breaks = 30)
  # feat_selected_p <- names(head(sort(pvalue), N))

  feat_p <- names(pvalue)[pvalue < alpha & !is.na(pvalue)]
  # At least one of the means have to be > EXPR
  feat_log2fc <- names(selected_median_logfc)[abs(selected_median_logfc) > LOGFC]
  print(sprintf("No. of features (p-value) = %d", length(feat_p)))
  print(sprintf("No. of features (log2-fc) = %d", length(feat_log2fc)))
  feat <- intersect(feat_p, feat_log2fc)
  return(feat)
}

# Factor to split data
subtypes_yeoh <- splitSubtype(data_yeoh, metadata_df)

X_subtypes <- subtypes_yeoh
X <- data_yeoh
normal_pid <- paste0("N0", c(1,2,4))
all_subtypes <- levels(metadata_df$subtype)
subtypes <- setdiff(all_subtypes, c("Hypodiploid", "Normal"))
pid_remission <- rownames(metadata_df)[metadata_df$label == 0]

# list_drug_genes <- list()
for (subtype in subtypes) {
  # print(c("Subtype:", subtype))
  
  # Select genes
  # X_subtype <- X_subtypes[[subtype]]
  #-- OPTION: Hyperdiploid
  col_pid <- pid_any_grp4[[1]]
  X_subtype <- X[, col_pid]
  pid <- c(col_pid, normal_pid)
  print(colnames(X_subtype))
  #--
  class_genes <- getLocalGenes(X_subtype, pid_remission)
  print(c("No. of selected genes = ", length(class_genes)))
  # list_drug_genes <- append(list_drug_genes, list(class_genes))
  
  selected_genes <- setdiff(class_genes, batch_genes)
  print(c("No. of final genes = ", length(selected_genes)))
  
  # # Subset pids in subtype
  # logi_idx <- rownames(metadata_df) %in% colnames(X) &
  #   metadata_df$subtype == subtype
  # subtype_pid <- rownames(metadata_df)[logi_idx]
  # subset_pid <- c(subtype_pid, normal_pid)

  # # Plot heatmap
  # X_class <- X[selected_genes, subset_pid]
  # pheatmap(X_class, col = brewer.pal(9, "Blues"),
  #          legend = T, border_color = "black", scale = "none",
  #          cluster_method = "ward.D2", cluster_rows = T, cluster_cols = T,
  #          show_colnames = T, show_rownames = F,
  #          annotation_col = metadata_df)
  # heatmap_class <- recordPlot()
  # HEATMAP_WPATH <- sprintf("~/Dropbox/temp/heatmap_drug-%s.pdf", subtype)
  # save_fig(heatmap_class, HEATMAP_WPATH,
  #          width = 10, height = 10)
  
  # Subtype and normal samples
  subset_yeoh <- X[selected_genes, pid] # OPTION!
  # subset_yeoh <- X[selected_genes, subset_pid] # OPTION!
  idx <- 1:(ncol(subset_yeoh)-3)
  response <- t(subset_yeoh)[idx,]
  normal <- t(subset_yeoh)[-idx,]
  print(colnames(subset_yeoh))
  print(rownames(response))

  # Collate MRD results as well
  results <- calcERM(response, normal)
  subset_features <- c("erm1", "erm1_ratio1", "erm1_ratio2",
                       "angle_d0d8_normal", "l2norm_d0_d8",
                       "l2norm_d0_normal", "l2norm_d8_normal",
                       "l2norm_ratio1", "l2norm_ratio2", "l2norm_ratio3")
  
  # # Plot heatmap of features
  # pheatmap(t(results[,subset_features]), col = brewer.pal(9, "Blues"),
  #          legend = T, border_color = "black", scale = "row",
  #          cluster_method = "ward.D2", cluster_rows = T, cluster_cols = T,
  #          show_colnames = F, show_rownames = T,
  #          annotation_col = metadata_df[, "label", drop = F])
  # heatmap_class <- recordPlot()
  # HEATMAP_WPATH <- sprintf("dump/heatmap_features-%s.pdf", subtype)
  # save_fig(heatmap_class, HEATMAP_WPATH,
  #          width = 10, height = 10)
  
  # features <- plotFeatures(results, metadata_df)
  # FEATURES_WPATH <- sprintf("~/Dropbox/temp/features_drug-%s.pdf", subtype)
  # ggsave(FEATURES_WPATH, features, width = 16, height = 10)
  
  # Plot
  prediction_parallel <- plotPrediction(results, metadata_df, yeoh_label)
  prediction_parallel
  
  # OPTION: Automatic filename
  # PREDICTION_WPATH <- sprintf("~/Dropbox/temp/prediction_top-%s.pdf", subtype)
  PREDICTION_WPATH <- "~/Dropbox/temp/prediction_grp4_HR_24.pdf"
  ggsave(PREDICTION_WPATH, prediction_parallel, width = 14, height = 7)
}

names(list_drug_genes) <- subtypes
saveRDS(list_drug_genes, "temp/list_drug_genes.rds")

hist(as.numeric(data_yeoh[,5]), breaks = 40)

## List of selected genes from each subtype
names(list_selected) <- subtypes
upset(fromList(list_selected),
      nsets = length(list_selected),
      nintersects = NA,
      order.by = "freq")
upset_plot <- recordPlot()
upset_plot
save_fig(upset_plot, "dump/upset-selected_genes.pdf",
         width = 10, height = 5)

# Subset of intersected genes
subset_selected <- list_selected[-c(1,4)]
intersect_genes <- Reduce(intersect, subset_selected)

# Hyperdiploid preprocessing ---------------------------------------------
## # Normalised: D0 data
## idx_d0 <- metadata_df[colnames(data_yeoh), "class_info"] == "D0"
## sum_d0 <- colSums(data_yeoh)[idx_d0]

## Plot: Sum of expression
# D <- data.frame(subtype = metadata_df[names(sum_d0), "subtype"],
#                 value = sum_d0)
# features_plot <- ggplot(D, aes(as.factor(subtype), value, colour = subtype)) +
#   geom_point(position = position_jitter(width=.1, height=0), cex = 2, show.legend = F) + # position = position_jitterdodge()
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_text(angle = 20, vjust = 0.5))
# ggsave("~/Dropbox/temp/colsum-scaled.pdf", features_plot,
#        width = 6, height = 5)

# # Raw: D0 data
# selected_raw <- removeProbesets(raw_yeoh)
# unnorm_raw <- log2_transform(filterProbesets(selected_raw, 0.7, metadata_df))
# sum_raw_d0 <- colSums(unnorm_raw)[idx_d0]
# 
# D1 <- data.frame(subtype = metadata_df[names(sum_raw_d0), "subtype"],
#                  value = sum_raw_d0)
# features_plot1 <- ggplot(D1, aes(as.factor(subtype), value, colour = subtype)) +
#   geom_point(position = position_jitter(width=.1, height=0), cex = 2, show.legend = F) + # position = position_jitterdodge()
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_text(angle = 20, vjust = 0.5))
# ggsave("~/Dropbox/temp/colsum-unnorm.pdf", features_plot1,
#        width = 6, height = 5)

##### SUBSETTING DATA ######
## Normal
idx_normal <- metadata_df[colnames(data_yeoh), "subtype"] == "Normal"
normal <- data_yeoh[,idx_normal]

## Hyperdiploid
idx_hyp <- metadata_df[colnames(data_yeoh), "subtype"] == "Hyperdiploid" &
  metadata_df[colnames(data_yeoh), "class_info"] == "D0"
hyperdiploid <- data_yeoh[,idx_hyp]
colnames(hyperdiploid) <- substring(colnames(hyperdiploid), 1, 4)

## TEL-AML1
idx_telaml1 <- metadata_df[colnames(data_yeoh), "subtype"] == "TEL-AML1" &
  metadata_df[colnames(data_yeoh), "class_info"] == "D0"
telaml1 <- data_yeoh[,idx_telaml1]
colnames(telaml1) <- substring(colnames(telaml1), 1, 4)

## MLL
idx_mll <- metadata_df[colnames(data_yeoh), "subtype"] == "MLL" &
  metadata_df[colnames(data_yeoh), "class_info"] == "D0"
mll <- data_yeoh[,idx_mll]
colnames(mll) <- substring(colnames(mll), 1, 4)

##### BATCH EFFECT CORRECTION #####
### ComBat
# Obtaining batch information of selected_yeoh df
# Rows of metadata are to be in same order as columns of edata
batch <- as.factor(metadata_df[colnames(data_yeoh), "batch_info"])
timepoint <- as.factor(metadata_df[colnames(data_yeoh), "class_info"])
# Covariates (subtype, label) are confounded
yeoh_metadata <- data.frame(batch, timepoint)

# Place adjustment/confounding variables in model.matrix (e.g. age)
# Do not put batch variables in model.matrix
## Put batch variables directly in combat function!
# OPTION: Include biological variable of interest as covariate
# model_combat <- model.matrix(~1, data = yeoh_metadata)
model_combat <- model.matrix(~timepoint, data = yeoh_metadata)
combat_yeoh <- ComBat(data.matrix(data_yeoh), batch, model_combat)
# Replacing negative values with 0
combat_yeoh[combat_yeoh < 0] <- 0

## Normal
idx_normal <- metadata_df[colnames(data_yeoh), "subtype"] == "Normal"
norm_combat <- combat_yeoh[,idx_normal]

## Hyperdiploid
idx_hyp <- metadata_df[colnames(data_yeoh), "subtype"] == "Hyperdiploid" &
  metadata_df[colnames(data_yeoh), "class_info"] == "D0"
hyp_combat <- combat_yeoh[,idx_hyp]
colnames(hyp_combat) <- substring(colnames(hyp_combat), 1, 4)

##### FILTER OUT ZEROS (ALL PATIENTS) #####
row_pct_nonzero <- rowSums(hyperdiploid != 0)/ncol(hyperdiploid)
PCT_THRESHOLD <- 0.7
ps_idx <- names(row_pct_nonzero)[row_pct_nonzero > PCT_THRESHOLD]

normal1 <- normal[ps_idx,]
hyperdiploid1 <- hyperdiploid[ps_idx,]
hyp_normal1 <- cbind(hyperdiploid1, normal1)

## Annotation: Chr location
ANNOT_RPATH <- "../info/microarray/HG-U133_Plus_2/affy/HG-U133_Plus_2.na35.annot.csv"
annot <- read.csv(ANNOT_RPATH,  row.names = 1, comment.char = "#")

DATA <- normal1 ## OPTION1
DATA <- normal ## OPTION2
ps_chrloc <- annot[rownames(DATA), "Chromosomal.Location"]
ps_chr <- sub("(chr.*?)(p|q|c).*", "\\1", ps_chrloc)
ps_chr[ps_chr == "---"] <- NA
names(ps_chr) <- rownames(DATA)

## FUNCTIONS
plot_chr_hyp <- function(X_subtype, X_norm, wpath1, wpath2) {
  X <- cbind(X_subtype, X_norm)
  long_chr <- melt(data.matrix(X), varnames = c("chr", "pid"))
  long_chr$chr <- factor(long_chr$chr,
                         levels = levels(long_chr$chr)[
                           c(1,12,16:22,2:11,13:15,23)])
  y_lim <- c(floor(min(X)), ceiling(max(X)))

  ## Create color map for chr
  g1_chr <- c("chr4", "chr6", "chr10", "chr14", "chr17",
              "chr18", "chr21", "chrX")
  g4_chr <- c("chr1", "chr7", "chr13", "chr15", "chr19", "chr20")
  col_chr <- setdiff(g1_chr, "chrX")
  all_chr <- paste0("chr", 1:22)
  nocol_chr <- setdiff(all_chr, col_chr)
  nocol_chr <- setdiff(all_chr, c(col_chr, g4_chr))
  col_map <- c(rep("darkolivegreen3", length(col_chr)),
               rep("tomato3", length(g4_chr)),
               rep("black", 22-length(col_chr)-length(g4_chr)))
  names(col_map) <- c(col_chr, g4_chr, nocol_chr)

  jitter1 <- ggplot(long_chr[1:(20*22),],
                        aes(chr, value, color = chr)) +
    geom_point(position = position_jitter(width=.1, height=0),
               cex = 2, show.legend = F) +
    facet_wrap(~pid, nrow = 4, ncol = 5,  scales = "free_x") +
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.x=element_text(angle = 90, vjust = 0.5)) +
    scale_color_manual(values = col_map) +
    ylim(y_lim[1], y_lim[2])

  jitter2 <- ggplot(long_chr[(20*22+1):902,], aes(chr, value, color = chr)) +
    geom_point(position = position_jitter(width=.1, height=0),
               cex = 2, show.legend = F) +
    facet_wrap(~pid, nrow = 4, ncol = 6,  scales = "free_x") +
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.x=element_text(angle = 90, vjust = 0.5)) +
    scale_color_manual(values = col_map) +
    ylim(y_lim[1], y_lim[2])

  ggsave(wpath1, jitter1, width = 16, height = 10)
  ggsave(wpath2, jitter2, width = 16, height = 10)
}


#' @param X_subtype dataframe of chr statistics containing only subtype patients
#' @param X_norm dataframe of chr statistics containing only subtype patients
plot_heatmap_batch <- function(X_subtype, X_norm, metadata, filename) {
  if (endsWith(colnames(X_subtype)[1], "D0"))
      stop("Colnames already suffixed with D0!")

  subset_metadata <- metadata[,c("batch_info", "subtype", "label"), drop = F]

  colnames(X_subtype) <- paste(colnames(X_subtype), "D0", sep = "_")
  ord_idx <- order(metadata[colnames(X_subtype), "batch_info"])
  X_subtype_ord <- X_subtype[,ord_idx]
  gaps <- cumsum(table(metadata[colnames(X_subtype), "batch_info"]))
  X_ord <- cbind(X_subtype_ord, X_norm)

  pheatmap(X_ord, col = brewer.pal(9, "Blues"),
           display_numbers = F, fontsize = 5.5,
           legend = T, border_color = "black", scale = "none",
           cluster_method = "ward.D2", cluster_rows = F, cluster_cols = F,
           show_colnames = T, show_rownames = T,
           annotation_col = subset_metadata,
           gaps_col = gaps, cellwidth = 8, cellheight = 10,
           filename = filename)
  cat("Heatmap saved!\n")
}

# ## Un-normalised data
# idx <- metadata_df[colnames(data_yeoh), "subtype"] %in% c("Hyperdiploid", "Normal") & 
#   metadata_df[colnames(data_yeoh), "class_info"] %in% c("D0", "N")
# hyp_raw <- unnorm_raw[,idx]
# ps_chrloc1 <- annot[rownames(hyp_raw), "Chromosomal.Location"]
# ps_chr1 <- sub("(chr.*?)(p|q|c).*", "\\1", ps_chrloc1)
# ps_chr1[ps_chr1 == "---"] <- NA
# list_chr_hypdip1 <- split.data.frame(hyp_raw, ps_chr1)
# hypdip_chr_mean1 <- t(sapply(list_chr_hypdip1, colMeans))
# hypdip_no_chrY1 <- hypdip_chr_mean1[1:23,]

## # Plot all genes by chr
## hyp_chr <- cbind(chr = ps_chr, hyperdiploid)
## chr_idx <- !(ps_chr %in% c("chrY", NA))

## # pid <- "P047_D0"
## for (pid in colnames(hyperdiploid)) {
##   plot_data <- droplevels(hyp_chr[chr_idx, c("chr", pid)])
##   colnames(plot_data) <- c("chr", "log2_expr")
##   # plot_data <- plot_data[plot_data$log2_expr != 0,]
  
##   groupby_chr <- ggplot(plot_data, aes(chr, log2_expr)) +
##     geom_point(position = position_jitter(width=.2, height=0)) +
##     geom_boxplot(alpha=.5) +
##     ggtitle(pid) +
##     theme(axis.title.x=element_blank(),
##           axis.text.x=element_text(angle = 30, vjust = 0.5),
##           plot.title = element_text(hjust = 0.5))
##   WPATH <- sprintf("~/Dropbox/temp/%s.pdf", pid)
##   ggsave(WPATH, groupby_chr, width = 12, height = 6)
## }

## MEDIAN
## # Hyperdiploid D0: Split into chr
## list_chr_hypdip <- split.data.frame(hyperdiploid, ps_chr)
## hypdip_chr_median <- t(sapply(list_chr_hypdip, apply, 2, median)) # median
## hypdip_median_no_chrY <- hypdip_chr_median[1:22,]
## ranked_chr_median <- apply(-hypdip_median_no_chrY, 2,
##                            function(x) names(sort(x)))
## # Values across samples may be affected by batch effects
## # Create relative values that remain constant
## # Relative to a basket of chr 1, 7, 9, 16
## median_ref <- colMeans(hypdip_median_no_chrY[c(1,20,22,8),])
## median_ratio_within <- sweep(hypdip_median_no_chrY, 2, median_ref, "/")

## # Plot all medians
## i <- 10
## plot(rep(0, 22), hypdip_median_no_chrY[,i],
##      main = colnames(hypdip_median_no_chrY)[i])
## text(rep(0, 22) + .5, hypdip_median_no_chrY[,i],
##      rownames(hypdip_median_no_chrY),
##      cex = .8)

# ## Ranks
# hypdip_rank <- apply(hyperdiploid, 2, rank, ties.method = "min")
# list_chr_rank <- split.data.frame(hypdip_rank, ps_chr)
# hypdip_chr_rank <- t(sapply(list_chr_rank, colMeans))
# hypdip_rank_no_chrY <- hypdip_chr_rank[1:22,]

## # Plot
## pheatmap(hypdip_no_chrY, col = brewer.pal(9, "Blues"),
##          legend = T, border_color = NA, scale = "none",
##          cluster_method = "ward.D2", cluster_rows = T, cluster_cols = T,
##          show_colnames = F, show_rownames = T,
##          annotation_col = metadata_df)
## heatmap <- recordPlot()
## save_fig(heatmap, "~/Dropbox/temp/heatmap_none-hypdip_rank_no_chrY.pdf",
##          width = 7, height = 6)

## # Plot batch effects
## long_batch <- cbind(long_no_chrY,
##                     batch = metadata_df[long_no_chrY$pid, "batch_info"])
## jitter_batch <- ggplot(long_batch, aes(pid, value, colour=pid)) +
##   geom_point(position = position_jitter(width=.1, height=0),
##              cex = 2, show.legend = F) +
##   facet_grid(~batch, scales = "free", space = "free") +
##   theme(axis.title.x=element_blank(),
##         axis.text.x=element_text(angle = 90, vjust = 0.5))
## ggsave("~/Dropbox/temp/jitter_batch-scaled.pdf", jitter_batch,
##        width = 16, height = 10)

## ## Plot: Raw data
## long_no_chrY1 <- melt(hypdip_no_chrY1, varnames = c("chr", "pid"))
## long_no_chrY1$chr <- factor(long_no_chrY1$chr,
##                             levels = levels(long_no_chrY1$chr)[
##                               c(1,12,16:22,2:11,13:15,23)])
## # Plot pid facet
## chr_jitter3 <- ggplot(long_no_chrY1[1:(20*23),], aes(chr, value)) +
##   geom_point(position = position_jitter(width=.1, height=0),
##              cex = 2, show.legend = F) +
##   facet_wrap(~pid, nrow = 4, ncol = 5,  scales = "free_x") +
##   theme(axis.title.x=element_blank(),
##         axis.text.x=element_text(angle = 90, vjust = 0.5))
## chr_jitter4 <- ggplot(long_no_chrY1[(20*23+1):943,], aes(chr, value)) +
##   geom_point(position = position_jitter(width=.1, height=0),
##              cex = 2, show.legend = F) +
##   facet_wrap(~pid, nrow = 4, ncol = 6,  scales = "free_x") +
##   theme(axis.title.x=element_blank(),
##         axis.text.x=element_text(angle = 90, vjust = 0.5))
## chr_jitter3
## ggsave("~/Dropbox/temp/jitter-chr2.pdf", chr_jitter2,
##        width = 16, height = 10)

## # Plot batch
## long_batch1 <- cbind(long_no_chrY1,
##                     batch = metadata_df[long_no_chrY$pid, "batch_info"])
## jitter_batch1 <- ggplot(long_batch1, aes(pid, value, colour=pid)) +
##   geom_point(position = position_jitter(width=.1, height=0),
##              cex = 2, show.legend = F) +
##   facet_grid(~batch, scales = "free", space = "free") +
##   theme(axis.title.x=element_blank(),
##         axis.text.x=element_text(angle = 90, vjust = 0.5))
## ggsave("~/Dropbox/temp/jitter_batch-raw.pdf", jitter_batch1,
##        width = 16, height = 10)

## ## Rank within patient
## rank_chr_mean <- apply(-hypdip_no_chrY, 2, rank)
## print(rank_chr_mean)

## # Rank: Mean and sd
## chr_rank_sd <- apply(rank_chr_mean, 1, sd)
## chr_rank_mean <- rowMeans(rank_chr_mean)
## chr_rank_mean
## chr_rank_sd
## plot(chr_rank_mean, chr_rank_sd,
##      xlim = c(0,25), ylim = c(0,6),
##      xlab = "Mean", ylab = "SD", main = "Chromosome ranks")
## text(chr_rank_mean+.8, chr_rank_sd+.15,
##      names(chr_rank_mean), cex = .8)
## rank_mean_sd <- recordPlot()
## save_fig(rank_mean_sd, "~/Dropbox/temp/rank_scatter-chr.pdf",
##          width = 7, height = 8)

## # Values across samples may be affected by batch effects
## # Create relative values that remain constant
## # Relative to a basket of chr 1, 7, 9, 16
## mean_ref <- colMeans(hypdip_no_chrY[c(1,20,22,8),])

## ratio_within <- sweep(hypdip_no_chrY, 2, mean_ref, "/")

## ### PLOT MLL
## list_chr_mll <- split.data.frame(mll, ps_chr)
## mll_nozero_median <- t(sapply(list_chr_mll, apply, 2,
##                             function(vec) median(vec[vec != 0])))
## mll_nozero_median1 <- mll_nozero_median[1:22,]

## std_chr_mll <- (mll_nozero_median1-ref_mean)/ref_sd

## long_standardised_chr <- melt(std_chr_mll, varnames = c("chr", "pid"))
## long_standardised_chr$chr <- factor(long_standardised_chr$chr,
##                                     levels = levels(long_standardised_chr$chr)[
##                                       c(1,12,16:22,2:11,13:15,23)])
## # Create color map for chr
## col_chr <- setdiff(lit_chr, "chrX")
## all_chr <- paste0("chr", 1:22)
## nocol_chr <- setdiff(all_chr, col_chr)
## nocol_chr <- setdiff(all_chr, c(col_chr, bad_chr))
## col_map <- c(rep("darkolivegreen3", length(col_chr)),
##              rep("tomato3", length(bad_chr)),
##              rep("black", 22-length(col_chr)-length(bad_chr)))
## names(col_map) <- c(col_chr, bad_chr, nocol_chr)

## mll_std <- ggplot(long_standardised_chr,
##                       aes(chr, value, color = chr)) +
##   geom_point(position = position_jitter(width=.1, height=0),
##              cex = 2, show.legend = F) +
##   facet_wrap(~pid, nrow = 2, ncol = 4,  scales = "free_x") +
##   theme(axis.title.x=element_blank(),
##         axis.title.y=element_blank(),
##         axis.text.x=element_text(angle = 90, vjust = 0.5)) +
##   scale_color_manual(values = col_map) +
##   ylim(-3.1, 8)

## ggsave("~/Dropbox/temp/mll_std-nozero_median.pdf", mll_std,
##        width = 16, height = 9)

### CLUSTERING HYPERDIPLOID ###
hyp_tsne <- Rtsne(t(hyperdiploid), perplexity = 5)$Y

### DENSITY PLOT ###
## Appending relevant info and converting to long format
t_hyp_normal <- t(hyp_normal1[
  !is.na(ps_chr) & !(ps_chr %in% c("chrX", "chrY")),]) # Filter out unassigned ps
batch_idx <- metadata_df[paste(rownames(t_hyp_normal), "D0", sep = "_"), "batch_info"]
batch_idx[is.na(batch_idx)] <- "Normal"
batch_hyp_normal <- data.frame(pid = rownames(t_hyp_normal),
                               batch = as.factor(batch_idx),
                               t_hyp_normal,
                               check.names = F)
long_hyp_normal <- melt(batch_hyp_normal,
                        id.vars = c("pid", "batch"),
                        variable.name = "probeset")
long_hyp_normal$chr <- ps_chr[as.character(long_hyp_normal$probeset)]
ps_density <- ggplot(long_hyp_normal) +
  geom_density(aes(x = value, group = pid, colour = batch)) +
  facet_wrap(~chr, nrow = 4, ncol = 6)
ggsave("~/Dropbox/temp/ps_density-fltr70.pdf", ps_density,
       width = 12, height = 8)

##### CHR SUMMARY STATISTICS #####
### NORMAL
list_chr_norm <- split.data.frame(normal, ps_chr)
## MEAN
norm_mean <- t(sapply(list_chr_norm, colMeans))[1:22,]
## MEDIAN
norm_median <- t(sapply(list_chr_norm, apply, 2, median))[1:22,]
## MEDIAN (NO ZERO)
norm_nozero_median <- t(sapply(list_chr_norm, apply, 2,
                                  function(vec) median(vec[vec != 0])))[1:22,]
## PCT ZERO
norm_pct_zero <- t(sapply(list_chr_norm, apply, 2,
                     function(vec) sum(vec == 0)/length(vec)))[1:22,]

### HYPERDIPLOID
list_chr_hyp <- split.data.frame(hyperdiploid, ps_chr)
## MEAN
hyp_mean <- t(sapply(list_chr_hyp, colMeans))[1:22,]
## MEDIAN
hyp_median <- t(sapply(list_chr_hyp, apply, 2, median))[1:22,]
## ## MEDIAN (NO ZERO)
hyp_nozero_median <- t(sapply(list_chr_hyp, apply, 2,
                            function(vec) median(vec[vec != 0])))[1:22,]
## PCT ZERO
hyp_pct_zero <- t(sapply(list_chr_hyp, apply, 2,
                         function(vec) sum(vec == 0)/length(vec)))[1:22,]

### NORMAL
list_chr_norm_combat <- split.data.frame(norm_combat, ps_chr)
## MEDIAN
norm_combat_median <- t(sapply(list_chr_norm_combat, apply, 2, median))[1:22,]
### HYPERDIPLOID
list_chr_hyp_combat <- split.data.frame(hyp_combat, ps_chr)
## MEDIAN
hyp_combat_median <- t(sapply(list_chr_hyp_combat, apply, 2, median))[1:22,]

## Plot
plot_heatmap_batch(hyp_mean, norm_mean, metadata_df,
                   "~/Dropbox/temp/fltr70-hyp_mean.pdf")
plot_heatmap_batch(hyp_median, norm_median, metadata_df,
                   "~/Dropbox/temp/fltr70-hyp_median.pdf")
plot_heatmap_batch(hyp_nozero_median, norm_nozero_median, metadata_df,
                   "~/Dropbox/temp/fltr70-hyp_nozero_median.pdf")
plot_heatmap_batch(hyp_pct_zero, norm_pct_zero, metadata_df,
                   "~/Dropbox/temp/fltr70-hyp_pct_zero.pdf")

plot_heatmap_batch(hyp_combat_median, norm_combat_median, metadata_df,
                   "~/Dropbox/temp/combat_median.pdf")

plot_chr_hyp(hyp_median, norm_median,
             "~/Dropbox/temp/median1.pdf",
             "~/Dropbox/temp/median2.pdf")

### TEL-AML1
list_chr_telaml1 <- split.data.frame(telaml1, ps_chr)
## MEAN
telaml1_mean <- t(sapply(list_chr_telaml1, colMeans))[1:22,]
## MEDIAN
telaml1_median <- t(sapply(list_chr_telaml1, apply, 2, median))[1:22,]
## MEDIAN (NO ZERO)
telaml1_nozero_median <- t(sapply(list_chr_telaml1, apply, 2,
                                  function(vec) median(vec[vec != 0])))[1:22,]
## PCT ZERO
telaml1_pct_zero <- t(sapply(list_chr_telaml1, apply, 2,
                     function(vec) sum(vec == 0)/length(vec)))[1:22,]

## ## Plot heatmap for batch effects
## plot_heatmap_batch(telaml1_mean, norm_mean, metadata_df,
##                    "~/Dropbox/temp/heatmap-telaml1_mean.pdf")
## plot_heatmap_batch(telaml1_nozero_median, norm_nozero_median, metadata_df,
##                    "~/Dropbox/temp/heatmap-telaml1_nozero_median.pdf")
## plot_heatmap_batch(telaml1_median, norm_median, metadata_df,
##                    "~/Dropbox/temp/heatmap-telaml1_median.pdf")
## plot_heatmap_batch(telaml1_pct_zero, norm_pct_zero, metadata_df,
##                    "~/Dropbox/temp/heatmap-telaml1_pct_zero.pdf")

### ESTIMATE BATCH EFFECTS FROM TEL-AML1
## No TEL-AML1 patients in batch 10
## Filtering out relapse patients in a batch would not generalise to unseen test batch
## TEL-AML1 (Ref batch: 2)
## Hyperdiploid is heterogeneous

#' @param X_subtype dataframe consisting of only subtype patients
calc_correction <- function(X_subtype, metadata, modify = F) {
  list_batch <- split.default(
    data.frame(X_subtype),
    metadata[paste(colnames(X_subtype), "D0", sep = "_"), "batch_info"]
  )
  batch_chr <- sapply(list_batch, apply, 1, median)
  correction <- batch_chr - batch_chr[,"2"]

  if (!modify) {
    return(correction)
  } else {
    list_corr <- lapply(names(list_batch),
                        function(batch) list_batch[[batch]] - correction[,batch])
    X_corr <- do.call(cbind, list_corr)
    return(X_corr[,colnames(X_subtype)])
  }
}

corr_tel_nozero_median <- calc_correction(telaml1_nozero_median, metadata_df)
corr_hyp_nozero_median <- calc_correction(hyp_nozero_median, metadata_df)

## Correlation between BE estimated from TEL-AML1 and hyperdiploid
for (batch in colnames(corr_tel_nozero_median)) {
  pdf(sprintf("~/Dropbox/temp/corr_hyp_tel-%s.pdf", batch))
  plot(corr_hyp_nozero_median[,batch],
       corr_tel_nozero_median[,batch],
       xlab = "Correction value (Hyperdiploid)",
       ylab = "Correction value (TEL-AML1)",
       main = sprintf("Batch %s", batch))
  dev.off()
}

## Perform correction using hyp and evaluate
X_corr_hyp <- calc_correction(hyp_median, metadata_df, modify = T)
plot_heatmap_batch(X_corr_hyp, norm_median, metadata_df,
                   "~/Dropbox/temp/fltr30-corr_hyp_median.pdf")

### REMOVE BATCH EFFECT GENES ###
fltr_tel <- telaml1[!(rownames(telaml1) %in% batch_genes),]
fltr_norm <- normal[!(rownames(normal) %in% batch_genes),]

## Annotation for fltr dataframe without batch genes
ps_chrloc1 <- annot[rownames(fltr_tel), "Chromosomal.Location"]
ps_chr1 <- sub("(chr.*?)(p|q|c).*", "\\1", ps_chrloc1)
ps_chr1[ps_chr1 == "---"] <- NA

list_fltr_tel <- split.data.frame(fltr_tel, ps_chr1)
list_fltr_norm <- split.data.frame(fltr_norm, ps_chr1)

## TEL-AML1 (FILTERED)
## MEAN
chr_mean <- t(sapply(list_chr_fltr_telaml1, colMeans))[1:22,]
## MEDIAN
chr_median <- t(sapply(list_chr_fltr_telaml1, apply, 2, median))[1:22,]
## MEDIAN (NO ZERO)
nozero_median <- t(sapply(list_chr_fltr_telaml1, apply, 2,
                          function(vec) median(vec[vec != 0])))[1:22,]

## ## Heatmap for batch effects
## plot_heatmap_batch(fltr_mean, metadata_df,
##                    "~/Dropbox/temp/heatmap-fltr_telaml1_mean.pdf")
## plot_heatmap_batch(fltr_median, metadata_df,
##                    "~/Dropbox/temp/heatmap-fltr_telaml1_median.pdf")
## plot_heatmap_batch(fltr_nozero_median, metadata_df,
##                    "~/Dropbox/temp/heatmap-fltr_telaml1_nozero_median.pdf")

## PCT ZERO
tel_pct_zero <- t(sapply(list_fltr_tel, apply, 2,
                         function(vec) sum(vec == 0)/length(vec)))[1:22,]
norm_pct_zero <- t(sapply(list_fltr_norm, apply, 2,
                          function(vec) sum(vec == 0)/length(vec)))[1:22,]
plot_heatmap_batch(tel_pct_zero, norm_pct_zero, metadata_df,
                   "~/Dropbox/temp/heatmap-fltr_telaml1_pct_zero.pdf")



#### TEST STATISTIC ####

### Z-SCORE (N01)
## Determine mean and sd of chr medians of normal
ref_normal <- norm_nozero_median[,"N01"]
ref_mean <- mean(ref_normal)
ref_sd <- sd(ref_normal)

std_hyp_nozero_median <- (hyp_nozero_median-ref_mean)/ref_sd
std_norm_nozero_median <- (norm_nozero_median-ref_mean)/ref_sd
plot_chr_hyp(std_hyp_median, std_norm_median,
             "~/Dropbox/temp/fltr30_median_std1.pdf",
             "~/Dropbox/temp/fltr30_median_std2.pdf")

### ComBat - Z-SCORE (N01)
## Determine mean and sd of chr medians of normal
ref_normal1 <- norm_combat_median[,"N01"]
ref_mean1 <- mean(ref_normal1)
ref_sd1 <- sd(ref_normal1)

hyp_combat_median_std <- (hyp_combat_median-ref_mean1)/ref_sd1
norm_combat_median_std <- (norm_combat_median-ref_mean1)/ref_sd1
plot_chr_hyp(hyp_combat_median_std, norm_combat_median_std,
             "~/Dropbox/temp/combat_median_std1.pdf",
             "~/Dropbox/temp/combat_median_std2.pdf")

## BATCH CORRECTED HYPERDIPLOID
std_corr_hyp_median <- (X_corr_hyp-ref_mean)/ref_sd
std_norm_median <- (norm_median-ref_mean)/ref_sd
plot_chr_hyp(std_corr_hyp_median, std_norm_median,
             "~/Dropbox/temp/fltr30-std_corr_median1.pdf",
             "~/Dropbox/temp/fltr30-std_corr_median2.pdf")q

### Z-SCORE (INDV CHR)
## Every chr has its own normal mu
normal_mu_chr <- rowMeans(norm_median)
normal_sigma_chr <- apply(norm_median, 1, sd)
avg_sigma_chr <- mean(normal_sigma_chr)
stdindv_hyp_median <- sweep(hyp_median, MARGIN = 1,
                            STATS = normal_mu_chr, FUN = `-`)/avg_sigma_chr
stdindv_norm_median <- sweep(norm_median, MARGIN = 1,
                             STATS = normal_mu_chr, FUN = `-`)/avg_sigma_chr

plot_chr_hyp(stdindv_hyp_median, stdindv_norm_median,
             "~/Dropbox/temp/fltr70_median-std_indv1.pdf",
             "~/Dropbox/temp/fltr70_median-std_indv2.pdf")

##### FEATURE ENGINEERING #####
#' @param chr_rank dataframe containing ranked chr of hyp patients only
#' @param n_rank number of ranks to consider
#' @return number of G4 chr in top n chr of patient
calc_n_g4 <- function(chr_rank, n_rank) {
  g4_chr <- c("chr1", "chr7", "chr13", "chr15", "chr19", "chr20")
  apply(chr_rank[1:n_rank,], 2, function(x) sum(x %in% g4_chr))
}

calc_n_g4_threshold <- function(X_chr, threshold, lower = T) {
  g4_chr <- c("chr1", "chr7", "chr13", "chr15", "chr19", "chr20")
  if (lower == T) {
    idx_logi <- X_chr[rownames(X_chr) %in% g4_chr,] < threshold
  }
  apply(idx_logi, 2, sum)
}

rank_hyp_nozero_median <- rank_chr(std_hyp_nozero_median)
n_g4 <- calc_n_g4(rank_hyp_nozero_median, 4)
n_g4_threshold <- calc_n_g4_threshold(std_hyp_nozero_median, -2)

## WEIXIN
calc_chr_grp <- function(x_chr) {
  g1_chr <- c("chr4", "chr6", "chr10", "chr14", "chr17", "chr18", "chr21", "chrX")
  g2_chr <- c("chr5", "chr8", "chr11", "chr12")
  g3_chr <- c("chr2", "chr3", "chr9", "chr16", "chr22")
  g4_chr <- c("chr1", "chr7", "chr13", "chr15", "chr19", "chr20")
  list_grp <- list(g1_chr, g2_chr, g3_chr, g4_chr)
  grp_mean <- sapply(list_grp, function(chr) mean(x_chr[names(x_chr) %in% chr]))
  names(grp_mean) <- c("g1_chr", "g2_chr", "g3_chr", "g4_chr")
  return(grp_mean)
}

norm_chr_mu <- apply(norm_median, 1, mean)
norm_chr_sigma <- apply(norm_median, 1, sd)

## pdf("~/Dropbox/temp/norm_median-mu_sigma.pdf")
## plot(norm_chr_mu, norm_chr_sigma)
## text(norm_chr_mu, norm_chr_sigma, names(norm_chr_mu))
## dev.off()

## Mean of grp1-4 for each patient
norm_grpchr_mean <- calc_chr_grp(norm_chr_mu)
hyp_grpchr_mean <- apply(hyp_median, 2, calc_chr_grp)
hyp_batch <- metadata_df[paste(colnames(hyperdiploid), "D0", sep = "_"), "batch_info"]
names(hyp_batch) <- colnames(hyperdiploid)

### ComBat
norm_chr_mu1 <- apply(norm_combat_median, 1, mean)
norm_chr_sigma1 <- apply(norm_combat_median, 1, sd)

## Mean of grp1-4 for each patient
norm_grpchr_mean1 <- calc_chr_grp(norm_chr_mu1)
hyp_grpchr_mean1 <- apply(hyp_combat_median, 2, calc_chr_grp)
hyp_batch <- metadata_df[paste(colnames(hyperdiploid), "D0", sep = "_"), "batch_info"]
names(hyp_batch) <- colnames(hyperdiploid)

## Z-SCORE (INDV CHR)
avg_sigma1 <- mean(norm_chr_sigma1)
hyp_combat_median_stdindv <- sweep(hyp_combat_median, MARGIN = 1,
                                   STATS = norm_chr_mu1, FUN = `-`)/avg_sigma1
norm_combat_median_stdindv <- sweep(norm_combat_median, MARGIN = 1,
                                    STATS = norm_chr_mu1, FUN = `-`)/avg_sigma1

plot_chr_hyp(hyp_combat_median_stdindv, norm_combat_median_stdindv,
             "~/Dropbox/temp/combat_median_stdindv1.pdf",
             "~/Dropbox/temp/combat_median_stdindv2.pdf")

## Hyperdiploid - Risk classification --------------------------------------
normal_pid <- paste0("N0", c(1,2,4))

#' @param X_chr Chromosomal summary of hyp and normal patients
#' @return Ranked chr of hyp patients only (only ranks autosomal chr)
rank_chr <- function(X_chr) {
  chr_rank <- apply(-X_chr[1:22,], 2, function(x) names(sort(x)))
  chr_rank <- chr_rank[,1:38]
  return(chr_rank)
}

#' @param chr_rank Ranked chr (does not include normal patients)
#' Returns pid of patients with top 4 chr in G4
get_pid_toprank <- function(chr_rank, n_rank) {
  g4_chr <- c("chr1", "chr7", "chr13", "chr15", "chr19", "chr20")
  top_chr <- chr_rank[1:n_rank,]
  list_topchr <- as.list(data.frame(top_chr))
  idx_topchr <- sapply(list_topchr, function(x) any(x %in% g4_chr))
  return(names(list_topchr[idx_topchr]))
}

#' @param X_subtype dataframe of chr summary containing only subtype patients
#' @param lower logical indicating if only lower threshold is to be used
#' @return patient IDs
get_pid_threshold <- function(X_subtype, threshold, lower) {
  g4_chr <- c("chr1", "chr7", "chr13", "chr15", "chr19", "chr20")
  ## g4_chr <- "chr13"
  ## Only group 4 chr
  ## X_fltr <- X_subtype[rownames(X_subtype) == g4_chr, , drop = F]
  X_fltr <- X_subtype[rownames(X_subtype) %in% g4_chr,]
  if (!lower) {
    idx_logi <- abs(X_fltr) > threshold
  } else {
    stopifnot(threshold < 0)
    idx_logi <- X_fltr < threshold
  }
  pid_idx_logi <- apply(idx_logi, 2, any)
  return(names(pid_idx_logi)[pid_idx_logi])
}

get_list_pid <- function(pid) {
  if (length(pid) == 0) stop("Empty list...")
  ## Assumption: Colnames of hyperdiploid has no D0/D8
  pid_not <- setdiff(colnames(hyperdiploid), pid)
  return(list(pid, pid_not))
}

#' Pastes D0 and D8 to pid
convert_list_pid <- function(list_pid) {
  lapply(list_pid, function(x) c(paste(x, "D0", sep = "_"),
                                 paste(x, "D8", sep = "_")))
}

#' @param list_pid list with first element containing high risk
#' patient ids and second element containing low risk patient ids
create_table <- function(list_pid) {
  relapse_pid <- c("P038", "P115", "P129", "P164", "P189")
  a <- sum(list_pid[[1]] %in% relapse_pid)
  b <- length(list_pid[[1]]) - a
  c <- sum(list_pid[[2]] %in% relapse_pid)
  d <- length(list_pid[[2]]) - c

  if (a < 4)
    warning("list_pid does not contain high risk patients in first element!")

  table_dimnames <- list(c("HR", "LR"),
                         c("Relapse", "Remission"))
  matrix(c(a, c, b, d), nrow = 2, dimnames = table_dimnames)
}

## ComBat - Median (filter indv) - Z-score
## HR: Top 4 in G4 OR Any G4 abs(chr) < threshold
pid_threshold <- get_pid_threshold(hyp_combat_median_std, -2.6, lower = T)
rank_hyp_median <- rank_chr(hyp_combat_median_std)
pid_top_hyp <- get_pid_toprank(rank_hyp_median, 4)
pid_union <- union(pid_threshold, pid_top_hyp)
list_pid <- get_list_pid(pid_union)
tab <- create_table(list_pid)
fisher <- fisher.test(tab)
chisq <- chisq.test(tab)
print(tab)
print(fisher)
print(chisq)

## Median (filter indv) - Z-score (N01)
## HR: Top 4 in G4 OR Any G4 abs(chr) < threshold
pid_threshold <- get_pid_threshold(std_hyp_nozero_median, -2.7, lower = T)
rank_hyp_median <- rank_chr(std_hyp_nozero_median)
pid_top_hyp <- get_pid_toprank(rank_hyp_median, 4)
pid_union <- union(pid_threshold, pid_top_hyp)
list_pid <- get_list_pid(pid_union)
tab <- create_table(list_pid)
fisher <- fisher.test(tab)
chisq <- chisq.test(tab)
print(tab)
print(fisher)
print(chisq)

### MEDIAN (FLTR 30%) - Z-SCORE (INDV)

## HR: Top 4 in G4 OR Any G4 abs(chr) < threshold
pid_threshold <- get_pid_threshold(stdindv_hyp_median, -8, lower = T)
rank_hyp_median <- rank_chr(stdindv_hyp_median)
pid_top_hyp <- get_pid_toprank(rank_hyp_median, 6)
pid_union <- union(pid_threshold, pid_top_hyp)
list_pid <- get_list_pid(pid_union)
tab <- create_table(list_pid)
fisher <- fisher.test(tab)
chisq <- chisq.test(tab)
print(tab)
print(fisher)
print(chisq)

### MEDIAN (FLTR 30%) - EST. CORRECTION - Z-SCORE (N01)
## HR: Top 4 in G4 OR Any G4 abs(chr) > threshold
## Most of the patients identified by overfitted threshold
X_std <- std_corr_hyp_median
pid_threshold <- get_pid_threshold(X_std, -2.5, lower = T)
rank_hyp_median <- rank_chr(X_std)
pid_top_hyp <- get_pid_toprank(rank_hyp_median, 3)
pid_union <- union(pid_threshold, pid_top_hyp)
list_pid <- get_list_pid(pid_union)
tab <- create_table(list_pid)
fisher <- fisher.test(tab)
chisq <- chisq.test(tab)
print(tab)
print(fisher)
print(chisq)

## MEDIAN (FLTR 30%) - Z-SCORE (N01)
X_std <- std_hyp_median
pid_threshold <- get_pid_threshold(X_std, -4, lower = T)
rank_hyp_median <- rank_chr(X_std)
pid_top_hyp <- get_pid_toprank(rank_hyp_median, 4)
pid_union <- union(pid_threshold, pid_top_hyp)
list_pid <- get_list_pid(pid_union)
tab <- create_table(list_pid)
fisher <- fisher.test(tab)
chisq <- chisq.test(tab)
print(tab)
print(fisher)
print(chisq)


# Investigate top ratios: No pattern
sorted_ratio <- apply(ratio_within, 2, sort, decreasing=TRUE)
sorted_ratio[1:3, d0_top3_1]
sorted_ratio[1:3, not_d0_top3_1]

# Plot PCA
plotPCA3DYeoh(hyperdiploid, metadata_df)
plotPCA3DYeoh(hypdip_no_chrY, metadata_df)
sort(table(ps_chr))

i <- 2
hist(list_chr_hypdip[[i]][,2], breaks = 30)

# Investigate hyperdiploid relapse!
metadata_df[colnames(hyperdiploid), "label", drop=F]
hyp_relapse <- c("P038_D0", "P115_D0", "P129_D0", "P164_D0", "P189_D0")
top5_relapse <- top_5[,hyp_relapse]
top5_remission <- top_5[, !(colnames(top_5) %in% hyp_relapse)]
xtable(t(top5_relapse))

# CNV data ----------------------------------------------------------------
# hyp_pid <- substring(colnames(hyperdiploid)[1:38], 1, 4)
# ID_RPATH <- "data/GSE67684/processed/metadata/lab_id.tsv"
# id_annot <- read.table(ID_RPATH, header = T, sep = "\t")
# hyp_annot <- id_annot[id_annot$pid %in% hyp_pid,]
# HYP_WPATH <- "~/Dropbox/temp/hyperdiploid_id.tsv"
# write.table(hyp_annot, HYP_WPATH, quote = F, sep = "\t", row.names = F)

CNV_RPATH <- "data/GSE67684/processed/hyperdiploid/hyperdiploid-cnv.txt"
raw_cnv <- read.table(CNV_RPATH, header = T, sep = "\t", row.names = 1,
                  strip.white = T)
raw_cnv[raw_cnv == "UPD"] <- 0 # Uniparental disomy
cnv <- trimws(as.matrix(raw_cnv))
class(cnv) <- "numeric"
colnames(cnv) <- substring(colnames(cnv), 4)

listCNV <- function(row) {
  two <- names(row)[row == 2]
  one <- names(row)[row == 1]
  zero <- names(row)[row == 0]
  neg_one <- names(row)[row == -1]
  list(two = two, one = one, zero = zero, neg_one = neg_one)
}

list_cnv <- apply(cnv, 1, listCNV)
paste_chr <- function(list) {
  lapply(list, function(vec) do.call(paste,
                                     c(as.list(vec), sep = ", ")))
}
list_concat <- lapply(list_cnv, paste_chr)
cnv_summary <- data.frame(sapply(list_concat, as.character))
pid_cnv <- paste(colnames(cnv_summary), "D0", sep = "_")
## Only show extra chr and hide P154
processed_cnv <- t(cnv_summary[1:2, colnames(cnv_summary) != "P154"])
print(processed_cnv)

## Cytogenetic data
raw_cyto_cnv <- merge(final_cyto, processed_cnv, by = "row.names")
rownames(raw_cyto_cnv) <- raw_cyto_cnv[,1]
cyto_cnv <- raw_cyto_cnv[,-1]
colnames(cyto_cnv) <- c("Cytogenetics", "CNV: Extra 2", "CNV: Extra 1")
print(xtable(cyto_cnv), type = "latex",
      file = "~/Dropbox/temp/cyto_cnv.txt")

## Mean
mean_top10 <- substring(ranked_mean[1:10, pid_cnv], 4) # remove char "chr"
colnames(mean_top10) <- substring(colnames(mean_top10), 1, 4)
cnv_mean_top10 <- rbind(cnv_summary, mean_top10)
cnv_mean_top10 <- t(cnv_mean_top10[
  c(1:2, 5:14), colnames(cnv_mean_top10) != "P154"])
colnames(cnv_mean_top10) <- c("Extra 2", "Extra 1", "1st", "2nd", "3rd",
                              paste0(4:10, "th"))
xtable(cnv_mean_top10)

## Median (no filtering)
subset_top_10_1 <- substring(ranked_chr_median[1:10, pid_cnv], 4) # remove char "chr"
colnames(subset_top_10_1) <- substring(colnames(subset_top_10_1), 1, 4)
cnv_top_1 <- rbind(cnv_summary, subset_top_10_1)
cnv_top10_1 <- t(cnv_top_1[c(1:2, 5:14), colnames(cnv_top_1) != "P154"])
colnames(cnv_top10_1) <- c("Extra 2", "Extra 1", "1st", "2nd", "3rd",
                           paste0(4:10, "th"))
xtable(cnv_top10_1)

## Median (filtered zeros)
ranked_nozero_median
# Remove char "chr"
nozero_median_top10 <- substring(ranked_nozero_median[1:10, pid_cnv], 4)
colnames(nozero_median_top10) <- substring(colnames(nozero_median_top10), 1, 4)
cnv_nozero_median <- rbind(cnv_summary, nozero_median_top10)
cnv_nozero_median <- t(cnv_nozero_median[
  c(1:2, 5:14), colnames(cnv_nozero_median) != "P154"])
colnames(cnv_nozero_median) <- c("Extra 2", "Extra 1", "1st", "2nd", "3rd",
                                 paste0(4:10, "th"))
xtable(cnv_nozero_median)

### No CNV results!
## Mean
not_cnv_mean <- substring(
  ranked_mean[1:10, !colnames(ranked_mean) %in% pid_cnv], 4) # remove char "chr"
not_cnv_mean <- t(not_cnv_mean)
colnames(not_cnv_mean) <- c("1st", "2nd", "3rd",
                               paste0(4:10, "th"))
rownames(not_cnv_mean) <- substring(rownames(not_cnv_mean), 1, 4)
xtable(not_cnv_mean)

## Median
median_not_cnv <- t(substring(
  ranked_nozero_median[1:10, !colnames(ranked_nozero_median) %in% pid_cnv],
  4)) # remove char "chr"
colnames(median_not_cnv) <- c("1st", "2nd", "3rd",
                              paste0(4:10, "th"))
rownames(median_not_cnv) <- substring(rownames(median_not_cnv), 1, 4)
xtable(median_not_cnv)

## CYTOGENETIC DATA -----
## ## Parse cytogenetic data
## CYTO_RPATH <- "data/GSE67684/raw/hyperdiploid/cytogenetics.tsv"
## raw_cyto <- read.table(CYTO_RPATH, header = T, sep = "\t", stringsAsFactors = F)
## list_raw <- strsplit(raw_cyto$Cytogenetics, "/")
## clone1 <- sapply(list_raw, function(x) x[1])
## clone2 <- sapply(list_raw, function(x) x[2])
## clone3 <- sapply(list_raw, function(x) x[3])

## list_clone1 <- strsplit(clone1, ",")
## n_chr1 <- sapply(list_clone1, function(x) x[1])
## sex_chr1 <- sapply(list_clone1, function(x) x[2])
## extra_chr1 <- sapply(list_clone1, function(x) x[-c(1,2)])
## cyto1 <- as.character(sapply(extra_chr1,
##                              function(x) do.call(paste, c(as.list(x), sep = ", "))))
## cyto1[cyto1 == "character(0)"] <- NA

## list_clone2 <- strsplit(clone2, ",")
## n_chr2 <- sapply(list_clone2, function(x) x[1])
## sex_chr2 <- sapply(list_clone2, function(x) x[2])
## extra_chr2 <- sapply(list_clone2, function(x) x[-c(1,2)])
## cyto2 <- as.character(sapply(extra_chr2,
##                              function(x) do.call(paste, c(as.list(x), sep = ", "))))
## cyto2[cyto2 == "character(0)"] <- NA

## list_clone3 <- strsplit(clone3, ",")
## n_chr3 <- sapply(list_clone3, function(x) x[1])
## sex_chr3 <- sapply(list_clone3, function(x) x[2])
## extra_chr3 <- sapply(list_clone3, function(x) x[-c(1,2)])
## cyto3 <- as.character(sapply(extra_chr3,
##                              function(x) do.call(paste, c(as.list(x), sep = ", "))))
## cyto3[cyto3 == "character(0)"] <- NA

## ## Write organised table
## cyto_clones <- cbind(raw_cyto[,1:2], n_chr1, sex_chr1, cyto1,
##                      n_chr2, sex_chr2, cyto2, n_chr3, sex_chr3, cyto3)
## rownames(cyto_clones) <- cyto_clones$pid
## cyto_clones <- cyto_clones[,-1]
## rownames(cyto_clones) <- raw_cyto$pid
## CYTO_WPATH <- "data/GSE67684/processed/hyperdiploid/cytogenetics_clones.tsv"
## write.table(cyto_clones, CYTO_WPATH, quote = F, sep = "\t")

## ## Filter out chromosomal aberrations
## ## Select char starting with +[0-9|X]
## curated_extra_chr1 <- lapply(extra_chr1,
##                              function(x) x[grepl("^\\+[0-9|X|Y]", x)])
## ## Delete square brackets at end of char
## deleteBracket <- function(vec) {
##   idx <- grepl("]$", vec)
##   if (!any(idx)) return(vec)
##   edited_vec <- sapply(vec[idx], gsub, pattern = "\\[.*\\]$", replacement = "")
##   vec[idx] <- edited_vec
##   return(vec)
##

## edited_extra_chr1 <- lapply(curated_extra_chr1, deleteBracket)
## ## Delete mar
## edited1_extra_chr1 <- lapply(edited_extra_chr1,
##                              function(x) x[!grepl("mar", x)])
## ## Named pid and get rid of patients with no info
## names(edited1_extra_chr1) <- raw_cyto$pid
## edited2_extra_chr1 <- Filter(function(x) length(x) != 0, edited1_extra_chr1)
## ## Get rid of duplicates and + sign
## edited3_extra_chr1 <- lapply(edited2_extra_chr1,
##                              function(x) substring(unique(x), 2))
## XCHR_WPATH <- "data/GSE67684/processed/hyperdiploid/cytogenetics_xchr.RDS"
## saveRDS(edited3_extra_chr1, XCHR_WPATH)

## Read cyto table and list of extra chr
CYTO_RPATH <- "data/GSE67684/processed/hyperdiploid/cytogenetics_clones.tsv"
cyto_clones <- read.table(CYTO_RPATH, header = T, sep = "\t")
XCHR_RPATH <- "data/GSE67684/processed/hyperdiploid/cytogenetics_xchr.RDS"
extra_chr <- readRDS(XCHR_RPATH)

## Cleaned cytogenetics table
processed_cyto <- sapply(extra_chr,
                         function(x) do.call(paste, c(as.list(x), sep = ", ")))
final_cyto <- data.frame(extra_chr = processed_cyto)
final_cyto1 <- cbind(cyto_clones[rownames(final_cyto), 1:2], final_cyto)
## print(xtable(final_cyto1), type = "latex",
##       file = "~/Dropbox/temp/final_cyto.txt")

## TABLE: CYTOGENETICS COMPARISON
## Median (filtered)
# Remove char "chr"
pid_cyto <- rownames(final_cyto1)
subset_nozero_median_top10 <- substring(ranked_nozero_median[1:10, pid_cyto], 4)
incyto_nozero_median <- sapply(colnames(subset_nozero_median_top10),
                               function(pid) sum(subset_nozero_median_top10[,pid] %in% extra_chr[[pid]]))
avg_incyto_nozero_median <- mean(incyto_nozero_median)
cyto_nozero_median_top10 <- cbind(final_cyto, t(subset_nozero_median_top10),
                                  incyto_nozero_median)
colnames(cyto_nozero_median_top10) <- c("Extra chr", "1st", "2nd", "3rd",
                                        paste0(4:10, "th"), "N")
print(xtable(cyto_nozero_median_top10), type = "latex",
      file = "~/Dropbox/temp/cyto_nozero_median.txt")

## Mean (filtered)
subset_mean_top10 <- substring(ranked_mean[1:10, pid_cyto], 4)
incyto_mean <- sapply(colnames(subset_mean_top10),
                      function(pid) sum(subset_mean_top10[,pid] %in% extra_chr[[pid]]))
avg_incyto_mean <- mean(incyto_mean)
cyto_mean_top10 <- cbind(final_cyto, t(subset_mean_top10), incyto_mean)
colnames(cyto_mean_top10) <- c("Extra chr", "1st", "2nd", "3rd",
                               paste0(4:10, "th"), "N")
print(xtable(cyto_mean_top10), type = "latex",
      file = "~/Dropbox/temp/cyto_mean.txt")

## Evaluate CNV using cyto
## Ignore sex chr
list_cnv_extra <- lapply(list_cnv,
                         function(l) c(l$two, l$one))
auto_extra_chr <- lapply(extra_chr,
                         function(x) setdiff(x, c("X", "Y")))
auto_cnv_extra <- lapply(list_cnv_extra,
                         function(x) setdiff(x, c("X", "Y")))
recall_cnv <- sapply(rownames(cyto_cnv),
                     function(pid) {
                       sum(auto_extra_chr[[pid]] %in% auto_cnv_extra[[pid]])/length(auto_extra_chr[[pid]])
                     })
precision_cnv <- sapply(rownames(cyto_cnv),
                        function(pid) {
                          sum(auto_cnv_extra[[pid]] %in% auto_extra_chr[[pid]])/length(auto_cnv_extra[[pid]])
                        })
mean(precision_cnv)

### INVESTIGATE BATCH EFFECTS
table(metadata_df$subtype, metadata_df$batch_info)
