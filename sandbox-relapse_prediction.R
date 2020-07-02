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
# library(Rtsne)
# library(dendextend)

## Custom
source("../functions.R")
source("bin/bcm.R")

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

# QPSP --------------------------------------------------------------------
# 1. Removes affymetrix ambiguous and control probesets
# 2. Map probesets to IDs
# Removes one-to-many probesets and probesets with no ID
# Selects maximum if two probesets match to same gene
# CHECK: What microarray platform is the data from?
ENTREZ_GPL570 <- "../info/microarray/HG-U133_Plus_2/annot_entrez-GPL570.tsv"
entrez_yeoh <- affy2id(data_yeoh, ENTREZ_GPL570)

## SYMBOL_GPL570 <- "../info/microarray/HG-U133_Plus_2/annot_genesymbol-GPL570.tsv"
## symbol_yeoh <- affy2id(data_yeoh, SYMBOL_GPL570)
## head(symbol_yeoh)

# Import NEA
NEA_RPATH <- paste0("../diff_expr/data/subnetwork/nea-hsa/",
                    "ovarian_cancer/geneset-nea_kegg_ovarian.tsv")
nea_df <- read.table(NEA_RPATH, sep = "\t", header = T, stringsAsFactors = F)
subnetwork_nea <- split(as.character(nea_df$gene_id), nea_df$subnetwork_id)

# Calculate QPSP profiles
gfs_yeoh <- normaliseGFS(entrez_yeoh, num_intervals = 4)
qpsp_yeoh <- calcQPSP(gfs_yeoh, subnetwork_nea)
X_qpsp <- calcQPSP(entrez_yeoh, subnetwork_nea) # No discrete-GFS normalisation
## dim(qpsp_yeoh)
## dim(subset_yeoh)

# No need for feature selection as features have been reduced
## OPTION 1
chosen_subtypes <- c("BCR-ABL", "T-ALL", "TEL-AML1")
## OPTION 2
all_subtypes <- levels(metadata_df$subtype)
## Some subtypes do not have same class samples in every batch (CANNOT USE BCM)
suitable_subtypes <- all_subtypes[-c(4:6)]
print(suitable_subtypes)
metadata_df[metadata_df$subtype == "Hypodiploid",]
table(metadata_df$subtype, metadata_df$batch_info)

for (subtype in suitable_subtypes) {
  print(c("Subtype:", subtype))
  normal_pid <- paste0("N0", c(1,2,4))
  
  logi_idx <- rownames(metadata_df) %in% colnames(qpsp_yeoh) &
    metadata_df$subtype == subtype
  subtype_pid <- rownames(metadata_df)[logi_idx]
  subset_pid <- c(subtype_pid, normal_pid)
  
  # Subtype and normal samples
  subtype_qpsp <- qpsp_yeoh[,subset_pid]
  print(subtype_qpsp[1:5,1:5])
  
  # BCM on QPSP transformed data
  bcm_qpsp <- correctGlobalBCM(subtype_qpsp, metadata_df)
  print(bcm_qpsp[1:5,1:5])
  
  # PCA
  pca_obj <- prcomp(t(bcm_qpsp))
  # PCA: Eigenvalues
  eigenvalues <- (pca_obj$sdev)^2
  var_pc <- eigenvalues/sum(eigenvalues)
  ## print(c("Variance:", var_pc))
  cum_var <- cumsum(var_pc)
  # Identify PC that has just above 70% variance
  pc_ind <- which.max(cum_var > 0.70)
  print(c("PC: ", pc_ind))
  
  # PCA: Coordinates
  pca_coord <- pca_obj$x[,1:pc_ind]
  # Response df and normal df
  n_idx <- nrow(pca_coord) - 3
  response_df <- pca_coord[1:n_idx,]
  normal_df <- pca_coord[-(1:n_idx),]
  print("Leukemia:")
  print(rownames(response_df))
  print("Normal:")
  print(rownames(normal_df))
  
  # Collate MRD results as well
  results_df <- calcERM(response_df, normal_df)
  SUBTYPE_WPATH <- sprintf("temp/pca_bcm_qpsp/pca_bcm_qpsp-%s.tsv", subtype)
  write.table(results_df, SUBTYPE_WPATH, sep = "\t", quote = F)
}

## ## Calculate likelihood ratio of feature
## ## Features: ERM1, l2norm_d0_d8, d0_normal_proj, angle_d0_d8
## features <- c("erm1", "l2norm_d0_d8", "d0_normal_proj", "angle_d0_d8")

## FEATURES_DIR <- "temp/pca_bcm_qpsp"
## list_rpaths <- list.files(FEATURES_DIR, full.names = TRUE)
## features_rpath <- list_rpaths[1]
## print(features_rpath)
## features_df <- read.table(features_rpath, sep = "\t",
##                           header = T, row.names = 1)
## print(features_df)
## X <- features_df[,features]
## y <- metadata_df[rownames(X), "label"]

## ## Investigate metadata_df
## metadata_df[metadata_df$subtype == "BCR-ABL",]

## # Gaussian Naive Bayes ----------------------------------------------------
## # Split according to labels for calculation of parameters
## list_X <- split.data.frame(X, y)
## # MOM estimate of mean
## X_mu <- sapply(list_X, colMeans)
## # Biased estimate of standard deviation
## X_sigma <- sapply(list_X, apply, 2, sd)

## SPLIT <- 0.3
## class_n <- c(10,6)
## # class_n <- sapply(list_X, nrow)
## test_n <- ceiling(class_n*SPLIT)
## test_id <- mapply(function(x,y) sample(class_n, test_n)
## print(test_id)

## # Calculate likelihood ratio for each feature of each sample
## #' @return vector of likelihood ratios for every feature
## calcLR <- function(x_vec) {
##   #' Calculates likelihood ratio for single feature
##   calcSingleLR <- function(idx) {
##     x <- x_vec[idx]
##     dnorm(x, X_mu[idx,2], X_sigma[idx,2])/
##       dnorm(x, X_mu[idx,1], X_sigma[idx,1])
##   }
##   feature_names <- names(x_vec)
##   sapply(feature_names, calcSingleLR, USE.NAMES = F)
## }

## # Predict likelihood ratios
## print(colnames(X_mu)) # Column name
## print(X)
## likelihood_ratios <- t(apply(X, 1, calcLR))
## X_y <- cbind(likelihood_ratios, y)
## X_y
## write.table(X_y,
##             "dump/sampled_Others-lr.tsv",
##             quote = F, sep = "\t")
## # # Plot
# par(mfrow=c(5,1))
# par(mar=rep(1,4))
# for(i in 1:5) {
#   plot(likelihood_ratios[,i], col = y+1)
# }
#
# # Experiment with averages of LR
# avg_lr <- rowMeans(likelihood_ratios[,1:2])
# plot(avg_lr, col = y+1)
# par(mfrow=c(1,1))
#
# plot(data.frame(likelihood_ratios[,1:4]), col = y+1)
#
# beta <- c(0.6,0.2,0.1,0.1)
# weighted_lr <- sweep(likelihood_ratios[,1:4], 2, beta, "*")
# plot(rowSums(weighted_lr), col = y+1)
#
# X_y[,5][X_y[,5] == 0] <- "remission"
# X_y[,5][X_y[,5] == 1] <- "relapse"
# model.matrix(X_y[,1:4])

# Logistic regression to assign betas

# ## PLOT
# pc_labels <- sprintf("PC%d (%.2f%%)", 1:3, (var_pc*100)[1:3])
# color1 <- rep("tomato3", 81)
# color1[c(10,49)] <- "darkolivegreen3"
# pch1 <- rep(21:23, c(39,39,3))
# plotPCA3D(pca_coord[,1:3], color1, pch1, pc_labels)

# par(mar=c(8,3,2,3))
# # Standardise each feature
# calibrated_features <- cbind(
#   results_df[,1:10],
#   results_df[,11:14]/60, # angle
#   -log10(results_df[,15])/2 # mrd
# )

# # Plotting parameters
# line_color <- ifelse(results_df[,16] == 1, "tomato3", "darkolivegreen3")
# feature_names <- c(colnames(results_df)[1:11], "ang_d0d8_normal",
#                    colnames(results_df)[13:14], "-log10(mrd)")
#
# # Plot parallel coordinates
# plot(0, xlim = c(1,15), ylim = c(-2,3),
#      type = "n", xaxt = "n", ann = F)
# axis(1, at = 1:15, feature_names, cex = 0.5)
# for (r in 1:nrow(calibrated_features)) {
#   lines(1:15, calibrated_features[r,], col = line_color[r])
# }
#
# subtype_parallel <- recordPlot()
# save_fig(subtype_parallel, sprintf("dump/parallel-%s.pdf", subtype),
#          width = 10, height = 7)

## ### DOWNSAMPLE HETEROGENEOUS subtypes ###
## # Sample size = 40 (33 remission, 7 relapse)
## subtype_metadata <- full_metadata_df[colnames(subtype_qpsp),]
## subtype_metadata1 <- subtype_metadata[endsWith(rownames(subtype_metadata), "_D0"),]
## pid_0 <- rownames(subtype_metadata1)[subtype_metadata1$label == 0]
## pid_1 <- rownames(subtype_metadata1)[subtype_metadata1$label == 1]

## sample_pid_0 <- sample(pid_0, 33)
## sample_pid_1 <- sample(pid_1, 7)
## sample_pid <- substring(c(sample_pid_0, sample_pid_1), 1, 4)

## logi_idx <- substring(colnames(subtype_qpsp), 1, 4) %in% sample_pid
## logi_idx[171:173] <- TRUE
## sampled_qpsp <- subtype_qpsp[,logi_idx]
## writeLines(colnames(sampled_qpsp), con = "dump/sampled_pid.txt")

## # Check that there are D0 samples in batch 2
## table(metadata_df[colnames(sampled_qpsp),])
## ###


# Visualisation -----------------------------------------------------------
## Quantile normalistion
quantile_yeoh <- normaliseQuantile(data_yeoh)

## CS-Quantile
colnames(data_yeoh)
quantile_d0 <- normaliseQuantile(subset_yeoh[, 1:201])
quantile_d8 <- normaliseQuantile(subset_yeoh[, 202:402])
quantile_normal <- normaliseQuantile(subset_yeoh[, 403:405])
csquantile_yeoh <- cbind(quantile_d0, quantile_d8, quantile_normal)

plotPCA3DYeoh(data_yeoh, metadata_df)
rgl.postscript("dump/pca_3d-sfl.pdf", "pdf")
ggsave("dump/jitter-sfl.pdf",
       plotJitterYeoh(data_yeoh, metadata_df),
       width = 10, height = 6)

plotPCA3DYeoh(csquantile_yeoh, metadata_df)
rgl.postscript("dump/pca_3d-cs_quantile.pdf", "pdf")
ggsave("dump/jitter-csquantile.pdf",
       plotJitterYeoh(csquantile_yeoh, metadata_df),
       width = 10, height = 6)

plotPCA3DYeoh(quantile_yeoh, metadata_df)
rgl.postscript("dump/pca_3d-quantile.pdf", "pdf")
ggsave("dump/jitter-quantile.pdf",
       plotJitterYeoh(quantile_yeoh, metadata_df),
       width = 10, height = 6)

plotPCA3DYeoh(gfs_yeoh, metadata_df)
rgl.postscript("dump/pca_3d-gfs.pdf", "pdf")
ggsave("dump/jitter-gfs.pdf",
       plotJitterYeoh(gfs_yeoh, metadata_df),
       width = 10, height = 6)

plotPCA3DYeoh(qpsp_yeoh*100, metadata_df)
rgl.postscript("dump/pca_3d-qpsp.pdf", "pdf")
jitter_qpsp <- plotJitterYeoh(qpsp_yeoh, metadata_df)
jitter_qpsp
ggsave("dump/jitter-qpsp.pdf", jitter_qpsp,
       width = 10, height = 6)

# Global QPSP -------------------------------------------------------------
# Discard PC2 and PC4
pca_obj <- prcomp(t(qpsp_yeoh))
eigenvalues <- (pca_obj$sdev)^2
var_pc <- eigenvalues/sum(eigenvalues)
## print(c("Variance:", var_pc))
cum_var <- cumsum(var_pc)
pc_ind <- which.max(cum_var > 0.70) # Identify PC theshold with >70% var
print(c("PC: ", pc_ind))

# PCA: Coordinates
pc_batch <- c(2,4) # Batch affected PCs!
selected_pc <- setdiff(1:10, pc_batch)
pca_coord <- pca_obj$x[,selected_pc]

# Response df and normal df
n_idx <- nrow(pca_coord) - 3
response_df <- pca_coord[1:n_idx,]
normal_df <- pca_coord[-(1:n_idx),]
print("Leukemia:")
print(rownames(response_df))
print("Normal:")
print(rownames(normal_df))

# Collate MRD results as well
results <- calcERM(response_df, normal_df)

SUBTYPE_WPATH <- sprintf("temp/remove_batch_pc/qpsp_pca-%s.tsv", subtype)
SUBTYPE_WPATH <- sprintf("temp/remove_batch_pc/global_qpsp.tsv")
write.table(results, SUBTYPE_WPATH, sep = "\t", quote = F)


globalqpsp_prediction <- plotPrediction(results, metadata_df, yeoh_label)
ggsave("dump/prediction-globalqpsp.pdf",
       globalqpsp_prediction, width = 12, height = 6)

# Subtype QPSP ------------------------------------------------------------
# Plot jitter
plotJitterYeoh <- function(X, metadata_df, n_pc = 10) {
  pca_obj <- prcomp(t(X))
  X_pca <- data.frame(pca_obj$x)
  batch <- as.factor(metadata_df[rownames(X_pca),"batch_info"])
  class <- as.factor(metadata_df[rownames(X_pca),"class_info"])
  X_meta <- cbind(batch, class, X_pca[,1:n_pc])
  X_long <- melt(X_meta, ids = c("batch", "class"), variable.name="PC")
  
  ax_batch <- ggplot(X_long, aes(x=PC, y=value)) +
    geom_boxplot(aes(fill=batch), alpha=0.3, outlier.shape=NA) +
    geom_point(aes(colour=batch, pch=class), position=position_jitterdodge(),
               size = 2, alpha = 1.0)
  
  ax_class <- ggplot(X_long, aes(x=PC, y=value)) +
    geom_boxplot(aes(fill=class), alpha=0.3, outlier.shape=NA) +
    geom_point(aes(colour=class), position=position_jitterdodge(),
               size = 2, alpha = 1.0)
  
  fig <- plot_grid(ax_batch, ax_class, nrow = 2)
  return(fig)  
}

normal_pid <- paste0("N0", c(1,2,4))
all_subtypes <- levels(metadata_df$subtype)
subtypes <- setdiff(all_subtypes, c("Hypodiploid", "Normal"))
# List of PCs to remove for each subtype
remove_pc <- list(1, c(1,2), numeric(), c(2,4), c(1,2), c(2,3), 2)
names(remove_pc) <- subtypes
print(remove_pc)
for (subtype in names(remove_pc)) {
  print(c("Subtype:", subtype))
  logi_idx <- rownames(metadata_df) %in% colnames(qpsp_yeoh) &
    metadata_df$subtype == subtype
  subtype_pid <- rownames(metadata_df)[logi_idx]
  subset_pid <- c(subtype_pid, normal_pid)
  
  # Subtype and normal samples
  subtype_qpsp <- qpsp_yeoh[,subset_pid]
  # print(subtype_qpsp[1:5,1:5])
  
  #   jitter_subtype <- plotJitterYeoh(subtype_qpsp, metadata_df)
  #   JITTER_WPATH <- sprintf("dump/jitter-%s.pdf", subtype)
  #   ggsave(JITTER_WPATH, jitter_subtype, width = 12, height = 6)
  #   jitter_subtype # Identify PCs with batch effects!
  
  pca_obj <- prcomp(t(subtype_qpsp))
  eigenvalues <- (pca_obj$sdev)^2
  var_pc <- eigenvalues/sum(eigenvalues)
  ## print(c("Variance:", var_pc))
  cum_var <- cumsum(var_pc)
  pc_ind <- which.max(cum_var > 0.70) # Identify PC theshold with >70% var
  # print(c("PC: ", pc_ind))
  
  # PCA: Coordinates
  pc_batch <- remove_pc[[subtype]] # Batch affected PCs!
  print(pc_batch)
  selected_pc <- setdiff(1:10, pc_batch)
  pca_coord <- pca_obj$x[,selected_pc]
  
  # Response df and normal df
  n_idx <- nrow(pca_coord) - 3
  response_df <- pca_coord[1:n_idx,]
  normal_df <- pca_coord[-(1:n_idx),]
  
  # # Collate MRD results as well
  results <- calcERM(response_df, normal_df)
  # SUBTYPE_WPATH <- sprintf("temp/remove_batch_pc/qpsp-%s.tsv", subtype)
  # write.table(results_df, SUBTYPE_WPATH, sep = "\t", quote = F)
  
  # Plot
  qpsp_prediction <- plotPrediction(results, metadata_df, yeoh_label)
  PREDICTION_WPATH <- sprintf("dump/remove_batchpc-%s.pdf", subtype)
  ggsave(PREDICTION_WPATH, qpsp_prediction, width = 12, height = 7)
}

# Feature selection (Wilson) ----------------------------------------------
#' Select top percentage genes with highest var
selectTopGenes <- function(X, f, percent) {
  list_group_df <- split.default(X, f, drop = F)
  group_means <- sapply(list_group_df, rowMeans)
  # Calculate group means for every gene and var of means
  group_means_var <- apply(group_means, 1, var)
  group_var_rank <- rank(-group_means_var) # rank desc
  top_group_genes <- names(group_var_rank)[
    group_var_rank < percent*length(group_var_rank)]
  return(top_group_genes)
}

# Only D0 samples
pid_d0 <- rownames(metadata_df)[metadata_df$class_info == "D0"]
pid_idx <- intersect(pid_d0, colnames(data_yeoh))
d0_yeoh <- data_yeoh[,pid_idx]

# Factor to split data
batch <- as.factor(metadata_df[colnames(data_yeoh), "batch_info"])
class <- as.factor(metadata_df[colnames(data_yeoh), "class_info"])
label <- as.factor(metadata_df[colnames(data_yeoh), "label"])
subtype <- as.factor(metadata_df[colnames(d0_yeoh), "label"]) # split D0 data only

print(table(batch, class))

PERCENT <- 0.1
batch_topgenes <- selectTopGenes(data_yeoh, batch, PERCENT)
class_topgenes <- selectTopGenes(data_yeoh, class, PERCENT)
label_topgenes <- selectTopGenes(data_yeoh, label, PERCENT)
subtype_topgenes <- selectTopGenes(d0_yeoh, subtype, PERCENT)

# No grouping of samples
gene_var <- apply(data_yeoh, 1, var)
var_rank <- rank(-gene_var) # rank desc
topgenes <- names(var_rank)[var_rank < PERCENT*length(var_rank)]

wilson_batch_genes <- setdiff(intersect(batch_topgenes, topgenes),
                              class_topgenes)
wilson_class_genes <- setdiff(intersect(class_topgenes, topgenes),
                              batch_topgenes)
wilson_label_genes <- setdiff(intersect(label_topgenes, topgenes),
                              batch_topgenes)
wilson_subtype_genes <- setdiff(intersect(subtype_topgenes, topgenes),
                                batch_topgenes)

class_yeoh <- data_yeoh[wilson_class_genes,]
batch_yeoh <- data_yeoh[wilson_batch_genes,]
label_yeoh <- data_yeoh[wilson_label_genes,]
subtype_yeoh <- data_yeoh[wilson_subtype_genes,]

## EVALUATION
subset_metadata <- metadata_df[,c("batch_info", "class_info")]
pheatmap(class_yeoh, col = brewer.pal(9, "Blues"),
         legend = T, border_color = NA, scale = "none",
         cluster_method = "ward.D2", cluster_rows = T, cluster_cols = T,
         show_colnames = F, show_rownames = F,
         annotation_col = subset_metadata)
heatmap <- recordPlot()
save_fig(heatmap, "dump/heatmap-class.pdf",
         width = 6, height = 6)

pheatmap(batch_yeoh, col = brewer.pal(9, "Blues"),
         legend = T, border_color = NA, scale = "none",
         cluster_method = "ward.D2", cluster_rows = T, cluster_cols = T,
         show_colnames = F, show_rownames = F,
         annotation_col = subset_metadata)
heatmap <- recordPlot()
save_fig(heatmap, "dump/heatmap-batch.pdf",
         width = 6, height = 6)

# Label genes
# subset_metadata <- metadata_df[,c("batch_info","class_info","label")]
pheatmap(label_yeoh, col = brewer.pal(9, "Blues"),
         legend = T, border_color = NA, scale = "none",
         cluster_method = "ward.D2", cluster_rows = T, cluster_cols = T,
         show_colnames = F, show_rownames = F,
         annotation_col = metadata_df)
heatmap <- recordPlot()
save_fig(heatmap, "dump/heatmap-label.pdf",
         width = 6, height = 6)

# Sort according to class info and subtype
subset_metadata <- metadata_df[colnames(subtype_yeoh),]
ord_pid <- rownames(subset_metadata)[order(subset_metadata$class_info,
                                           subset_metadata$subtype)]
ordered_subtype_yeoh <- subtype_yeoh[,ord_pid]

# Clustering using pearson correlation
# Pairwise correlation between rows (genes)
rows_cor <- cor(t(ordered_subtype_yeoh),
                use = "pairwise.complete.obs",
                method = "pearson")

pheatmap(ordered_subtype_yeoh, col = brewer.pal(9, "Blues"),
         legend = T, border_color = NA, scale = "none",
         cluster_rows = T, cluster_cols = F,
         show_colnames = F, show_rownames = F,
         clustering_distance_rows = as.dist(1 - rows_cor),
         annotation_col = metadata_df)
heatmap <- recordPlot()
save_fig(heatmap, "dump/heatmap-subtype_pearson.pdf",
         width = 6, height = 6)

## Evaluate batch effects
ggsave("dump/jitter-class.pdf",
       plotJitterYeoh(class_yeoh, metadata_df),
       width = 10, height = 6)
ggsave("dump/jitter-batch.pdf",
       plotJitterYeoh(batch_yeoh, metadata_df),
       width = 10, height = 6)

# Calculate features (without PCA)
# Class genes
idx <- 1:(ncol(class_yeoh)-3)
classgenes_response <- t(class_yeoh)[idx,]
classgenes_normal <- t(class_yeoh)[-idx,]
features1 <- calcERM(classgenes_response, classgenes_normal)
global_prediction_class <- plotPrediction(features1, metadata_df, yeoh_label)
ggsave("dump/prediction-global_classgenes.pdf", global_prediction_class,
       width = 8, height = 6)

# Batch genes
idx <- 1:(ncol(batch_yeoh)-3)
batchgenes_response <- t(batch_yeoh)[idx,]
batchgenes_normal <- t(batch_yeoh)[-idx,]
features2 <- calcERM(batchgenes_response, batchgenes_normal)
global_prediction_batch <- plotPrediction(features2, metadata_df, yeoh_label)
ggsave("dump/prediction-global_batchgenes.pdf", global_prediction_batch,
       width = 8, height = 6)

# Label genes
idx <- 1:(ncol(label_yeoh)-3)
labelgenes_response <- t(label_yeoh)[idx,]
labelgenes_normal <- t(label_yeoh)[-idx,]
features2 <- calcERM(labelgenes_response, labelgenes_normal)
global_prediction_label <- plotPrediction(features2, metadata_df, yeoh_label)
ggsave("dump/prediction-global_labelgenes.pdf", global_prediction_label,
       width = 8, height = 6)

# Subtype genes
idx <- 1:(ncol(subtype_yeoh)-3)
subtypegenes_response <- t(subtype_yeoh)[idx,]
subtypegenes_normal <- t(subtype_yeoh)[-idx,]
features2 <- calcERM(subtypegenes_response, subtypegenes_normal)
global_prediction_subtype <- plotPrediction(features2, metadata_df, yeoh_label)
ggsave("dump/prediction-global_subtypegenes.pdf", global_prediction_subtype,
       width = 8, height = 6)

## Subtype specific
normal_pid <- paste0("N0", c(1,2,4))
all_subtypes <- levels(metadata_df$subtype)
subtypes <- setdiff(all_subtypes, c("Hypodiploid", "Normal"))
X <- subtype_yeoh
for (subtype in subtypes) {
  print(c("Subtype:", subtype))
  logi_idx <- rownames(metadata_df) %in% colnames(X) &
    metadata_df$subtype == subtype
  subtype_pid <- rownames(metadata_df)[logi_idx]
  subset_pid <- c(subtype_pid, normal_pid)
  
  # Subtype and normal samples
  subset_yeoh <- X[,subset_pid]
  idx <- 1:(ncol(subset_yeoh)-3)
  response <- t(subset_yeoh)[idx,]
  normal <- t(subset_yeoh)[-idx,]
  print(colnames(subset_yeoh))
  print(rownames(response))
  
  # # Collate MRD results as well
  results <- calcERM(response, normal)
  # SUBTYPE_WPATH <- sprintf("temp/remove_batch_pc/qpsp-%s.tsv", subtype)
  # write.table(results_df, SUBTYPE_WPATH, sep = "\t", quote = F)
  
  # Plot
  prediction_parallel <- plotPrediction(results, metadata_df, yeoh_label)
  PREDICTION_WPATH <- sprintf("dump/prediction-subtypegenes_%s.pdf", subtype)
  ggsave(PREDICTION_WPATH, prediction_parallel, width = 12, height = 7)
}

# Subtype specific (PCA) to pick up genes with relevant
# signal (hopefully higher var genes are with signal)
normal_pid <- paste0("N0", c(1,2,4))
all_subtypes <- levels(metadata_df$subtype)
subtypes <- setdiff(all_subtypes, c("Hypodiploid", "Normal"))
X <- subtype_yeoh
for (subtype in subtypes) {
  print(c("Subtype:", subtype))
  logi_idx <- rownames(metadata_df) %in% colnames(X) &
    metadata_df$subtype == subtype
  subtype_pid <- rownames(metadata_df)[logi_idx]
  subset_pid <- c(subtype_pid, normal_pid)
  
  # Subtype and normal samples
  subset_yeoh <- X[,subset_pid]
  JITTER_WPATH <- sprintf("dump/jitter_subtypegenes-%s.pdf", subtype)
  ggsave(JITTER_WPATH,
         plotJitterYeoh(subset_yeoh, metadata_df),
         width = 10, height = 6)
  
  pca_obj <- prcomp(t(subset_yeoh))
  eigenvalues <- (pca_obj$sdev)^2
  var_pc <- eigenvalues/sum(eigenvalues)
  ## print(c("Variance:", var_pc))
  cum_var <- cumsum(var_pc)
  # pc_ind <- which.max(cum_var > 0.70) # Identify PC theshold with >70% var
  
  pca_coord <- pca_obj$x[,1:3]
  
  # Response df and normal df
  n_idx <- nrow(pca_coord) - 3
  response_df <- pca_coord[1:n_idx,]
  normal_df <- pca_coord[-(1:n_idx),]
  
  # # Collate MRD results as well
  results <- calcERM(response_df, normal_df)
  # SUBTYPE_WPATH <- sprintf("temp/remove_batch_pc/qpsp-%s.tsv", subtype)
  # write.table(results_df, SUBTYPE_WPATH, sep = "\t", quote = F)
  
  # Plot
  prediction <- plotPrediction(results, metadata_df, yeoh_label)
  PREDICTION_WPATH <- sprintf("dump/prediction-subtypegenes_pca_%s.pdf", subtype)
  ggsave(PREDICTION_WPATH, prediction, width = 12, height = 7)
}

GENES_DIR <- "data/yeoh_2002/README/chi_square_probesets"
filenames <- list.files(GENES_DIR, full.names = T)

# Prediction (Label genes) ------------------------------------------------
X <- data_yeoh
normal_pid <- paste0("N0", c(1,2,4))
all_subtypes <- levels(metadata_df$subtype)
subtypes <- setdiff(all_subtypes, c("Hypodiploid", "Normal"))
for (subtype in subtypes) {
  print(c("Subtype:", subtype))
  # Select genes
  d8_subtype <- subtypes_d8[[subtype]]
  d8_subtype_label <- metadata_df[colnames(d8_subtype), "label"]
  ord_d8_subtype <- d8_subtype[, order(d8_subtype_label)]
  size_remission <- sum(d8_subtype_label == 0)
  label_pvalue <- calc_ttest(ord_d8_subtype, size_remission)
  hist(label_pvalue, breaks = 20, main = subtype)
  hist_label <- recordPlot()
  HIST_WPATH <- sprintf("dump/hist_label-%s.pdf", subtype)
  save_fig(hist_label, HIST_WPATH,
           width = 6, height = 6)
  
  ALPHA <- 0.05
  label_genes <- names(label_pvalue)[label_pvalue < ALPHA & !is.na(label_pvalue)]
  print(c("No. of label genes = ", length(label_genes)))
  
  d8_label <- d8_subtype[selected_genes,]
  pheatmap(d8_label, col = brewer.pal(9, "Blues"),
           legend = T, border_color = "black", scale = "none",
           cluster_method = "ward.D2", cluster_rows = T, cluster_cols = T,
           show_colnames = F, show_rownames = F,
           annotation_col = metadata_df)
  heatmap_label <- recordPlot()
  HEATMAP_WPATH <- sprintf("dump/heatmap_label-%s.pdf", subtype)
  save_fig(heatmap_label, HEATMAP_WPATH,
           width = 10, height = 10)
  selected_genes <- setdiff(label_genes, batch_genes)
  print(c("No. of selected genes = ", length(selected_genes)))
  
  # Subset pids in subtype
  logi_idx <- rownames(metadata_df) %in% colnames(X) &
    metadata_df$subtype == subtype
  subtype_pid <- rownames(metadata_df)[logi_idx]
  subset_pid <- c(subtype_pid, normal_pid)
  
  # Subtype and normal samples
  subset_yeoh <- X[selected_genes, subset_pid]
  idx <- 1:(ncol(subset_yeoh)-3)
  response <- t(subset_yeoh)[idx,]
  normal <- t(subset_yeoh)[-idx,]
  print(colnames(subset_yeoh))
  print(rownames(response))
  
  # # Collate MRD results as well
  results <- calcERM(response, normal)
  # SUBTYPE_WPATH <- sprintf("temp/remove_batch_pc/qpsp-%s.tsv", subtype)
  # write.table(results_df, SUBTYPE_WPATH, sep = "\t", quote = F)
  
  # Plot
  prediction_parallel <- plotPrediction(results, metadata_df, yeoh_label)
  PREDICTION_WPATH <- sprintf("dump/prediction-label_%s.pdf", subtype)
  ggsave(PREDICTION_WPATH, prediction_parallel, width = 12, height = 7)
}


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
  print(c("Subtype:", subtype))
  
  # Select genes
  # X_subtype <- X_subtypes[[subtype]]
  #-- OPTION: Hyperdiploid
  X_subtype <- X[, not_pid_top3_2]
  pid <- c(not_pid_top3_2, normal_pid)
  print(colnames(X_subtype))
  #--
  class_genes <- getLocalGenes(X_subtype, pid_remission)
  print(c("No. of selected genes = ", length(class_genes)))
  # list_drug_genes <- append(list_drug_genes, list(class_genes))
  
  selected_genes <- setdiff(class_genes, batch_genes)
  print(c("No. of final genes = ", length(selected_genes)))
  
  # Subset pids in subtype
  logi_idx <- rownames(metadata_df) %in% colnames(X) &
    metadata_df$subtype == subtype
  subtype_pid <- rownames(metadata_df)[logi_idx]
  subset_pid <- c(subtype_pid, normal_pid)

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
  # TODO: Investigate
  metadata_df[rownames(results), "label", drop=F]
  results[order(results$l2norm_ratio2), "l2norm_ratio2", drop=F]
  
  # OPTION: Automatic filename
  # PREDICTION_WPATH <- sprintf("~/Dropbox/temp/prediction_top-%s.pdf", subtype)
  PREDICTION_WPATH <- "~/Dropbox/temp/prediction_not_top3_2-Hyperdiploid.pdf"
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

# Prediction (Subtype genes - Chi2) --------------------------------------------
pid_d0 <- rownames(metadata_df)[metadata_df$class_info == "D0"]
pid_remission <- rownames(metadata_df)[metadata_df$label == 0]
# Recursive intersect
pid_idx <- intersect(
  intersect(pid_remission, pid_d0),
  colnames(data_yeoh)
)
X_d0_remission <- data_yeoh[,pid_idx]
X_gfs <- normaliseGFS(X_d0_remission, upper = .05,
                      lower = .55, num_intervals = 4)
X_factor <- data.frame(apply(X_gfs, 2, factor))

X <- X_factor
# Factor to split data
subtype_factor <- as.factor(metadata_df[colnames(data_yeoh), "subtype"])
subtypes_yeoh <- split.default(data_yeoh, subtype_factor, drop = F) # Split by subtype
length(batch_genes)
normal_pid <- paste0("N0", c(1,2,4))
all_subtypes <- levels(metadata_df$subtype)
subtypes <- setdiff(all_subtypes, c("Hypodiploid", "Normal"))
ALPHA <- 0.05
list_subtype_genes <- list()
for (subtype in subtypes) {
  print(c("Subtype:", subtype))
  
  # Assigning subtype labels (one vs rest)
  subtype_info <- metadata_df[colnames(X), "subtype"]
  levels(subtype_info) <- c(levels(subtype_info), "Rest")
  subtype_info[subtype_info != subtype] <- "Rest"
  subtype_info <- droplevels(subtype_info)
  
  X_t <- t(X)
  gene_p <- numeric()
  # Loop through each gene
  for (i in 1:ncol(X_t)) {
    # Create contingency table for each gene
    single_gene <- data.frame(subtype = subtype_info,
                              value = X_t[,i])
    count_tab <- table(single_gene$value, single_gene$subtype)
    
    if (nrow(count_tab) == 1) {
      gene_p <- c(gene_p, NaN)
    } else {
      chisq_obj <- chisq.test(count_tab,
                              simulate.p.value = T)
      gene_p <- c(gene_p, chisq_obj$p.value)
    }
  }
  names(gene_p) <- colnames(X_t)
  
  # hist(gene_p, breaks = 20, main = subtype)
  # hist_p <- recordPlot()
  # HIST_WPATH <- sprintf("dump/hist_chi2_p-%s.pdf", subtype)
  # save_fig(hist_p, HIST_WPATH)
  
  subtype_genes <- names(gene_p)[gene_p < ALPHA & !is.na(gene_p)]
  print(sprintf("No. of NaNs = %d", sum(is.na(gene_p))))
  print(sprintf("No. of selected genes = %d", length(subtype_genes)))
  list_subtype_genes <- append(list_subtype_genes, list(subtype_genes))

  X_subtype <- X_d0_remission[subtype_genes,]
  heatmap_subtype <- plotHeatmapSubtype(X_subtype, metadata_df)
  HEATMAP_WPATH <- sprintf("dump/heatmap_subtype1-%s.pdf", subtype)
  save_fig(heatmap_subtype, HEATMAP_WPATH,
           width = 6, height = 6)
  
  # selected_genes <- setdiff(class_genes, batch_genes)
  # print(c("No. of final genes = ", length(selected_genes)))
  
  # # Subset pids in subtype
  # logi_idx <- rownames(metadata_df) %in% colnames(X) &
  #   metadata_df$subtype == subtype
  # subtype_pid <- rownames(metadata_df)[logi_idx]
  # subset_pid <- c(subtype_pid, normal_pid)
  # 
  # # Subtype and normal samples
  # # subset_yeoh <- X[class_genes, subset_pid] # TODO
  # subset_yeoh <- X[selected_genes, subset_pid]
  # idx <- 1:(ncol(subset_yeoh)-3)
  # response <- t(subset_yeoh)[idx,]
  # normal <- t(subset_yeoh)[-idx,]
  # print(colnames(subset_yeoh))
  # print(rownames(response))
  # 
  # # Collate MRD results as well
  # results <- calcERM(response, normal)
  # 
  # # Plot
  # prediction_parallel <- plotPrediction(results, metadata_df, yeoh_label)
  # PREDICTION_WPATH <- sprintf("dump/prediction_s5b2-%s.pdf", subtype)
  # ggsave(PREDICTION_WPATH, prediction_parallel, width = 12, height = 7)
}
names(list_subtype_genes) <- subtypes
saveRDS(list_subtype_genes, "temp/list_subtype_genes.rds")

plotHeatmapSubtype <- function(X, metadata_df) {
  subtype_factor <- metadata_df[colnames(X), "subtype"]
  set3_pal <- brewer.pal(9, "Set3")
  subtype_col <- set3_pal[subtype_factor]
  
  par(mar = c(1,1,1,1))
  heatmap(data.matrix(X),
          col = brewer.pal(9, "Blues"),
          ColSideColors = subtype_col,
          scale = "none",
          labRow = NA, labCol = NA)
  
  legend(x = "topleft", legend = levels(subtype_factor),
         col = set3_pal[factor(levels(subtype_factor))],
         pch = 15, cex = .6)
  heatmap_subtype <- recordPlot()
  par(mar = c(5.1, 4.1, 4.1, 2.1)) # Reset to default
  return(heatmap_subtype)
}

# Prediction (Subtype genes - T-test) --------------------------------------
pid_d0 <- rownames(metadata_df)[metadata_df$class_info == "D0"]
pid_remission <- rownames(metadata_df)[metadata_df$label == 0]
# Recursive intersect
pid_idx <- intersect(
  intersect(pid_remission, pid_d0),
  colnames(data_yeoh)
)
X_d0_remission <- data_yeoh[,pid_idx]
X <- data_yeoh
# Factor to split data
length(batch_genes)
normal_pid <- paste0("N0", c(1,2,4))
all_subtypes <- levels(metadata_df$subtype)
subtypes <- setdiff(all_subtypes, c("Hypodiploid", "Normal"))
ALPHA <- 0.05
THRESHOLD <- 1
list_subtype_genes <- list()
for (subtype in subtypes) {
  subtype <- subtypes[[5]]
  print(c("Subtype:", subtype))
  
  # Assigning subtype labels (one vs rest)
  subtype_info <- metadata_df[colnames(X_d0_remission), "subtype"]
  X_subtype <- X_d0_remission[,subtype_info == subtype]
  X_rest <- X_d0_remission[,subtype_info != subtype]
  p <- calc_ttest(cbind(X_subtype, X_rest), ncol(X_subtype))
  feat_p <- names(p)[p < ALPHA & !is.na(p)]
  logfc <- calc_logfc(X_subtype, X_rest)
  feat_logfc <- names(logfc)[abs(logfc) > THRESHOLD]
  feat_intersect <- intersect(feat_p, feat_logfc)
  print(sprintf("No. of features (p-value) = %d", length(feat_p)))
  print(sprintf("No. of features (log2-fc) = %d", length(feat_logfc)))
  print(sprintf("No. of features (intersect) = %d", length(feat_intersect)))
  # list_subtype_genes <- append(list_subtype_genes, list(feat_intersect))
  
  # hist(gene_p, breaks = 20, main = subtype)
  # hist_p <- recordPlot()
  # HIST_WPATH <- sprintf("dump/hist_chi2_p-%s.pdf", subtype)
  # save_fig(hist_p, HIST_WPATH)
  
  # # Plot heatmap
  # X_subtype_ttest <- X_d0_remission[feat_intersect,]
  # pheatmap(X_subtype_ttest, col = brewer.pal(9, "Blues"),
  #          legend = T, border_color = "black", scale = "none",
  #          cluster_method = "ward.D2", cluster_rows = T, cluster_cols = T,
  #          show_colnames = F, show_rownames = F,
  #          annotation_col = metadata_df)
  # heatmap_class <- recordPlot()
  # HEATMAP_WPATH <- sprintf("dump/heatmap_subtype_ttest1-%s.pdf", subtype)
  # save_fig(heatmap_class, HEATMAP_WPATH,
  #          width = 10, height = 10)
  
  # X_subtype_ttest <- X_d0_remission[feat_intersect,]  
  # heatmap_subtype <- plotHeatmapSubtype(X_subtype_ttest, metadata_df)
  # HEATMAP_WPATH <- sprintf("dump/heatmap_subtype_ttest-%s.pdf", subtype)
  # save_fig(heatmap_subtype, HEATMAP_WPATH,
  #          width = 6, height = 6)
  
  selected_genes <- setdiff(feat_intersect, batch_genes)
  print(c("No. of final genes = ", length(selected_genes)))

  # Subset pids in subtype
  logi_idx <- rownames(metadata_df) %in% colnames(X) &
    metadata_df$subtype == subtype
  subtype_pid <- rownames(metadata_df)[logi_idx]
  subset_pid <- c(subtype_pid, normal_pid)

  # Subtype and normal samples
  # subset_yeoh <- X[class_genes, subset_pid] # TODO
  subset_yeoh <- X[selected_genes, subset_pid]
  idx <- 1:(ncol(subset_yeoh)-3)
  response <- t(subset_yeoh)[idx,]
  normal <- t(subset_yeoh)[-idx,]
  print(colnames(subset_yeoh))
  print(rownames(response))

  # Collate MRD results as well
  results <- calcERM(response, normal)

  # Plot
  prediction_parallel <- plotPrediction(results, metadata_df, yeoh_label)
  PREDICTION_WPATH <- sprintf("dump/prediction_subtype_ttest-%s.pdf", subtype)
  ggsave(PREDICTION_WPATH, prediction_parallel, width = 12, height = 7)
}
str(list_subtype_genes)
names(list_subtype_genes) <- subtypes
saveRDS(list_subtype_genes, "temp/list_subtype_genes_ttest.rds")

plotHeatmapSubtype <- function(X, metadata_df) {
  subtype_factor <- metadata_df[colnames(X), "subtype"]
  set3_pal <- brewer.pal(9, "Set3")
  subtype_col <- set3_pal[subtype_factor]
  
  par(mar = c(1,1,1,1))
  heatmap(data.matrix(X),
          col = brewer.pal(9, "Blues"),
          ColSideColors = subtype_col,
          scale = "none",
          labRow = NA, labCol = NA)
  
  legend(x = "topleft", legend = levels(subtype_factor),
         col = set3_pal[factor(levels(subtype_factor))],
         pch = 15, cex = .6)
  heatmap_subtype <- recordPlot()
  par(mar = c(5.1, 4.1, 4.1, 2.1)) # Reset to default
  return(heatmap_subtype)
}

# Prediction (Drug AND Subtype genes) ------------------------------------
plotVenn <- function(vec1, vec2) {
  # Generate overlap list
  overlap_list <- calculate.overlap(list(vec1,vec2))
  # Calculate venndiagram areas
  venn_area <- sapply(overlap_list, length)
  grid.newpage()
  venn_plot <- draw.pairwise.venn(venn_area[1], venn_area[2], venn_area[3],
                                  category = c("Drug", "Subtype"),
                                  cex = 1, fontfamily = "sans",
                                  cat.cex = 1, cat.fontfamily = "sans",
                                  margin = 0.1)
  union <- (venn_area[1] + venn_area[2] - venn_area[3])
  
  print(sprintf("Drug = %d", venn_area[1]))
  print(sprintf("Subtype = %d", venn_area[2]))
  print(sprintf("Intersect = %d", venn_area[3]))
  print(sprintf("Jaccard coefficient = %.4f",
                unname(venn_area[3]/union)))
  print(strrep("=", 20))
  
  return(overlap_list[[3]])
}
str(list_drug_genes)
str(list_subtype_genes)
list_intersect <- mapply(plotVenn, list_drug_genes, list_subtype_genes)
str(list_intersect)

subtype_factor <- as.factor(metadata_df[colnames(data_yeoh), "subtype"])
subtypes_yeoh <- split.default(data_yeoh, subtype_factor, drop = F) # Split by subtype
X <- data_yeoh
for (i in 2:7) {
  genes <- list_intersect[[i]]
  subtype <- names(list_intersect)[[i]]
  X_intersect <- data_yeoh[genes,]
  # pheatmap(X_intersect, col = brewer.pal(9, "Blues"),
  #          legend = T, border_color = "black", scale = "none",
  #          cluster_method = "ward.D2", cluster_rows = T, cluster_cols = T,
  #          show_colnames = F, show_rownames = F,
  #          annotation_col = metadata_df)
  # 
  # heatmap_class <- recordPlot()
  # HEATMAP_WPATH <- sprintf("dump/heatmap_intersect-%s.pdf", subtype)
  # save_fig(heatmap_class, HEATMAP_WPATH,
  #          width = 10, height = 10)
  
  selected_genes <- setdiff(genes, batch_genes)
  print(c("No. of final genes = ", length(selected_genes)))

  # Subset pids in subtype
  logi_idx <- rownames(metadata_df) %in% colnames(X) &
    metadata_df$subtype == subtype
  subtype_pid <- rownames(metadata_df)[logi_idx]
  subset_pid <- c(subtype_pid, normal_pid)

  # Subtype and normal samples
  subset_yeoh <- X[selected_genes, subset_pid]
  idx <- 1:(ncol(subset_yeoh)-3)
  response <- t(subset_yeoh)[idx,]
  normal <- t(subset_yeoh)[-idx,]
  print(colnames(subset_yeoh))
  print(rownames(response))

  # Collate MRD results as well
  results <- calcERM(response, normal)

  # Plot
  prediction_parallel <- plotPrediction(results, metadata_df, yeoh_label)
  PREDICTION_WPATH <- sprintf("dump/prediction_intersect_ttest-%s.pdf", subtype)
  ggsave(PREDICTION_WPATH, prediction_parallel, width = 12, height = 7)
}

# Prediction (Batch genes) QPSP --------------------------------------------
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
d0_telaml1 <- X_qpsp[,pid_idx]
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

hist(aov_pvalue, breaks = 20)
# hist_pvalue <- recordPlot()
# save_fig(hist_pvalue, "dump/hist_pvalue-batch.pdf",
#          width = 6, height = 6)

# Selecting by pvalue threshold
batch_genes <- names(aov_pvalue)[aov_pvalue < 0.05 & !is.na(aov_pvalue)]
# welch_genes <- names(welch_pvalue)[welch_pvalue < 0.05 & !is.na(welch_pvalue)]
# kruskal_genes <- names(kruskal_pvalue)[kruskal_pvalue < 0.05 & !is.na(kruskal_pvalue)]
length(batch_genes)

X_batch <- data_yeoh[batch_genes,]
pheatmap(X_batch, col = brewer.pal(9, "Blues"),
         legend = T, border_color = "black", scale = "none",
         cluster_method = "ward.D2", cluster_rows = T, cluster_cols = T,
         show_colnames = F, show_rownames = F,
         annotation_col = metadata_df)
heatmap_batch <- recordPlot()
save_fig(heatmap_batch, "dump/heatmap-batch_2565.pdf",
         width = 10, height = 10)

# Prediction (Drug genes) QPSP --------------------------------------------
getLocalGenes <- function(X_subtype, pid_remission,
                          alpha = 0.05, EXPR = 6, N = 50, LOGFC = 1) {
  pid_idx <- intersect(pid_remission, colnames(X_subtype))
  print(pid_idx)
  X_subtype_remission <- X_subtype[,pid_idx, drop = F]
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

X <- X_qpsp
X_subtypes <- splitSubtype(X_qpsp, metadata_df)
length(batch_genes)
normal_pid <- paste0("N0", c(1,2,4))
all_subtypes <- levels(metadata_df$subtype)
subtypes <- setdiff(all_subtypes, c("Hypodiploid", "Normal"))
pid_remission <- rownames(metadata_df)[metadata_df$label == 0]
for (subtype in subtypes) {
  print(c("Subtype:", subtype))
  
  # Select genes
  X_subtype <- X_subtypes[[subtype]]
  class_genes <- getLocalGenes(X_subtype, pid_remission)
  print(c("No. of selected genes = ", length(class_genes)))
  # list_drug_genes <- append(list_drug_genes, list(class_genes))
  
  selected_genes <- setdiff(class_genes, batch_genes)
  print(c("No. of final genes = ", length(selected_genes)))
  
  # Subset pids in subtype
  logi_idx <- rownames(metadata_df) %in% colnames(X) &
    metadata_df$subtype == subtype
  subtype_pid <- rownames(metadata_df)[logi_idx]
  subset_pid <- c(subtype_pid, normal_pid)

  # Plot heatmap
  X_class <- X[selected_genes, subset_pid]
  pheatmap(X_class, col = brewer.pal(9, "Blues"),
           legend = T, border_color = "black", scale = "none",
           cluster_method = "ward.D2", cluster_rows = T, cluster_cols = T,
           show_colnames = F, show_rownames = F,
           annotation_col = metadata_df)
  heatmap_class <- recordPlot()
  HEATMAP_WPATH <- sprintf("dump/heatmap_qpsp_drug_wobatch-%s.pdf", subtype)
  save_fig(heatmap_class, HEATMAP_WPATH,
           width = 10, height = 10)

  # Subtype and normal samples
  subset_yeoh <- X[selected_genes, subset_pid] # OPTION
  idx <- 1:(ncol(subset_yeoh)-3)
  response <- t(subset_yeoh)[idx,]
  normal <- t(subset_yeoh)[-idx,]
  print(colnames(subset_yeoh))
  print(rownames(response))

  # Collate MRD results as well
  results <- calcERM(response, normal)

  # Plot
  prediction_parallel <- plotPrediction(results, metadata_df, yeoh_label)
  PREDICTION_WPATH <- sprintf("dump/prediction_qpsp_drug_wobatch-%s.pdf", subtype)
  ggsave(PREDICTION_WPATH, prediction_parallel, width = 12, height = 7)
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

# Prediction (Batch genes) BCM --------------------------------------------
batch <- metadata_df[colnames(data_yeoh), "batch_info"]
metadata <- metadata_df[colnames(data_yeoh), 2:3]
data_bcm <- correctMultiBCM(data_yeoh, batch, metadata, ref_batch = "2")

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
d0_telaml1 <- data_bcm[,pid_idx]
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

# ### PLOTS ###
# X_batch <- data_bcm[batch_genes,]
# pheatmap(X_batch, col = brewer.pal(9, "Blues"),
#          legend = T, border_color = "black", scale = "none",
#          cluster_method = "ward.D2", cluster_rows = T, cluster_cols = T,
#          show_colnames = F, show_rownames = F,
#          annotation_col = metadata_df)
# heatmap_batch <- recordPlot()
# save_fig(heatmap_batch, "dump/heatmap-batch_2565.pdf",
#          width = 10, height = 10)
# 
# table(metadata_df$batch_info, metadata_df$subtype)
# 
# list_batch_genes <- list(anova = batch_genes, welch = welch_genes,
#                          kruskal = kruskal_genes)
# 
# upset(fromList(list_batch_genes),
#       nsets = length(list_selected),
#       nintersects = NA,
#       order.by = "freq")
# upset_plot <- recordPlot()
# upset_plot
# save_fig(upset_plot, "dump/upset-batch_genes.pdf",
#          width = 8, height = 8)
# 
# ## Label genes
# pid_d8 <- rownames(metadata_df)[metadata_df$class_info == "D8"]
# pid_idx <- intersect(pid_d8, colnames(data_yeoh))
# d8_yeoh <- data_yeoh[,pid_idx]
# 
# subtype_factor1 <- as.factor(metadata_df[colnames(d8_yeoh), "subtype"])
# subtypes_d8 <- split.default(d8_yeoh, subtype_factor1, drop = F) # Split by subtype
# 
# plotPCA3DYeoh(data_bcm, metadata_df)
# plotPCA3DYeoh(data_yeoh, metadata_df)

# Prediction (Drug genes) BCM ---------------------------------------------
# Nonlocals: pid_remission
getLocalGenes <- function(X_subtype, pid_remission,
                          alpha = 0.05, EXPR = 6, N = 50, LOGFC = 1) {
  pid_idx <- intersect(pid_remission, colnames(X_subtype))
  print(pid_idx)
  X_subtype_remission <- X_subtype[,pid_idx, drop = F]
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
subtypes_yeoh <- splitSubtype(data_bcm, metadata_df)

length(batch_genes)
X_subtypes <- subtypes_yeoh
X <- data_bcm
normal_pid <- paste0("N0", c(1,2,4))
all_subtypes <- levels(metadata_df$subtype)
subtypes <- setdiff(all_subtypes, c("Hypodiploid", "Normal"))
pid_remission <- rownames(metadata_df)[metadata_df$label == 0]
# list_drug_genes <- list()
for (subtype in subtypes) {
  print(c("Subtype:", subtype))
  
  # Select genes
  X_subtype <- X_subtypes[[subtype]]
  class_genes <- getLocalGenes(X_subtype, pid_remission)
  print(c("No. of selected genes = ", length(class_genes)))
  # list_drug_genes <- append(list_drug_genes, list(class_genes))
  
  selected_genes <- setdiff(class_genes, batch_genes)
  print(c("No. of final genes = ", length(selected_genes)))
  
  # Subset pids in subtype
  logi_idx <- rownames(metadata_df) %in% colnames(X) &
    metadata_df$subtype == subtype
  subtype_pid <- rownames(metadata_df)[logi_idx]
  subset_pid <- c(subtype_pid, normal_pid)
  
  # X_class <- X[selected_genes, subset_pid]
  # plotPCA3DYeoh(X_class, metadata_df)
  
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
  subset_yeoh <- X[selected_genes, subset_pid] # OPTION!
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
  
  # # Plot heatmap
  # pheatmap(t(results[,subset_features]), col = brewer.pal(9, "Blues"),
  #          legend = T, border_color = "black", scale = "row",
  #          cluster_method = "ward.D2", cluster_rows = T, cluster_cols = T,
  #          show_colnames = F, show_rownames = T,
  #          annotation_col = metadata_df[,  "label", drop = F])
  # heatmap_class <- recordPlot()
  # HEATMAP_WPATH <- sprintf("dump/heatmap_features-%s.pdf", subtype)
  # save_fig(heatmap_class, HEATMAP_WPATH,
  #          width = 10, height = 10)
  
  # features <- plotFeatures(results, metadata_df)
  # FEATURES_WPATH <- sprintf("~/Dropbox/temp/features_drug-%s.pdf", subtype)
  # ggsave(FEATURES_WPATH, features, width = 16, height = 10)

  # Plot
  prediction_parallel <- plotPrediction(results, metadata_df, yeoh_label)
  PREDICTION_WPATH <- sprintf("~/Dropbox/temp/prediction_bcm_drug-%s.pdf", subtype)
  ggsave(PREDICTION_WPATH, prediction_parallel, width = 12, height = 7)
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

# Hyperdiploid classification ---------------------------------------------
## Plot: Sum of expression
# Normalised: D0 data
idx_d0 <- metadata_df[colnames(data_yeoh), "class_info"] == "D0"
sum_d0 <- colSums(data_yeoh)[idx_d0]

D <- data.frame(subtype = metadata_df[names(sum_d0), "subtype"],
                value = sum_d0)
features_plot <- ggplot(D, aes(as.factor(subtype), value, colour = subtype)) +
  geom_point(position = position_jitter(width=.1, height=0), cex = 2, show.legend = F) + # position = position_jitterdodge()
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(angle = 20, vjust = 0.5))
ggsave("~/Dropbox/temp/colsum-scaled.pdf", features_plot,
       width = 6, height = 5)

# Raw: D0 data
selected_raw <- removeProbesets(raw_yeoh)
unnorm_raw <- log2_transform(filterProbesets(selected_raw, 0.7, metadata_df))
sum_raw_d0 <- colSums(unnorm_raw)[idx_d0]

D1 <- data.frame(subtype = metadata_df[names(sum_raw_d0), "subtype"],
                 value = sum_raw_d0)
features_plot1 <- ggplot(D1, aes(as.factor(subtype), value, colour = subtype)) +
  geom_point(position = position_jitter(width=.1, height=0), cex = 2, show.legend = F) + # position = position_jitterdodge()
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(angle = 20, vjust = 0.5))
ggsave("~/Dropbox/temp/colsum-unnorm.pdf", features_plot1,
       width = 6, height = 5)

## Annotation: Chr location
ANNOT_RPATH <- "../info/microarray/HG-U133_Plus_2/affy/HG-U133_Plus_2.na35.annot.csv"
annot <- read.csv(ANNOT_RPATH,  row.names = 1, comment.char = "#")

idx_hyp <- metadata_df[colnames(data_yeoh), "subtype"] %in% c("Hyperdiploid", "Normal") & 
  metadata_df[colnames(data_yeoh), "class_info"] %in% c("D0", "N")
hyperdiploid <- data_yeoh[,idx_hyp]

ps_chrloc <- annot[rownames(data_yeoh), "Chromosomal.Location"]
ps_chr <- sub("(chr.*?)(p|q|c).*", "\\1", ps_chrloc)
ps_chr[ps_chr == "---"] <- NA

## Un-normalised data
idx <- metadata_df[colnames(data_yeoh), "subtype"] %in% c("Hyperdiploid", "Normal") & 
  metadata_df[colnames(data_yeoh), "class_info"] %in% c("D0", "N")
hyp_raw <- unnorm_raw[,idx]
ps_chrloc1 <- annot[rownames(hyp_raw), "Chromosomal.Location"]
ps_chr1 <- sub("(chr.*?)(p|q|c).*", "\\1", ps_chrloc1)
ps_chr1[ps_chr1 == "---"] <- NA
list_chr_hypdip1 <- split.data.frame(hyp_raw, ps_chr1)
hypdip_chr_mean1 <- t(sapply(list_chr_hypdip1, colMeans))
hypdip_no_chrY1 <- hypdip_chr_mean1[1:23,]

## Hyperdiploid D0: Split into chr
list_chr_hypdip <- split.data.frame(hyperdiploid, ps_chr)
hypdip_chr_mean <- t(sapply(list_chr_hypdip, colMeans))
hypdip_no_chrY <- hypdip_chr_mean[1:23,]

## Ranks
hypdip_rank <- apply(hyperdiploid, 2, rank, ties.method = "min")
list_chr_rank <- split.data.frame(hypdip_rank, ps_chr)
hypdip_chr_rank <- t(sapply(list_chr_rank, colMeans))
hypdip_rank_no_chrY <- hypdip_chr_rank[1:23,]

# Entire dataset
list_chr <- split.data.frame(data_yeoh, ps_chr)
chr_mean <- t(sapply(list_chr, colMeans))

# Plot
pheatmap(hypdip_no_chrY, col = brewer.pal(9, "Blues"),
         legend = T, border_color = NA, scale = "none",
         cluster_method = "ward.D2", cluster_rows = T, cluster_cols = T,
         show_colnames = F, show_rownames = T,
         annotation_col = metadata_df)
heatmap <- recordPlot()
save_fig(heatmap, "~/Dropbox/temp/heatmap_none-hypdip_rank_no_chrY.pdf",
         width = 7, height = 6)

## Plot scaled data
long_no_chrY <- melt(hypdip_no_chrY, varnames = c("chr", "pid"))
long_no_chrY$chr <- factor(long_no_chrY$chr,
                           levels = levels(long_no_chrY$chr)[
                             c(1,12,16:22,2:11,13:15,23)])
# Plot pid facet
chr_jitter1 <- ggplot(long_no_chrY[1:(20*23),], aes(chr, value)) +
  geom_point(position = position_jitter(width=.1, height=0),
             cex = 2, show.legend = F) +
  facet_wrap(~pid, nrow = 4, ncol = 5,  scales = "free_x") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(angle = 90, vjust = 0.5))
chr_jitter2 <- ggplot(long_no_chrY[(20*23+1):943,], aes(chr, value)) +
  geom_point(position = position_jitter(width=.1, height=0),
             cex = 2, show.legend = F) +
  facet_wrap(~pid, nrow = 4, ncol = 6,  scales = "free_x") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(angle = 90, vjust = 0.5))
ggsave("~/Dropbox/temp/jitter-chr2.pdf", chr_jitter2,
       width = 16, height = 10)

# Plot batch effects
long_batch <- cbind(long_no_chrY,
                    batch = metadata_df[long_no_chrY$pid, "batch_info"])
jitter_batch <- ggplot(long_batch, aes(pid, value, colour=pid)) +
  geom_point(position = position_jitter(width=.1, height=0),
             cex = 2, show.legend = F) +
  facet_grid(~batch, scales = "free", space = "free") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(angle = 90, vjust = 0.5))
ggsave("~/Dropbox/temp/jitter_batch-scaled.pdf", jitter_batch,
       width = 16, height = 10)

## Plot: Raw data
long_no_chrY1 <- melt(hypdip_no_chrY1, varnames = c("chr", "pid"))
long_no_chrY1$chr <- factor(long_no_chrY1$chr,
                            levels = levels(long_no_chrY1$chr)[
                              c(1,12,16:22,2:11,13:15,23)])
# Plot pid facet
chr_jitter3 <- ggplot(long_no_chrY1[1:(20*23),], aes(chr, value)) +
  geom_point(position = position_jitter(width=.1, height=0),
             cex = 2, show.legend = F) +
  facet_wrap(~pid, nrow = 4, ncol = 5,  scales = "free_x") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(angle = 90, vjust = 0.5))
chr_jitter4 <- ggplot(long_no_chrY1[(20*23+1):943,], aes(chr, value)) +
  geom_point(position = position_jitter(width=.1, height=0),
             cex = 2, show.legend = F) +
  facet_wrap(~pid, nrow = 4, ncol = 6,  scales = "free_x") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(angle = 90, vjust = 0.5))
chr_jitter3
ggsave("~/Dropbox/temp/jitter-chr2.pdf", chr_jitter2,
       width = 16, height = 10)

# Plot batch
long_batch1 <- cbind(long_no_chrY1,
                    batch = metadata_df[long_no_chrY$pid, "batch_info"])
jitter_batch1 <- ggplot(long_batch1, aes(pid, value, colour=pid)) +
  geom_point(position = position_jitter(width=.1, height=0),
             cex = 2, show.legend = F) +
  facet_grid(~batch, scales = "free", space = "free") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(angle = 90, vjust = 0.5))
ggsave("~/Dropbox/temp/jitter_batch-raw.pdf", jitter_batch1,
       width = 16, height = 10)

## TODO: Other methods of consolidating chr info aside from colMean?

## Rank within patient
rank_chr_mean <- apply(-hypdip_no_chrY, 2, rank)
print(rank_chr_mean)

# Rank: Mean and sd
chr_rank_sd <- apply(rank_chr_mean, 1, sd)
chr_rank_mean <- rowMeans(rank_chr_mean)
chr_rank_mean
chr_rank_sd
plot(chr_rank_mean, chr_rank_sd,
     xlim = c(0,25), ylim = c(0,6),
     xlab = "Mean", ylab = "SD", main = "Chromosome ranks")
text(chr_rank_mean+.8, chr_rank_sd+.15,
     names(chr_rank_mean), cex = .8)
rank_mean_sd <- recordPlot()
save_fig(rank_mean_sd, "~/Dropbox/temp/rank_scatter-chr.pdf",
         width = 7, height = 8)

# Values across samples may be affected by batch effects
# Create relative values that remain constant
# Relative to a basket of chr 1, 7, 9, 16
mean_ref <- colMeans(hypdip_no_chrY[c(1,20,22,8),])

ratio_within <- sweep(hypdip_no_chrY, 2, mean_ref, "/")
long_ratio <- melt(ratio_within, varnames = c("chr", "pid"))
long_ratio$chr <- factor(long_ratio$chr,
                         levels = levels(long_ratio$chr)[
                           c(1,12,16:22,2:11,13:15,23)])
ratio_jitter1 <- ggplot(long_ratio[1:(20*23),], aes(chr, value)) +
  geom_point(position = position_jitter(width=.1, height=0),
             cex = 2, show.legend = F) +
  facet_wrap(~pid, nrow = 4, ncol = 5,  scales = "free_x") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(angle = 90, vjust = 0.5))
ratio_jitter2 <- ggplot(long_ratio[(20*23+1):943,], aes(chr, value)) +
  geom_point(position = position_jitter(width=.1, height=0),
             cex = 2, show.legend = F) +
  facet_wrap(~pid, nrow = 4, ncol = 6,  scales = "free_x") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(angle = 90, vjust = 0.5))
ggsave("~/Dropbox/temp/jitter-ratio1.pdf", ratio_jitter1,
       width = 16, height = 10)
ratio_jitter1

# Manual identification
# Rank without chrX
ranked_chr <- apply(-ratio_within[1:22,], 2, function(x) names(sort(x)))
top_5 <- ranked_chr[1:5,]

# Subset patients
idx <- colSums(ratio_within > 1.1) >= 5 # pid
list_top5 <- as.list(data.frame(top_5[,idx]))
lit_chr <- c("chr4", "chr6", "chr10", "chr14", "chr17",
             "chr18", "chr20", "chr21", "chrX")
idx1 <- sapply(list_top5, function(x) all(x %in% lit_chr))
selected_list_top5 <- list_top5[idx1]
sort_top5 <- lapply(selected_list_top5, sort)
paste_chr5 <- sapply(sort_top5, function(vec) do.call(paste0, as.list(vec)))
table(paste_chr5)
d0_top5 <- names(paste_chr5)[paste_chr5 == "chr14chr17chr18chr21chr6"]
d8_top5 <- paste0(substring(d0_top5, 1, 6), "8")
pid_top5 <- c(d0_top5, d8_top5)

# Unique top 3 chr
top_3 <- ranked_chr[1:3,]
list_top3 <- as.list(data.frame(top_3))
idx_lit3 <- sapply(list_top3, function(x) all(x %in% lit_chr))
selected_top3 <- list_top3[idx_lit3]
sort_top3 <- lapply(selected_top3, sort)
paste_chr3 <- sapply(sort_top3, function(vec) do.call(paste0, as.list(vec)))
table(paste_chr3)
d0_top3 <- names(paste_chr3)[paste_chr3 == "chr14chr18chr21"]
d8_top3 <- paste0(substring(d0_top3, 1, 6), "8")
pid_top3 <- c(d0_top3, d8_top3)
pid_top3

# All top 3 chr have to be in literature
d0_top3_1 <- names(selected_top3)
d8_top3_1 <- paste0(substring(d0_top3_1, 1, 6), "8")
pid_top3_1 <- c(d0_top3_1, d8_top3_1)
not_d0_top3_1 <- setdiff(colnames(hyperdiploid),
                         c(d0_top3_1, normal_pid))
not_d8_top3_1 <- paste0(substring(not_d0_top3_1, 1, 6), "8")
not_pid_top3_1 <- c(not_d0_top3_1, not_d8_top3_1)

# Top 3 have to be in special set
set_popular <- paste0("chr", c(6,14,17,18,21))
logi_idx <- apply(top_3, 2, function(x) sum(x %in% set_popular) == 3)
d0_top3_2 <- colnames(top_3)[logi_idx]
d8_top3_2 <- paste0(substring(d0_top3_2, 1, 6), "8")
pid_top3_2 <- c(d0_top3_2, d8_top3_2)
not_d0_top3_2 <- setdiff(colnames(hyperdiploid),
                         c(d0_top3_2, normal_pid))
not_d8_top3_2 <- paste0(substring(not_d0_top3_2, 1, 6), "8")
not_pid_top3_2 <- c(not_d0_top3_2, not_d8_top3_2)

# Special set of chr {6,14,17,18,21} if 3 are present in top 5
set_popular <- paste0("chr", c(6,14,17,18,21))
logi_idx <- apply(top_5, 2, function(x) sum(x %in% set_popular) >= 4)
pid_popular <- colnames(top_5)[logi_idx]
pid_unpopular <- setdiff(colnames(top_5), pid_popular)
pid_unpopular
top_5[,pid_unpopular]

# Investigate top ratios: No pattern
sorted_ratio <- apply(ratio_within, 2, sort, decreasing=TRUE)
sorted_ratio[1:3, d0_top3_1]
sorted_ratio[1:3, not_d0_top3_1]
top_5

colnames(hyperdiploid)

# Unique top 2 chr
top_2 <- ranked_chr[1:2,]
list_top2 <- as.list(data.frame(top_2))
idx_lit2 <- sapply(list_top2, function(x) all(x %in% lit_chr))
selected_top2 <- list_top2[idx_lit2]
sort_top2 <- lapply(selected_top2, sort)
paste_chr2 <- sapply(sort_top2, function(vec) do.call(paste0, as.list(vec)))
table(paste_chr2)
d0_top2 <- names(paste_chr2)[paste_chr2 == "chr14chr18chr21"]
d8_top2 <- paste0(substring(d0_top2, 1, 6), "8")
pid_top2 <- c(d0_top2, d8_top2)
pid_top2

# Patients with at least 3 chr >1.1
idx1 <- colSums(ratio_within > 1.1) >= 3
d0_small <- colnames(ratio_within)[idx1]
d8_small <- paste0(substring(d0_small, 1, 6), "8")
pid_small <- c(d0_small, d8_small)

# Patients with at least 1 chr >1.2
idx2 <- colSums(ratio_within > 1.2) >= 1
d0_big <- colnames(ratio_within)[idx2]
d8_big <- paste0(substring(d0_big, 1, 6), "8")
pid_big <- c(d0_big, d8_big)

# Batch effect correction

# Plot PCA
plotPCA3DYeoh(hyperdiploid, metadata_df)
plotPCA3DYeoh(hypdip_no_chrY, metadata_df)
sort(table(ps_chr))

i <- 2
hist(list_chr_hypdip[[i]][,2], breaks = 30)
