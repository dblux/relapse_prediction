library(reshape2)
library(dplyr)

## Plotting
library(ggplot2)
library(cowplot)
library(rgl)
library(RColorBrewer)
library(pheatmap)
# library(Rtsne)
# library(dendextend)

## Custom
source("../functions.R")
source("bin/bcm.R")

# theme_set(theme_dark())
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
  ## ERM1 / projection of D0-N on L-N
  erm1_ratio <- erm1/d0_normal_proj
  
  d8_normal_vstack <- normal_centroid - t(d8_df)
  ### D8-Normal projection ###
  d8_normal_proj <- colSums(d8_normal_vstack * unit_leuk_normal)
  
  stopifnot(identical(names(erm1), names(erm1_ratio)))
  
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
  
  ### Concatenate all features ###
  features_df <- cbind(erm1, erm1_ratio, erm2, erm2_ratio, erm3, erm3_ratio,
                       d0_normal_proj, d8_normal_proj,
                       l2norm_d0_d8, diff_l2norm, angle_d0_d8,
                       angle_d0d8_normal, angle_d0_normal, angle_d8_normal)
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
