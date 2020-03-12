library(reshape2)
library(ggplot2)
library(cowplot)
library(rgl)

source("../functions.R")

theme_set(theme_cowplot())

# IMPORT DATA -------------------------------------------------------------
MAQC_RPATH <- "../diff_expr/data/MAQC-I/processed/mas5_original-ref.tsv"
raw_maqc <- read.table(MAQC_RPATH, sep = "\t", header = T, row.names = 1)

# MAQC metadata
batch_info <- as.factor(rep(1:6, each = 10))
class_info <- rep(rep(LETTERS[1:2], each = 5), 6)
metadata_df <- data.frame(batch_info, class_info)
rownames(metadata_df) <- colnames(raw_maqc)

# SCALE->REMOVE->FILTER->LOG
scaled_maqc <- removeProbesets(normaliseMeanScaling(raw_maqc))
filtered_maqc <- filterProbesets(scaled_maqc, 0.7, metadata_df)
log_maqc <- log2_transform(filtered_maqc)

# Subset 2 batches
subset_maqc <- log_maqc[,1:20]
filtered_subset <- filterProbesets(subset_maqc, 0.5, metadata_df)

# INITIAL TESTS -----------------------------------------------------------
row_means <- rowMeans(filtered_subset)
row_sd <- apply(filtered_subset, 1, sd)
row_cv <- row_means/row_sd

row_idx1 <- names(sort(row_cv)[1:10])
subset_low_mean <- filtered_subset[row_idx1,]

plotMulti <- function(X) {
  cat(dim(X))
  par(mfrow=c(4,3))
  n_feat <- nrow(X)
  n <- ncol(X)
  for (i in 1:n_feat) {
    plot(1:n, X[i,], main = rownames(X)[i])
  }
  par(mfrow=c(1,1))
}

plotMulti(subset_low_mean)
head(subset_low_mean)
