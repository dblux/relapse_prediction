library(reshape2)
library(rgl)
library(ggplot2)
library(cowplot)
library(reshape2)
library(RColorBrewer)
library(sva)

source("../functions.R")
source("bcm.R")

theme_set(theme_cowplot())

# FUNCTIONS ---------------------------------------------------------------
plotPCA2D <- function(df1, metadata_df) {
  # PCA transformation
  pca_obj <- prcomp(t(df1))
  pca_df <- data.frame(pca_obj$x[,1:5])
  eigenvalues <- (pca_obj$sdev)^2
  var_pc <- eigenvalues[1:5]/sum(eigenvalues)
  pc_labels <- sprintf("PC%d (%.2f%%)", 1:5, var_pc*100)
  
  # Plotting annotations
  batch_factor <- as.factor(metadata_df[colnames(df1),"batch"])
  class_factor <- metadata_df[colnames(df1),"class"]

  pc1_pc2 <- ggplot(pca_df, aes(x = PC1, y = PC2,
                                shape = class_factor,
                                color = batch_factor)) +
    geom_vline(xintercept = 0, color = "black", alpha = 0.5) +
    geom_hline(yintercept = 0, color = "black", alpha = 0.5) +
    geom_point(size = 3, show.legend = F) +
    # geom_point(size = 3, fill = colour_vec, shape = shape_vec, show.legend = F) +
    labs(x = pc_labels[1], y = pc_labels[2])
    # theme(plot.title = element_text(hjust = 0.5))
  return(pc1_pc2)
}

# Corrects for each class
#' @param cov logical: include class as covariate?
correctComBat <- function(df1, metadata_df, cov) {
  # Place adjustment/confounding variables in model.matrix (e.g. age)
  # Do not put batch variables in model.matrix
  # Put batch variables directly in combat function
  subset_metadata <- metadata_df[colnames(df1), "class", drop = FALSE]
  
  ### INCLUDE CLASS AS COVARIATE OR NOT
  mod_mat <- ifelse(cov,
                    model.matrix(~class, data = subset_metadata),
                    model.matrix(~1, data = subset_metadata))
  print(head(mod_mat))
  
  # Batch information
  subset_batch <- metadata_df[colnames(df1), "batch"]
  print(head(subset_batch))
  
  # ComBat batch correction
  combat_df <- ComBat(data.matrix(df1), subset_batch, mod_mat)
  # Replacing negative values with 0
  combat_df[combat_df < 0] <- 0
  return(combat_df)  
}

#' Batch correction using class-specific ComBat
correctCSComBat <- function(df1, metadata_df) {
  class <- metadata_df[colnames(df1), "class"]
  # Splits df according to class
  list_class_df <- split.default(df1, class)
  
  # Performs CS-Combat
  list_combat_df <- lapply(list_class_df,
                           correctComBat, metadata_df, cov = FALSE)
  # Reorders columns
  return(do.call(cbind, list_combat_df)[,colnames(df1)])
}

# IMPORT DATA -------------------------------------------------------------
rma_zhou <- read.table("data/GSE22544_GSE8977/rma_original.txt",
                       sep = "\t", header = T, row.names = 1)

zhou_metadata <- read.table("data/GSE22544_GSE8977/metadata.tsv",
                            sep = "\t", header = T, row.names = 1)
# Convert numeric to factor
zhou_metadata <- data.frame(apply(zhou_metadata, 2, factor))

# LOG & SCALE & FILTER ----------------------------------------------------
# rma data is already log2
# Trimmed-mean scaling
scaled_zhou <- normaliseMeanScaling(rma_zhou)
# Filter probesets
filtered_zhou <- filterProbesets(scaled_zhou, 0.3, zhou_metadata)

# Remove outlier: S20
subset_zhou <- filtered_zhou[,-20]
plotPCA3DBatchEffects(subset_zhou, zhou_metadata)
rgl.postscript("dump/pca_3d-zhou_refbcm.pdf", "pdf")

# SFL
sfl_zhou <- log2_transform(
  filterProbesets(normaliseMeanScaling(2^rma_zhou, 500), 0.3, zhou_metadata)
)
sfl_zhou1 <- sfl_zhou[,-20]
plot_sflzhou <- plotExplore(sfl_zhou1, zhou_metadata)
ggsave("dump/eval-sfl_zhou.pdf", plot_sflzhou,
       width = 8, height = 10)

# BCM ---------------------------------------------------------------------
refbcm_zhou <- correctRefBCM(subset_zhou, zhou_metadata, ref_batch = 1)
svdbcm_obj <- correctSVDBCM(subset_zhou, zhou_metadata, ref_batch = 1)

plotExplore(refbcm_zhou, zhou_metadata)
plotPCA3DBatchEffects(svdbcm_zhou, zhou_metadata)
rgl.postscript("dump/pca_3d-zhou_refbcm.pdf", "pdf")

class(unlist(svdbcm_obj$plot))
plot_svdbcm <- svdbcm_obj$plot[[1]]
class(plot_svdbcm)
ggsave("dump/pca_pairwise-corrected_zhou.pdf", plot_svdbcm)
plotPCA3DBatchEffects(svdbcm_obj$data, zhou_metadata)
rgl.postscript("dump/pca_3d-zhou_svdbcm.pdf", "pdf")


plot_svdbcm <- svdbcm_obj$plot[[1]]
plot_svdbcm

# ComBat ------------------------------------------------------------------
combat_zhou <- correctComBat(subset_zhou, zhou_metadata, FALSE)
plotPCA3DBatchEffects(combat_zhou, zhou_metadata)
rgl.postscript("dump/pca_3d-zhou_combat_cov.pdf", "pdf")

covcombat_zhou <- correctComBat(subset_zhou, zhou_metadata, TRUE)
# CS-ComBat ---------------------------------------------------------------
cscombat_zhou <- correctCSComBat(subset_zhou, zhou_metadata)
plotPCA3DBatchEffects(cscombat_zhou, zhou_metadata)
rgl.postscript("dump/pca_3d-zhou_cscombat.pdf", "pdf")
