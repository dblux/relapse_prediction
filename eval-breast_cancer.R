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
# Corrects for each class
#' @param cov logical: include class as covariate?
correctComBat <- function(df1, metadata_df, cov = FALSE) {
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
# SUBSET DATA: Remove outlier
nolog_zhou <- 2^rma_zhou[,-20]
sfl_zhou <- log2_transform(
  filterProbesets(normaliseMeanScaling(nolog_zhou), 0.3,zhou_metadata)
)

plot_sflzhou <- plotExplore(sfl_zhou, zhou_metadata)
plot_sflzhou
ggsave("dump/eval-sfl_zhou.pdf", plot_sflzhou,
       width = 8, height = 10)

# BCM ---------------------------------------------------------------------
refbcm_zhou <- correctRefBCM(subset_zhou, zhou_metadata, ref_batch = 1)
svdbcm_obj <- correctSVDBCM(sfl_zhou, zhou_metadata, ref_batch = 1)

pca_plot <- svdbcm_obj$plot[[1]]
pca_plot <- svdbcm_obj$corrected_plot[[1]]
ggsave("dump/pca-corrected_zhou.pdf", pca_plot)


dim(sfl_zhou)
plotPCA2D(sfl_zhou[,1:19], zhou_metadata)
plotPCA2D(sfl_zhou[,20:41], zhou_metadata)

plotExplore(refbcm_zhou, zhou_metadata)
plotPCA3DBatchEffects(svdbcm_zhou, zhou_metadata)
rgl.postscript("dump/pca_3d-zhou_refbcm.pdf", "pdf")

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
