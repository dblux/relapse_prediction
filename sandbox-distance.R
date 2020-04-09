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

# library(xtable)

source("../functions.R")
source("bin/bcm.R")

# theme_set(theme_dark())
# theme_set(theme_gray())
theme_set(theme_cowplot())

# FUNCTIONS ---------------------------------------------------------------
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

##### SANDBOX: Same or different #####
# # Cosine normalise data
# normalised_yeoh <- normaliseCosine(data_yeoh)

### Subset all D0 patients
data_d0 <- data_yeoh[,endsWith(colnames(data_yeoh), "0")]
d0_metadata <- metadata_df[colnames(data_d0),]
## Order metadata according to batch
d0_metadata_ord <- d0_metadata[order(d0_metadata$batch_info),]
data_d0_ordered <- data_d0[,rownames(d0_metadata_ord)]

# Choosing the batches
pid_1 <- rownames(d0_metadata_ord)[d0_metadata_ord$batch_info == 8]
pid_2 <- rownames(d0_metadata_ord)[d0_metadata_ord$batch_info == 9]
all_pid <- c(pid_1, pid_2)

### FEATURE SELECTION (CHI2) ###
# Add class labels to each individual sample (row)
## Run the chi2 algorithm n times (one vs rest)
t_d0 <- t(data_d0)
subtype <- d0_metadata[rownames(t_d0), "subtype"] -> original_subtype
subtype <- original_subtype

chi2_features <- character()
for (i in 1:9) {
  subtype <- original_subtype
  levels(subtype)[-i] <- "Rest"
  weka_d0 <- data.frame(t_d0, subtype)
  selected_d0 <- chi2.algorithm(weka_d0, numeric(), threshold = 0)
  # Selected features
  chi2_features <- c(chi2_features, colnames(selected_d0[[1]]))
}

writeLines(substring(chi2_features, 2), con = "dump/chi2_features.txt")

X <- data_d0[substring(chi2_features, 2),]
y <- metadata_df[colnames(X),]
X <- X[,order(y$subtype)]

subtype_df <- metadata_df[colnames(data_d0), "subtype", drop = F]
pheatmap(X, col = brewer.pal(9, "Blues"),
         legend = T, border_color = NA, scale = "none",
         cluster_rows = T, cluster_cols = F,
         show_colnames = F, show_rownames = F,
         annotation_col = subtype_df)
heatmap <- recordPlot()
save_fig(heatmap, "dump/heatmap-chi2.pdf")]

# ## Save only D0 patients for WEKA
# # Add class labels to each individual sample (row)
# t_d0 <- t(data_d0)
# subtype <- d0_metadata[rownames(t_d0), "subtype"]
# weka_d0 <- cbind(t_d0, subtype)
# 
# write.table(weka_d0, "temp/weka_d0.csv", sep = ",", quote = F,
#             row.names = F, col.names = T)

## Batch-subtype configuration
table(d0_metadata_ord$batch_info, d0_metadata_ord$subtype)

## Entire batch 1 D0
pw_dist <- dist(t(data_d0_ordered), method = "euclidean")
dist_mat <- as.matrix(pw_dist)

i <- 1
x_i <- dist_mat[i,]
x_i1 <- x_i[x_i != 0]
hist(x_i1, breaks = 30)
x_i1

# Batch colours
batch_ordered <- d0_metadata[colnames(data_d0_ordered), "batch_info"]
batch_factor <- as.factor(batch_ordered)
levels(batch_factor) <- brewer.pal(10, "Set3")
batch_col <- as.character(batch_factor)

heatmap(dist_mat, Rowv = NA, Colv = "Rowv", scale = "none",
        labRow = NA, labCol = NA, col = brewer.pal(9, "Blues"),
        ColSideColors = batch_col, RowSideColors = batch_col)

# heatmap_plot <- recordPlot()
# save_fig(heatmap_plot, "dump/heatmap-dist_all_batches.pdf",
#          width = 8, height = 8)

# Choosing the batches
pid_1 <- rownames(d0_metadata_ord)[d0_metadata_ord$batch_info == 8]
pid_2 <- rownames(d0_metadata_ord)[d0_metadata_ord$batch_info == 9]
all_pid <- c(pid_1, pid_2)

# Batch 1 vs Batch 2
pair_dist <- dist_mat[pid_1, pid_2]
pheatmap(pair_dist, col = brewer.pal(9, "Blues"),
         legend = F, border_color = NA, scale = "none",
         cluster_rows = F, cluster_cols = F)

heatmap_b8b9 <- recordPlot()
save_fig(heatmap_b8b9, "dump/heatmap-dist_b8b9.pdf",
         width = 12, height = 8)

# Subtype breakdown for pair of batches
pair_metadata <- d0_metadata_ord[c(pid_1, pid_2),]
table(pair_metadata$subtype, pair_metadata$batch_info)

getNN <- function(pairwise_mat, flag) {
  # Distance: Choose smallest distance
  # Correlation: Choose largest correlation
  if(flag == "dist") FUNC <- which.min
  else FUNC <- which.max
  
  # 1-NN for batch 1 samples
  for(p1 in rownames(pairwise_mat)) {
    p1_subtype <- as.character(d0_metadata[p1,"subtype"])
    
    p2 <- names(FUNC(pairwise_mat[p1,]))
    p2_subtype <- as.character(d0_metadata[p2,"subtype"])
    
    msg <- sprintf("%s (%s): %s (%s)\n", p1, p1_subtype, p2, p2_subtype)
    cat(msg)
    
    # dist_incr <- sort(pairwise_mat[p1,])
    # names(dist_incr) <- metadata_df[names(dist_incr), "subtype"]
    # print(dist_incr)
    
    # hist(pairwise_mat[p1,], breaks = 10, main = p1)
    # abline(v = min(pairwise_mat[p1,]), col = "red")
  }
  
  cat("\n")
  
  # 1-NN for batch 1 samples
  for(p2 in colnames(pairwise_mat)) {
    p2_subtype <- as.character(d0_metadata[p2,"subtype"])
    
    p1 <- names(FUNC(pairwise_mat[,p2]))
    p1_subtype <- as.character(d0_metadata[p1,"subtype"])
    
    msg <- sprintf("%s (%s): %s (%s)\n", p2, p2_subtype, p1, p1_subtype)
    cat(msg)
    
    # dist_incr <- sort(pairwise_mat[,p2])
    # names(dist_incr) <- metadata_df[names(dist_incr), "subtype"]
    # print(dist_incr)
    
    # hist(pairwise_mat[,p2], breaks = 10, main = p2)
    # abline(v = min(pairwise_mat[,p2]), col = "red")
  }  
}

getNN(pair_dist, "dist")
getNN(pair_cor, "cor")
getNN(pair_cor1, "cor")

pair_cor1[1,]



plotPCA3DYeoh1(data_yeoh[,all_pid], metadata_df)
# rgl.postscript("dump/pca_3d-b7b8_d0.pdf", "pdf")

### Chi-square selected genes from yeoh_2002 ###
# Neither probesets nor gene symbols correspond with current technology
CHI_GENES <- "data/yeoh_2002/README/chi_square_probesets/combined.tsv"
chi_genes <- read.delim(CHI_GENES)
chi_symbols <- levels(chi_genes$gene_symbol)[-1]

SYMBOL_GPL570 <- "../info/microarray/HG-U133_Plus_2/annot_genesymbol-GPL570.tsv"
affy_yeoh <- affy2id(data_yeoh, SYMBOL_GPL570)
intersect(chi_symbols, rownames(affy_yeoh))

plotPCA3DYeoh1(data_yeoh[chi_genes,all_pid], metadata_df)

### SELECTING GENE SIGNATURES OF SUBTYPES ###
# Using only D0 samples
# T-test based on the ranks
getTTestFeatures <- function(data, subtype) {
  col_idx <- which(metadata_df[colnames(data), "subtype"] == subtype)
  ordered_data <- cbind(data[,col_idx], data[,-col_idx])
  pvalue <- calc_ttest(ordered_data, length(col_idx), is_paired = F)
  names(head(sort(pvalue), 20))
}

filtered_d0 <- filterProbesets(data_d0, 0.95, metadata_df)

list_probesets <- sapply(levels(metadata_df$subtype),
                    getTTestFeatures, data = filtered_d0)
# Hypodiploid and Normal are not included
probesets <- unique(unlist(list_probesets, use.names = F))

# Ordering data_d0 according to subtype
subtype_info <- metadata_df[colnames(data_d0), "subtype"]
col_idx <- order(subtype_info)

subtype_df <- metadata_df[colnames(data_d0), "subtype", drop = F]
ordered_d0 <- data_d0[probesets, col_idx]
ps_median <- apply(ordered_d0, 1, median)
ps_mad <- apply(ordered_d0, 1, mad)
znorm_d0 <- (ordered_d0 - ps_median)/ ps_mad
znorm_d0[1:5,1:5]
cust_pal <- c(rep("#67001F", 45),
              brewer.pal(11, "RdBu"),
              rep("#053061", 5))

# FOR PROGNOSTIC MARKERS. USUALLY UNIQUELY HIGH EXPRESSION..
# ..FOR THAT SUBTYPE, REST ARE ZERO VALUES

pheatmap(znorm_d0, color = cust_pal,
         legend = T, border_color = NA, scale = "none",
         cluster_rows = F, cluster_cols = F,
         annotation_col = subtype_df)

##### CORRELATION AS A DISTANCE METRIC #####

# Correlation ignoring all zeros
# Correlated data
a <- rnorm(100)
b <- rnorm(100) + a
plot(a,a*3)

# All zeroes
# Both zeroes does not affect correlation
c <- rnorm(100)
d <- c(rep(0,100), c)
e <- c(rep(0,100), rnorm(100))
# cbind(d,e)
plot(d,e)
cor(d,e)

# Zero vs RV
# More zero-signal values decrease correlation
f <- rnorm(100) + 5
g <- rep(f, 21)
h <- c(rep(0, 2000), f)
plot(g, h, xlim = c(0,10))
cor(g,h)

# Pearson correlation: Zero-inflated data
cor(a, b, method = "pearson")
cor(a, b, method = "spearman")
cor(a, b, method = "kendall")

## TEST
# Batch 8 and 9 metadata
d0_metadata_ord[all_pid,]
p1 <- data_d0_ordered[,"P103_D0"]
p2 <- data_d0_ordered[,"P151_D0"] # CORRECT
p3 <- data_d0_ordered[,"P137_D0"] # WRONG

calcCor(p1, p2)
calcCor(p1, p3)

# corr_p1 <- sapply(pid_1, function(pid) cor(p1, data_d0_ordered[,pid]))
# which.max(corr_p1)

fig1 <- ggplot(data.frame(p1, p2), aes(x = p1, y = p2)) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon")

fig2 <- ggplot(data.frame(p1, p3), aes(x = p1, y = p3)) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon")

density_plot <- plot_grid(fig1, fig2, ncol = 2)
density_plot
 
par(mfrow=c(1,2))
plot(p1, p2) # Correct
plot(p1, p3) # Wrong
cor_plot <- recordPlot()
save_fig(cor_plot, "dump/scatter-cor.pdf", 12, 6)

# Investigating patterns in zero values
# No pattern in dispersion of values that accompany zeros
logi1 <- p1 == 0
logi2 <- p2 == 0
logi3 <- p3 == 0
xor_logi1 <- xor(logi1, logi2)
xor_logi2 <- xor(logi1, logi3)

sum(xor_logi1)
sum(logi1 & logi2)

sum(xor_logi2)
sum(logi1 & logi3)
# P3 has disproportionately little zeros? Try other patientss

mat1 <- cbind(p1, p2)[xor_logi1,]
mat2 <- cbind(p1, p3)[xor_logi2,]

hist(mat1[mat1[,1] == 0, 2], breaks = 15,
     xlab = NULL, main = "Zero-signal values [P1 (0) vs P2]")
hist(mat2[mat2[,1] == 0, 2], breaks = 16,
     xlab = NULL, main = "Zero-signal values [P1 (0) vs P3]")
hist_plot <- recordPlot()
save_fig(hist_plot, "dump/hist-xzero_ysignal.pdf", 12, 6)

cor(p1, p2)
cor(p2, p1)

#' @param X df with samples as columns and features as rows
calcPairwise <- function(X, pid_grp1, pid_grp2, FUNC, ...) {
  sapply(pid_grp2, function(p2) {
    sapply(pid_grp1, function(p1) FUNC(X[,p2], X[,p1], ...))
  })
}

calcCor <- function(v1, v2) {
  logi_idx <- v1 != 0 & v2 != 0
  cor(v1[logi_idx], v2[logi_idx])
}

pair_cor <- calcPairwise(data_d0_ordered, pid_1, pid_2, cor)
pair_cor1 <- calcPairwise(data_d0_ordered, pid_1, pid_2, calcCor)

pheatmap(pair_cor, col = brewer.pal(9, "Blues"),
         legend = F, border_color = NA, scale = "none",
         cluster_rows = F, cluster_cols = F)

heatmap_b8b9_cor <- recordPlot()
save_fig(heatmap_b8b9_cor, "dump/heatmap-cor_b8b9.pdf",
         width = 12, height = 8)
