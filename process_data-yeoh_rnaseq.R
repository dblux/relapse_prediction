# LABID_RPATH <- "data/GSE67684/processed/metadata/lab_id.tsv"
# labid_df <- read.table(LABID_RPATH, sep = "\t", header = T, row.names = 1)
# 
# RNA_RPATH <- "data/GSE67684/processed/rna_seq/count.matrix.txt"
# raw_rna <- read.table(RNA_RPATH, sep = "\t")
# 
# ### SUBSETTING PATIENTS AND CONVERTING PID ###
# labid_data <- unlist(lapply(strsplit(colnames(raw_rna), "_"),
#                             function(vec) vec[2]))
# # Only select patients present in microarray data
# logi_idx <- labid_data %in% rownames(labid_df)
# subset_rna <- raw_rna[,logi_idx]
# 
# # Replace labid of subsetted patients with pid
# colnames(subset_rna) <- labid_df[labid_data[logi_idx],]
# # Add in D0 and D8 labels
# colnames(subset_rna)[1:67] <- paste(colnames(subset_rna)[1:67],
#                                     "D0", sep = "_")
# colnames(subset_rna)[68:134] <- paste(colnames(subset_rna)[68:134],
#                                       "D8", sep = "_")

# ### DISCARDING GENES & REPLACING SYNONYMS ###
# DISCARD_RPATH <- "temp/discarded_genes.tsv"
# discard_synonyms <- readLines(DISCARD_RPATH)
# 
# SYNONYM_RPATH <- "temp/gene_synonym-lookup.tsv"
# synonym_lookup <- read.table(SYNONYM_RPATH, sep = "\t", header = T,
#                              row.names = 1, stringsAsFactors = F)
# 
# # Discard genes
# discard_rna <- subset_rna[!(rownames(raw_rna) %in% discard_synonyms),]
# 
# # Replace genes
# # Synonym table contains all synonyms with no corresponding length info
# logi_idx <- rownames(discard_rna) %in% rownames(synonym_lookup)
# # Synonyms to be replaced
# replace_synonyms <- rownames(discard_rna)[logi_idx]
# replacement_genes <- synonym_lookup[replace_synonyms,]
# names(replacement_genes) <- replace_synonyms
# 
# # Second discard: Synonyms with replacement gene names that are already present
# logi_idx2 <- replacement_genes %in% rownames(discard_rna)
# discard_synonyms2 <- names(replacement_genes[logi_idx2])
# discard_rna2 <- discard_rna[!(rownames(discard_rna) %in% discard_synonyms2),]
# 
# replacement_genes2 <- replacement_genes[!logi_idx2]
# 
# # Remove replacement genes that occur more than once
# # Multiple synonyms map to same gene
# ambiguous_genes <- names(table(replacement_genes2))[
#   table(replacement_genes2) > 1]
# # Third discard: Multiple synonyms that map to same gene
# logi_idx3 <- replacement_genes2 %in% ambiguous_genes
# ambiguous_synonyms <- names(replacement_genes2)[logi_idx3]
# discard_rna3 <- discard_rna2[!(rownames(discard_rna2) %in% ambiguous_synonyms),]
# # Filter out ambiguous genes
# replacement_genes3 <- replacement_genes2[!logi_idx3]
# 
# # Find idx of gene synonyms to be replaced
# int_idx <- match(names(replacement_genes3), rownames(discard_rna3))
# rownames(discard_rna3)[int_idx] <- replacement_genes3
# 
# write.table(discard_rna3, "data/GSE67684/processed/subpid_mapFeat-rna.tsv",
#             sep = "\t", quote = F)

##### RANK VS RANK #####
source("../functions.R")
library(rgl)
library(ggplot2)

# FUNCTIONS ---------------------------------------------------------------
plotScatter <- function(pid, X, Y,
                        xlab, ylab, flag = "log") {
  if (flag == "log") {
    plot(log2_transform(X[,pid]),
         log2_transform(Y[,pid]),
         xlab = xlab, ylab = ylab, main = pid)
  } else {
    plot(X[,pid],
         Y[,pid],
         xlab = xlab, ylab = ylab, main = pid)
  }
}

plotScatter1 <- function(pid1, pid2, X, Y,
                        xlab, ylab, flag = "log") {
  if (flag == "log") {
    plot(log2_transform(X[,pid1]),
         log2_transform(Y[,pid2]),
         xlab = xlab, ylab = ylab, main = paste(pid1, pid2))
  } else {
    plot(X[,pid1],
         Y[,pid2],
         xlab = xlab, ylab = ylab)
  }
}

calcTPM <- function(X, annot, length = "mean_length", c = 1e+6) {
  feature_length <- annot[rownames(X), length]
  # Stop if any of the rownames are not present in annot
  if (sum(is.na(feature_length)) > 0)
    stop("Rownames with no corresponding annotation exist..")
  # Divide by feature length in kilobases
  X_length <- X/(feature_length/1e+3)
  # Divide by sample sum
  tpm <- sweep(X_length, 2, colSums(X_length), "/") * c
  return(tpm)
}

calcRPKM <- function(X, annot, length = "mean_length", c = 1e+6) {
  feature_length <- annot[rownames(X), length]
  # Divide by sample sum
  rpm <- sweep(X, 2, colSums(X)/c, "/")
  # Divide by feature length in kilobases
  rpkm <- rpm/(feature_length/1e+3)
  return(rpkm)
}

# Subset of patients that overlap
# Susbet of features that are present in ncbi length annotation
RNA_RPATH <- "data/GSE67684/processed/subpid_mapFeat-rna.tsv"
rna_df <- read.table(RNA_RPATH, sep = "\t")

LENGTH_DF <- "../info/ref_genome/GRCh38/ncbi_GRCh38_p13-gene_lengths.tsv"
ncbi_length <- read.table(LENGTH_DF, sep = "\t",
                          header = T, row.names = 1, quote = "")

## Subset of original data
SUBSET_RPATH <- "data/GSE67684/processed/subset_yeoh.tsv"
raw_yeoh <- read.table(SUBSET_RPATH, sep = "\t")

## Metadata
# Preprocessed metadata
METADATA_RPATH <- "data/GSE67684/processed/metadata/full_metadata.tsv"
metadata_df <- read.table(METADATA_RPATH, sep = "\t")

# SCALE->REMOVE->FILTER->LOG
scaled_yeoh <- removeProbesets(normaliseMeanScaling(raw_yeoh))
SYMBOL_GPL570 <- "../info/microarray/HG-U133_Plus_2/annot_genesymbol-GPL570.tsv"
symbol_yeoh <- affy2id(scaled_yeoh, SYMBOL_GPL570)

# nrow(symbol_yeoh) # ~8000
# mean(colSums(symbol_yeoh == 0)) # ~4000
# nrow(rna_tpm) # ~20000
# mean(colSums(rna_tpm == 0)) # ~7000

# NORMALISE
# FILTER
# LOG (DOES NOT AFFECT RANK)

# Patients present in both data types
common_pid <- intersect(colnames(rna_df), colnames(symbol_yeoh))

# Subset of affy
subset_affy <- symbol_yeoh[,common_pid]
subset_rna <- rna_df[,common_pid]

# FILTER PROBESETS
filtered_affy <- filterProbesets(subset_affy, 0.7, metadata_df)
filtered_rna <- filterProbesets(subset_rna, 0.7, metadata_df)

# Genes present in both affy and RNA
genes <- intersect(rownames(filtered_affy), rownames(filtered_rna))
length(genes)

# Remove probesets that have 0 in affy
intersect_affy <- filtered_affy[genes,]
intersect_rna <- filtered_rna[genes,]

# Normalise affy by total intensity
total_intensity <- colSums(intersect_affy)
percent_affy <- sweep(intersect_affy, 2, total_intensity, "/") * 1e+6

# TPM
rna_tpm_max <- calcTPM(intersect_rna, ncbi_length, "max_length")

##### PID TO INVESTIGATE #####
pid <- common_pid[6]

# Re-rank after removing probesets with zero
keep_ps <- rownames(percent_affy)[-which(percent_affy[,pid] == 0)]

nozero_affy <- percent_affy[keep_ps,]
nozero_rna <- rna_tpm_max[keep_ps,]

plotScatter(pid, nozero_affy, nozero_rna,
            xlab = "log2(norm_intensity)", ylab = "log2(tpm_max)")
scatter_count <- recordPlot()
save_fig(scatter_count, "dump/same-nozero.pdf", 7, 7)

# Rank
nozero_affy <- percent_affy[keep_ps,]
nozero_rna <- rna_tpm_max[keep_ps,]

plotScatter(pid, nozero_affy, nozero_rna,
            xlab = "log2(norm_intensity)", ylab = "log2(tpm_max)")
scatter_count <- recordPlot()
save_fig(scatter_count, "dump/same-nozero.pdf", 7, 7)

rank_affy <- apply(nozero_affy, 2, rank, ties.method = "min")
rank_rna <- apply(nozero_rna, 2, rank, ties.method = "min")
plotScatter(pid, rank_affy, rank_rna,
            xlab = "rank(norm_intensity)", ylab = "rank(tpm_max)", flag = "nolog")
scatter_count <- recordPlot()
save_fig(scatter_count, "dump/same-rank_nozero.pdf", 7, 7)

# GFS to data
gfs_affy <- normaliseGFS(percent_affy, upper = 0.05, lower = 0.50, num_intervals = 4)
gfs_rna <- normaliseGFS(rna_tpm_max, upper = 0.05, lower = 0.50, num_intervals = 4)

# Plot intensity vs counts
plotScatter(pid, gfs_affy, gfs_rna, flag = "nolog",
            "gfs(norm_itensity, 0.05, 0.50)", "gfs(tpm_max0.05, 0.50)")
scatter_count <- recordPlot()
# save_fig(scatter_count, "dump/same-gfs2.pdf", 7, 7)

# 3D plots
plot_3d <- MASS::kde2d(data.matrix(gfs_affy[,pid]), data.matrix(gfs_rna[,pid]))
persp3d(plot_3d, col = "lightblue",
        xlab = "gfs(norm_itensity, 0.05, 0.50)",
        ylab = "gfs(tpm_max0.05, 0.50)", zlab = "",
        aspect = c(2,2,1.5))
persp3d(plot_3d, front = "lines", back = "lines",
        aspect = c(2,2,1.5), add = T)
rgl.postscript("dump/same-gfs3.pdf", "pdf")

# Get rid of the (0,0) points
row_idx <- which(gfs_affy[,pid] != 0 & gfs_rna[,pid] != 0)
plot_3d_2 <- MASS::kde2d(data.matrix(gfs_affy[row_idx, pid]),
                         data.matrix(gfs_rna[row_idx, pid]))
persp3d(plot_3d_2, col = "lightblue",
        xlab = "gfs(norm_itensity, 0.05, 0.50)",
        ylab = "gfs(tpm_max0.05, 0.50)", zlab = "",
        aspect = c(2,2,1))
persp3d(plot_3d_2, front = "lines", back = "lines",
        aspect = c(2,2,1), add = T)
rgl.postscript("dump/same-gfs4.pdf", "pdf")


X <- data.frame(gfs_affy[,pid], gfs_rna[,pid])
colnames(X) <- c("affy", "rna")

# Plot intensity vs counts
plotScatter(pid, intersect_affy, intersect_rna,
            "log2(itensity)", "log2(count)")
scatter_count <- recordPlot()
save_fig(scatter_count, "dump/same-count_intensity.pdf", 7, 7)

rna_tpm_mean <- calcTPM(intersect_rna, ncbi_length, "mean_length")
rna_tpm_max <- calcTPM(intersect_rna, ncbi_length, "max_length")

plotScatter(pid, affy_tpm_max, rna_tpm_max,
            "log2(tpm[max]_itensity)", "log2(tpm[max])")
scatter_tpm <- recordPlot()
save_fig(scatter_tpm, "dump/same-tpm_maxVStpm_max_intensity.pdf", 7, 7)

# intersect_affy
# percent_affy
# 
# rna_tpm_mean
# rna_tpm_max

plotScatter(pid, intersect_affy, rna_tpm_max,
            "log2(itensity)", "log2(tpm[max])")
scatter_tpm <- recordPlot()
save_fig(scatter_tpm, "dump/same-tpm_maxVSintensity.pdf", 7, 7)

plotScatter(pid, percent_affy, rna_tpm_max,
            "log2(norm_itensity)", "log2(tpm[max])")
scatter_tpm <- recordPlot()
save_fig(scatter_tpm, "dump/same-tpm_maxVSnorm_intensity.pdf", 7, 7)

pid1 <- common_pid[6]
pid2 <- common_pid[11]
metadata_df[pid2, "subtype"]

plotScatter1(pid1, pid2, percent_affy, rna_tpm_max,
            "log2(norm_itensity)", "log2(tpm[max])")
scatter_tpm <- recordPlot()
save_fig(scatter_tpm, "dump/diff-tpm_maxVSnorm_intensity.pdf", 7, 7)

plotScatter(pid, intersect_affy, rna_tpm_mean,
            "log2(itensity)", "log2(tpm[mean])")
scatter_tpm <- recordPlot()
save_fig(scatter_tpm, "dump/same-tpm_meanVSintensity.pdf", 7, 7)

plotScatter(pid, percent_affy, rna_tpm_mean,
            "log2(norm_itensity)", "log2(tpm[mean])")
scatter_tpm <- recordPlot()
save_fig(scatter_tpm, "dump/same-tpm_meanVSnorm_intensity.pdf", 7, 7)

# RANK
rank_affy <- apply(percent_affy, 2, rank, ties.method = "min")
rank_rna <- apply(rna_tpm_max, 2, rank, ties.method = "min")

plotScatter(pid, rank_affy, rank_rna,
            "rank(norm_itensity)", "rank(tpm_max)", "nolog")
scatter_rank <- recordPlot()
save_fig(scatter_rank, "dump/same-rank.pdf", 7, 7)

nrow(filtered_affy)
nrow(rna_tpm)
sum(rna_tpm[,pid] == 0)
sum(filtered_affy[,pid] == 0)

rna_rpkm <- calcRPKM(filtered_rna, ncbi_length)
plotScatter(pid, filtered_affy, rna_rpkm,
            "log2(itensity)", "log2(rpkm)")
scatter_rpkm <- recordPlot()
save_fig(scatter_rpkm, "dump/scatter-rpkm.pdf", 7, 7)

rna_rpkm[1:10,1:5]
rna_tpm[1:10,1:5]
