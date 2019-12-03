#!/usr/bin/env Rscript
library(dplyr)

# Cosine normalisation / L2-norm normalisation
normaliseCosine <- function(df1) {
  l2norm_vec <- apply(df1, 2, calcL2Norm)
  return(mapply(`/`, df1, l2norm_vec))
}

# Min-max scaling function
# Returns: Scaled vector with range [0,1]
normaliseMinmax <- function(vec) {(vec-min(vec))/(max(vec)-min(vec))}

# Trimmed mean scaling
# Non-log values
# Trimmed mean_scaling does not remove all tied values
normaliseMeanScaling <- function(df, target_mean = 500, trim = 0.02) {
  trimmed_mean <- apply(df, 2, mean, trim = trim)
  scaling_factor <- target_mean/trimmed_mean
  print(scaling_factor)
  scaled_df <- as.data.frame(mapply(function(a,b) a*b, df, scaling_factor))
  rownames(scaled_df) <- rownames(df)
  return(scaled_df)
}

# Quantile normalisation: 0 values are assigned 0 automatically
# Takes in df where columns are samples and rows are genes
normaliseQuantile <- function(df) {
  zero_filter <- df == 0
  sort_arr <- apply(df, 2, sort)
  # Creates reference distribution
  ref_distr <- apply(sort_arr, 1, mean)
  rank_arr <- apply(df, 2, rank, ties.method = "min")
  qnorm_arr <- apply(rank_arr, c(1,2), function(x) ref_distr[x])
  rownames(qnorm_arr) <- rownames(df)
  qnorm_df <- as.data.frame(qnorm_arr)
  qnorm_df[zero_filter] <- 0
  return(qnorm_df)
}

# Problem: When there are too many zeros and fewer values are assigned to be 0
# .. No zeroes will be assigned
# Gene Fuzzy Scoring function transforms gene expression values
# Wilson Goh's paper
# Dense rank is used
normaliseGFS <- function(A, upper = 0.05, lower = 0.15, num_intervals = 0) {
  # Bins score with range [0,1] into intervals
  # E.g. 4 Intervals: Binned into 0.2, 0.4, 0.6, 0.8
  bin <- function(score, num_intervals) {
    for (i in 1:num_intervals) {
      if (score <= i/num_intervals) {
        return (i/(num_intervals+1))
      }
    }
  }
  cat(sprintf("Top %.2f of expressed genes are assigned GFS scores of 1\n", upper))
  cat(sprintf("Genes below the top %.2f of expressed genes are assigned GFS scores of 0\n", lower))
  # Rank function ranks largest value as 1 [-A is used]
  # Handle NaN?
  ranked_A <- apply(-A, 2, dense_rank)
  rownames(ranked_A) <- rownames(A)
  # Returns [1,] = upper, [2,] = lower
  qtile <- apply(ranked_A, 2, quantile, probs=c(upper, lower), names=F)
  
  if (num_intervals <= 0) {
    for (c in 1:ncol(ranked_A)) {
      # Calculate qtile range
      q_range <- qtile[2,c] - qtile[1,c]
      for (r in 1:nrow(ranked_A)) {
        if (ranked_A[r,c] <= qtile[1,c]) {
          # Assign 1s
          ranked_A[r,c] <- 1
        } else if (ranked_A[r,c] > qtile[2,c]){
          # Assign 0s
          ranked_A[r,c] <- 0
        } else {
          # Assign score
          score <- (qtile[2,c] - ranked_A[r,c]) / q_range
          ranked_A[r,c] <- score
        }
      }
    }
  } else {
    # Discrete intervals
    for (c in 1:ncol(ranked_A)) {
      # Calculate qtile range
      q_range <- qtile[2,c] - qtile[1,c]
      for (r in 1:nrow(ranked_A)) {
        if (ranked_A[r,c] <= qtile[1,c]) {
          # Assign 1s
          ranked_A[r,c] <- 1
        } else if (ranked_A[r,c] > qtile[2,c]){
          # Assign 0s
          ranked_A[r,c] <- 0
        } else {
          # Assign score
          score <- (qtile[2,c] - ranked_A[r,c]) / q_range
          # Round off score
          ranked_A[r,c] <- bin(score, num_intervals)
        }
      }
    }
  }
  return (as.data.frame(ranked_A))
}

normaliseCDF <- function(df) {
  for (c in 1:ncol(df)) {
    notzero <- df[,c] != 0
    df[,c][notzero] <- rank(df[,c][notzero])
    df[,c] <- df[,c]/sum(notzero)
  }
  return(df)
}

# Filter probes with too many zeros
#' @param df dataframe
#' @param percent_threshold percentage threshold of non-zero values
#' @param metadata_df df containing class labels of samples
#' @param logical_func a function that is either "all" or "any". (Either all or
#' just one class have to pass the threshold)
#' @return vector containing rownames of rows that meet threshold of non-zero
#' values
filterProbesets <- function(df1, percent_threshold, metadata_df = NULL,
                            logical_func = any) {
  if (is.null(metadata_df)) {
    logical_df <- df1 != 0
    selected_logvec <- rowSums(logical_df) > percent_threshold * ncol(df1)
    print(paste("No. of probesets removed =",
                nrow(df1) - sum(selected_logvec)))
    return(df1[selected_logvec,])
  } else {
    class_factor <- metadata_df[colnames(df1),"class"]
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

#' Removes ambiguous and AFFY probesets from dataframe
#' Rowname of affymetrix probesets
removeProbesets <- function(df) {
  logical_vec <- grepl("[0-9]_at", rownames(df)) & !startsWith(rownames(df),
                                                               "AFFX")
  print(paste0("No. of ambiguous and AFFY probesets removed: ",
               nrow(df) - sum(logical_vec)))
  return(df[logical_vec, , drop=F])
}

# Function that calculates the matrix of mean differences from the patient and control matrix
# Arguments: a <- larger matrix, b <- smaller matrix
# CLT: As n -> infinity, approximates more accurately a normal distribution.
# If underlying distribution is normal, then no need for n to be large
# Returns: mean_diff matrix with dimension of smaller matrix
old_calc_diff <- function(a, b) {
  a_colsum <- apply(a, 1, sum)
  a_size <- ncol(a)
  mean_diff <- array(numeric(), dim(b))
  rownames(mean_diff) <- rownames(b)
  
  for (r in 1:nrow(b)) {
    for (c in 1:ncol(b)) {
      mean_diff[r,c] <- (a_colsum[r]/a_size) - b[r,c]
    }
  }
  return (mean_diff)
}

# Function that calculates the matrix of mean differences from the patient and control matrix
# Arguments: a <- larger matrix, b <- smaller matrix
# CLT: As n -> infinity, approximates more accurately a normal distribution.
# If underlying distribution is normal, then no need for n to be large
# Returns: mean_diff matrix with dimension of smaller matrix
# pct= 0.75 smallest=4. pct=0.8 smallest=5
calc_diff <- function(a, b, sample_pct) {
  a_size <- ncol(a)
  sample_size <- floor(sample_pct * a_size)
  if (sample_size == a_size) {print("Sample size too small..")}
  
  mean_diff <- array(numeric(), dim(b))
  rownames(mean_diff) <- rownames(b)
  
  for (r in 1:nrow(b)) {
    for (c in 1:ncol(b)) {
      subsample <- unlist(sample(a[r,], sample_size, replace = F))
      mean_diff[r,c] <- mean(subsample) - b[r,c]
    }
  }
  return (mean_diff)
}

all_diff <- function(a, b) {
  all_diff <- array(numeric(), c(nrow(a), ncol(a)*ncol(b)))
  rownames(all_diff) <- rownames(a)
  for (r in 1:nrow(a)) {
    print(r)
    for (col_a in 1:ncol(a)) {
      for (col_b in 1:ncol(b)) {
        c <- col_b + (col_a - 1) * ncol(b)
        all_diff[r,c] <- b[r,col_b] - a[r,col_a]
      }
    }
  }
  return (all_diff)
}

# Calculate one-sample t-test statistic
# Using sample sd as estimator of population sd
# Use t-distribution to calculate p value
ttest_onesample <- function(vector, mu) {
  n <- length(vector)
  t_stat <- (mean(vector) - mu) / (sd(vector)/sqrt(n))
  if (is.infinite(t_stat)|is.na(t_stat)) {
    return (NaN)
  }
  if (t_stat < 0) {
    p_value <- pt(t_stat, n-1) * 2
  } else if (t_stat > 0) {
    p_value <- (1-pt(t_stat, n-1)) * 2
  } else {
    # t_stat == 0
    p_value <- 1
  }
  return (p_value)
}

# Naive row-wise two-sample t-test for every probe
# Does a t-test between every row of matrices a and b
# Arguments: Dataframe with both clases, size of class A, ...
calc_ttest <- function(df, size_a, flag = "pvalue", is_paired = F) {
  row_pvalue <- function(row_vec) {
    return(t.test(row_vec[1:size_a],
                  row_vec[-(1:size_a)],
                  paired = is_paired)$p.value)
  }
  row_tstat <- function(row_vec) {
    return(t.test(row_vec[1:size_a],
                  row_vec[-(1:size_a)],
                  paired = is_paired)$statistic)
  }
  
  if (flag == "pvalue") {
    ttest_vec <- apply(df, 1, function(row) tryCatch(row_pvalue(row), error = function(e) NA))
  } else if (flag == "tstat") {
    ttest_vec <- apply(df, 1, function(row) tryCatch(row_tstat(row), error = function(e) NA))
  } else {
    stop("Flag not in options pvalue or tstat..")
  }
  return(ttest_vec)
}

# Arguments: 2 dataframes that are not log-transformed
# Log-fold change (class1/class2)
calc_logfc <- function(df1, df2, func = mean) {
  # Minimum value of both df besides 0 chosen as prior value
  prior_value <- min(c(df1[df1 != 0], df2[df2 != 0]))
  print(paste("Prior value:", prior_value))
  vec1 <- apply(df1, 1, func)
  vec2 <- apply(df2, 1, func)
  # log2(0) = -Inf; log2(Inf) = -Inf; log2(0/0) = NaN
  # Reassigns 0s with prior_expr
  vec1[vec1 == 0] <- prior_value
  vec2[vec2 == 0] <- prior_value
  fc <- vec1/vec2
  return(log2(fc))
}

#' Calculates QPSP profiles of samples
#' Rownames of df1 has to be the same annotation type as list of network vectors
#' @param df1 dataframe of discrete GFS-transformed expression data with
#' features in rows and samples in columns
#' @param list_complex list of vectors that contain components of networks
#' @return dataframe of QPSP profiles
calcQPSP <- function(df1, list_complex) {
  # All components in complex have to be measured in df1
  df1_proteins <- rownames(df1)
  # Logvec of whether entire complex is present
  complex_logvec <- sapply(list_complex,
                           function(vec) all(vec %in% df1_proteins))
  print(sum(complex_logvec))
  # Error if any NA in complex_logvec
  if (anyNA(complex_logvec)) stop("NA present in logvec")
  sublist_complex <- list_complex[complex_logvec]
  
  # Driver function that calculates QPSP profile for single sample
  # @param col_vec Column vector representing single sample
  calcSingleQPSP <- function(col_vec) {
    sapply(sublist_complex, function(vec) mean(col_vec[vec]))
  }
  return(data.frame(apply(df1, 2, calcSingleQPSP)))
}

# Argument: Recorded plot
# Save figure as file format indicated
save_fig <- function(recorded_plot, fpath, width = 8, height = 5) {
  if (endsWith(fpath, ".eps")) {
    setEPS()
    postscript(fpath, width = width, height = height)
    replayPlot(recorded_plot)
    dev.off()
  } else if (endsWith(fpath, ".pdf")) {
    pdf(fpath, width = width, height = height)
    replayPlot(recorded_plot)
    dev.off()
  } else if (endsWith(fpath, ".png")) {
    png(fpath, width = width, height = height)
    replayPlot(recorded_plot)
    dev.off()
  } else if (endsWith(fpath, ".jpg")) {
    jpeg(fpath, width = width, height = height)
    replayPlot(recorded_plot)
    dev.off()
  } else {
    stop("File extension not supported...")
  }
}

# Arguments: Dataframe, probeset annotation filepath
# Maps affy probesets to ID
# Removes ambiguous probesets and probesets with no ID
# Selects maximum if two probesets match to same gene
affy2id <- function(df, annot_fpath) {
  probeset_annot <- read.table(annot_fpath,
                               sep="\t", header=T, row.names=1,
                               stringsAsFactors=F, strip.white = T)
  # Filters out ambiguous and AFFY probesets from annot
  fltr_annot <- probeset_annot[grepl("[0-9]_at", rownames(probeset_annot))
                               & !startsWith(rownames(probeset_annot), "A"), , drop=F]
  # Returns entrez ID for all probe sets
  id <- unname(sapply(rownames(df), function(x) probeset_annot[x,]))
  
  # Indices of ambiguous probe sets and probe sets with no corresponding entrez ID to be deleted
  list_del <- which(grepl("///", id) | id == "")
  print(paste0("No. of probesets mapping to multiple IDs removed: ", sum(grepl("///", id))))
  print(paste0("No. of probesets with no ID removed: ", sum(id == "")))
  # Identifies genes that have multiple probesets mapping to it
  freq_gene <- table(id)
  dup_genes <- names(freq_gene[freq_gene > 1])
  for (i in dup_genes) {
    # Rows of dataframe with the same entrez ID
    same_rows <- df[id == i,]
    # Assign indices as rownames
    rownames(same_rows) <- which(id == i)
    # Rows that do not have the maximum sum are deleted
    row_del <- as.integer(rownames(same_rows[-which.max(apply(same_rows,1,sum)),]))
    # Concat with existing list of indices to be deleted
    list_del <- c(list_del, row_del)
  }
  # Rows are deleted
  df_genes <- df[-list_del,]
  fltr_id <- id[-list_del]
  # Assigning entrez ID to df
  rownames(df_genes) <- fltr_id
  # # CONCEPT CHECK: Deleted rows
  # df_genes_del <- df[list_del,]
  # entrez_del <- entrez[list_del]
  print(paste0("Total no. of probesets removed (incl. probesets mapping to same gene): ",
               length(list_del)))
  return(df_genes)
}

# Evaluates DE analysis method
# Arguments: 2 logical vectors - predicted and truth labels
# Returns: Accuracy, Sensitivity, Precision, Specificity, NPR
evaluation_report <- function(predict_vec, label_vec) {
  accuracy <- sum(predict_vec == label_vec)/length(predict_vec)
  TP <- sum(predict_vec & label_vec)
  TN <- sum(!predict_vec & !label_vec)
  sensitivity_1 <- TP/sum(label_vec)
  precision_1 <- TP/sum(predict_vec)
  specificity_1 <- TN/sum(!label_vec)
  negative_predictive_rate_1 <- TN/sum(!predict_vec)
  Sensitivity <- c(sensitivity_1, specificity_1)
  Precision <- c(precision_1, negative_predictive_rate_1)
  df <- cbind(Sensitivity, Precision)
  rownames(df) <- c("DE", "Not DE")
  print(df)
  print(paste("Accuracy:", accuracy))
  metrics <- c(accuracy, sensitivity_1, precision_1,
               specificity_1, negative_predictive_rate_1)
  names(metrics) <- c("Accuracy", "Sensitivity", "Precision",
                      "Specificity", "NPR")
  return(metrics)
}

# Log2 transforms data and handles -Inf values
log2_transform <- function(df) {
  log2_df <- log2(df)
  logical_df <- is.infinite(data.matrix(log2_df))
  log2_df[logical_df] <- 0
  return(log2_df)
}

# Substring without n last characters
substring_head <- function(string, n) substring(string, 1, nchar(string)-n)

# Substring tail starting from length - n characters
substring_tail <- function(string, n) substring(string, nchar(string)-n+1)

# Plot venn diagram
jacc_coeff <- function(vec1, vec2) {
  # Generate overlap list
  overlap_list <- calculate.overlap(list(vec1,vec2))
  # Calculate venndiagram areas
  venn_area <- sapply(overlap_list, length)
  grid.newpage()
  venn_plot <- draw.pairwise.venn(venn_area[1], venn_area[2], venn_area[3],
                                  category = c("D1", "D2"),
                                  cex = 3, fontfamily = "sans",
                                  cat.cex = 3, cat.fontfamily = "sans",
                                  margin = 0.1)
  union <- (venn_area[1] + venn_area[2] - venn_area[3])
  print(unname(venn_area[3]/union))
  return(venn_plot)
}

# Side effects: Writes KEGG pathway dataframes (Entrez)
kegg_df <- function(kegg_fpath) {
  pathway_id <- substring(kegg_fpath, 29, 36)
  wpath <- sprintf("../info/KEGG/hsa_df/%s.tsv", pathway_id)
  # Parses file into dataframe
  df_hsa <- parseKGML2DataFrame(kegg_fpath)
  # Convert KEGG ID to Entrez ID
  df_hsa$from <- substring(df_hsa$from, 5)
  df_hsa$to <- substring(df_hsa$to, 5)
  write.table(df_hsa, wpath,
              quote = F, sep = "\t", row.names = F)
}

# Assumes that dataframe has been log-transformed
plotExplore <- function(df1, metadata_df) {
  # Obtaining batch and class annotations
  batch_factor <- as.factor(metadata_df[colnames(df1),"batch"])
  class_factor <- metadata_df[colnames(df1),"class"]
  
  # Melt dataframe
  melt_df <- melt(df1, measure.vars = colnames(df1), variable.name = "ID")
  melt_df$batch <- as.factor(metadata_df[melt_df$ID,"batch"])
  melt_df$class <- metadata_df[melt_df$ID,"class"]
  
  # Plot means
  mean_tibble <- melt_df %>% group_by(ID) %>% summarise(mean = mean(value))
  mean_scatter <- ggplot(mean_tibble, aes(x = ID, y = mean,
                                          col = batch_factor,
                                          pch = class_factor)) +
    geom_point(size = 3, show.legend = F) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  # Plot boxplot
  boxplot <- ggplot(melt_df, aes(x = ID, y = value, col = ID)) + 
    geom_boxplot(show.legend = F) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  # Plot density curve
  print(head(melt_df))
  print(tail(melt_df))
  pdf <- ggplot(melt_df, aes(x = value, group = ID, col = batch)) +
    geom_density(show.legend = F, alpha = 0.3) +
    facet_wrap(~class) +
    scale_color_viridis_d()
  
  # Plot PCA
  # Filters out rows with all zero values
  nonzero_logvec <- rowSums(df1) != 0 & apply(df1, 1, var) != 0
  pca_obj <- prcomp(t(df1[nonzero_logvec,]),
                    center = T, scale. = F)
  pca_df <- data.frame(pca_obj$x[,1:4])
  eigenvalues <- (pca_obj$sdev)^2
  var_pc <- eigenvalues[1:5]/sum(eigenvalues)
  pc_labels <- sprintf("PC%d (%.2f%%)", 1:5, var_pc*100)
  
  pc1_pc2 <- ggplot(pca_df, aes(x = PC1, y = PC2, col = batch_factor,
                                pch = class_factor)) +
    geom_point(size = 3, show.legend = F) +
    labs(x = pc_labels[1], y = pc_labels[2])
  pc1_pc3 <- ggplot(pca_df, aes(x = PC1, y = PC3, col = batch_factor,
                                pch = class_factor)) +
    geom_point(size = 3, show.legend = F) +
    labs(x = pc_labels[1], y = pc_labels[3])
  
  # Plot all graphs
  pca <- plot_grid(pc1_pc2, pc1_pc3)
  multiplot <- plot_grid(mean_scatter, pdf, pca, nrow = 3)
  return(multiplot)
}

# Returns: ggplot2 of PCA plot
plot_pca <- function(df, batch_info) {
  # Principal component analysis
  # Removes columns with all zeroes
  col_logical <- apply(t(df), 2, var) != 0
  pca_df <- t(df)[, col_logical]
  pca_obj <- prcomp(pca_df, center = T, scale. = T)
  top_pc <- as.data.frame(pca_obj$x[,1:2])
  pc1_pc2 <- ggplot(top_pc, aes(x = PC1, y = PC2, col = factor(batch_info))) +
    geom_point(size = 2, show.legend = F) +
    geom_label(label = rownames(pca_obj$x),
               nudge_x = 1, nudge_y = 2, size = 4,
               show.legend = F)
  return(pc1_pc2)
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

# Plot PCA 3D: Batch effects
# Plot batches in different colours and classes in different shapes
plotPCA3DBatchEffects <- function(df1, metadata_df) {
  # Batch and class annotations
  batch_factor <- metadata_df[colnames(df1), "batch"]
  batch_palette <- generateGgplotColours(length(unique(batch_factor)))
  batch_colour <- batch_palette[batch_factor]
  
  class_factor <- metadata_df[colnames(df1), "class"]
  all_pch <- 21:25
  # Error if there are more classes than pch symbols (> 5)
  stopifnot(length(unique(class_factor)) <= 5)
  class_pch <- all_pch[class_factor]
  plotPCA3D(df1, batch_colour, class_pch)
}

# Plots ROC and calculates AUC in a primitive fashion (i.e. ROC is step function)
# Does not resolve ties in the score
# Assumption: Lower score will be labelled preferentially as 1, ROC is step function
# Assumption that score vec and label vecs are corresponding
plotROC <- function(score_list, label_vec,
                    name_vec = NULL, plot_title = NULL,
                    is_bigger_better_vec = rep(F, length(score_list))) {
  # Function to plot a single ROC curve and calculate AUC
  plotSingleROC <- function(score_vec, is_bigger_better, color) {
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
    # Plot single ROC curve
    lines(roc_df$FPR, roc_df$TPR,
          col = color, lwd = 3)
    return(AUC)
  }
  # Initialise plot
  colour_palette <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A",
                      "#66A61E", "#E6AB02", "#A6761D", "#666666")
  color_index <- colour_palette[1:length(score_list)]
  plot(NULL, main = plot_title,
       xlab = "1 - Specificity", ylab = "Sensitivity",
       xlim = c(0,1), ylim = c(0,1),
       xaxs = "i", yaxs = "i",
       cex.main = 1.8, cex.lab = 1.7, cex.axis = 1.6)
  abline(0, 1, lty = 5, lwd = 2)
  auc_vec <- mapply(plotSingleROC, score_list, is_bigger_better_vec,
                    color_index, SIMPLIFY = T)
  # If name_vec is not NULL display legend
  if (!is.null(name_vec)) {
    format_text <- function(name, auc) sprintf("%s (%.3f)", name, auc)
    legend_text <- mapply(format_text, name_vec, auc_vec)
    legend("bottomright", inset = 0.03, lty = 1, lwd = 2,
           cex = 1.2, bty = "o", text.width = 0.35,
           legend = legend_text, col = color_index)
  }
  return(auc_vec)
}

eval_batch_effects <- function(df, batch_info, class_info) {
  # Tranpose df to n x p
  pca_obj <- prcomp(t(df), center = T, scale. = F)
  # Eigenvalues
  pca_df <- data.frame(pca_obj$x)
  eig_value <- (pca_obj$sdev)^2
  # Percentage variance of each PC
  var_pct <- eig_value/sum(eig_value)
  
  # Libraries: cluster, gpca
  # Argument: PC coordinates for single PC
  calc_pc_var <- function(vec) {
    # Argument: ANOVA attributes; Calculates percentage of variance
    # SS_between/(SS_between + SS_within)
    calc_var_percentage <- function(vec) unname(vec[3]/(vec[3] + vec[4]))
    pc_metadata <- data.frame(pc = vec,
                              batch = as.factor(batch_info),
                              class = as.factor(class_info))
    batch_anova_attr <- unlist(summary(aov(pc~batch, data = pc_metadata)))
    class_anova_attr <- unlist(summary(aov(pc~class, data = pc_metadata)))
    return(c(calc_var_percentage(batch_anova_attr),
             calc_var_percentage(class_anova_attr)))
  }
  var_composition <- sapply(pca_df, calc.pc_var)
  # Dataframe showing prop. var and pi_BE and pi_BV
  pca_var_df <- data.frame(var_pct,
                           batch_pct = var_composition[1,],
                           class_pct = var_composition[2,])
  total_batch_pct <- sum(pca_var_df$var_pct * pca_var_df$batch_pct)
  total_class_pct <- sum(pca_var_df$var_pct * pca_var_df$class_pct)
  print(head(pca_var_df))
  
  library(cluster)
  # Calculating average silhouette scores
  dist_pca <- dist(pca_df, "euclidean")
  batch_df <- silhouette(batch_info, dist_pca)
  class_df <- silhouette(class_info, dist_pca)
  
  # gPCA
  # gpca_obj <- gPCA.batchdetect(t(df), batch_info)
  
  # # Collating results
  # metrics <- c(total_batch_pct*100, total_class_pct*100,
  #              mean(batch_df[,3]), mean(class_df[,3]),
  #              gpca_obj$delta, gpca_obj$p.val)
  # Collating results
  metrics <- c(total_batch_pct*100, total_class_pct*100,
               mean(batch_df[,3]), mean(class_df[,3]))
  
  # names(metrics) <- c("pi_BE", "pi_BV",
  #                     "silhouette_BE", "silhouette_BV",
  #                     "gPCA_delta", "gPCA_pvalue")[]
  print(xtable(t(matrix(metrics)), digits = 3))
  cat(metrics, sep = "\t")
  return(metrics)
}

calc_var_preservation <- function(df_bef, df_aft) {
  return(sum(apply(df_aft, 1, var))/sum(apply(df_bef, 1, var)))
}

#' Calculates proportion of variance in dataframe due to batch effects and biological variable of interest
#' 
#' @param df p x n gene expression dataframe (p: no. of probes, n: no. of samples)
#' @param batch_info vector containing batch labels of samples (ordered the same way as in the dataframe)
#' @param class_info vector containing class labels of samples (ordered the same way as in the dataframe)
#' @return vector containing total proportion of variance in dataframe due to batch effects and biological variable of interest, respectively
calc_var_prop <- function(df, batch_info, class_info) {
  # Tranpose df to n x p
  pca_obj <- prcomp(t(df), center = T, scale. = F)
  # Eigenvalues
  pca_df <- data.frame(pca_obj$x)
  eig_value <- (pca_obj$sdev)^2
  # Percentage variance of each PC
  var_pct <- eig_value/sum(eig_value)
  
  # Libraries: cluster, gpca
  # Argument: PC coordinates for single PC
  calc_pc_var <- function(vec) {
    # Argument: ANOVA attributes; Calculates percentage of variance
    # SS_between/SS_between + SS_within
    calc_var_percentage <- function(vec) unname(vec[3]/(vec[3] + vec[4]))
    pc_metadata <- data.frame(pc = vec,
                              batch = as.factor(batch_info),
                              class = as.factor(class_info))
    batch_anova_attr <- unlist(summary(aov(pc~batch, data = pc_metadata)))
    class_anova_attr <- unlist(summary(aov(pc~class, data = pc_metadata)))
    return(c(calc_var_percentage(batch_anova_attr),
             calc_var_percentage(class_anova_attr)))
  }
  var_composition <- sapply(pca_df, calc_pc_var)
  # Dataframe showing prop. var and pi_BE and pi_BV
  pca_var_df <- data.frame(var_pct,
                           batch_pct = var_composition[1,],
                           class_pct = var_composition[2,])
  total_batch_pct <- sum(pca_var_df$var_pct * pca_var_df$batch_pct)
  total_class_pct <- sum(pca_var_df$var_pct * pca_var_df$class_pct)
  metrics <- c(total_batch_pct, total_class_pct)
  names(metrics) <- c("var_batch", "var_class")
  return(metrics)
}

#' Calculates l2-norm of vector
calcL2Norm <- function(vec) sqrt(sum(vec^2))

#' Calculates cosine similarity between two vectors
#' @return Scalar cos(theta)
calcCosineSim <- function(vec1, vec2) {
  sum(vec1*vec2)/(calcL2Norm(vec1)*calcL2Norm(vec2))
}

#' Converts radians to degrees
rad2degree <- function(rad) rad/pi * 180

#' Calculates angle between two vectors (in degrees)
calcAngleVectors <- function(vec1, vec2) {
  # Prevents error when number is incorrectly calculated..
  # to be slightly above 1
  cosine_sim <- calcCosineSim(vec1, vec2)
  if (cosine_sim > 1) {
    print(sprintf("Cosine similarity: %.20f -> Rounded off!", cosine_sim))
    cosine_sim <- 1
  }
  return(rad2degree(acos(cosine_sim)))
}

#' @param rad rotation angle (counter-clockwise) in radians
#' @return rotation matrix for 2D vector
calcRotationMatrix <- function(rad) matrix(c(cos(rad), -sin(rad),
                                             sin(rad), cos(rad)),
                                           2, 2, byrow = T)

# Generate default ggplot colours
generateGgplotColours <- function(n) {
  hues = seq(15, 375, length = n + 1)
  return(hcl(h = hues, c = 100, l = 65)[1:n])
}
