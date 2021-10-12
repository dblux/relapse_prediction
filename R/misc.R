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
                              batch_info = as.factor(batch_info),
                              class_info = as.factor(class_info))
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