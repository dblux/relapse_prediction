#' Calculates proportion of variance in dataframe due to batch effects and biological variable of interest
#' 
#' @param df p x n gene expression dataframe (p: no. of probes, n: no. of samples)
#' @param batch_info vector containing batch labels of samples (ordered the same way as in the dataframe)
#' @param class_info vector containing class labels of samples (ordered the same way as in the dataframe)
#' @return Vector containing total proportion of variance in dataframe due to batch effects and biological variable of interest, respectively
calc.explanatory_var <- function(df, batch_info, class_info) {
  # Tranpose df to n x p
  pca_obj <- prcomp(t(df), center = T, scale. = F)
  # Eigenvalues
  pca_df <- data.frame(pca_obj$x)
  eig_value <- (pca_obj$sdev)^2
  # Percentage variance of each PC
  var_pct <- eig_value/sum(eig_value)
  
  # Libraries: cluster, gpca
  # Argument: PC coordinates for single PC
  calc.pc_var <- function(vec) {
    # Argument: ANOVA attributes; Calculates percentage of variance
    # SS_between/SS_between + SS_within
    calc.var_percentage <- function(vec) unname(vec[3]/(vec[3] + vec[4]))
    pc_metadata <- data.frame(pc = vec,
                              batch = as.factor(batch_info),
                              class = as.factor(class_info))
    batch_anova_attr <- unlist(summary(aov(pc~batch, data = pc_metadata)))
    class_anova_attr <- unlist(summary(aov(pc~class, data = pc_metadata)))
    return(c(calc.var_percentage(batch_anova_attr),
             calc.var_percentage(class_anova_attr)))
  }
  var_composition <- sapply(pca_df, calc.pc_var)
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
