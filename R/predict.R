#' Identify features from dataframe that are afflicted with batch effects
#' @param X dataframe with features x samples
identify_batch_features <- function(
  X, annot, method = c('welch', 'aov', 'kruskal'), alpha = .05
) {
  #' @param X matrix with samples as rows and features as columns
  test_row_anova <- function(X, batch, method) {
    #' @param x vector of feature values
    #' @param batch vector indicating the batch sample values belong to
    test_anova <- function(x, batch, method) {
      X <- data.frame(gene = x, batch = as.factor(batch))
      
      if (method == "welch") {
        return(oneway.test(gene ~ batch, X)$p.value)
      } else if (method == "aov") {
        return(unname(unlist(summary(aov(gene ~ batch, X)))[9]))
      } else if (method == "kruskal") {
        return(kruskal.test(gene ~ batch, X)$p.value)
      }
    }
    
    apply(X, 1, test_anova, batch, method)
  }
  
  method <- match.arg(method)
  # Subsetting D0 TEL-AML1 remission samples
  annot <- annot[colnames(X), ] # Rearrange annot
  sid_d0 <- rownames(annot)[annot$class_info == "D0"]
  sid_telaml1 <- rownames(annot)[annot$subtype == "TEL-AML1"]
  sid_remission <- rownames(annot)[annot$label == 0]
  sid <- Reduce(intersect, list(sid_d0, sid_telaml1, sid_remission))
  d0_telaml1 <- X[, sid]
  d0_telaml1 <- remove_rows(d0_telaml1, var(row) == 0)
  
  batch <- annot[sid, "batch_info"]
  pvalues <- test_row_anova(d0_telaml1, batch, method)
  
  n_nan <- sum(sapply(pvalues, is.na))
  print(sprintf("No. of NaNs = %d", n_nan))
  
  # Thresholding by p-value
  names(pvalues)[pvalues < alpha & !is.na(pvalues)]
}


# Prediction (Drug genes) --------------------------------------------
## Drug responsive genes
#' @param X_subtype df of patients from a specific subtype (D0 followed by D8)
getLocalGenes <- function(X_subtype, sid_remission,
                          alpha = 0.05, EXPR = 6, N = 50, LOGFC = 1) {
  if (!is_paired(X_subtype))
    stop("Patient IDs are not paired..")
  
  sid_idx <- intersect(sid_remission, colnames(X_subtype))
  X_subtype_remission <- X_subtype[, sid_idx, drop = F]
  n_pairs <- ncol(X_subtype_remission) / 2
  
  # P-value
  pvalue <- calc_ttest(X_subtype_remission, n_pairs, is_paired = T) # nan values!

  # # Q-value
  # calc_qvalue <- function(p) length(p)*p/rank(p)
  # qvalue <- calc_qvalue(pvalue) # FDR threshold
  # hist(qvalue, breaks =20)
  D0 <- X_subtype_remission[, 1:n_pairs, drop = F]
  D8 <- X_subtype_remission[, -(1:n_pairs), drop = F]
  
  # Median paired log-FC
  paired_logfc <- D8 - D0
  median_logfc <- apply(paired_logfc, 1, median)
  print(sprintf("No. of NaN values in log-fc = %d",
                 sum(is.na(median_logfc))))
  median_logfc1 <- median_logfc[!is.na(median_logfc)]
  
  d0_mu <- rowMeans(D0)
  d8_mu <- rowMeans(D8)
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


 #' @param response_df dataframe with samples x features
#' @param normal_df dataframe with samples x features
# D0 centroid used to define D0-Normal vector
compute_features <- function(
  response_df, normal_df,
  sid_train, sid_remission
) {
  # Split response df into D0 and D8 df
  n <- nrow(response_df)/2
  d0_df <- response_df[1:n, , drop = F]
  d8_df <- response_df[-(1:n), , drop = F]
  
  if (!is_paired(t(response_df)))
    stop("Patient IDs are not paired..")
  
  # Calculate centroids
  # Only use remission patients in training set to calculate centroid
  sid_leuk <- Reduce(
    intersect,
    list(rownames(d0_df), sid_train, sid_remission)
  )
  
  print(sprintf("NO. OF SAMPLES IN CENTROID: %d", length(sid_leuk)))
  leuk_centroid <- apply(d0_df[sid_leuk, , drop = F], 2, median)
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
  ## Along a chosen PC that represents timepoint
  PC <- 1
  # Be careful of direction of D0-N (may be negative)
  # If negative, a larger shift will lead to a smaller ERM3
  dir <- sign(median(normal_df[,PC]) - median(d0_df[,PC]))
  erm3 <- (d8_df[,PC] - d0_df[,PC]) * dir # direction is normalised
  # Divide by D0-Normal along PC
  erm3_ratio <- erm3/(median(normal_df[,PC]) - d0_df[,PC])
  
  stopifnot(identical(names(erm3), names(erm3_ratio)))
  
  ### l2norm ###
  l2norm_d0_d8 <- apply(d0_d8_hstack, 1, calcL2Norm)
  l2norm_d0 <- apply(d0_df, 1, calcL2Norm)
  l2norm_d8 <- apply(d8_df, 1, calcL2Norm)
  diff_l2norm <- l2norm_d8 - l2norm_d0
  
  ### Angle between D0-D8 and Leuk-Normal ###
  angle_d0d8_normal <- apply(
    d0_d8_hstack, 1, function(row_vec) calcAngleVectors(row_vec, leuk_normal)
  )
  
  ### Angle between D0-D8 and D0-Normal ###
  angle_d0d8_d0normal <- mapply(calcAngleVectors,
                                data.frame(t(d0_d8_hstack)),
                                data.frame(d0_normal_vstack))
  
  ### Angle between D0 and D8 ###
  angle_d0_d8 <- mapply(calcAngleVectors,
                        data.frame(t(d0_df)), data.frame(t(d8_df)))
  
  ### Angle between D0 and normal ###
  angle_d0_normal <- apply(
    d0_df, 1, function(row_vec) calcAngleVectors(row_vec, normal_centroid)
  )
  
  ### Angle between D8 and Normal ###
  angle_d8_normal <- apply(
    d8_df, 1, function(row_vec) calcAngleVectors(row_vec, normal_centroid)
  )
  
  ### Angle between N-D0 and N-D8 ###
  # Equivalent to angle between D0-N and D8-N
  angle_nd0_nd8 <- mapply(calcAngleVectors,
                          data.frame(d0_normal_vstack),
                          data.frame(d8_normal_vstack))
  
  ### Angle between N-centroid(D0) N-D8 ###
  # Equivalent to angle between centroid(D0)-N and D8-N
  angle_nl_nd8 <- sapply(data.frame(d8_normal_vstack),
                         function(x, y) calcAngleVectors(x, y),
                         leuk_normal)
  
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
    diff_l2norm, angle_d0_d8, angle_nd0_nd8, angle_nl_nd8,
    angle_d0d8_normal, angle_d0d8_d0normal,
    angle_d0_normal, angle_d8_normal,
    l2norm_d0_normal, l2norm_d8_normal,
    l2norm_ratio1, l2norm_ratio2,
    l2norm_diff, l2norm_diff_ratio,
    erm1_ratio2, erm1_ratio3
  )
  
  rownames(features_df) <- substring(rownames(features_df), 1, 4)
  return(features_df)

}


#' Calculate probability of remission as percentage of remission cases with
#' scores that are worse than or equal to the current score
#' @param X_train dataframe of training set (incl. MRD) with patients x features
#' @param Y dataframe of metadata with samples x info
#' @param bigpos_names vector of feature names where bigger is positive
#' @param smallpos_names vector of feature names where smaller is positive
#' @param X_predict dataframe containing samples to be predicted (incl. MRD)
calc_p_remission_x <- function(X_train, Y,
                             bigpos_names,
                             smallpos_names,
                             X_predict) {
  #' Pct of remission cases with worse than or equal to the current score
  #' Bigger values indicate it being worse
  #' @param x vector of feature scores from diff samples
  #' @param x_remission vector of feature scores from relapse samples
  calc_pct_remission_xi <- function(x, x_remission) {
    sapply(
      x,
      function(x_i) sum(x_i <= x_remission) / length(x_remission)
    )
  }

  #' Standardises features so that bigger values indicate positive label
  select_orientate_features <- function(X, bigpos_names, smallpos_names) {
    # Assumption: Either smallpos_names or bigpos_names will not be NULL
    if (is.null(smallpos_names)) {
      return(X[, bigpos_names, drop =  F])
    } else if (is.null(bigpos_names)) {
      # reverse order - bigger values now indicate relapse
      return(-X[, smallpos_names, drop =  F])
    } else {
      return(cbind(-X[, smallpos_names, drop =  F],
                   X[, bigpos_names, drop =  F]))
    }
  }
  
  X1_train <- select_orientate_features(X_train, bigpos_names, smallpos_names)
  idx <- paste0(rownames(X1_train), "_D0") # in order to access metadata
  Y1_train <- Y[idx, , drop = F]
  X1_train_remission <- X1_train[Y1_train$label == 0,]  # Only remission

  X1_predict <- select_orientate_features(X_predict, bigpos_names, smallpos_names)

  pct_remission <- mapply(calc_pct_remission_xi,
                          data.frame(X1_predict),
                          data.frame(X1_train_remission),
                          SIMPLIFY = F)
  pct_remission <- data.frame(pct_remission)
  rownames(pct_remission) <- rownames(X1_predict)
  
  # Without MRD
  pct_wo_mrd <- pct_remission[, colnames(pct_remission) != "mrd"]
  p_wo_mrd <- apply(pct_wo_mrd, 1, mean, na.rm = T)
  
  # Arithmetic average
  p_avg <- apply(pct_remission, 1, mean, na.rm = T)
  
  # Geometric average
  p_geomavg <- apply(pct_remission, 1, function(x) exp(mean(log(x))))
                     
  label <- as.factor(
    Y[paste0(rownames(X1_predict), "_D0"), "label"]
  )

  data.frame(
    pid = rownames(pct_remission),
    label = label,
    p = p_avg,
    p_wo_mrd = p_wo_mrd,
    pct_remission
  )
}


#' ASSUMPTION: X_train is filtered of NA MRD values and contains all features!
#' @param X_train dataframe of training set (incl. MRD) with patients x features
#' @param Y dataframe of metadata with samples x info
#' @param bigpos_names vector of feature names where bigger is positive
#' @param smallpos_names vector of feature names where smaller is positive
#' @param X_test dataframe of test set (incl. MRD)
predict_plot <- function(X_train, Y,
                         bigpos_names,
                         smallpos_names,
                         X_test = NULL) {  
  # If test set is present, predict test set
  if (is.null(X_test)) {
    X_predict <- X_train
  } else {
    X_predict <- X_test
  }
  
  p_remission_x <- calc_p_remission_x(
    X_train, Y,
    bigpos_names,
    smallpos_names,
    X_predict
  )
  proba <- p_remission_x # OPTION!
  
  # Select features
  X_fltr_train <- X_train[, c(bigpos_names, smallpos_names)]
  X_fltr_predict <- X_predict[, c(bigpos_names, smallpos_names)]
  
  # Select p(remission|x)
  p <- proba[, "p", drop = F] # OPTION!
  colnames(p) <- "p_rem"
  # Concatenate features and probabilities
  X_y <- cbind(
    X_fltr_predict,
    p,
    label = proba$label
  )
  
  X_y$mrd <- log10(X_y$mrd) # log-transform mrd
  colnames(X_y)[colnames(X_y) == "mrd"] <- "log_mrd"

  long_X_y <- melt(X_y, id = "label", variable.name = "feature")             
  
  FEAT_ORDER <- c(
    "erm1_ratio2", "l2norm_ratio2",
    "angle_d0d8_d0normal", "log_mrd", "p_rem"
  )
  FEAT_LABS <- c(
    "'ERM Ratio'", "'ARM Ratio'", "theta",
    "log[10](MRD)", "paste('P(Remission|', bold(x), ')')"
  )
  long_X_y$feature <- factor(
    long_X_y$feature,
    levels = FEAT_ORDER,
    labels = FEAT_LABS
  ) # Reorder levels

  ##### PLOTS #####
  ax_jitter <- ggplot(
    long_X_y,
    aes(x = feature, y = value, colour = label)
  ) +
    geom_boxplot(alpha = 0, show.legend = F) +
    geom_point(position = position_jitterdodge(),
               cex = 2, show.legend = F) +
    scale_color_manual(values = COL_LABEL) +
    facet_wrap(
      ~feature,
      nrow = 1, scales = "free",
      labeller = label_parsed
    ) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "none"
    )
  
  # Jitter plot: p-value label
  # Both group sizes must be > 1
  if (length(table(X_y$label)) > 1 && min(table(X_y$label)) > 1) {
    list_p_rem <- split(X_y$p_rem, X_y$label)
    
    try({
      ttest <- t.test(list_p_rem[[1]], list_p_rem[[2]])
      p_lab <- sprintf("p = %.3f", ttest$p.value)
      
      ann_text <- data.frame(
        feature = factor(
          FEAT_ORDER[5],
          levels = FEAT_ORDER,
          labels = FEAT_LABS
        ),
        value = Inf,
        lab =  p_lab
      )

      ax_jitter <- ax_jitter +
        geom_text(data = ann_text,
                  aes(x = feature, y = value, label = lab),
                  size = 3, colour = "black",
                  vjust = 3, hjust = -0.1)
      })
  }
  
  ## Plot: Parallel coordinates - Pct
  proba1 <- proba[, !(colnames(proba) == "p_wo_mrd")]
  long_proba <- melt(proba1, id = c("pid", "label"),
                    variable.name = "feature")
             
  ax_parallel <- ggplot(long_proba,
                        aes(feature, value, colour = label, group = pid)) +
    geom_line(show.legend = F) +
    scale_color_manual(values = COL_LABEL)
  
  ## PLOT: CDF
  emp_cdf <- ggplot(proba, aes(x = p, colour = label)) +
    stat_ecdf(show.legend = F) +
    scale_color_manual(values = COL_LABEL)
  
  ## PLOT: RELATIVE RISK & ODDS RATIO
  p_sorted <- proba[order(proba$p),]
  p_sorted$label <- as.numeric(as.character(p_sorted$label))
  p_sorted$total_le <- rank(p_sorted$p, ties.method = "max")
  p_sorted$total_g <- nrow(p_sorted) - p_sorted$total_le
  p_sorted$relapse_le <- sapply(p_sorted$total_le,
                                function(i) sum(p_sorted$label[1:i]))
  p_sorted$relapse_g <- sum(p_sorted$label) - p_sorted$relapse_le
  
  p_sorted <- within(
    p_sorted,
    relative_risk <- (relapse_le/total_le) / (relapse_g/total_g)
  )
  
  p_sorted <- within(
    p_sorted,
    odds_ratio <- (relapse_le/(total_le-relapse_le)) / (relapse_g/(total_g-relapse_g))
  )
                                 
  ax_rr_or <- ggplot(p_sorted) +
    geom_step(aes(p, relative_risk, colour = "RR"), direction = "hv") + 
    geom_step(aes(p, odds_ratio, colour = "OR"), direction = "hv") +
    scale_color_manual("",
                       breaks = c("RR", "OR"),
                       values = c("RR" = "orange", "OR" = "steelblue3")) +
    theme(axis.title.y = element_blank())
  
  ## Plot: ROC
  # ERM1 evaluated is not from global GSS model
  proba_x <- cbind(proba, erm = X_predict$erm1, d33_mrd = X_predict$mrd) # subset mrd
                                
  x_names <- c("p", "erm", "d33_mrd")
  # WARNING: Change bigger.positive according to features!
  bigger.positive <- c(F, T, F) # bigger means relapse
  
  # ROC can only be plotted when there are both positive and negative samples
  if (all(table(proba_x$label) != 0)) {
    ax_roc <- plot_roc(proba_x, "label", x_names)
    # Able to plot ROC
    ax2 <- plot_grid(ax_parallel, ax_roc,
                     ncol = 2, rel_widths = c(1.8, 1))
  } else{
    ax2 <- ax_parallel # unable to plot ROC
  }
  
  # Plot: MRD v.s. Risk of relapse
  mrd_p <- ggplot(proba_x) +
    geom_point(aes(p, log10(d33_mrd), colour = label),
               cex = 3, show.legend = F) +
    scale_color_manual(values = COL_LABEL)
                                
  ax1 <- plot_grid(ax_jitter, mrd_p,
                   ncol = 2, rel_widths = c(2.8, 1))
  
  fig <- plot_grid(ax1, ax2, nrow = 2)
  
  list(
    p_rem = p,
    P = proba,
    X_y = X_y,
    plot = fig
  )
}


#' Does not perform PCA transform on data
#' Used to predict relapse for all subtypes
#' X df containing all subtypes of patients in arg: pid and normal patients
#' @param pid vector of pid belonging to both D0 and D8 patients (identically ordered)
#' @return list containing prediction plot and vector coordinates
predict_pipeline <- function(X_subtype, X_normal,
                             metadata, metadata_mrd,
                             batch_genes = NULL) {
  sid_remission <- colnames(X_subtype)[
    metadata[colnames(X_subtype), "label"] == 0
  ]
  
  class_genes <- getLocalGenes(X_subtype, sid_remission)
  
  if (is.null(batch_genes)) {
    selected_genes <- class_genes
  } else {
    selected_genes <- setdiff(class_genes, batch_genes)
  }
  
  print(c("No. of selected genes = ", length(class_genes)))
  print(c("No. of final genes = ", length(selected_genes)))
  
  # Subtype and normal samples
  response <- t(X_subtype[selected_genes, ])
  normal <- t(X_normal[selected_genes, ])
  
  # Collate MRD results as well
  V <- compute_features(response, normal, colnames(X_subtype), sid_remission)
  V$mrd <- metadata_mrd[rownames(V), "d33_mrd"]
  
  prediction_obj <- predict_plot(
    V, metadata,                             
    bigpos_names = "angle_d0d8_d0normal",
    smallpos_names = c("erm1_ratio2", "l2norm_ratio2", "mrd")
  )
  
  return(prediction_obj)
}