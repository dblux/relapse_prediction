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
      
      if (method == 'welch') {
        return(oneway.test(gene ~ batch, X)$p.value)
      } else if (method == 'aov') {
        return(unname(unlist(summary(aov(gene ~ batch, X)))[9]))
      } else if (method == 'kruskal') {
        return(kruskal.test(gene ~ batch, X)$p.value)
      }
    }
    
    apply(X, 1, test_anova, batch, method)
  }
  
  method <- match.arg(method)
  # Subsetting D0 TEL-AML1 remission samples
  annot <- annot[colnames(X), ] # Rearrange annot
  sid_d0 <- rownames(annot)[annot$class_info == 'D0']
  sid_telaml1 <- rownames(annot)[annot$subtype == 'TEL-AML1']
  sid_remission <- rownames(annot)[annot$label == 'Remission']
  sid <- Reduce(intersect, list(sid_d0, sid_telaml1, sid_remission))
  d0_telaml1 <- X[, sid]
  d0_telaml1 <- remove_rows(d0_telaml1, var(row) == 0)
  
  batch <- annot[sid, 'batch_info']
  pvalues <- test_row_anova(d0_telaml1, batch, method)
  
  n_nan <- sum(sapply(pvalues, is.na))
  print(sprintf('No. of NaNs = %d', n_nan))
  
  # Thresholding by p-value
  names(pvalues)[pvalues < alpha & !is.na(pvalues)]
}


# Prediction (Drug genes) --------------------------------------------
## Drug responsive genes
#' @param X_subtype df of patients from a specific subtype (D0 followed by D8)
identify_DE <- function(
  X_subtype, sid_remission, alpha = 0.05, EXPR = 6, N = 50, LOGFC = 1
) {
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
  cat(sprintf(
    "No. of NaN values in log-fc = %d\n",
    sum(is.na(median_logfc))
  ))
  median_logfc1 <- median_logfc[!is.na(median_logfc)]
  
  d0_mu <- rowMeans(D0)
  d8_mu <- rowMeans(D8)
  selected_median_logfc <- median_logfc1[d0_mu > EXPR | d8_mu > EXPR]
  cat(sprintf(
    "No. of probesets excluded by expr threshold = %d\n",
    length(median_logfc1) - length(selected_median_logfc)
  ))
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
  cat(sprintf("No. of features (p-value) = %d\n", length(feat_p)))
  cat(sprintf("No. of features (log2-fc) = %d\n", length(feat_log2fc)))

  intersect(feat_p, feat_log2fc)
}


#' @param response_df dataframe with samples x features
#' @param normal_df dataframe with samples x features
# D0 centroid used to define D0-Normal vector
compute_features <- function(
  response_df, normal_df,
  sid_train, sid_remission
) {
  if (!is_paired(t(response_df)))
    stop("Patient IDs are not paired..")
  
  # Split response df into D0 and D8 df
  n <- nrow(response_df)/2
  d0_df <- response_df[1:n, , drop = F]
  d8_df <- response_df[-(1:n), , drop = F]
  
  # Calculate centroids
  # Only use remission patients in training set to calculate centroid
  sid_leuk <- Reduce(
    intersect,
    list(rownames(d0_df), sid_train, sid_remission)
  )
  
  cat(sprintf("No. of samples in centroid = %d\n", length(sid_leuk)))
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
  erm2_ratio2 <- erm2 / (l2norm_d0_normal - erm2)
  
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
  
  L_D0 <- d0_df - leuk_centroid
  L_D8 <- d8_df - leuk_centroid
  
  # Angle between LD0 and LD8
  angle_LD0_LD8 <- mapply(
    calcAngleVectors,
    data.frame(t(L_D0)),
    data.frame(t(L_D8))
  ) 

  # Angle between LD0 and LN
  angle_LD0_LN <- mapply(
    calcAngleVectors,
    data.frame(t(L_D0)),
    data.frame(matrix(leuk_normal))
  ) 

  # Angle between LD8 and LN
  angle_LD8_LN <- mapply(
    calcAngleVectors,
    data.frame(t(L_D8)),
    data.frame(matrix(leuk_normal))
  )

  angle_LD0_LD8_ratio1 <- angle_LD0_LD8 / angle_LD0_LN
  angle_LD0_LD8_ratio2 <- angle_LD0_LD8 / angle_LD8_LN

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
  features <- data.frame(
    erm1, erm1_ratio1, erm1_ratio2, erm1_ratio3,
    erm2, erm2_ratio, erm2_ratio2,
    erm3, erm3_ratio,
    d0_normal_proj, d8_normal_proj, l2norm_d0_d8, diff_l2norm,
    angle_d0_d8, angle_nd0_nd8, angle_nl_nd8,
    angle_d0d8_normal, angle_d0d8_d0normal,
    angle_d0_normal, angle_d8_normal,
    angle_LD0_LD8, angle_LD0_LN, angle_LD8_LN,
    angle_LD0_LD8_ratio1, angle_LD0_LD8_ratio2, 
    l2norm_d0_normal, l2norm_d8_normal,
    l2norm_ratio1, l2norm_ratio2,
    l2norm_diff, l2norm_diff_ratio
  )
  rownames(features) <- substring(rownames(features), 1, 4)
  
  features
}

#' Calculate probability of remission as percentage of remission cases with
#' scores that are worse than or equal to the current score. Dataframes are
#' assumed to only have selected features.
#'
#' @param X_train dataframe of training set (incl. MRD) with patients x features
#' @param X_predict dataframe containing samples to be predicted (incl. MRD)
#' @param metadata_pid dataframe of metadata with patient x info
#' @param direction character vector containing "<", ">". "<" indicates that a 
#' larger feature value indicates a higher probability of remission.
#' @param samples numeric indicating no. of samples to augment. Defaults to NA.
calc_p_remission_x <- function(
  X_train, X_predict, metadata_pid, direction, samples, include_tp2
) {
  #' Helper function that calculates probability of remission as the number of
  #' of remission cases with a score worse than or equal to the current score
  #' for each feature.
  #'
  #' @param x vector of feature scores from diff samples
  #' @param x_train vector of feature scores from remission samples training set
  #' @param direction character of either "<" or ">"
  calc_p_remission_xi <- function(x, x_train, direction) {
    p_remission <- if (direction == "<") {
      sapply(x, function(x_i) sum(x_i >= x_train) / length(x_train))
    } else if (direction == ">") {
      sapply(x, function(x_i) sum(x_i <= x_train) / length(x_train))
    }
    p_remission
  }

  # assert that X_train does not have null values
  if (sum(is.na(X_train)) != 0 | sum(is.na(X_predict)) != 0)
    stop("Missing values present in either the training or test set.")
  
  X_remission <- X_train[
    metadata_pid[rownames(X_train), "label"] == 'Remission', , drop = F
  ]
  n <- nrow(X_remission)
  cat(sprintf("No. of remission samples in training set = %d\n", n))
  
  # Augment data by simulating samples 
  if (!is.null(samples)) {
    # Estimate parameters of P(x_i|s, y) ~ Normal
    # assert: X_remission is dataframe
    mu_vec <- sapply(X_remission, median)
    sigma_vec <- sapply(X_remission, sd)
    simulated_data <- mapply(
      function(mu, sigma) rnorm(samples, mu, sigma),
      mu_vec, sigma_vec
    )
    X_remission <- rbind(X_remission, simulated_data)
    print(sprintf("Simulated %d samples.", samples))
  }  

  p_remission_xi <- data.frame(mapply(
    calc_p_remission_xi,
    data.frame(X_predict),
    data.frame(X_remission),
    as.list(direction),
    SIMPLIFY = F
  ))
  rownames(p_remission_xi) <- rownames(X_predict)
  colnames(p_remission_xi) <- paste0("p_", colnames(p_remission_xi))
  
  # Without MRD
  # ASSERT: Columns 5, 6 are D33, TP2 MRD
  p_d8 <- rowMeans(p_remission_xi[1:3], na.rm = T)
  p_d33 <- rowMeans(p_remission_xi[1:4], na.rm = T)
  
  label <- as.factor(metadata_pid[rownames(X_predict), "label"])
  
  if (include_tp2) {
    p_tp2 <- rowMeans(p_remission_xi[1:5], na.rm = T)
    proba <- data.frame(p_d8, p_d33, p_tp2, label)
  } else {
    proba <- data.frame(p_d8, p_d33, label)
  }
  # Geometric average
  # p_geomavg <- apply(p_remission_xi, 1, function(x) exp(mean(log(x))))

  list(
    p_remission_xi = p_remission_xi,
    p = proba
  )
}


#' Does not perform PCA transform on data
#' Used to predict relapse for all subtypes
#' X df containing all subtypes of patients in arg: pid and normal patients
#' @param pid vector of pid belonging to both D0 and D8 patients (identically ordered)
#' @param sid_train_test list of length 2 in the form of (sid_train, sid_test)
#' @return list containing prediction plot and vector coordinates
predict_pipeline <- function(
  X_subtype,
  X_normal,
  metadata_sid,
  metadata_pid,
  batch_genes = NULL,
  class_genes = NULL,
  samples = NULL,
  include_tp2 = FALSE,
  sid_train_test = NULL,
  return_features = FALSE,
  features = c(
    "erm1_ratio2", "l2norm_ratio2",
    "angle_LD0_LD8_ratio2", "log_mrd_d33"
  ),
  direction = c("<", "<", "<", ">")
) {
  if (!is.null(sid_train_test)) {
    stopifnot(length(sid_train_test) == 2)
    sid_train <- intersect(sid_train_test[[1]], colnames(X_subtype))
    sid_test <- intersect(sid_train_test[[2]], colnames(X_subtype))
    # assert that D0 and D8 samples match
    stopifnot(is_paired(sid_train)) 
    stopifnot(is_paired(sid_test))
  }
  sid_remission <- colnames(X_subtype)[
    metadata_sid[colnames(X_subtype), "label"] == 'Remission' 
  ]
  # Feature selection 
  if (is.null(class_genes)) { 
    # Identify DE features between D0 and D8 samples
    if (!is.null(sid_train_test)) {
      class_genes <- identify_DE(
        X_subtype[, sid_train, drop = FALSE],
        sid_remission
      )
    } else {
      class_genes <- identify_DE(X_subtype, sid_remission)
    }
    # # Save drug response genes
    # subtype <- unique(metadata_sid[colnames(X_subtype), "subtype"])
    # writeLines(class_genes, sprintf("tmp/response-%s.txt", subtype))
  }
  
  if (is.null(batch_genes)) {
    selected_genes <- class_genes
  } else {
    selected_genes <- setdiff(class_genes, batch_genes)
  }
  cat(sprintf("No. of DE features = %d\n", length(class_genes)))
  cat(sprintf("No. of final features = %d\n", length(selected_genes)))
  
  # Subtype and normal samples
  if (is.null(sid_train_test)) {
    response <- t(X_subtype[selected_genes, ])
    normal <- t(X_normal[selected_genes, ])
    V <- compute_features(response, normal, colnames(X_subtype), sid_remission)
  } else {
    sid <- sort_sid(c(sid_train, sid_test))
    response <- t(X_subtype[selected_genes, sid])
    normal <- t(X_normal[selected_genes, ])
    V <- compute_features(response, normal, sid_train, sid_remission)
  }
  
  # Collate MRD results
  V$log_mrd_d33 <- log10(metadata_pid[rownames(V), "d33_mrd"])
  if (include_tp2) {
    # Include MRD TP2
    V$log_mrd_tp2 <- log10(metadata_pid[rownames(V), "wk12_mrd"]) 
    features <- c(features, "log_mrd_tp2")
    direction <- c(direction, ">")  
  }
  
  V_sub <- V[features] # select features and specify order
  V_sub <- na.omit(V_sub) # Removes patients that have NA MRD values
  pid_omitted <- setdiff(rownames(V), rownames(V_sub))
  if (length(pid_omitted) > 0) {
    sid_omitted <- c(
      paste0(pid_omitted, '_D0'),
      paste0(pid_omitted, '_D8')
    )
    sid_train <- setdiff(sid_train, sid_omitted)
    sid_test <- setdiff(sid_test, sid_omitted)
    cat(sprintf("Omitted patients: %s!\n", pid_omitted))
  }
  stopifnot(!any(is.na(V_sub))) # assert no NA values 
  
  if (return_features) {
    return(na.omit(V))
  }

  # If sid_train_test is supplied, returns predictions for both train and test
  if (!is.null(sid_train_test)) {
    X_train <- V_sub[unique(substring(sid_train, 1, 4)), , drop = F]
    X_test <- V_sub[unique(substring(sid_test, 1, 4)), , drop = F]
    prediction_train <- calc_p_remission_x(
      X_train, X_train, metadata_pid, direction, samples, include_tp2
    )
    prediction_test <- calc_p_remission_x(
      X_train, X_test, metadata_pid, direction, samples, include_tp2
    )
    train_results <- list(
      p_remission_xi = prediction_train$p_remission_xi,
      X_y = cbind(X_train, prediction_train$p)
    )
    test_results <- list(
      p_remission_xi = prediction_test$p_remission_xi,
      X_y = cbind(X_test, prediction_test$p)
    )

    cat("Prediction complete!\n\n")
    return(list(train = train_results, test = test_results))
  } else { 
    prediction <- calc_p_remission_x(
      V_sub, V_sub, metadata_pid, direction, samples, include_tp2
    )
    # Concatenate features and probabilities
    X_y <- cbind(V_sub, prediction$p)

    cat("Prediction complete!\n\n")
    return(list(
      p_remission_xi = prediction$p_remission_xi,
      X_y = X_y
    ))
  }
}
