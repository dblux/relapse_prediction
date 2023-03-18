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
  # D0 samples
  OD0 <- response_df[1:n, , drop = F]
  # D8 samples
  OD8 <- response_df[-(1:n), , drop = F]
  
  # Calculate centroids
  # Only use remission patients in training set to calculate centroid
  sid_leuk <- Reduce(intersect,
    list(rownames(OD0), sid_train, sid_remission)
  )
  cat(sprintf("No. of samples in centroid = %d\n", length(sid_leuk)))
  # Leukemia centroid
  OL <- apply(OD0[sid_leuk, , drop = F], 2, median)
  # Normal centroid
  ON <- apply(normal_df, 2, median)
  
  LN <- ON - OL
  unit_LN <- LN / l2norm(LN)
  D0D8 <- OD8 - OD0
  D0N_T <- ON - t(OD0)
  D8N_T <- ON - t(OD8)
  LD0 <- sweep(OD0, 2, OL, FUN = `-`)
  LD8 <- sweep(OD8, 2, OL, FUN = `-`)
  
  # Scalar projection of D0N on LN
  comp_LN_D0N <- colSums(D0N_T * unit_LN)
  # Scalar projection of D8N on LN
  comp_LN_D8N <- colSums(D8N_T * unit_LN)
  
  ### l2norm ###
  l2norm_D0D8 <- apply(D0D8, 1, l2norm)
  l2norm_OD0 <- apply(OD0, 1, l2norm)
  l2norm_OD8 <- apply(OD8, 1, l2norm)
  diff_l2norm <- l2norm_OD8 - l2norm_OD0
  l2norm_D0N <- apply(D0N_T, 2, l2norm)
  l2norm_D8N <- apply(D8N_T, 2, l2norm)
  # l2norm of rejection of ND0 from NL
  l2norm_rej_NL_ND0 <- sqrt(l2norm_D0N ^ 2 - comp_LN_D0N ^ 2) 
  # l2norm of rejection of ND8 from NL
  l2norm_rej_NL_ND8 <- sqrt(l2norm_D8N ^ 2 - comp_LN_D8N ^ 2) 
  
  # Unit vectors
  # Calculate vstack of unit D0-Normal vectors
  unit_D0N_T <- sweep(D0N_T, 2, l2norm_D0N, "/")
  
  # l2norm ratios
  l2norm_ratio1 <- l2norm_D0D8 / l2norm_D0N
  l2norm_ratio2 <- l2norm_D0D8 / l2norm_D8N
  l2norm_diff <- l2norm_D0N - l2norm_D8N
  l2norm_diff_ratio <- l2norm_diff / l2norm_D0D8
  
  ### ERM1 ###
  # Calculate scalar projection by dot product of a and unit b
  erm1 <- colSums(t(D0D8) * unit_LN)
  erm1_ratio1 <- erm1 / abs(comp_LN_D0N)
  erm1_ratio2 <- erm1 / abs(comp_LN_D8N)
  erm1_ratio3 <- erm1 / l2norm_D0D8
  stopifnot(identical(names(erm1), names(erm1_ratio1)))
  # Patients whose D8 have exceeded N along LN
  if (any(comp_LN_D8N < 0)) {
    pid_exceed <- names(comp_LN_D8N)[comp_LN_D8N < 0]
    cat(sprintf("%s has a negative comp_LN_D8N!\n", pid_exceed))
  }
  
  ### ERM2 ###
  # Projection of D0-D8 on D0-N
  erm2 <- colSums(t(D0D8) * unit_D0N_T)
  erm2_ratio1 <- erm2 / l2norm_D0N
  erm2_ratio2 <- erm2 / (l2norm_D0N - erm2)
  stopifnot(identical(names(erm2), names(erm2_ratio1)))
  
  ### ERM3 ###
  ## Along a chosen PC that represents timepoint
  PC <- 1
  # Be careful of direction of D0-N (may be negative)
  # If negative, a larger shift will lead to a smaller ERM3
  dir <- sign(median(normal_df[,PC]) - median(OD0[,PC]))
  erm3 <- (OD8[,PC] - OD0[,PC]) * dir # direction is normalised
  # Divide by D0-Normal along PC
  erm3_ratio <- erm3 / (median(normal_df[,PC]) - OD0[,PC])
  stopifnot(identical(names(erm3), names(erm3_ratio)))
  
  # Angle between D0-D8 and Leuk-Normal
  angle_D0D8_LN <- apply(
    D0D8, 1, function(row_vec) calcAngleVectors(row_vec, LN)
  )
  # Angle between D0-D8 and D0-Normal
  angle_D0D8_D0N <- mapply(
    calcAngleVectors,
    data.frame(t(D0D8)),
    data.frame(D0N_T)
  )
  # Angle between O-D0 and O-D8
  angle_OD0_OD8 <- mapply(
    calcAngleVectors,
    data.frame(t(OD0)),
    data.frame(t(OD8))
  )
  # Angle between O-D0 and O-N
  angle_OD0_ON <- apply(
    OD0, 1, function(row_vec) calcAngleVectors(row_vec, ON)
  )
  # Angle between O-D8 and O-N
  angle_OD8_ON <- apply(
    OD8, 1, function(row_vec) calcAngleVectors(row_vec, ON)
  )
  # Angle between N-D0 and N-D8
  # Equivalent to angle between D0-N and D8-N
  angle_ND0_ND8 <- mapply(
    calcAngleVectors,
    data.frame(D0N_T),
    data.frame(D8N_T)
  )
  # Angle between NL and ND0
  # Equivalent to angle between LN and D0-N
  angle_NL_ND0 <- sapply(data.frame(D0N_T), function(x) calcAngleVectors(x, LN))
  angle_NL_ND8 <- sapply(data.frame(D8N_T), function(x) calcAngleVectors(x, LN))
  angle_LD0_LD8 <- mapply(
    calcAngleVectors,
    data.frame(t(LD0)),
    data.frame(t(LD8))
  ) 
  angle_LD0_LN <- mapply(
    calcAngleVectors,
    data.frame(t(LD0)),
    data.frame(matrix(LN))
  ) 
  angle_LD8_LN <- mapply(
    calcAngleVectors,
    data.frame(t(LD8)),
    data.frame(matrix(LN))
  )

  angle_LD0_LD8_ratio1 <- angle_LD0_LD8 / angle_LD0_LN
  angle_LD0_LD8_ratio2 <- angle_LD0_LD8 / angle_LD8_LN
  

  features <- data.frame(
    erm1, erm1_ratio1, erm1_ratio2, erm1_ratio3,
    erm2, erm2_ratio1, erm2_ratio2, erm3, erm3_ratio,
    l2norm_D0N, l2norm_D8N,
    l2norm_ratio1, l2norm_ratio2,
    l2norm_diff, l2norm_diff_ratio, l2norm_D0D8, diff_l2norm,
    l2norm_rej_NL_ND0, l2norm_rej_NL_ND8, 
    comp_LN_D0N, comp_LN_D8N,
    angle_NL_ND0, angle_NL_ND8,
    angle_D0D8_LN, angle_D0D8_D0N,
    angle_LD0_LD8, angle_LD0_LN, angle_LD8_LN,
    angle_OD0_ON, angle_OD8_ON, angle_OD0_OD8, angle_ND0_ND8,
    angle_LD0_LD8_ratio1, angle_LD0_LD8_ratio2 
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
  return_genes = FALSE,
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
  
  if (return_genes)
    return(selected_genes)

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
  
  if (return_features)
    return(V_sub) # returns patients that have NA MRD values
  
  V_sub <- na.omit(V_sub) # Removes patients that have NA MRD values
  pid_omitted <- setdiff(rownames(V), rownames(V_sub))
  print(pid_omitted)
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


#' Constructs linear models between all pairwise combinations of samples
#'
#' @param X dataframe of dim (n_features, n_samples)
#' @param metadata dataframe of dim (n_samples, :)
#' @param pos_ctrl positive control features 
#' @param batch numeric of batch number to be imputed
#' @return object of class pairwise, which is a list containing:
#'   \item{rms} matrix of residual mean squares with dim (n_samples, n_samples)
#'   \item{intercept} matrix of intercepts of all linear models with dim (n_samples, n_samples)
#'   \item{residuals} array of residuals with dim (n_samples, n_samples, n_features)
pairwise_lm <- function(
  X, metadata, pos_ctrl, batch = NULL 
) {
  if (is.null(batch)) {
    ids <- colnames(X)
  } else {
    is_batch <- metadata[colnames(X), "batch_info"] == batch
    ids <- colnames(X)[is_batch]
  }

  rms <- intercepts <- matrix(
    NA, length(ids), ncol(X),
    dimnames = list(ids, colnames(X))
  )
  residuals <- array(
    NA, dim = c(length(ids), ncol(X), length(pos_ctrl)),
    dimnames = list(ids, colnames(X), pos_ctrl) 
  )
  for (i in ids) {
    for (j in colnames(X)) {
      if (i != j) {
        # print(sprintf("%s ~ %s", i, j))
        # only use positive control features to fit linear model
        pair <- X[pos_ctrl, c(i, j)]
        # only keep features where both values are not missing 
        pair <- log2(pair[rowSums(pair != 0) == 2, ])
        colnames(pair) <- c("y", "x")
        # slope of lm is fixed at 1
        model <- lm(y - x ~ 1, data = pair)
        residuals[i, j, rownames(pair)] <- model$residuals
        # used residual mean square (rms) instead of rss with df of n - 1
        # different pairs have different number of present features
        rms[i, j] <- sum(model$residuals^2) / (nrow(pair) - 1)
        intercepts[i, j] <- model$coefficients
      }
    }
  }
  pairwise <- list(
    rms = rms, intercepts = intercepts, residuals = residuals
  )
  class(pairwise) <- "pairwise"

  pairwise
}


#' @param features character of features to test for consistent bias. features
#' has to be be in pairwise$residuals
# Assumption: residuals is a symmetric matrix
bias.pairwise <- function(
  pairwise, metadata, group, features = NULL,
  missing_threshold = 0.5, min_sample = 5
) {
  if (is.null(features))
    features <- dimnames(residuals)[[3]]

  residuals <- pairwise$residuals
  ids <- colnames(residuals)
  batches <- as.character(sort(unique(metadata[ids, group])))
  pvalues <- array(
    NA, dim = c(length(batches), length(batches), length(features)),
    dimnames = list(batches, batches, features)
  )
  means <- array(
    NA, dim = c(length(batches), length(batches), length(features)),
    dimnames = list(batches, batches, features)
  )
  for (feature in features) {
    residuals_feature <- as.numeric(residuals[, , feature])
    pct_missing <- sum(is.na(residuals_feature)) / length(residuals_feature)
    if (pct_missing > missing_threshold) {
      cat(sprintf(
        "%s has %.1f%% missing residuals.\n", feature, pct_missing * 100
      ))
      next
    }
    for (i in seq_along(batches)) {
      batch_i <- ids[metadata[ids, group]  == batches[i]]
      for (j in seq(i)) {
        batch_j <- ids[metadata[ids, group]  == batches[j]]
        residuals_vec <- if (i != j) {
          as.numeric(residuals[batch_i, batch_j, feature])
        } else {
          # dealing with a symmetric matrix
          residuals_mat <- residuals[batch_i, batch_j, feature]
          residuals_mat[lower.tri(residuals_mat)]
        }
        # t.test handles na according to na.action (default is na.omit)
        n_residuals <- length(na.omit(residuals_vec))
        if (n_residuals < min_sample) {
          cat(sprintf(
            "Small sample size of %d between %s and %s. T-test not performed.\n",
            n_residuals, batches[[i]], batches[[j]]
          ))
          next
        }
        ttest <- t.test(residuals_vec)
        pvalues[batches[i], batches[j], feature] <- ttest$p.value
        means[batches[i], batches[j], feature] <- mean(residuals_vec, na.rm = TRUE)
      }
    }
  }

  list(pvalues = pvalues, means = means)
}


#' Impute missing values using pairwise linear models
#'
#' rownames(pairwise$rms) are the samples that will be imputed
#' @param features character of features to be imputed
impute.pairwise <- function(pairwise, X, features, k = 5) {
  rms <- pairwise$rms
  intercepts <- pairwise$intercepts
  samples <- rownames(rms)  # samples to be imputed
  is_present <- X[features, ] != 0 
  n_present <- colSums(is_present)
  
  # # filter out samples that cannot predict a single feature
  # # shift up to pairwise function
  # rest_n_present <- n_present[colnames(r_sq)]
  # ids_ms <- names(rest_n_present)[rest_n_present == 0]
  # if (length(ids_ms) != 0)
  #   stop('TODO: Filter out samples that cannot be used.')
  
  # each id has a predictions matrix of shape (n_missing_features, k)
  # different ids have different missing features
  knn_predictions <- kneighbours <- vector("list", length = nrow(rms))
  names(knn_predictions) <- samples 
  names(kneighbours) <- samples 
  for (id in samples) {
    missing_features <- features[!is_present[, id]]
    nearest_neighbours <- colnames(rms)[
      order(as.numeric(rms[id, ]), decreasing = TRUE)]
    
    predictions <- kns <- matrix(
      NA, length(missing_features), k,
      dimnames = list(missing_features, paste0("nn", seq(k))) 
    )
    # if all features can be represented by k nearest neighbours
    knns <- nearest_neighbours[seq(k)]
    knns_npresent <- colSums(is_present[, knns])
    if (all(knns_npresent == length(missing_features))) {
      cat("knns are the same across all missing features!\n")
      for (i in seq(k)) {
        nn <- nearest_neighbours[i]
        b_0 <- intercepts[id, nn]
        for (feature in missing_features) {
          predictions[feature, i] <- 2^(b_0 + log2(X[feature, nn]))
        }
      }
      knn_predictions[[id]] <- predictions
      kneighbours[[id]] <- t(replicate(
        length(missing_features), nearest_neighbours[seq(k)]
      ))
    } else{
      # select the best knns for each feature
      # check whether all features have enough samples with non-zero values
      if (any(rowSums(is_present) < k)) {
        stop(paste(
          "Insufficient samples with non-missing values to impute some",
          "features. Try decreasing k."
        ))
      }
      # count to keep track of no. of values for each feature
      cnt <- rep(1, length(missing_features))
      names(cnt) <- missing_features
      i <- 1
      while (any(cnt <= k)) {
        nn <- nearest_neighbours[i]
        b_0 <- intercepts[id, nn]
        for (feature in missing_features) {
          if (is_present[feature, nn] && cnt[feature] <= k) {
            predictions[feature, cnt[feature]] <-
              2^(b_0 + log2(X[feature, nn]))
            kns[feature, cnt[feature]] <- nn
            cnt[feature] <- cnt[feature] + 1
          }
        }
        i <- i + 1
      }
      knn_predictions[[id]] <- predictions
      kneighbours[[id]] <- kns
    }
  }
  
  # consolidate predictions and execute prediction
  for (id in names(knn_predictions)) {
    predictions <- knn_predictions[[id]]
    values <- apply(predictions, 1, median)
    for (feature in names(values)) {
      X[feature, id] <- values[feature]
    }
  }

  list(X = X, knn_predictions = knn_predictions, kneighbours = kneighbours)
}


#' Scales affymetrix probesets with batch scale factor
#'
#' @param X dataframe of shape (n_features, n) with raw expression values
#' @param metadata dataframe of shape (n, n_info)
#' @param affy_pos_ctrl character containing positive control affymetrix probesets
#' @param batches numeric consisting of target batches to scale
#' @param @return dataframe of shape (n_affy, n)
scale_batch_affy <- function(
  X, metadata, affy_pos_ctrl, trim = 0.02,
  batches = NULL, plot = FALSE
) {
  if (is.null(batches))
    batches <- sort(unique(metadata[colnames(X), "batch_info"]))
  batches <- as.character(batches)
  classes <- as.character(sort(unique(metadata[colnames(X), "class_info"])))

  is_affy <- startsWith(rownames(X), 'AFFX')
  X_affy <- X[affy_pos_ctrl, ]
  X_notaffy <- X[!is_affy, ]
  affy_ratio <- data.frame(
    # mean of non-missing values only!
    mean_affy = sapply(X_affy, mean, trim = trim),
    mean_notaffy = sapply(X, mean, trim = trim),
    metadata[colnames(X), ]
  )
  
  for (cls in classes) {
    slopes <- sapply(batches, function(batch) {
      grp_ratio <- subset(affy_ratio, batch_info == batch & class_info == cls)
      # linear model with no intercept
      unname(coef(lm(mean_notaffy ~ mean_affy + 0, data = grp_ratio)))
    })
    target_slope <- median(slopes)
    cat(sprintf("Class %s (median slope = %.3f):\n", cls, target_slope))
    
    for (batch in names(slopes)) {
      slope <- slopes[batch]
      scale_factor <- target_slope / slope
      samples <- rownames(subset(
        affy_ratio, batch_info == batch & class_info == cls
      ))
      X_notaffy[, samples] <- X_notaffy[, samples] * scale_factor
      cat(sprintf(
        "Scaled batch %s with original slope of %.3f to %.3f\n",
        batch, slope, target_slope
      ))
    }
    cat("\n")
  }  
  X_scaled <- rbind(X_notaffy, X[is_affy, ])
  if (!identical(rownames(X), rownames(X_scaled)))
    warning("Order of rows have not been preserved!")
  
  X_scaled
}
