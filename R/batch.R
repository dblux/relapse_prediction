# obtain beta matrix for each batch i by g
# quality check on matrix to make sure var is tight
# obtain beta vector

#' Prepends if string values start with 0-9
prepend_ifnumeric <- function(x, prefix) {
  x <- as.character(x)
  has_numeric_prefix <- grepl("[0-9]", x)
  x[has_numeric_prefix] <- paste0(prefix, x[has_numeric_prefix])
  x
}


# TODO: compute beta: using median vs all samples
#' Compute batch scaling factor for each feature
#' If batch scaling factors are unable to be computed due to missing values,
#' NA will be assigned as the value of beta 
estimate_parameters <- function(X, metadata) {
  # probeset i, sample j, batch k, class g
  #' If vector has less than 3 non-zero values: return NA
  #' Median of vector with all zero values results in NA
  median_nonzero <- function(x, ...) {
    if (length(x[x != 0]) < 3) {
      cnt <- cnt + 1
      return(NA)
    }
    median(x[x != 0], ...)
  }
  mean_nonzero <- function(x, ...) {
    if (length(x[x != 0]) < 3) {
      cnt <- cnt + 1
      return(NA)
    }
    mean(x[x != 0], ...)
  }
  cnt <- 0
  metadata <- metadata[colnames(X), , drop = FALSE] # metadata only of samples in X
  sample_classes <- prepend_ifnumeric(metadata$class_info, "C")
  sample_batches <- prepend_ifnumeric(metadata$batch_info, "B")
  n_classes <- length(unique(sample_classes))
  n_batches <- length(unique(sample_batches))
  # initialise arrays
  beta <- mu <- array(
    NA, dim = c(nrow(X), n_batches, n_classes),
    dimnames = list(rownames(X), unique(sample_batches), unique(sample_classes))
  )
  ref_mu <- array(
    NA, dim = c(nrow(X), n_classes),
    dimnames = list(rownames(X), unique(sample_classes))
  )
  X_classes <- split.default(X, sample_classes)
  
  # Estimation of beta
  for (g in names(X_classes)) {
    X_g <- X_classes[[g]]
    batch_g <- prepend_ifnumeric(metadata[colnames(X_g), "batch_info"], "B")
    X_batches <- split.default(X_g, batch_g)
    mu_g <- sapply(X_batches, apply, 1, mean_nonzero, na.rm = TRUE)
    ref_mu_g <- apply(mu_g, 1, median_nonzero, na.rm = TRUE)
    if (any(is.na(ref_mu_g)))
      warning(sprintf(
        "Class %s: %d out of %d features have median of batch medians with value zero.",
        as.character(g), sum(is.na(ref_mu_g)), nrow(X)
      ))
    beta_g <- data.frame(mu_g / ref_mu_g)
    for (k in colnames(mu_g)) {
      mu[, k, g] <- mu_g[, k]
      beta[, k, g] <- beta_g[, k]
    }
    ref_mu[, g] <- ref_mu_g
  }
  
  # beta_arr <- abind::abind(list_betas, rev.along = 0) # beta for all classes
  beta_hat <- apply(beta, c(1, 2), mean, na.rm = FALSE) #
  beta_hat[is.nan(beta_hat)] <- NA # mean of vector with all NA values returns NaN
  cat(sprintf(
    "No. of NAs in beta_hat: %d/%d\n",
    sum(is.na(beta_hat)), length(beta_hat)
  ))
  beta_hat[is.na(beta_hat)] <- 1 # replace NA with 1 (i.e. do not correct if beta = NA)
  beta_sigma2 <- apply(beta, c(1, 2), var, na.rm = FALSE)
  
  ## Estimating outlier gamma
  # Compute class pair ratios
  n_pairs <- 1
  # Assume: 2 classes have been chosen for beta
  pair_classes <- c("A", "B")
  rho <- mu[, , pair_classes[1]] / mu[, , pair_classes[2]]
  rho_mu <- apply(rho, 1, median_nonzero, na.rm = TRUE)
  cat(sprintf("No. of features with missing rho means = %d\n", sum(is.na(rho_mu))))
  rho_sigma <- apply(rho, 1, sd, na.rm = TRUE)
  # 1: estimate gamma from rho
  gamma_1 <- rho / rho_mu
  
  cat(sprintf(
    "Median was computed for %d non-zero vectors with length < 3.\n", cnt
  ))
  
  list(
    X = X,
    beta = beta,
    beta_hat = beta_hat, # mean of beta
    beta_sigma2 = beta_sigma2,
    mu = mu, # mean of each (batch, class)
    ref_mu = ref_mu,
    rho = rho,
    gamma_1 = gamma_1,
    sample_batches = sample_batches,
    sample_classes = sample_classes
  )
}

correct_batch_effects <- function(
  obj,
  beta.var.threshold = 0.05,
  gamma.threshold = 1.7
) {
  median_nonzero <- function(x, ...) {
    if (length(x[x != 0]) < 3) {
      cnt <- cnt + 1
      return(NA)
    }
    median(x[x != 0], ...)
  }
  mean_nonzero <- function(x, ...) {
    if (length(x[x != 0]) < 3) {
      cnt <- cnt + 1
      return(NA)
    }
    mean(x[x != 0], ...)
  }
  X <- obj$X
  beta <- obj$beta
  beta_hat <- obj$beta_hat # mean of beta
  beta_sigma2 <- obj$beta_sigma2
  mu <- obj$mu # mean of each (batch, class)
  ref_mu <- obj$ref_mu
  rho <- obj$rho
  gamma_1 <- obj$gamma_1
  sample_batches <- obj$sample_batches
  sample_classes <- obj$sample_classes
  
  # TODO: beware of all rhos for genes is NA
  # Identifying outlier (batch, class, feature) that requires gamma correction
  # Identification using gamma_1 value only
  # gamma_1 cannot be zero or NA
  is_outlier <- gamma_1 != 0 &
    (gamma_1 > gamma.threshold | gamma_1 < (1 / gamma.threshold))
  outlier_indices <- which(is_outlier, arr.ind = TRUE)
  pair_classes <- c("A", "B") # assuming there only is one pair
  beta_1 <- beta[, , pair_classes[1]]
  beta_2 <- beta[, , pair_classes[2]]
  outlier_beta <- cbind(beta_1[outlier_indices], beta_2[outlier_indices])
  # class with beta closer to 1 is reference class
  outlier_class <- apply(abs(log(outlier_beta)), 1, which.max)
  # 2: estimate gamma from beta
  outliers <- cbind(
    outlier_indices,
    rho = rho[outlier_indices],
    beta_1 = beta_1[outlier_indices],
    beta_2 = beta_2[outlier_indices],
    gamma_1 = gamma_1[outlier_indices],
    gamma_2 = beta_1[outlier_indices] / beta_2[outlier_indices],
    outlier_class = outlier_class
  )
  # use gamma_1 as gamma_2 may be biased by possible errors in estimating ref in classes
  
  # Correcting gamma_1 -> Update X and beta
  combinations <- split.data.frame(
    outliers, list(outliers[, "col"], outliers[, "outlier_class"])
  )
  cat(sprintf("No. of gamma outliers = %d\n", nrow(outliers)))
  # iterating through all (batch, outlier_class) combinations
  for (outlier_kg in combinations) {
    if (dim(outlier_kg)[1] == 0) next
    k <- unique(sample_batches)[outlier_kg[1, "col"]]
    class_idx <- outlier_kg[1, "outlier_class"]
    g <- pair_classes[class_idx]
    # subset target patients
    sids <- colnames(X)[sample_batches == k & sample_classes == g]
    gamma_kg <- outlier_kg[, "gamma_1"]
    stopifnot(length(gamma_kg) == length(rownames(outlier_kg)))
    if (class_idx == 2) {
      # if outlier_class == 2 -> multiply by gamma_1
      X[rownames(outlier_kg), sids] <- X[rownames(outlier_kg), sids] * gamma_kg
    } else if (class_idx == 1) {
      # if outlier_class == 1 -> divide by gamma_1
      # gamma defined as "1"/"2"
      X[rownames(outlier_kg), sids] <- X[rownames(outlier_kg), sids] / gamma_kg
    } else {
      warning("Gamma was not corrected!")
    }
  }
  # TODO: Do we update mu to recalculate some betas, and then identify betas to correct?
  # TODO: Can directly scale mus

#   # Correction using beta_hat
#   X_batches <- split.default(X, batch)
#   if (sum(is.na(beta_hat)) > 0 | sum(beta_hat == 0, na.rm = TRUE) > 0)
#     stop("Beta_hat matrix contains zeros or NAs.")
#   for (k in colnames(beta_hat)) {
#     X_batches[[k]] <- X_batches[[k]] / beta_hat[, k]
#   }
#   X1 <- do.call(cbind, unname(X_batches))
#   X1 <- X1[, colnames(X)]

  list(X = X, outliers = outliers)
}
