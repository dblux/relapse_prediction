#' Calculates proportion of variance in dataframe due to batch effects
#' 
#' @param X dataframe with dim (n_features, n_samples) 
#' @param batch vector containing batch labels of samples (ordered the same way as in X)
#' @param class vector containing class labels of samples (ordered the same way as in X)
#' @param ret.obj logical indicating whether to return object or proportion of variance
#' @param pca logical indicating whether to perform PCA on X
#' @return numeric containing total proportion of variance in dataframe due to batch effects
rvp <- function(X, batch, class = NULL, ret.obj = FALSE) {
  X[is.na(X)] <- 0
  Z <- t(X)
  batch <- as.character(batch)
  if (is.list(class)) {
    print('list')
  } else if (is.vector(class) || is.factor(class)) {
    class <- as.character(class)
  }
  if (length(unique(batch)) == 1) {
    cat("All samples are from the same batch!\n")
    # Use NA as is.na works on lists
    return(list(percentage = 0, sum.squares = NA))
  }
  if (nrow(Z) != length(batch))
    stop("Length of batch vector does not match no. of samples in X")
  
  if (is.null(class)) {
    feature_means <- colMeans(Z)
    ss_total_class <- colSums(sweep(Z, 2, feature_means, `-`) ^ 2)
    Z_batches <- split.data.frame(Z, batch)
    batch_means <- sapply(Z_batches, function(Z) colMeans(Z))
    nperbatches <- sapply(Z_batches, nrow)
    squares <- (batch_means - feature_means) ^ 2
    ss_between_batch <- rowSums(sweep(squares, 2, nperbatches, `*`))

    stopifnot(length(ss_between_batch) == ncol(Z))
    pi_total <- sum(ss_between_batch) / sum(ss_total_class) 
    if (ret.obj) {
      return(list(
        percentage = pi_total,
        sum.squares = data.frame(
          ss_between = ss_between_batch,
          ss_total = ss_total_class
        )
      ))
    } else {
      return(pi_total)
    }
  } else {
    Z_classes <- split.data.frame(Z, class)
    Z_classes <- Filter(function(X) nrow(X) != 0, Z_classes)
    classes_string <- do.call(paste, as.list(names(Z_classes)))
    cat(sprintf("Split into classes: %s\n", classes_string))
    Z_t_classes <- lapply(Z_classes, t)
    batch_classes <- split(batch, class)
    batch_classes <- Filter(function(x) length(x) != 0, batch_classes)
    # Warning: Recursive call
    objs <- mapply(
      rvp, Z_t_classes, batch_classes,
      MoreArgs = list(class = NULL, ret.obj = TRUE),
      SIMPLIFY = FALSE
    )
    ss_total <- colSums(sweep(Z, 2, colMeans(Z), `-`) ^ 2)
    sumsquares_classes <- lapply(objs, function(obj) obj$sum.squares)
    # filters out obj$features == NA
    sumsquares_classes <- sumsquares_classes[!is.na(sumsquares_classes)]
    ss_batch_classes <- lapply(sumsquares_classes, function(X) X$ss_between)
    stopifnot(is.list(ss_batch_classes))
    ss_batch <- Reduce(`+`, ss_batch_classes)
    stopifnot(length(ss_batch) == length(ss_total))
    pi_total <- sum(ss_batch) / sum(ss_total)
    if (ret.obj) {
      return(list(
        percentage = pi_total,
        percentage_classes = sapply(objs, function(obj) obj$percentage),
        sum.squares = data.frame(ss_batch, ss_total)
      ))
    } else {
      return(pi_total)
    }
  }
}
