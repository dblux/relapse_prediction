#' Calculates percentage of variance in data due to batch effects
#' 
#' @param X dataframe with dim (n_samples, n_features)
#' @param batch vector containing batch labels of samples
#' @param class vector or list of vectors containing class labels of samples
#' @param ret.obj logical indicating whether to return object or percentage of variance
#' @return numeric containing total percentage of variance in data due to batch effects
RVP <- function(X, batch, class = NULL, ret.obj = FALSE) {
  X[is.na(X)] <- 0
  batch <- as.character(batch)
  
  if (is.vector(class) || is.factor(class))
    class <- as.character(class)
  if (length(unique(batch)) == 1) {
    message("All samples are from the same batch!")
    # use NA as is.na works on lists
    return(list(percentage = 0, sum_squares = NA))
  }
  if (nrow(X) != length(batch))
    stop("Length of batch vector does not match no. of samples in X")

  # COMPUTE RVP 
  if (is.null(class)) {
    feature_means <- colMeans(X)
    ss_total <- colSums(sweep(X, 2, feature_means, `-`) ^ 2)
    X_batches <- split.data.frame(X, batch)
    # rm(X)
    batch_means <- sapply(X_batches, function(X) colMeans(X))
    nperbatches <- sapply(X_batches, nrow)
    squares <- (batch_means - feature_means) ^ 2
    ss_batch <- rowSums(sweep(squares, 2, nperbatches, `*`))

    stopifnot(length(ss_batch) == ncol(X))
    total_percentage <- sum(ss_batch) / sum(ss_total) 
    if (ret.obj) {
      return(list(
        percentage = total_percentage,
        sum_squares = data.frame(ss_batch, ss_total)
      ))
    } else {
      return(total_percentage)
    }
  } else {
    ss_total <- colSums(sweep(X, 2, colMeans(X), `-`) ^ 2)
    X_classes <- split.data.frame(X, class)
    # rm(X)
    X_classes <- Filter(function(X) nrow(X) != 0, X_classes)
    classes_string <- do.call(paste, as.list(names(X_classes)))
    # message(sprintf("Split into classes: %s", classes_string))
    batch_classes <- split(batch, class)
    batch_classes <- Filter(function(x) length(x) != 0, batch_classes)
    # warning: recursive call
    objs <- mapply(
      RVP, X_classes, batch_classes,
      MoreArgs = list(class = NULL, ret.obj = TRUE),
      SIMPLIFY = FALSE
    )
    sumsquares_classes <- lapply(objs, function(obj) obj$sum_squares)
    # filters out obj$sum_squares == NA
    sumsquares_classes <- sumsquares_classes[!is.na(sumsquares_classes)]
    ss_batch_classes <- lapply(sumsquares_classes, function(X) X$ss_batch)
    stopifnot(is.list(ss_batch_classes))
    ss_batch <- Reduce(`+`, ss_batch_classes)
    stopifnot(length(ss_batch) == length(ss_total))
    total_percentage <- sum(ss_batch) / sum(ss_total)
    if (ret.obj) {
      return(list(
        percentage = total_percentage,
        percentage_classes = sapply(objs, function(obj) obj$percentage),
        sum_squares = data.frame(ss_batch, ss_total)
      ))
    } else {
      return(total_percentage)
    }
  }
}
