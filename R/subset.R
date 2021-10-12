#' Subset columns according to features in annotation
#'
#' Colnames of X has to match with rownames of annot.
#'
#' @example subset_cols(X, annot, subtype == 'BCR-ABL' & class_info == 'D0')
subset_cols <- function(data, annot, ...) {
  X[, rownames(subset(annot[colnames(X), ], ...))]
}


#' Remove rows of dataframe which meet the condition
#'
#' Special syntax: 'row' indicates row of dataframe `X`
#'
#' @examples
#' remove_rows(X, sum(row) == 0)
#' remove_rows(X, var(row) == 0)
remove_rows <- function(X, condition) {
  idx <- vector()
  for (i in seq_len(nrow(X))) {
    row <- as.numeric(X[i, ])
    # eval env defaults to calling env
    # i.e. row is automatically used for expression
    is_cond <- eval(substitute(condition))
    idx <- c(idx, is_cond)
  }
  
  X[!idx, ]
}


#' Select top n highly variable genes
#' @param X matrix or dataframe with rownames
select_hvg <- function(X, n, ret.genes = FALSE) {
  if (is.null(rownames(X)))
    stop("X does not have rownames")
  
  row_var <- apply(X, 1, var)
  sorted_var <- sort(row_var, decreasing = TRUE)
  hvg <- names(sorted_var)[1:n]
  
  if (ret.genes)
    return(hvg)
  
  X[hvg, ]
}


# Filter probes with too many zeros
#' @param df dataframe
#' @param percent_threshold percentage threshold of non-zero values
#' @param metadata_df df containing class labels of samples
#' @param logical_func a function that is either "all" or "any". (Either all or
#' just one class have to pass the threshold)
#' @return dataframe containing rows that meet threshold of non-zero
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
    class_factor <- metadata_df[colnames(df1),"class_info"]
    print(head(class_factor))
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