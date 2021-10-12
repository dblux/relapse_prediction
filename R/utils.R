#' Annotates affymetrix probesets to gene ID of choice
#'
#' Annotates affymetrix probesets according to platform and gene
#' naming convention provided in the annotation file. Option of
#' returning vector of annotaions with probesets as names
#'
#' Removes probesets with no matching ID. If multiple probesets
#' map to a single ID, the probeset with the max sum is used.
#' @param X data.frame with affy probesets as rownames
#' @param file name of the annotation file
#' @param ret.annot logical indicating whether to return vector of annotations
affy2id <- function(X, file, ret.annot = F) { 
  # probesets are rownames of dataframe
  annot_table <- read.table(
    file, sep = "\t", header = T, row.names = 1,
    stringsAsFactors = F, strip.white = T
  )
#   # Filters out ambiguous and AFFY probesets from annot 
#   fltr_annot <- annot_table[
#     grepl("[0-9]_at", rownames(annot_table)) 
#     & !startsWith(rownames(annot_table), "A"),
#     , drop = F
#   ]
  orig_rownames <- annot_table[rownames(X), ]
                      
  msg_no_id <- sprintf( 
    "No. of probesets with no ID removed: %d\n", sum(orig_rownames == "")
  )
  cat(msg_no_id)
   
  # Indices of probe sets with no corresponding ID to be deleted 
  idx_del <- which(orig_rownames == "")
   
  # Identifies genes that have multiple probesets mapping to it 
  freq_gene <- table(orig_rownames)
  dup_genes <- names(freq_gene[freq_gene > 1])
  dup_genes <- setdiff(dup_genes, "")
  for (gene in dup_genes) {
    # subset rows of dataframe with the same id
    same_rows <- X[orig_rownames == gene, , drop = F] 
    # assign numeric indices as rownames 
    rownames(same_rows) <- which(orig_rownames == gene)
    # rows that do not have the maximum sum are deleted 
    row_del <- as.numeric(
      rownames(same_rows)[-which.max(rowSums(same_rows))]
    )
    # concat with existing list of indices to be deleted 
    idx_del <- c(idx_del, row_del)
  }
  
  msg_total <- sprintf(
    "Total no. of probesets removed (incl. probesets mapping to same gene): %d\n", 
    length(idx_del)
  ) 
  cat(msg_total)
  
  # Rows are deleted 
  X_genes <- X[-idx_del, ]
  fltr_rownames <- orig_rownames[-idx_del]
  names(fltr_rownames) <- rownames(X)[-idx_del]
  # Assigning id to X
  rownames(X_genes) <- fltr_rownames
  
  if (ret.annot) {
    orig_rownames[idx_del] <- NA
    names(orig_rownames) <- rownames(X)
    return(orig_rownames)
  }
  
  X_genes
}


#' Sorts sample IDs so that they are paired
sort_sid <- function(x) {
  pid <- sapply(x, substring, 1, 4)
  time <- sapply(x, substring, 6, 7)
  time_pid <- mapply(paste0, time, pid)
  return(x[order(time_pid)])
}


#' Boolean function checking that dataframe provided has matching pair names
#' @param X dataframe with paired samples
#' @return logical indicating if it is paired or not
is_paired <- function(x) {
  #' @param character of sample ids
  #' @return logical indicating if it is paired or not
  .is_paired <- function(sid) {
    n <- length(sid)
    pid <- substring(sid, 1, 4)
    all(pid[1:(n / 2)] == pid[(n / 2 + 1):n])
  }
  
  if (is.data.frame(x) | is.matrix(x)) {
    return(.is_paired(colnames(x)))
  } else if (is.character(x)) {
    return(.is_paired(x))
  } else{
    stop("x is not a dataframe, matrix or character.")
  }
}

#' @param y_true numeric vector of true labels with 1 as positive
calc_recall <- function(y_true, y_pred) {
  sum(y_pred[y_true == 1] == 1) / sum(y_true == 1)
}


#' @param y_true numeric vector of true labels with 1 as positive
calc_specificity <- function(y_true, y_pred) {
  sum(y_pred[y_true == 0] == 0) / sum(y_true == 0)
}

# Substring without n last characters
substring_head <- function(string, n) substring(string, 1, nchar(string)-n)

# Substring tail starting from length - n characters
substring_tail <- function(string, n) substring(string, nchar(string)-n+1)