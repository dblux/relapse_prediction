# Arguments: combined df, batch, class and order info
# Order info: Left most is anchor data (Recursively corrected)
# Returns: Combined df of corrected data
norm.CBC <- function(df, batch_info, class_info, order_batch, correction_wpath = "dump/correction_vectors.tsv") {
  # Arguments: df1 is anchor batch, df2 is second batch
  # Returns: Combined df of df2 mapped to df1
  pairwise_correction <- function(df1, df2) {
    class_info1 <- metadata_df[colnames(df1), "class_factor"]
    class_info2 <- metadata_df[colnames(df2), "class_factor"]
    # Find classes that appear in both batches (Impt to sort!!!)
    intersect_classes <- sort(intersect(class_info1, class_info2))
    
    # Create list of df split by class for both batches
    # Only intersect_classes appear in list_class_df1_intersect
    print("intersect_classes"); print(intersect_classes)
    list_class_df1_intersect <- lapply(intersect_classes,
                                       function(class_id) df1[, class_info1 == class_id, drop = F])
    names(list_class_df1_intersect) <- intersect_classes
    # All classes in df2 is present in list
    # Sort to avoid error later on
    list_class_df2 <- lapply(sort(unique(class_info2)),
                             function(class_id) df2[, class_info2 == class_id, drop = F])
    names(list_class_df2) <- sort(unique(class_info2))
    
    # Calculate centroid for classes that appear in both batches for df1
    list_df1_centroid_intersect <- lapply(list_class_df1_intersect, rowMeans)
    # Calculate centroid for all classes for df2
    list_df2_centroid <- lapply(list_class_df2, rowMeans)
    list_df2_centroid_intersect <- list_df2_centroid[intersect_classes]
    
    # Calculate correction vectors for each class
    # Anchor: df1, hence df2_centroid - df1_centroid (order matters!!!)
    correction_vectors <- data.frame(mapply(`-`,
                                            list_df2_centroid_intersect,
                                            list_df1_centroid_intersect))
    
    # Identify classes in df2 that are not present in df1
    extra_classes <- setdiff(unique(class_info2), intersect_classes)
    print(paste("Extra classes:", extra_classes))
    print(paste("Extra classes:", length(extra_classes)))
    if(length(extra_classes) != 0) {
      # Arguments: Class string (one of the extra classes)
      # Returns df of correction vector for class
      # Correction vector is replicated from correction vector of nearest class
      calc_correction_vec <- function(class) {
        # # Logical vector selecting current class
        # class_logvec <- names(list_df2_centroid) == class
        # centroid_vec <- unlist(list_df2_centroid[class_logvec])
        # list_other_centroids <- list_df2_centroid[!class_logvec]
        centroid_vec <- list_df2_centroid[[class]]
        # Calculate euclidean distance of current centroid to intersection centroids
        distance_vec <- sapply(list_df2_centroid_intersect,
                               function(vec) l2_norm(vec-centroid_vec))
        nearest_class <- names(sort(distance_vec))[1]
        # Use correction vector of nearest class
        correction_df <- correction_vectors[,nearest_class, drop = F]
        colnames(correction_df) <- class
        return(correction_df)
      }
      # Returns list of single column df for all extra_classes
      list_correction <- lapply(extra_classes, calc_correction_vec)
      extra_correction_df <- do.call(cbind, list_correction)
      # Append to correction_vectors
      all_correction_df <- cbind(correction_vectors, extra_correction_df)
      correction_vectors <- all_correction_df[,sort(colnames(all_correction_df))]
    }
    # TODO: Give different filenames for different pairwise correction
    write.table(correction_vectors, correction_wpath,
                quote = F, sep = "\t")
    # Ensure that correction_vectors is ordered the same way as list_class_df2
    print(names(list_class_df2))
    print(colnames(correction_vectors))
    stopifnot(identical(names(list_class_df2), colnames(correction_vectors)))
    # Apply correction for each class in df2 (list_class_df2)
    list_corrected_df2 <- unname(mapply(function(df,vec) df-vec,
                                        list_class_df2, correction_vectors,
                                        SIMPLIFY = F))
    corrected_df2 <- do.call(cbind, list_corrected_df2)
    # Replace all negative values with 0
    corrected_df2[corrected_df2 < 0] <- 0
    combined_df <- cbind(df1, corrected_df2)
    return(combined_df)
  }
  
  col_order <- colnames(df)
  order_factor <- as.factor(order_batch)
  class_factor <- as.factor(class_info)
  batch_factor <- as.factor(batch_info)
  batch_levels <- levels(batch_factor)
  # Create metadata df
  metadata_df <- data.frame(batch_factor, class_factor)
  rownames(metadata_df) <- colnames(df)
  
  # Split df by batches into list
  list_batch_df <- lapply(batch_levels, function(batch_id) df[, batch_factor == batch_id])
  names(list_batch_df) <- batch_levels
  # Order list of split dfs
  ordered_list_batch_df <- list_batch_df[order_factor]
  str(ordered_list_batch_df)
  
  # Fold left on ordered_list_batch_df
  # Sorts columns according to initial df order
  return(Reduce(pairwise_correction, ordered_list_batch_df)[,col_order])
  # TODO: Check if all batches in batch info are provided in the order (IDENTICAL)
  # Print if less batches in order, the batch not included will not be shown!
  
  # Instead of using centroids...
  # Within the class, calculate a nearest neighbour for each sample (both ways)
  
  # TODO: How to determine the order!!!
  # ERROR IF A BATCH HAS ONLY ONE CLASS I THINK!!!
  # ERROR: CLASS HAS TO BE GIVEN IN ALPHABETS!!!
}
