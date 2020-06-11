# Arguments: combined df, batch, class and order info
# Order info: Left most is anchor data (Recursively corrected)
# Returns: Combined df of corrected data
correctBCM <- function(df, batch_info, class_info, order_batch,
                       correction_wpath = "dump/correction_vectors.tsv") {
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


correctRefBCM <- function(df1, metadata_df, ref_batch = 1,
                          correction_wpath = "dump/correction_vectors.tsv") {
  # Obtaining batch and class annotations
  batch_factor <- as.factor(metadata_df[colnames(df1),"batch_info"])
  class_factor <- metadata_df[colnames(df1),"class_info"]
  batch_levels <- levels(batch_factor)
  
  # Split df by batches into list
  list_batch_df <- split.default(df1, batch_factor)
  # Define first reference batch and list of other df
  ref_df <- list_batch_df[[ref_batch]]
  list_df <- list_batch_df[-ref_batch]
  
  # # Recursive function that returns ref_df
  # # Arguments: ref_df is ref batch, list_df is list of other batches
  # # Returns: list(ref_df, list_df)
  recursive_correction <- function(ref_df, list_df) {
    if (length(list_df) == 0) return(ref_df)
    else {
      # BATCH EFFECT CORRECTION
      ref_class_info <- as.character(metadata_df[colnames(ref_df),
                                                 "class_info"])
      ref_levels <- unique(ref_class_info)
      
      list_other_classes <- lapply(
        list_df, function(nonref_df) metadata_df[colnames(nonref_df),
                                                 "class_info"]
      )
      
      # Logical vector selecting batches that have at least one intersecting
      # class
      is_intersect <- sapply(list_other_classes,
                             function(vec) any(unique(vec) %in% ref_levels))
      list_intersect_batches <- list_df[is_intersect]
      # To be returned
      list_remaining <- list_df[!is_intersect]
      
      # CALCULATE CENTROIDS FOR EACH CLASS IN REF_DF
      # Split into list of class df
      list_ref_df <- split.default(ref_df, ref_class_info)
      # Order list according to name
      list_ref_df <- list_ref_df[order(names(list_ref_df))]
      print(str(ref_df))
      print(str(list_ref_df))
      list_ref_centroids <- lapply(list_ref_df, apply, 1, mean, trim = 0.2)

      # Corrects to ref_df
      # Argument: Uncorrected df with intersecting class
      # Return: Corrected df of a single batch
      # Global: Variables belonging to ref_batch
      single_correction <- function(other_df) {
        other_class_info <- as.character(metadata_df[colnames(other_df),
                                                     "class_info"])
        print("other_class_info"); print(other_class_info)
        other_levels <- unique(other_class_info)
        # Find classes that appear in both batches (Impt to sort!!!)
        intersect_classes <- sort(intersect(ref_levels, other_levels))
        
        # Split other_df into classes
        list_other_df <- split.default(other_df, other_class_info)
        # Order list according to name
        list_other_df <- list_other_df[order(names(list_other_df))]
        print(str(other_df))
        print(str(list_other_df))
        list_other_centroids <- lapply(list_other_df,
                                       apply, 1, mean, trim = 0.2)
        str(list_other_centroids)
        
        # Calculate correction vectors for intersect classes
        # Anchor: df1, hence df2_centroid - df1_centroid (order matters!!!)
        calc.correction <-  function(class_char) {
          list_other_centroids[[class_char]] - list_ref_centroids[[class_char]]
        }
          # Convert to df and name according to classes
        correction_vectors <- data.frame(lapply(intersect_classes,
                                                calc.correction))
        colnames(correction_vectors) <- intersect_classes
        
        # Identify classes in other_df that are not present in intersect_classes
        extra_classes <- setdiff(other_levels, intersect_classes)
        print(paste("Extra classes:", extra_classes))
        print(paste("Extra classes:", length(extra_classes)))
        if(length(extra_classes) != 0) {
          # Arguments: Class string (one of the extra classes)
          # Returns df of correction vector for class
          # Correction vector is replicated from correction vector of nearest class
          calc.correction_vec <- function(class) {
            # # Logical vector selecting current class
            # class_logvec <- names(list_df2_centroid) == class
            # centroid_vec <- unlist(list_df2_centroid[class_logvec])
            # list_other_centroids <- list_df2_centroid[!class_logvec]
            centroid_vec <- list_other_centroids[[class]]
            # Calculate euclidean distance of current centroid to intersection centroids
            distance_vec <- sapply(list_other_centroids[intersect_classes],
                                   function(vec) calcL2Norm(vec-centroid_vec))
            nearest_class <- names(sort(distance_vec))[1]
            # Use correction vector of nearest class
            correction_df <- correction_vectors[,nearest_class, drop = F]
            colnames(correction_df) <- class
            return(correction_df)
          }
          # Returns list of single column df for all extra_classes
          list_correction <- lapply(extra_classes, calc.correction_vec)
          extra_correction_df <- do.call(cbind, list_correction)
          # Append to correction_vectors
          correction_vectors <- cbind(correction_vectors, extra_correction_df)
        }
        # Sort correction vectors
        correction_vectors <- correction_vectors[,order(colnames(correction_vectors))]
        # TODO: Give different filenames for different pairwise correction
        write.table(correction_vectors, correction_wpath,
                    quote = F, sep = "\t")
        # Ensure that correction_vectors is ordered the same way as list_class_df2
        print(names(list_other_df))
        print(colnames(correction_vectors))
        stopifnot(identical(names(list_other_df), colnames(correction_vectors)))
        # Apply correction for each class in df2 (list_class_df2)
        list_other_corrected <- unname(mapply(function(df, vec) df-vec,
                                              list_other_df, correction_vectors,
                                              SIMPLIFY = F))
        corrected_df <- do.call(cbind, list_other_corrected)
        # Replace all negative values with 0
        corrected_df[corrected_df < 0] <- 0
        return(corrected_df)
      }
      list_corrected_df <- unname(lapply(list_intersect_batches,
                                         single_correction))
      # Corrected df contains all samples
      final_df <- cbind(ref_df, do.call(cbind, list_corrected_df))
      # Replace all negative values with 0
      final_df[final_df < 0] <- 0
      # Recursive call
      recursive_correction(final_df, list_remaining)
    }
  }
  
  final_df <- recursive_correction(ref_df, list_df)
  # Sort the corrected df
  return(final_df[,colnames(df1)])
}

# Global BCM (Using D0 correction vector)
# Batch 2 hard coded to be ref batch
correctGlobalBCM <- function(df1, metadata_df) {
  df1_batch <- metadata_df[colnames(df1), "batch_info"]
  # TODO: Error if df1 is array instead of df
  stopifnot(class(df1) == "data.frame")
  # Split df by batches into list
  list_batch_df <- split.default(df1, df1_batch)
  # Subset only D0 patients in each batch
  list_subset_d0 <- lapply(list_batch_df,
                           function(df1) df1[,endsWith(colnames(df1), "D0"),
                                             drop = F])
  print(str(list_subset_d0))
  # TODO: Throw error when any member of list has an empty dataframe!!!
  # Calculate trimmed mean centroids
  list_d0_centroids <- lapply(list_subset_d0,
                              function(df1) apply(df1, 1, mean, trim = 0.2))
  # Use batch 2 as reference batch
  list_correction_vec <- lapply(list_d0_centroids,
                                function(vec) vec - list_d0_centroids[[2]])
  stopifnot(identical(names(list_batch_df), names(list_correction_vec)))
  list_corrected_df <- mapply(function(df1, vec) df1 - vec,
                              list_batch_df,
                              list_correction_vec)
  corrected_df <- do.call(cbind, unname(list_corrected_df))
  corrected_df[corrected_df < 0] <- 0
  # Re-ordering of columns
  bcm_df <- corrected_df[,colnames(df1)]
  return(bcm_df)
}

# All columns of the metadata have to be factors
# Have to have column names "batch" and "class"
### PAIRWISE PCA ###
correctSVDBCM <- function(df1, metadata_df, ref_batch) {
  # TODO: DELETE
  df1 <- subset_maqc
  metadata_df
  ref_batch <- 1
  
  all_batch_info <- metadata_df[colnames(df1), "batch_info"]
  # Split df by batches into list
  list_batch_df <- split.default(df1, all_batch_info)
  # Initialise empty lists
  list_plot <- list()
  list_corrected_plot <- list()
  list_corrected_df <- list()
  
  j <- 1
  ref_df <- list_batch_df[[ref_batch]]
  # Convert factor into numeric before unique
  batch_ids <- unique(as.numeric(as.character(all_batch_info)))
  nonref_ids <- batch_ids[batch_ids != ref_batch]
  print(nonref_ids)
  for (i in nonref_ids) {
    # i <- 1 # TODO
    other_df <- list_batch_df[[i]]
    pair_batch <- cbind(ref_df, other_df)
    pair_prcomp <- prcomp(t(pair_batch), scale. = F)
    # Subset samples from ref_df
    ref_pca <- pair_prcomp$x[1:ncol(ref_df),]
    other_pca <- pair_prcomp$x[-(1:ncol(ref_df)),]
    
    # Linear regression to get gradient
    ref_coef <- coef(lm(ref_pca[,2] ~ ref_pca[,1]))
    other_coef <- coef(lm(other_pca[,2] ~ other_pca[,1]))
    mean_gradient <- mean(c(ref_coef[2], other_coef[2]))
    # Biological vector in PCA space (unit norm)
    bio_vec_pca <- c(1, mean_gradient) / calcL2Norm(c(1, mean_gradient))
    rotation_90 <- calcRotationMatrix(pi/2)
    # Batch effect vector in PCA space (unit norm)
    batch_vec_pca <- rotation_90 %*% bio_vec_pca
    print(mean_gradient)
    print(as.vector(batch_vec_pca))
    
    ## Plotting parameters
    plot_title <- sprintf("B2 vs. B%d", i)
    eigenvalues <- (pair_prcomp$sdev)^2
    var_pct <- eigenvalues[1:5]/sum(eigenvalues)
    pc_labels <- sprintf("PC%d (%.2f%%)", 1:5, var_pct*100)
    
    batch_pairwise <- metadata_df[colnames(pair_batch), "batch_info"]
    
    class_pairwise <- metadata_df[colnames(pair_batch), "class_info"]
    all_pch <- 21:25
    class_pch <- all_pch[class_pairwise]
    
    # Plot BEFORE CORRECTION
    list_plot[[j]] <- ggplot(data.frame(rbind(ref_pca, other_pca)),
                             aes(x = PC1, y = PC2)) +
      geom_point(aes(fill = batch_pairwise), pch = class_pch,
                 size = 3, show.legend = F) + 
      geom_vline(xintercept = 0, color = "black", alpha = 0.5) +
      geom_hline(yintercept = 0, color = "black", alpha = 0.5) +
      geom_abline(slope = ref_coef[2], intercept = ref_coef[1],
                  color = "blue", alpha = 0.5) +
      geom_abline(slope = other_coef[2], intercept = other_coef[1],
                  color = "blue", alpha = 0.5) +
      geom_abline(slope = mean_gradient,
                  color = "orange", alpha = 0.5) +
      geom_abline(slope = batch_vec_pca[2]/batch_vec_pca[1],
                  color = "orange", alpha = 0.5) +
      labs(x = pc_labels[1], y = pc_labels[2], title = plot_title) +
      theme(plot.title = element_text(hjust = 0.5)) +
      coord_fixed(ratio = 1)
    
    # BCM: In PC1 & PC2 & PC3 subspace
    # other_ref_vec <- colMeans(ref_pca[,1:2]) - colMeans(other_pca[,1:2]) # TODO
    # TODO: Identify intersecting classes
    ref_class <- metadata_df[rownames(ref_pca), "class_info"]
    list_ref_pca <- split.data.frame(ref_pca[,1:3], ref_class)
    list_ref_centroid <- lapply(list_ref_pca, apply, 2, mean, trim = 0.2)
    
    other_class <- metadata_df[rownames(other_pca), "class_info"]
    list_other_pca <- split.data.frame(other_pca[,1:3], other_class)
    list_other_centroid <- lapply(list_other_pca, apply, 2, mean, trim = 0.2)
    
    # Correction vectors for each intersected class in other batch
    list_correction_vec <- mapply(`-`, list_ref_centroid, list_other_centroid,
                                  SIMPLIFY = FALSE)
    
    # Adding correction vectors to other batch
    list_corrected_other <- mapply(sweep,
                                   x = list_other_pca,
                                   STATS = list_correction_vec,
                                   MoreArgs = list(MARGIN = 2, FUN = "+"),
                                   SIMPLIFY = FALSE)
    
    corrected_df <- do.call(rbind, list_corrected_other)[rownames(other_pca),]
    other_pca[,1:3] <- corrected_df
    
    # Edit: Trimmed mean BCM
    # ref_centroid <- apply(ref_pca[,1:2], 2, mean, trim = 0.2)
    # other_centroid <- apply(other_pca[,1:2], 2, mean, trim = 0.2)
    # other_ref_vec <- ref_centroid - other_centroid
    # other_pca[,1:2] <- sweep(other_pca[,1:2], 2, other_ref_vec, `+`)
    
    # Transform corrected non-reference batch to original space
    rotated_other <- other_pca %*% t(pair_prcomp$rotation)
    # Return corrected non-reference df
    list_corrected_df[[j]] <- t(sweep(rotated_other, 2,
                                      pair_prcomp$center, "+"))
    
    # Plot single pairwise
    list_corrected_plot[[j]] <- ggplot(data.frame(rbind(ref_pca, other_pca)),
                                       aes(x = PC1, y = PC2)) +
      geom_point(aes(fill = batch_pairwise), pch = class_pch,
                 size = 3, show.legend = F) + 
      geom_vline(xintercept = 0, color = "black", alpha = 0.5) +
      geom_hline(yintercept = 0, color = "black", alpha = 0.5) +
      geom_abline(slope = ref_coef[2], intercept = ref_coef[1],
                  color = "blue", alpha = 0.5) +
      geom_abline(slope = other_coef[2], intercept = other_coef[1],
                  color = "blue", alpha = 0.5) +
      geom_abline(slope = mean_gradient,
                  color = "orange", alpha = 0.5) +
      geom_abline(slope = batch_vec_pca[2]/batch_vec_pca[1],
                  color = "orange", alpha = 0.5) +
      labs(x = pc_labels[1], y = pc_labels[2], title = plot_title) +
      theme(plot.title = element_text(hjust = 0.5)) +
      coord_fixed(ratio = 1)
    
    j <- j + 1
  }
  
  # Combine reference and corrected non-references
  # Reorder columns according to initial df
  corrected_nonref_df <- do.call(cbind, list_corrected_df)
  corrected_df <- cbind(ref_df, corrected_nonref_df)[,colnames(df1)]
  return(list(data = corrected_df, plot = list_plot,
              corrected_plot = list_corrected_plot))
}


# SANDBOX -----------------------------------------------------------------
source("../functions.R")

## Subset of original data
# Removed outliers, patients with timepoints from different batches and batch 5
SUBSET_RPATH <- "data/GSE67684/processed/subset_yeoh.tsv"
subset_yeoh <- read.table(SUBSET_RPATH, sep = "\t")

## Metadata
# Preprocessed metadata
METADATA_RPATH <- "data/GSE67684/processed/metadata/full_metadata.tsv"
metadata_df <- read.table(METADATA_RPATH, sep = "\t")

# BATCH_RPATH <- "data/GSE67684/processed/metadata/metadata-batch.tsv"
# LABEL_RPATH <- "data/GSE67684/processed/metadata/metadata-label_mrd_subtype.tsv"
# yeoh_batch <- read.table(BATCH_RPATH, sep = "\t", header = T, row.names = 1)
# yeoh_label <- read.table(LABEL_RPATH, sep = "\t", header = T, row.names = 1)

# SCALE->REMOVE->FILTER->LOG
scaled_yeoh <- normaliseMeanScaling(subset_yeoh)
subset_yeoh <- removeProbesets(scaled_yeoh)
data_yeoh <- log2_transform(filterProbesets(subset_yeoh, 0.7, metadata_df))

# table(metadata_df$batch_info, metadata_df$subtype)

library(Rgraphviz)
library(igraph)

#' @param X dataframe of gene expression data
#' @param batch vector indicating batch
#' @param metadata dataframe with columns encoding relevant class information
correctMultiBCM <- function(X, batch, metadata, ref_batch, size = 3,
                            n_outlier = 2)
{
  ## Assert no negative values
  if (sum(X < 0) > 0) stop("Negative values present in X..")
  
  # Coercion of argument datatypes
  batch <- as.character(batch)
  ref_batch <- as.character(ref_batch)
  complete_metadata <- data.frame(batch,
                                  lapply(metadata, as.character),
                                  row.names = rownames(metadata),
                                  stringsAsFactors = F)
  metadata <- complete_metadata[,2:ncol(complete_metadata)]
  
  ## Check that batches share at least one category (graph is fully connected)
  ## Create G = (V: Batches, E: Common categories)
  # Split according to all possible combinations of categories
  list_combi <- split(complete_metadata, metadata,
                      drop = T, sep = "_")
  # No. of samples in each batch for each combination
  list_combi_table <- lapply(list_combi, function(X) table(X$batch))
  # Filter batches with samples > size!
  list_combi_batches <- lapply(list_combi_table,
                               function(x) names(x)[x > size])
  nbatches_combi <- sapply(list_combi_batches, length)
  # All pairwise combinations between batches in each combination
  list_pairwise <- lapply(list_combi_batches[nbatches_combi > 1], # More than 1 batch
                          function(x) t(combn(x, 2)))
  list_edges <- do.call(rbind, list_pairwise) # Edge list for graph
  n_pairwise <- sapply(list_pairwise, nrow)
  edge_names <- rep(names(n_pairwise), n_pairwise)
  # Assumption: Metadata does not contain "_"
  edge_cols <- do.call(rbind, strsplit(edge_names, "_"))
  edge_list <- data.frame(cbind(list_edges, edge_cols)) # With edge attributes
  colnames(edge_list) <- c("from", "to", colnames(metadata))
  batch_df <- data.frame(batch = as.character(sort(unique(batch))))
  G <- graph_from_data_frame(edge_list, directed = F,
                             vertices = batch_df)
  simple_G <- simplify(G)
  plot(simple_G)
  plot(G)
  stopifnot(is_connected(G)) # There is a batch which has no common class
  # TODO: Decide reference batch (Simplify graph and one with highest degree?
  # Also need to consider the magnitude of correction?)
  
  ## Split df by batches into list
  list_batches <- split.default(X, batch)
  # Define first reference batch and list of other df
  X_ref <- list_batches[[ref_batch]]
  list_X_other <- list_batches[names(list_batches) != ref_batch]
  
  #' Recursive function that selects batches with shared class with 
  #' reference batch and corrects accordingly returning corrected data
  #' @param ref_df dataframe of reference batch
  #' @param list_df list of dataframes of all other batches aside
  #' from reference batch
  #' @return batch corrected data
  #' Global metadata, ref_batch
  recursive_correction <- function(ref_df, list_df, graph) {
    if (length(list_df) == 0) return(ref_df) # base case
    else {
      # Filter batches that have at least one intersecting class
      adjacent_batches <- names(neighbors(graph, ref_batch))
      list_intersect_batches <- list_df[adjacent_batches]
      # USED IN RECURSIVE CALL!
      list_remaining <- list_df[!names(list_df) %in% adjacent_batches]
      
      # Split ref_df into all possible combinations
      list_ref <- split.default(ref_df, metadata[colnames(ref_df),],
                                drop = T, sep = "_")
      # Preliminary filtering based on size
      list_nsamples <- sapply(list_ref, ncol)
      sublist_ref <- list_ref[list_nsamples > size]
      
      #' Converts gene expression df to contingency df with features:
      #' mean_positive, n_positive, n_zero
      tierExpression <- function(X) {
        tierFeature <- function(x) {
          x_positive <- x[x > 0]
          mean_positive <- mean(x_positive) # TODO: Median?
          n_positive <- length(x_positive)
          # ASSUMPTION: No negative values
          c(mean_positive, n_positive, length(x)-n_positive)
        }
        features <- data.frame(t(apply(X, 1, tierFeature)))
        colnames(features) <- c("mean_positive", "n_positive", "n_zero")
        return(features)
      }
      
      sublist_ref_tiered <- lapply(sublist_ref, tierExpression)
      
      #' Pairwise batch correction with reference batch
      #' @param batch_name name of batch to be corrected
      #' @param list_batches list of df from different batches
      #' @return corrected df
      # Global: Variables belonging to ref_batch, G, metadata
      single_correction <- function(batch_name, list_batches) {
        print(sprintf("Batch %s", batch_name))
        X <- list_batches[[batch_name]]
        # Split by combinations
        list_combi <- split.default(X, metadata[colnames(X),],
                                    drop = T, sep = "_")
        # Filter intersecting classes with reference batch
        intersecting_edges <- E(G)[ref_batch %--% batch_name]
        all_edges_attr <- edge_attr(G)
        intersecting_edges_attr <-
          lapply(all_edges_attr,
                 function(attr, edge_idx) attr[as.numeric(edge_idx)],
                 intersecting_edges)
        intersecting_combi <- do.call(paste,
                                      c(intersecting_edges_attr, sep = "_"))
        sublist_combi <- list_combi[intersecting_combi]
        # Transform
        sublist_combi_tiered <- lapply(sublist_combi, tierExpression)

        #' Calculate correction vectors for a single intersecting combination
        #' Global: sublist_ref_tiered
        calcCorrection <- function(combi_name, other_tiered, ref_tiered) {
          X_ref <- ref_tiered[[combi_name]]
          X_other <- other_tiered[[combi_name]]
          
          # Initialise correction vector
          correction_vector <- rep(0, nrow(X_ref))
          names(correction_vector) <- rownames(X_ref)
          
          ## Calculate zero - non-zero correction
          ref_allzero <- X_ref$n_positive == 0 # No positive samples
          ref_allpositive <- X_ref$n_zero == 0 # Non zero samples
          other_allzero <- X_other$n_positive == 0 # No positive samples
          other_allpositive <- X_other$n_zero == 0 # Non zero samples
          
          correction_vector[ref_allzero & other_allpositive] <- 
            -X_other$mean_positive[ref_allzero & other_allpositive]
          correction_vector[other_allzero & ref_allpositive] <- 
            X_ref$mean_positive[other_allzero & ref_allpositive]
          big_correction <-
            (other_allzero & ref_allpositive)|(ref_allzero & other_allpositive)
          print(sprintf("No. of features with big correction = %d",
                        sum(big_correction)))
          
          ## Calculate non-zero - non-zero correction
          ref_positive <- X_ref$n_positive > n_outlier
          other_positive <- X_other$n_positive > n_outlier
          logi_idx <- ref_positive & other_positive
          correction_vector[logi_idx] <-
            X_other$mean_positive[logi_idx] - X_ref$mean_positive[logi_idx]
          
          # Check if both logi_idx overlap
          stopifnot(sum(which(big_correction) %in% which(logi_idx)) == 0)
          return(correction_vector)
        }
        # Calculates correction vectors for all combinations
        correction_df <- data.frame(
          sapply(names(sublist_combi_tiered),
                 calcCorrection,
                 sublist_combi_tiered,
                 sublist_ref_tiered)
        )
        CORRECTION_WPATH <- sprintf("temp/correction-%s.tsv", batch_name)
        write.table(correction_df, CORRECTION_WPATH, sep = "\t")
        
        correction_vector <- apply(correction_df, 1, median)
        X_corrected <- X - correction_vector
        # Replace all negative values with 0
        X_corrected[X_corrected < 0] <- 0
        return(X_corrected)
      }
      
      list_X_corrected <- lapply(names(list_intersect_batches),
                                 single_correction,
                                 list_intersect_batches)
      # Corrected df contains all samples
      X_concat <- cbind(ref_df, do.call(cbind, list_X_corrected))
      # Replace all negative values with 0
      stopifnot(sum(X_concat < 0) == 0)
       
      ## Merge vertices of corrected batches in graph
      V_adj_idx <- as.numeric(neighbors(graph, ref_batch))
      V_ref_idx <- as.numeric(V(graph)[ref_batch])
      V_mapping <- 1:vcount(graph) # Initialise mapping
      V_mapping[V_adj_idx] <- V_ref_idx
      print(V_mapping)
      graph <- contract(graph, V_mapping)
      V_names <- V(graph)$name
      nbatch_names <- sapply(V_names, length)
      # Rename merged vertex
      V(graph)[nbatch_names > 1]$name <- ref_batch
      # Delete empty vertices
      graph <- delete.vertices(graph, which(nbatch_names == 0))
      graph <- simplify(graph)
      
      # Recursive call
      recursive_correction(X_concat, list_remaining, graph)
    }
  }
  
  X_final <- recursive_correction(X_ref, list_X_other, simple_G)
  return(X_final[,colnames(X)]) # Sort columns according to X
}

batch <- metadata_df[colnames(data_yeoh), "batch_info"]
metadata <- metadata_df[colnames(data_yeoh), 2:3]
data_bcm <- correctMultiBCM(data_yeoh, batch, metadata, ref_batch = "2")

correction_1 <- read.table("temp/correction-1.tsv", sep = "\t")
head(correction_1, 30)
hist(correction_1[,4], breaks = 30)
logi_df <- abs(correction_1) > 4
idx <- rownames(logi_df)[apply(logi_df, 1, any)]
correction_1[idx,]

