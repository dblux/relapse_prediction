#' Plots PCA plot using ggplot2
#'
#' @import ggplot2
#' @param X dataframe with features as rows and samples as columns
#' @param annot dataframe of annotation with samples as rows
#' @param newdata dataframe of data to be predicted by prcomp object
#' @param ... optional arguments are passed to aes_string in ggplot
ggplot_pca <- function(X, annot, cex = 2, newdata = NULL, ...) {
  # PCA
  pca_obj <- prcomp(t(X))
  Z <- data.frame(pca_obj$x[, 1:3])
  eigenvalues <- (pca_obj$sdev)^2
  var_pc <- eigenvalues[1:3]/sum(eigenvalues)
  pc_labels <- sprintf("PC%d (%.2f%%)", 1:3, var_pc*100)

  # Projects newdata into PCA space
  if (!is.null(newdata)) {
    Z_new <- predict(pca_obj, newdata = t(newdata))[, 1:3]
    # Remove duplicate rows
    Z_new <- Z_new[!(rownames(Z_new) %in% rownames(Z)), ]
    Z <- rbind(Z, Z_new)
  }

  # concat with annotations
  annot_cols <- unlist(list(...))
  annot1 <- annot[rownames(Z), annot_cols, drop = F]
  Z_annot <- cbind(Z, annot1)

  ggplot(Z_annot, aes_string(x = "PC1", y = "PC2", ...)) +
    geom_point(cex = cex) +
    labs(x = pc_labels[1], y = pc_labels[2]) +
    geom_vline(xintercept = 0, color = "black", alpha = 0.5) +
    geom_hline(yintercept = 0, color = "black", alpha = 0.5)
}


#' Plots UMAP plot using ggplot2
#'
#' @import ggplot2
#' @param obj UMAP object from umap function
#' @param ... optional arguments are passed to aes_string in ggplot
ggplot_umap <- function(obj, annot, cex = 2, alpha = 1, ...) {
  Z <- data.frame(obj$layout)
  colnames(Z) <- c("UMAP1", "UMAP2")
  annot_cols <- unlist(list(...))
  annot1 <- annot[rownames(Z), annot_cols, drop = F]
  Z_annot <- cbind(Z, annot1)
  
  ggplot(Z_annot, aes_string(x = "UMAP1", y = "UMAP2" , ...)) +
    geom_point(cex = cex, alpha = alpha)
}


#' Plots ROC curve using ggplot2
#' @import ggplot2
#' @param X dataframe containing scores / predictions and labels
#' @param response character containing name of column in X with labels
#' @param predictor character containing name of column/s in X of predictors
#' @param pauc.limits numeric vector of length 2 indicating limits
ggplot_roc <- function(
  X, response, predictor,
  pauc.limits = FALSE,
  pauc.axis = c("specificity", "sensitivity"),
  pauc.correct = TRUE,
  return.auc = FALSE,
  lwd = 1
) {
  if (!is.logical(pauc.limits)) {
    pauc.limits <- sort(pauc.limits, decreasing = TRUE)
  }
  
  if (length(predictor) == 1) {
    str_formula <- paste(response, predictor, sep = " ~ ")
  } else if (length(predictor) > 1) {
    predictors <- do.call(paste, c(as.list(predictor), sep = " + "))
    str_formula <- paste(response, predictors, sep = " ~ ")
  } else {
    stop("arg predictor is of non-positive length.")
  }
  
  roc_objs <- roc(
    as.formula(str_formula),
    data = X,
    partial.auc = pauc.limits,
    partial.auc.focus = match.arg(pauc.axis),
    partial.auc.correct = pauc.correct
  )

  if (is(roc_objs, "roc")) {
    aucs <- roc_objs$auc
    d <- data.frame(
      FPR = 1 - roc_objs$specificities,
      TPR = roc_objs$sensitivities
    )
    d <- d[nrow(d):1, ]
    d <- cbind(names = predictor, d)
  } else {
    aucs <- sapply(roc_objs, function(obj) obj$auc)
    list_d <- lapply(
      roc_objs,
      function(obj) data.frame(
        FPR = 1 - obj$specificities,
        TPR = obj$sensitivities
      )
    )
    list_d <- lapply(list_d, function(d) d[nrow(d):1, ])
    d <- do.call(rbind, list_d)
    n_rows <- sapply(list_d, nrow)
    predictors_col <- rep(predictor, n_rows)
    d <- cbind(names = predictors_col, d)
  }
            
  # label colors are assigned by alphabetical order
  labels <- sort(sprintf("%s (%.3f)", predictor, aucs))
                     
  ax_roc <- ggplot() +
    geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1),
                 inherit.aes = FALSE,
                 lty = "dashed", lwd = lwd,
                 colour = "black", alpha = .4) +
    geom_step(data = d,  # to avoid mutiple plotting of geom_segment
              aes(x = FPR, y = TPR, col = names, lty = names),
              direction = "hv", lwd = lwd) +
    scale_color_discrete(name = element_blank(), label = labels) +
    theme_bw() +
    labs(x = "FPR", y = "TPR") +
    theme(
      legend.position = c(.95, .05),
      legend.justification = c("right", "bottom"),
      legend.background = element_rect(fill = NA)
    )
  
  if (is.logical(pauc.limits)) {
    # if no pauc.limits is provided
    ax_roc <- ax_roc +
      coord_cartesian(xlim = c(0, 1)) +
      coord_cartesian(ylim = c(0, 1))
  } else if (pauc.axis == "specificity") {
    fpr_limits <- 1 - pauc.limits
    d_rect <- data.frame(
      xmin = fpr_limits[1], xmax = fpr_limits[2],
      ymin = 0, ymax = 1
    )
    ax_roc <- ax_roc +
#       coord_cartesian(xlim = fpr_limits) +
      geom_rect(
        d_rect, mapping = aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        fill = "blue", color = "black", alpha = 0.3
      )
  } else if (pauc.axis == "sensitivity") {
    d_rect <- data.frame(
      xmin = 0, xmax = 1,
      ymin = pauc.limits[2], ymax = pauc.limits[1]
    )
    ax_roc <- ax_roc +
#       coord_cartesian(ylim = pauc.limits) +
      geom_rect(
        data = d_rect,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        fill = "blue", color = "black", alpha = 0.3
      )
  }
  
  if (return.auc) {
    return(list(plot = ax_roc, auc = aucs))
  }
  ax_roc
}


# Plots ROC and calculates AUC in a primitive fashion (i.e. ROC is step function)
# Does not resolve ties in the score
# Assumption: Lower score will be labelled preferentially as 1, ROC is step function
# Assumption that score vec and label vecs are corresponding
plot_roc <- function(score_list, label_vec,
                    name_vec = NULL, plot_title = NULL,
                    is_bigger_better_vec = rep(F, length(score_list))) {
  # Function to plot a single ROC curve and calculate AUC
  plotSingleROC <- function(score_vec, is_bigger_better, color) {
    # Sort label vector according to score vector in ascending order
    sort_label <- label_vec[order(score_vec, decreasing = is_bigger_better)]
    # Dataframe of TPR and FPR
    df <- data.frame(TPR = cumsum(sort_label)/sum(sort_label),
                     FPR = cumsum(!sort_label)/sum(!sort_label))
    # Insert 0, 0 at the first row
    roc_df <- rbind(c(0,0), df)
    # Calculates change in FPR at each step
    dFPR <- c(0, diff(roc_df$FPR))
    # Sum of individual rectangle steps
    AUC <- sum(roc_df$TPR * dFPR)
    # Plot single ROC curve
    lines(roc_df$FPR, roc_df$TPR,
          col = color, lwd = 3)
    return(AUC)
  }
  # Initialise plot
  colour_palette <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A",
                      "#66A61E", "#E6AB02", "#A6761D", "#666666")
  color_index <- colour_palette[1:length(score_list)]
  plot(NULL, main = plot_title,
       xlab = "1 - Specificity", ylab = "Sensitivity",
       xlim = c(0,1), ylim = c(0,1),
       xaxs = "i", yaxs = "i",
       cex.main = 1.8, cex.lab = 1.7, cex.axis = 1.6)
  abline(0, 1, lty = 5, lwd = 2)
  auc_vec <- mapply(plotSingleROC, score_list, is_bigger_better_vec,
                    color_index, SIMPLIFY = T)
  # If name_vec is not NULL display legend
  if (!is.null(name_vec)) {
    format_text <- function(name, auc) sprintf("%s (%.3f)", name, auc)
    legend_text <- mapply(format_text, name_vec, auc_vec)
    legend("bottomright", inset = 0.03, lty = 1, lwd = 2,
           cex = 1.2, bty = "o", text.width = 0.35,
           legend = legend_text, col = color_index)
  }
  return(auc_vec)
}


#' Performs PCA and plots vectors. Assume: 3 normal samples
#' @param X dataframe with features x samples and containing
#' D0, D8 and 3 Normal patients (ordered correctly)
plot_vectors <- function(
  X, metadata_df, pca = T, cex = 3, main = NULL
) {
  # PCA
  if (pca) {
    pca_obj <- prcomp(t(X), scale = T)
    pca_df <- data.frame(pca_obj$x[,1:2])
    n <- nrow(pca_df)
    d0_pca <- pca_df[1:((n-3) / 2),]
    d8_pca <- pca_df[((n-1) / 2):(n-3),]
    
    stopifnot(all(substring(rownames(d0_pca), 1, 4) == substring(rownames(d8_pca), 1, 4)))
    
    subtype_pca <- cbind(d0_pca, d8_pca)
    colnames(subtype_pca) <- c("start_x", "start_y", "end_x", "end_y")
    norm_pca <- pca_df[(n-2):n,]
    
    # Obtaining batch and class annotations
    label <- as.numeric(metadata_df[rownames(d0_pca), "label"]) + 1
    batch <- as.factor(metadata_df[rownames(d0_pca), "batch_info"])
    print(label)
    
    # Axis labels
    eigenvalues <- (pca_obj$sdev)^2
    var_pc <- eigenvalues[1:4]/sum(eigenvalues)
    pc_labels <- sprintf("PC%d (%.2f%%)", 1:4, var_pc*100)
  } else {
    print("No PCA performed!")
    pca_df <- data.frame(X)
  }
  
  scatter_pca <- ggplot(data = subtype_pca) +
    geom_point(aes(x = start_x, y = start_y, colour = batch), shape = 15,
               size = cex, show.legend = T) +
    geom_point(aes(x = end_x, y = end_y, colour = batch), shape = 16,
               size = cex, show.legend = F) +
    geom_segment(aes(x = start_x, y = start_y, xend = end_x, yend = end_y),
                 arrow = arrow(length = unit(0.3, "cm")),
                 alpha = 0.5, colour = label) +
#     geom_text(aes(x = PC1_A, y = PC2_A,
#                   label = rownames(subtype_pca)),
#               position = position_nudge(x = 4, y = 2), size = 2.5) +
    geom_point(data = norm_pca, aes(x = PC1, y = PC2),
               size = cex, shape = 17)
  
  if (pca) {
    scatter_pca <- scatter_pca +
      xlab(pc_labels[1]) + ylab(pc_labels[2])
  }
  
  if (!is.null(main)) {
    scatter_pca <- scatter_pca + labs(title = main)
  }
  
  return(scatter_pca)
}


#' Plots boxplot of features, p_rem WITHOUT global ERM
#' Provides p-values from wilcoxon rank-sum test
plot_boxplots <- function(X_y) {
  FEAT_ORDER <- c(
    "erm1_ratio2", "l2norm_ratio2", "angle_d0d8_d0normal",
    "log_mrd", "p"
  )
  FEAT_LABS <- c(
    "'ERM Ratio'", "'ARM Ratio'", "theta",
    "log[10](MRD)", "paste('P(Remission|', bold(x), ', s)')"
  )
  
  long_X_y <- melt(X_y, id = "label", variable.name = "feature")
  # Reorder levels and label features
  long_X_y$feature <- factor(
    long_X_y$feature,
    levels = FEAT_ORDER,
    labels = FEAT_LABS
  )
  levels(long_X_y$label) <- c("Remission", "Relapse")
  
  ax_jitter <- ggplot(
    long_X_y,
    aes(x = label, y = value)
  ) +
    geom_boxplot(
      aes(group = label),
      col = "black", alpha = 0,
      show.legend = F
    ) +
    geom_point(
      aes(fill = label),
      position = position_jitterdodge(jitter.width = 1),
      cex = 2.5, pch = 21, alpha = 1,
      col = "black", show.legend = F
    ) +
    scale_fill_manual(values = COL_LABEL) +
    facet_wrap(
      ~feature,
      nrow = 1, scales = "free",
      labeller = label_parsed
    ) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 15),
      axis.title.y = element_blank(),
      legend.position = "none"
    )
  
  # Jitter plot: p-value label
  # Both group sizes must be > 1
  if (length(table(X_y$label)) > 1 && min(table(X_y$label)) > 1) {
    p_labs <- sapply(
      FEAT_ORDER,
      function(idx) {
        list_feat <- split(X_y[, idx], X_y$label)
        tryCatch(
          {
            htest <- wilcox.test(
              list_feat[[1]], list_feat[[2]], exact  = T
            )
            sprintf("p = %.3f", htest$p.value)
          },
          error = function(err) {
            print(err)
            return(err)
          }
        )
      }
    )
    
    ann_text <- data.frame(
      feature = factor(
        FEAT_ORDER,
        levels = FEAT_ORDER,
        labels = FEAT_LABS
      ),
      value = c(
        max(X_y[FEAT_ORDER[1]]),
        max(X_y[FEAT_ORDER[2]]),
        min(X_y[FEAT_ORDER[3]]),        
        min(X_y[FEAT_ORDER[4]]),
        max(X_y[FEAT_ORDER[5]])
      ),
      lab = p_labs
    )
    ann_text <- na.omit(ann_text)

    ax_jitter <- ax_jitter +
      geom_text(
        data = ann_text,
        aes(x = 1.5, y = value, label = lab),
        size = 3, # colour = "black",
        hjust = 0
      )
  }
  
  ax_jitter
}


# 3D PCA plot
rglplot_scatter <- function(
  df, colour = "blue", pch = 21, pc_labels = NULL,
  ratio_list = list(2,1,1)
) {
  # RGL plot parameters
  rgl.open()
  rgl.bg(color="white")
  rgl.viewpoint(zoom = 0.8)
  # rgl.viewpoint(theta = 110, phi = 5, zoom = 0.8)
  par3d(windowRect = c(50, 20, 500, 500))
  pch3d(df[,1], df[,2], df[,3], bg = colour,
        pch = pch, cex = 0.5, lwd = 1.5)
  box3d(col = "black")
  # title3d(xlab = pc_labels[1], ylab = pc_labels[2],
  #         zlab = pc_labels[3], col = "black")
  # Plot aspect ratios of axis according to variance
  do.call(aspect3d, ratio_list)
}


# 3D PCA plot
plotPCA3D <- function(df, colour, pch, pc_labels = NULL,
                      ratio_list = list(2,1,1)) {
  if (is.null(pc_labels)) {
    print("PCA performed!")
    pca_obj <- prcomp(t(df), center = T, scale. = F)
    pca_df <- as.data.frame(pca_obj$x[,1:3])
    eigenvalues <- (pca_obj$sdev)^2
    var_pc <- eigenvalues[1:3]/sum(eigenvalues)
    print(var_pc)
    pc_labels <- sprintf("PC%d (%.2f%%)", 1:3, var_pc*100)
  } else {
    print("No PCA performed!")
    pca_df <- as.data.frame(df)
  }
  
  # RGL plot parameters
  rgl.open()
  rgl.bg(color="white")
  rgl.viewpoint(zoom = 0.8)
  # rgl.viewpoint(theta = 110, phi = 5, zoom = 0.8)
  par3d(windowRect = c(50, 20, 500, 500))
  with(pca_df, pch3d(PC1, PC2, PC3, bg = colour,
                     pch = pch, cex = 0.5, lwd = 1.5))
  box3d(col = "black")
  title3d(xlab = pc_labels[1], ylab = pc_labels[2],
          zlab = pc_labels[3], col = "black")
  # Plot aspect ratios of axis according to variance
  do.call(aspect3d, ratio_list)
}


# Plot PCA 3D: Batch effects
# Plot batches in different colours and classes in different shapes
plotPCA3DBatchEffects <- function(df1, metadata_df) {
  # Batch and class annotations
  batch_factor <- metadata_df[colnames(df1), "batch_info"]
  batch_palette <- generateGgplotColours(length(unique(batch_factor)))
  batch_colour <- batch_palette[batch_factor]
  
  class_factor <- metadata_df[colnames(df1), "class_info"]
  all_pch <- 21:25
  # Error if there are more classes than pch symbols (> 5)
  stopifnot(length(unique(class_factor)) <= 5)
  class_pch <- all_pch[class_factor]
  plotPCA3D(df1, batch_colour, class_pch)
}


# 3D PCA plot
plotPCA3D <- function(df, colour, pch, pc_labels = NULL,
                      ratio_list = list(2,1,1)) {
  if (is.null(pc_labels)) {
    print("PCA performed!")
    pca_obj <- prcomp(t(df), center = T, scale. = F)
    pca_df <- as.data.frame(pca_obj$x[,1:3])
    eigenvalues <- (pca_obj$sdev)^2
    var_pc <- eigenvalues[1:3]/sum(eigenvalues)
    print(var_pc)
    pc_labels <- sprintf("PC%d (%.2f%%)", 1:3, var_pc*100)
  } else {
    print("No PCA performed!")
    pca_df <- as.data.frame(df)
  }
  
  # RGL plot parameters
  rgl.open()
  rgl.bg(color="white")
  rgl.viewpoint(zoom = 0.8)
  # rgl.viewpoint(theta = 110, phi = 5, zoom = 0.8)
  par3d(windowRect = c(50, 20, 500, 500))
  with(pca_df, pch3d(PC1, PC2, PC3, bg = colour,
                     pch = pch, cex = 0.5, lwd = 1.5))
  box3d(col = "black")
  title3d(xlab = pc_labels[1], ylab = pc_labels[2],
          zlab = pc_labels[3], col = "black")
  # Plot aspect ratios of axis according to variance
  do.call(aspect3d, ratio_list)
}


# Plot PCA before selecting features
# Batch information of all the timepoints
plotPCA3DYeoh <- function(df1, metadata_df) {
  batch_info <- metadata_df[colnames(df1), "batch_info"]
  generate_colour <- colorRampPalette(c("lightblue", "darkblue"))
  batch_palette <- generate_colour(10)
  # batch_palette <- brewer.pal(10, "Set3")
  batch_colour <- batch_palette[batch_info]
  # Shape of all timepoints
  class_info <- metadata_df[colnames(df1), "class_info"]
  print(levels(class_info))
  levels(class_info) <- 21:23
  timepoint_shape <- as.numeric(as.character(class_info))
  plotPCA3D(df1, batch_colour, timepoint_shape)
}


# Plot PCA before selecting features
# Batch information of all the timepoints
plotPCA3DYeoh1 <- function(df1, metadata_df) {
  batch_info <- metadata_df[colnames(df1), "batch_info"]
  batch_factor <- droplevels(as.factor(batch_info))
  print(batch_factor)
  print(levels(batch_factor))
  levels(batch_factor) <- 21:22
  pch <- as.numeric(as.character(batch_factor))
  # generate_colour <- colorRampPalette(c("lightblue", "darkblue"))
  # batch_palette <- generate_colour(10)
  
  # Shape of all timepoints
  class_info <- metadata_df[colnames(df1), "subtype"]
  palette <- brewer.pal(10, "Set3")
  col <- palette[class_info]
  
  plotPCA3D(df1, col, pch)
}


plotHeatmapSubtype <- function(X, metadata_df) {
  subtype_factor <- metadata_df[colnames(X), "subtype"]
  set3_pal <- brewer.pal(9, "Set3")
  subtype_col <- set3_pal[subtype_factor]
  
  par(mar = c(1,1,1,1))
  heatmap(data.matrix(X),
          col = brewer.pal(9, "Blues"),
          ColSideColors = subtype_col,
          scale = "none",
          labRow = NA, labCol = NA)
  
  legend(x = -.04, y = 1350, legend = levels(subtype_factor),
         col = set3_pal[factor(levels(subtype_factor))],
         pch = 15, cex = .7)
  heatmap_subtype <- recordPlot()
  par(mar = c(5.1, 4.1, 4.1, 2.1)) # Reset to defaults
  return(heatmap_subtype)
}


plot_mean <- function(df, batch_vec1) {
  # Melt dataframe
  melt_df <- melt(df, variable.name = "ID")
  print(head(melt_df))
  # Trimmed mean probe intensities for each chip
  mean_tibble <- melt_df %>% group_by(ID) %>%
    summarise(mean = mean(value))
  mean_batch_tibble <- cbind(mean_tibble,
                             batch_vec1 = batch_vec1[mean_tibble$ID])
  
  mean_scatter <- ggplot(mean_batch_tibble, aes(x = ID, y = mean)) +
    geom_point(aes(col = factor(batch_vec1)),
               show.legend = F, size = 3) +
    facet_wrap(factor(batch_vec1), scales = "free_x") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  return(mean_scatter)
}


# Assumes that dataframe has been log-transformed
plotExplore <- function(df1, metadata_df) {
  # Obtaining batch and class annotations
  batch_factor <- as.factor(metadata_df[colnames(df1),"batch_info"])
  class_factor <- metadata_df[colnames(df1),"class_info"]
  
  # Melt dataframe
  melt_df <- melt(df1, measure.vars = colnames(df1), variable.name = "ID")
  melt_df$batch <- as.factor(metadata_df[melt_df$ID,"batch_info"])
  melt_df$class <- metadata_df[melt_df$ID,"class_info"]
  
  # Plot means
  mean_tibble <- melt_df %>% group_by(ID) %>% summarise(mean = mean(value))
  mean_scatter <- ggplot(mean_tibble, aes(x = ID, y = mean,
                                          col = batch_factor,
                                          pch = class_factor)) +
    geom_point(size = 3, show.legend = F) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  # Plot boxplot
  boxplot <- ggplot(melt_df, aes(x = ID, y = value, col = ID)) + 
    geom_boxplot(show.legend = F) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  # Plot density curve
  print(head(melt_df))
  print(tail(melt_df))
  pdf <- ggplot(melt_df, aes(x = value, group = ID, col = batch_info)) +
    geom_density(show.legend = F, alpha = 0.3) +
    facet_wrap(~class_info) +
    scale_color_viridis_d()
  
  # Plot PCA
  # Filters out rows with all zero values
  nonzero_logvec <- rowSums(df1) != 0 & apply(df1, 1, var) != 0
  pca_obj <- prcomp(t(df1[nonzero_logvec,]),
                    center = T, scale. = F)
  pca_df <- data.frame(pca_obj$x[,1:4])
  eigenvalues <- (pca_obj$sdev)^2
  var_pc <- eigenvalues[1:5]/sum(eigenvalues)
  pc_labels <- sprintf("PC%d (%.2f%%)", 1:5, var_pc*100)
  
  pc1_pc2 <- ggplot(pca_df, aes(x = PC1, y = PC2, col = batch_factor,
                                pch = class_factor)) +
    geom_point(size = 3, show.legend = F) +
    labs(x = pc_labels[1], y = pc_labels[2])
  pc1_pc3 <- ggplot(pca_df, aes(x = PC1, y = PC3, col = batch_factor,
                                pch = class_factor)) +
    geom_point(size = 3, show.legend = F) +
    labs(x = pc_labels[1], y = pc_labels[3])
  
  # Plot all graphs
  pca <- plot_grid(pc1_pc2, pc1_pc3)
  multiplot <- plot_grid(mean_scatter, pdf, pca, nrow = 3)
  return(multiplot)
}


# Plot venn diagram
jacc_coeff <- function(vec1, vec2) {
  # Generate overlap list
  overlap_list <- calculate.overlap(list(vec1,vec2))
  # Calculate venndiagram areas
  venn_area <- sapply(overlap_list, length)
  grid.newpage()
  venn_plot <- draw.pairwise.venn(venn_area[1], venn_area[2], venn_area[3],
                                  category = c("D1", "D2"),
                                  cex = 3, fontfamily = "sans",
                                  cat.cex = 3, cat.fontfamily = "sans",
                                  margin = 0.1)
  union <- (venn_area[1] + venn_area[2] - venn_area[3])
  print(unname(venn_area[3]/union))
  return(venn_plot)
}


#' Get column names in same order as pheatmap plot
#' @param obj pheatmap object
get_colnames <- function(obj) {
  obj$tree_col$labels[obj$tree_col$order]
}


#' Get row names in same order as pheatmap plot
#' @param obj pheatmap object
get_rownames <- function(obj) {
  obj$tree_row$labels[obj$tree_row$order]
}


#' Generate default ggplot colours
ggplot_palette <- function(n) {
  hues = seq(15, 375, length = n + 1)
  return(hcl(h = hues, c = 100, l = 65)[1:n])
}


#  ## Plot: Parallel coordinates - Pct
#  proba1 <- proba[, !(colnames(proba) == "p_wo_mrd")]
#  long_proba <- melt(proba1, id = c("pid", "label"),
#                    variable.name = "feature")
#             
#  ax_parallel <- ggplot(long_proba,
#                        aes(feature, value, colour = label, group = pid)) +
#    geom_line(show.legend = F) +
#    scale_color_manual(values = COL_LABEL)
#  
#  ## PLOT: CDF
#  emp_cdf <- ggplot(proba, aes(x = p, colour = label)) +
#    stat_ecdf(show.legend = F) +
#    scale_color_manual(values = COL_LABEL)
#  
#  ## PLOT: RELATIVE RISK & ODDS RATIO
#  p_sorted <- proba[order(proba$p),]
#  p_sorted$label <- as.numeric(as.character(p_sorted$label))
#  p_sorted$total_le <- rank(p_sorted$p, ties.method = "max")
#  p_sorted$total_g <- nrow(p_sorted) - p_sorted$total_le
#  p_sorted$relapse_le <- sapply(p_sorted$total_le,
#                                function(i) sum(p_sorted$label[1:i]))
#  p_sorted$relapse_g <- sum(p_sorted$label) - p_sorted$relapse_le
#  
#  p_sorted <- within(
#    p_sorted,
#    relative_risk <- (relapse_le/total_le) / (relapse_g/total_g)
#  )
#  
#  p_sorted <- within(
#    p_sorted,
#    odds_ratio <- (relapse_le/(total_le-relapse_le)) / (relapse_g/(total_g-relapse_g))
#  )
#                                 
#  ax_rr_or <- ggplot(p_sorted) +
#    geom_step(aes(p, relative_risk, colour = "RR"), direction = "hv") + 
#    geom_step(aes(p, odds_ratio, colour = "OR"), direction = "hv") +
#    scale_color_manual("",
#                       breaks = c("RR", "OR"),
#                       values = c("RR" = "orange", "OR" = "steelblue3")) +
#    theme(axis.title.y = element_blank())
#  
#  ## Plot: ROC
#  # ERM1 evaluated is not from global GSS model
#  proba_x <- cbind(proba, erm = X_predict$erm1, d33_mrd = X_predict$mrd) # subset mrd
#                                
#  x_names <- c("p", "erm", "d33_mrd")
#  # WARNING: Change bigger.positive according to features!
#  bigger.positive <- c(F, T, F) # bigger means relapse
#  
#  # ROC can only be plotted when there are both positive and negative samples
#  if (all(table(proba_x$label) != 0)) {
#    ax_roc <- plot_roc(proba_x, "label", x_names)
#    # Able to plot ROC
#    ax2 <- plot_grid(ax_parallel, ax_roc,
#                     ncol = 2, rel_widths = c(1.8, 1))
#  } else{
#    ax2 <- ax_parallel # unable to plot ROC
#  }
#  
#  # Plot: MRD v.s. Risk of relapse
#  mrd_p <- ggplot(proba_x) +
#    geom_point(aes(p, log10(d33_mrd), colour = label),
#               cex = 3, show.legend = F) +
#    scale_color_manual(values = COL_LABEL)
#                                
#  ax1 <- plot_grid(ax_jitter, mrd_p,
#                   ncol = 2, rel_widths = c(2.8, 1))
#  
#  fig <- plot_grid(ax1, ax2, nrow = 2)
