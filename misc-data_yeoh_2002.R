library(RColorBrewer)
library(ggplot2)
library(dendextend)
library(cowplot)
source("../functions.R")
theme_set(theme_cowplot())

import_data <- function(fpath, alphabet) {
  df <- read.table(fpath, sep = "\t", header = T, row.names = 1, comment.char = "")
  original_names <- unlist(strsplit(readLines(fpath, n = 1), "\t"))
  new_names <- sprintf("%s%02d", alphabet, 1:length(original_names))
  annot_arr <- cbind(new_names, original_names)
  colnames(df) <- new_names
  return(list(df, annot_arr))
}

fpath <- list.files("data/yeoh_2002/processed", full.names = T)
# Ordered in order to match subset of class alphabets with MILE dataset
ordered_fpath <- fpath[c(6,3,2,9,1,4,7,5,8)]
alphabet_vec <- c("N", LETTERS[1:8])

# Mapply concatenates all the lists returned creating one-layer list
list_alternate <- mapply(import_data, ordered_fpath, alphabet_vec)
# Returns all odd index of the list
list_df <- list_alternate[seq(1,length(list_alternate),2)]
all_df <- do.call(cbind, list_df)
# Returns all even index of the list
list_annot <- list_alternate[seq(2,length(list_alternate),2)]
annot_df <- do.call(rbind, list_annot)

col_logvec <- startsWith(colnames(all_df), "E") | startsWith(colnames(all_df), "F")
reordered_df <- cbind(all_df[,!col_logvec], all_df[,col_logvec])
write.table(reordered_df, "data/yeoh_2002/processed/mas5_original.tsv",
            sep = "\t", quote = F)

# Scan dates
scan_dates <- read.table("data/yeoh_2002/README/scan_dates.tsv",
                         sep = "\t", header = F, comment.char = "", stringsAsFactors = F)
scan_dates[,1] <- substring_head(scan_dates[,1], 4)
scan_dates[,2] <- substring(scan_dates[,2], 1, 8)
scan_dates[,2] <- as.character(as.Date(strptime(scan_dates[,2], format = "%m/%d/%y")))
colnames(scan_dates) <- c("patient_names","scan_date")

metadata_df <- merge(data.frame(annot_df), scan_dates, by.x = "original_names", by.y = "patient_names")[,c(2,1,3)]
colnames(metadata_df)[1] <- "pid"
ordered_metadata <- metadata_df[order(metadata_df$pid),]
write.table(ordered_metadata, "data/yeoh_2002/processed/metadata.tsv",
            quote = F, sep = "\t", row.names = F)

ordered_metadata
# yeoh_metadata$month <- as.Date(cut(yeoh_metadata$scan_date, breaks = "month"))

# Plot count of scans by date
batch_date <- ggplot(scan_dates, aes(x = as.Date(scan_date))) +
  geom_bar(show.legend = F) +
  # xlim(as.Date(c("2003-01-01", "2016-01-01"))) +
  scale_x_date(date_breaks = "1 weeks", date_minor_breaks = "1 days") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("dump/scan_dates-yeoh_2002.pdf", batch_date,
       width = 12, height = 6)

# Hierachical clustering
rownames(ordered_metadata) <- ordered_metadata$pid
mapped_dates <- ordered_metadata[colnames(reordered_df), 3]
colnames(reordered_df) <- sprintf("%s (%s)", colnames(reordered_df), mapped_dates)

# Normalisation: Scaling
scaled_df <- norm_mean_scaling(reordered_df)
selected_probesets <- filter_probesets(scaled_df, 0.1)
filtered_df <- scaled_df[selected_probesets,]

# Creation of dendrogram object
pairwise_dist <- dist(t(filtered_df))
hcluster <- hclust(pairwise_dist)
dendo_obj <- as.dendrogram(hcluster)
# sample_id <- labels(dendo_obj)
# nodePar <- list(lab.cex = 0.3, pch = c(NA, NA),
#                 cex = 0.5, col = "blue")

# Subtype labels
subtype_info <- as.numeric(as.factor(substring(colnames(filtered_df), 1, 1)))
subtype_palette <- brewer.pal(9, "Set1")
subtype_colour <- subtype_palette[subtype_info]
# Date labels
date_numeric <- as.numeric(as.Date(mapped_dates))
generate_colour <- colorRampPalette(c("lightblue", "darkblue"))
date_palette <- generate_colour(length(unique(date_numeric)))
names(date_palette) <- sort(unique(date_numeric))
date_colour <- date_palette[as.character(date_numeric)]

# Settings of dendogram
dendo_obj <- set(dendo_obj, "labels_cex", 0.2)
plot(dendo_obj, horiz = F)
colored_bars(cbind(subtype_colour, date_colour), dendo_obj,
             rowLabels = c("Subtype  ", "Date  "), y_shift = -3e5,
             sort_by_labels_order = T)
dendogram <- recordPlot()
save_fig(dendogram, "dump/hclust-yeoh_2002_raw.pdf",
         width = 12, height = 6)

# PCA
pca_plot <- function(untransformed_df, colour_vec, shape_vec = 19) {
  pca_obj <- prcomp(t(untransformed_df))
  df <- data.frame(pca_obj$x[,1:5])
  eig_value <- (pca_obj$sdev)^2
  var_pc <- eig_value[1:5]/sum(eig_value)
  pc_labels <- sprintf("PC%d (%.2f%%)", 1:5, var_pc*100)
  # Scatter plots
  pc1_pc2 <- ggplot(df, aes(x = PC1, y = PC2)) +
    geom_point(size = 2, col = colour_vec, shape = shape_vec, show.legend = F) +
    xlab(pc_labels[1]) + ylab(pc_labels[2])
  pc2_pc3 <- ggplot(df, aes(x = PC2, y = PC3)) +
    geom_point(size = 2, col = colour_vec, shape = shape_vec, show.legend = F) +
    xlab(pc_labels[2]) + ylab(pc_labels[3])
  pc1_pc3 <- ggplot(df, aes(x = PC1, y = PC3)) +
    geom_point(size = 2, col = colour_vec, shape = shape_vec, show.legend = F) +
    xlab(pc_labels[1]) + ylab(pc_labels[3])
  pc3_pc4 <- ggplot(df, aes(x = PC3, y = PC4)) +
    geom_point(size = 2, col = colour_vec, shape = shape_vec, show.legend = F) +
    xlab(pc_labels[3]) + ylab(pc_labels[4])
  multiplot <- plot_grid(pc1_pc2, pc2_pc3, pc1_pc3, pc3_pc4,
                         ncol = 2, nrow = 2)
  return(multiplot)
}
subtype_shape <- as.numeric(as.factor(subtype_colour)) - 1
pca_yeoh <- pca_plot(filtered_df, date_colour, subtype_shape)
ggsave("dump/pca-yeoh_2002_scaled.pdf", pca_yeoh,
       width = 9, height = 9)
