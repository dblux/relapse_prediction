# Perform an inner join on both dataframes to obtain probesets found in both
GPL570 <- read.table("data/GSE67684/processed/mas5-GPL570.tsv",
                     sep = "\t", header = T, row.names = 1)
GPL96 <- read.table("data/GSE67684/processed/mas5-GPL96.tsv",
                     sep = "\t", header = T, row.names = 1)
sum(rownames(GPL96) %in% rownames(GPL570))
yeoh_data <- merge(GPL570, GPL96, by = "row.names")
rownames(yeoh_data) <- yeoh_data[,1]
yeoh_data <- yeoh_data[,-1]
colnames(yeoh_data)[1:381] <- substring_head(colnames(yeoh_data)[1:381], 7)
colnames(yeoh_data) <- substring(colnames(yeoh_data), 12)
write.table(yeoh_data, "data/GSE67684/processed/mas5_unordered.tsv",
            quote = F, sep = "\t", row.names = T, col.names = T)

yeoh_data <- read.table("data/GSE67684/processed/mas5_unordered.tsv",
                        sep = "\t", header = T, row.names = 1)
# Rename column names
yeoh_col <- colnames(yeoh_data)
pid <- sapply(substring(yeoh_col, 2),
              function(string) {
                sprintf("%03d_%s",
                        as.numeric(unlist(strsplit(string, "_"))[1]),
                        unlist(strsplit(string, "_"))[2]) 
              })
colnames(yeoh_data) <- paste0("P", pid)

# Sort columns
ordered_data0 <- yeoh_data[, order(colnames(yeoh_data))]
colnames(ordered_data0)
ordered_index <- c(which(sapply(colnames(ordered_data0), endsWith, "0")),
                   which(sapply(colnames(ordered_data0), endsWith, "8")))
ordered_data1 <- ordered_data0[, ordered_index]

write.table(ordered_data1, "data/GSE67684/processed/mas5_ordered.tsv",
            quote = F, sep = "\t")

##### SUBSETTING DATA #####
yeoh_d0d8 <- read.table("data/GSE67684/processed/mas5_ordered.tsv",
                        sep = "\t", header = T, row.names = 1)
yeoh_normal <- read.table("data/leuk_normal/processed/mas5_filtered.tsv",
                          sep = "\t", header = T, row.names = 1)
yeoh_d33 <- read.table("data/leuk_D33/processed/mas5_filtered.tsv",
                       sep = "\t", header = T, row.names = 1)

## Removal of outlier samples and their associated pairs
# Patients "P198_D8", "P186_D0" are associated pairs
outlier_samples <- c("P198_D0", "P186_D8", "N03", "P198_D8", "P186_D0")

## Identify patients with D0 and D8 profiles from different batches
# for (i in seq(1,420,2)) {
#   if (yeoh_batch[i,"batch"] != yeoh_batch[i+1,"batch"]) {
#     print(yeoh_batch[i:(i+1),])
#   }
diffbatch_pid <- c("P107", "P110", "P112", "P113", "P114", "P118", "P168")
outlier_pid <- c("P198", "P186", "N03")
remove_pid <- c(diffbatch_pid, outlier_pid)

yeoh_all <- cbind(yeoh_d0d8, yeoh_normal)
logi_idx <- !(substring(colnames(yeoh_all), 1, 4) %in% remove_pid)
subset_yeoh <- yeoh_all[,logi_idx]
SUBSET_WPATH <- "data/GSE67684/processed/subset_yeoh.tsv"
write.table(subset_yeoh, SUBSET_WPATH,
            quote = F, sep = "\t")

# # Selection of d33 remission samples
# d33_label <- yeoh_label[substring(colnames(yeoh_d33), 1, 4), "label", drop = F]
# # D33 patients that experience remission
# d33_remission <- rownames(d33_label)[d33_label == 0 & !is.na(d33_label)]
# yeoh_remission <- yeoh_d33[, paste0(d33_remission, "_D33")]

# # Added D33 samples and then removed them as batch 5 samples!!!
# yeoh_combined <- cbind(yeoh_d0d8, yeoh_remission, yeoh_normal)
## Remove samples with D0 and D8 from different batches and batch 5

# Batch 5 contains all D33 samples and P107_D*=8
# As a result P107_D0 is removed as well

### PREPROCESSING: SCALE->REMOVE->FILTER->LOG
scaled_yeoh <- removeProbesets(normaliseMeanScaling(subset_yeoh))
data_yeoh <- log2_transform(filterProbesets(scaled_yeoh, 0.7, metadata_df))
