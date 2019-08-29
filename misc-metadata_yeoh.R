library(tidyr)
library(ggplot2)
library(cowplot)
source("../functions.R")

dates_df1 <-read.table("data/GSE67684/processed/scan_dates-GPL570.tsv",
                       sep = "\t", header = F, stringsAsFactors = F)
dates_df2 <-read.table("data/GSE67684/processed/scan_dates-GPL96.tsv",
                       sep = "\t", header = F, stringsAsFactors = F)
dates_df1 <- cbind(dates_df1, platform = "GPL570")
dates_df2 <- cbind(dates_df2, platform = "GPL96")
metadata_df <- rbind(dates_df1, dates_df2)
metadata_df1 <- separate(metadata_df,
                         col = 1,
                         into = c("geo_accession", "pid", "time_point"),
                         sep = "_")

metadata_df1[,3] <- substring(metadata_df1[,3], 1, 2)
metadata_df1[146:381,4] <- as.character(as.Date(strptime(metadata_df1[146:381,4], format = "%Y-%m-%d")))
metadata_df1[-(146:381), 4] <- as.character(as.Date(strptime(metadata_df1[-(146:381), 4], format = "%m/%d/%y")))
colnames(metadata_df1)[4] <- "scan_date"

write.table(metadata_df1, "data/GSE67684/processed/metadata.tsv",
            quote = F, sep = "\t", row.names = F, col.names = T)

# Further processing
metadata_df1 <- read.table("data/GSE67684/processed/metadata.tsv",
                           sep = "\t", header = T)
rownames(metadata_df1) <- metadata_df1[,1]

# metadata_df1 <- metadata_df1[order(metadata_df1[,1]),]
info_gpl570 <- read.table("data/GSE67684/README/GPL570.txt",
                          sep = "\t", header = F)
info_gpl570 <- t(info_gpl570)[-1,]
# Mistake in original annotation file
info_gpl570[324:381, 3:4] <- info_gpl570[324:381, 2:3]
info_gpl570[,4] <- substring(info_gpl570[,4], 10)

metadata_df1[info_gpl570[,1], 6] <- info_gpl570[,4]

info_gpl96 <- read.table("data/GSE67684/README/GPL96.txt",
                          sep = "\t", header = F)
info_gpl96 <- t(info_gpl96)[-1,]
info_gpl96[,4] <- substring(info_gpl96[,4], 10)
metadata_df1[info_gpl96[,1], 6] <- info_gpl96[,4]

labels <- read.table("data/GSE67684/README/labels_MRD.tsv",
                     sep = "\t", header = F, skip = 1)
metadata_df2 <- merge(metadata_df1, labels[, 2:6], by.x = "pid", by.y = "V2")
colnames(metadata_df2)[6:10] <- c("training_test", "d33_mrd",	"d33_mrd_risk", "status", "event_code")
metadata_df3 <- metadata_df2[,c(2,1,3:10)]

write.table(metadata_df3, "data/GSE67684/processed/metadata_all.tsv",
            quote = F, sep = "\t", row.names = F, col.names = T)
table(paste(as.character(metadata_df3$scan_date), metadata_df3$training_test))

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


# BATCH ASSIGNMENT --------------------------------------------------------
yeoh_metadata <- read.table("data/GSE67684/processed/metadata_batch.tsv",
                            sep = "\t", header = T, row.names = 1)
# Convert to type date
yeoh_metadata$month <- as.Date(yeoh_metadata$month)
# # Convert dates to the same month
# yeoh_metadata$month <- as.Date(cut(yeoh_metadata$scan_date, breaks = "month"))

batch_date <- ggplot(yeoh_metadata, aes(x = month, fill = factor(batch))) +
  geom_bar(width = 25, show.legend = F) +
  xlim(as.Date(c("2003-01-01", "2016-01-01"))) +
  scale_x_date(date_breaks = "12 months", date_minor_breaks = "1 months") +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
batch_date

ggsave("dump/batch_date.pdf", batch_date,
       width = 12, height = 6)
x <- yeoh_metadata[yeoh_metadata$batch == 9,]

# Batch 1 has all GPL96
yeoh_metadata[yeoh_metadata$scan_date >= "2006-01-01" & yeoh_metadata$scan_date <= "2007-01-01", 11] <- 2
yeoh_metadata[yeoh_metadata$scan_date >= "2007-01-01" & yeoh_metadata$scan_date <= "2009-01-01", 11] <- 3
yeoh_metadata[yeoh_metadata$training_test == "Test set" & yeoh_metadata$scan_date <= "2012-01-01", 11] <- 4
yeoh_metadata[yeoh_metadata$scan_date == "2011-06-24", 11] <- 5
yeoh_metadata[yeoh_metadata$scan_date >= "2007-01-01" & yeoh_metadata$scan_date <= "2009-01-01", 11] <- 6
yeoh_metadata[yeoh_metadata$scan_date >= "2013-01-01" & yeoh_metadata$scan_date <= "2013-06-01", 11] <- 7
yeoh_metadata[yeoh_metadata$scan_date >= "2013-07-01" & yeoh_metadata$scan_date <= "2014-01-01", 11] <- 8
yeoh_metadata[yeoh_metadata$scan_date >= "2015-01-01" & yeoh_metadata$scan_date <= "2016-01-01", 11] <- 9
# Inconsistencies
yeoh_metadata[yeoh_metadata$pid == 107 & yeoh_metadata$time_point == "D0", 11] <- 3

yeoh_metadata1 <- yeoh_metadata[,c(1:4,12,11,5:10)]
write.table(yeoh_metadata1, "data/GSE67684/processed/metadata_batch.tsv",
            quote = F, sep = "\t", row.names = F)

table(yeoh_metadata1$batch)

# My assignment of batches
yeoh_metadata <- read.table("data/GSE67684/processed/metadata_batch.tsv",
                            sep = "\t", header = T, row.names = 1)
colnames(yeoh_metadata)

ans_metadata <- read.table("data/leuk_D33/README/batch.info.txt",
                            sep = "\t", header = T, row.names = 1)
old_pid <- as.character(ans_metadata$ID_GEO)
patient_num <- as.numeric(sapply(strsplit(old_pid, "_"), `[[`, 1))
time_point <- sapply(strsplit(old_pid, "_"), `[[`, 2)
pid <- sprintf("P%03d_%s", patient_num, time_point)
# Subsets based on order in ans_metadata
compare <- cbind(yeoh_metadata[pid,"batch"], ans_metadata$Batch_No)
head(yeoh_metadata)

subset_metadata <- yeoh_metadata[,c(1,4,7)]
head(subset_metadata)

# Add D33 samples
d33_dates <- read.table("data/leuk_D33/processed/scan_dates.tsv",
                        sep = "\t", header = F, row.names = 1)

scan_date <- as.character(as.Date(strptime(d33_dates[,1], format = "%Y-%m-%d")))
d33_metadata <- data.frame(geo_accession = NA, scan_date, platform = "GPL570")
patient_num <- as.numeric(sapply(strsplit(substring_head(rownames(d33_dates), 4), "_"), `[[`, 1))
d33_pid <- sprintf("P%03d_D33", patient_num)
rownames(d33_metadata) <- d33_pid
d33_metadata

# Add normal samples
normal_dates <- read.table("data/leuk_normal/processed/scan_dates.tsv",
                           sep = "\t", header = F, row.names = 1)
scan_date <- as.character(as.Date(strptime(normal_dates[,1], format = "%m/%d/%y")))
normal_metadata <- data.frame(geo_accession = NA, scan_date, platform = "GPL570")
rownames(normal_metadata) <- sprintf("N%02d", 1:4)

# Combine all samples
combined_metadata <- rbind(subset_metadata, d33_metadata, normal_metadata)
View(combined_metadata)

# Obtain batch information
rownames(ans_metadata) <- ans_metadata[,1]
batch_info <- ans_metadata[rownames(combined_metadata)[1:479], 2]
batch_metadata <- cbind(combined_metadata, batch = c(batch_info, rep(2, 4)))

write.table(batch_metadata, "data/GSE67684/processed/metadata_combined-batch.tsv",
            quote = F, sep = "\t")

# Appending training set information and labels to MRD info
mrd_df <- read.table("data/GSE67684/README/labels_MRD.tsv",
                     sep = "\t", header = T, comment.char = "")
rownames(mrd_df) <- sprintf("P%03d", mrd_df$geo_patient_.)
training_test <- yeoh_metadata[paste0(rownames(mrd_df), "_D0"), "training_test"]
colnames(mrd_df)
patient_metadata <- cbind(mrd_df[,3:6], label = mrd_df$event_code, training_test)
patient_metadata$label[patient_metadata$label == 2] <- 1
write.table(patient_metadata, "data/GSE67684/processed/metadata_combined-label.tsv",
            quote = F, sep = "\t")

# Appending subtype information to metadata labels
subtype_df <- read.table("data/GSE67684/README/subtype_erm_data.txt",
                         sep = "\t", header = T)
mrd_df <- read.table("data/GSE67684/README/labels_MRD.tsv",
                     sep = "\t", header = T, comment.char = "")
merged_df <- merge(subtype_df, mrd_df[,1:2], by.x = "lab_id", by.y = "labid")
rownames(merged_df) <- sprintf("P%03d", merged_df[,3])

metadata_label <- read.table("data/GSE67684/processed/metadata_combined-label.tsv",
                             sep = "\t", header = T)
head(metadata_label1)
subtype <- merged_df[rownames(metadata_label), "subtype"]
metadata_label1 <- cbind(metadata_label, subtype)

write.table(metadata_label1, "data/GSE67684/processed/metadata_combined-label_subtype.tsv",
            quote = F, sep = "\t")
