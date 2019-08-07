library(ggplot2)

subtype_annot <- read.table("data/GSE13204/README/subtype_annot.txt",
                            header = F, row.names = 1)
subtype_df <- t(subtype_annot)
subtype_levels <- unique(subtype_df[,2])
normal_id <- subtype_df[subtype_df[,2] == "leukemia class: Non-leukemia and healthy bone marrow", 1]
length(normal_id)

# Write lines to text file
normal_fileconn <- file("normal_id.txt", open = "w")
writeLines(normal_id, normal_fileconn)
close(normal_fileconn)

lymph_leukemia_id <- subtype_df[subtype_df[,2] %in% subtype_levels[1:8], 1]
# Write lines to text file
lymph_leukemia_fileconn <- file("lymph_leukemia_id.txt", open = "w")
writeLines(lymph_leukemia_id, lymph_leukemia_fileconn)
close(lymph_leukemia_fileconn)

# Create annotation file for subtypes
normal_df <- subtype_df[subtype_df[,2] == "leukemia class: Non-leukemia and healthy bone marrow",]
lymph_leukemia_df <- subtype_df[subtype_df[,2] %in% subtype_levels[1:8],]
normal_ALL_df <- rbind(lymph_leukemia_df, normal_df)
write.table(normal_ALL_df, "data/GSE13204/README/normal_ALL_subtype_annot.tsv",
            quote = F, sep = "\t", row.names = F, col.names = T)

# Create metadata file incorporating scan dates
subtype_df <- read.table("data/GSE13204/processed/normal_ALL_subtype_annot.tsv",
                         sep = "\t", header = T)
rownames(subtype_df) <- subtype_df[,1]

dates_rpaths <- c("data/GSE13204/processed/scan_dates-ALL_1.tsv",
                  "data/GSE13204/processed/scan_dates-ALL_2.tsv",
                  "data/GSE13204/processed/scan_dates-ALL_3.tsv",
                  "data/GSE13204/processed/scan_dates-normal.tsv")

dates_list <- lapply(dates_rpaths, read.table, sep = "\t", header = F, stringsAsFactors = F)
dates_df <- do.call(rbind, dates_list)

dates_df[,1] <- substring(dates_df[,1], 1, 9)
dates_df[,2] <- as.Date(strptime(dates_df[,2], "%m/%d/%y"))
metadata_df <- merge(subtype_df, dates_df, by.x = "X.Sample_geo_accession", by.y = "V1")

sample_type <- read.table("data/GSE13204/README/sample_type.txt",
                          sep = "\t", header = F, row.names = 1)
sample_type <- t(sample_type)
metadata_df1 <- merge(metadata_df, sample_type, by.x = 1, by.y = 1)
colnames(metadata_df1) <- c("geo_accession", "subtype", "scan_date", "sample_type")
write.table(metadata_df1, "data/GSE13204/processed/metadata.tsv",
            quote = F, sep = "\t", row.names = F)
metadata_df2 <- with(metadata_df1, metadata_df1[order(subtype, scan_date),])

subtype_factors <- metadata_df2[,2]
levels(subtype_factors) <- c(LETTERS[1:6], "N", LETTERS[7:8])
num_subtype <- table(subtype_factors)
pid <- unname(unlist(mapply(function(num, subtype) sprintf("%s%03d", subtype, 1:num),
                            num_subtype, names(num_subtype))))
metadata_df3 <- cbind(pid, metadata_df2)
metadata_df3 <- metadata_df3[order(metadata_df3$pid),]
write.table(metadata_df3, "data/GSE13204/processed/metadata.tsv",
            quote = F, sep = "\t", row.names = F)

# Consolidate all expression data
mas5_rpaths <- c("data/GSE13204/processed/mas5-ALL_1.tsv",
                 "data/GSE13204/processed/mas5-ALL_2.tsv",
                 "data/GSE13204/processed/mas5-ALL_3.tsv",
                 "data/GSE13204/processed/mas5-normal.tsv")

mas5_list <- lapply(mas5_rpaths, read.table, sep = "\t", header = T, row.names = 1)
mas5_df <- do.call(cbind, mas5_list)

# Select probes common with other dataset
yeoh_df <- read.table("data/GSE67684/processed/mas5_unordered.tsv",
                      sep = "\t", header = T, row.names = 1)
mas5_commonprobes <- mas5_df[rownames(mas5_df) %in% rownames(yeoh_df),]

# Change colnames from GEO accession to patient ID
rownames(metadata_df3) <- metadata_df3$geo_accession
colnames(mas5_commonprobes) <- metadata_df3[colnames(mas5_commonprobes), 1]
mas5_ordered <- mas5_commonprobes[,order(colnames(mas5_commonprobes))]
write.table(mas5_ordered, "data/GSE13204/processed/mas5_ordered.tsv",
            quote = F, sep = "\t")
saveRDS(mas5_ordered, "data/GSE13204/processed/mas5_ordered.rds")


# Batch information -------------------------------------------------------
mile_metadata <- read.table("data/GSE13204/processed/metadata.tsv",
                            sep = "\t", header = T, row.names = 1)

mile_metadata[,3] <- as.Date(strptime(mile_metadata[, 3], format = "%Y-%m-%d"))
# # Convert dates to the same month
# yeoh_metadata$month <- as.Date(cut(yeoh_metadata$scan_date, breaks = "month"))

plot_date <- ggplot(mile_metadata, aes(x = scan_date)) +
  geom_bar(show.legend = F) +
  # xlim(as.Date(c("2003-01-01", "2016-01-01"))) +
  # scale_x_date(date_breaks = "12 months", date_minor_breaks = "1 months") +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
