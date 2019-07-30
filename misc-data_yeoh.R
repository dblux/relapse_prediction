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
