raw_data <- read.table("data/D33/processed/mas5-original.tsv",
                       sep = "\t", header = T, row.names = 1)
yeoh_data <- read.table("data/GSE67684/processed/mas5_ordered.tsv",
                        sep = "\t", header = T, row.names = 1)

select_logvec <- rownames(raw_data) %in% rownames(yeoh_data)
write.table(raw_data[select_logvec,], "data/D33/processed/mas5_filtered.tsv",
            quote = F, sep = "\t")
