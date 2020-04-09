RNASEQ_RPATH <- "data/GSE67684/README/rnaseq-metadata.tsv"
AFFY_RPATH <- "data/GSE67684/README/microarray-metadata.tsv"

rnaseq_meta <- read.table(RNASEQ_RPATH, header = T, row.names = 1,
                          sep = "\t", comment.char = "")
affy_meta <- read.table(AFFY_RPATH, header = T, row.names = 1,
                        sep = "\t", comment.char = "")
head(rnaseq_meta)
head(affy_meta)

common_pid <- intersect(rownames(rnaseq_meta), rownames(affy_meta))

a <- merge(affy_meta, rnaseq_meta, by = "row.names")
names(rnaseq_meta)
