source("../functions.R")

# IMPORT MICROARRAY DATA -------------------------------------------------------------
## Subset of original data
SUBSET_RPATH <- "data/GSE67684/processed/subset_yeoh.tsv"
raw_yeoh <- read.table(SUBSET_RPATH, sep = "\t")

## Metadata
# Preprocessed metadata
METADATA_RPATH <- "data/GSE67684/processed/metadata/full_metadata.tsv"
metadata_df <- read.table(METADATA_RPATH, sep = "\t")

# SCALE->REMOVE->FILTER->LOG
scaled_yeoh <- removeProbesets(normaliseMeanScaling(raw_yeoh))
data_yeoh <- log2_transform(filterProbesets(scaled_yeoh, 0.7, metadata_df))

# Map probesets to gene symbols
SYMBOL_GPL570 <- "../info/microarray/HG-U133_Plus_2/annot_genesymbol-GPL570.tsv"
symbol_yeoh <- affy2id(data_yeoh, SYMBOL_GPL570)

## RNA-Seq
# Temporary lab id lookup table
# Has to be updated to include new patients
LABID_RPATH <- "data/GSE67684/processed/metadata/lab_id.tsv"
labid_df <- read.table(LABID_RPATH, sep = "\t", header = T, row.names = 1)

RNA_RPATH <- "data/GSE67684/processed/rna_seq/count.matrix.txt"
raw_rna <- read.table(RNA_RPATH, sep = "\t")

labid_data <- unlist(lapply(strsplit(colnames(raw_rna), "_"),
                            function(vec) vec[2]))
# Only select patients present in microarray data
logi_idx <- labid_data %in% rownames(labid_df)
subset_rna <- raw_rna[,logi_idx]
# Replace labid of subsetted patients with pid
colnames(subset_rna) <- labid_df[labid_data[logi_idx],]
colnames(subset_rna)[1:67] <- paste(colnames(subset_rna)[1:67], "D0", sep = "_")
colnames(subset_rna)[68:134] <- paste(colnames(subset_rna)[68:134], "D8", sep = "_")

# Check whether the feature names intersect
intersect_genes <- intersect(rownames(symbol_yeoh), rownames(subset_rna))
# Select features that are present in both
data_rna <- subset_rna[intersect_genes,]
data_yeoh <- symbol_yeoh[intersect_genes,]

# TEMP DATA FILES
RNA_WPATH <- "data/GSE67684/processed/rna_seq/temp_rna.tsv"
AFFY_WPATH <- "data/GSE67684/processed/rna_seq/temp_affy.tsv"
# write.table(data_rna, RNA_WPATH, quote = F, sep = "\t")
# write.table(data_yeoh, AFFY_WPATH, quote = F, sep = "\t")
data_rna <- read.table(RNA_WPATH, sep = "\t")
data_yeoh <- read.table(AFFY_WPATH, sep = "\t")

LENGTH_RPATH <- "../info/ref_genome/GRCh38/gene_length.tsv"
length_df <- read.table(LENGTH_RPATH, sep = "\t", header = T)
head(length_df)
# MAIN --------------------------------------------------------------------
rownames(subset_rna)

# Log2 normalised microarray vs log2 raw counts
# Same patients are not present in microarray data
intersect_pid <- intersect(colnames(data_rna), colnames(data_yeoh)) 

pid <- intersect_pid[1]
print(pid)
x <- log2_transform(data_rna[,pid])
colnames(data_yeoh)
r <- cor(data_yeoh[,pid], x)
XLAB <- "Normalised microarray (log2)"
YLAB <- "RNA-Seq raw counts (log2)"
MAIN <- sprintf("%s (r = %.3f)", pid, r)
plot(data_yeoh[,pid], x,
     xlab = XLAB, ylab = YLAB, main = MAIN)
log_rna <- recordPlot()
save_fig(log_rna, "dump/scatter-log_rna.pdf", 8, 8)



hist(x, breaks = 100)
