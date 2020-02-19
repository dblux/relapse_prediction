RPATH1 <- "data/GSE67684/README/batch.info.txt"
RPATH2 <- "data/GSE67684/processed/rna_seq/rna_seq-pid.txt"

batch <- read.table(RPATH1, sep = "\t", header = TRUE)
existing_pid <- sub("_.*", "", batch$ID)

new_pid <- scan(RPATH2, what = "character", nlines = 1, sep = "\t")
new_pid1 <- sub("_.*", "", substring(new_pid[1:234], 4))

new_pid
new_pid[!(new_pid1 %in% existing_pid)]


##### NOT ALL GENE SYMOBLS HAVE CORRESPONDING GENE LENGTHS #####
### USE NCBI ANNOTATION TO CONVERT SYNONYMS
# Use latest version of NCBI annotation
RNA_RPATH <- "data/GSE67684/processed/rna_seq/count.matrix.txt"
raw_rna <- read.table(RNA_RPATH, sep = "\t")

# Table of gene lengths
LENGTH_RPATH <- "../info/ref_genome/GRCh38/GRCh38_99-gene_lengths.tsv"
length_df <- read.table(LENGTH_RPATH, sep = "\t", header = T, row.names = 1)

## Identify (no length) genes 
no_length <- rownames(raw_rna)[!(rownames(raw_rna) %in% rownames(length_df))]

### Compile genes that are represented in NCBI ###
NCBI_RPATH <- "../info/NCBI/refseq-Homo_sapiens.gene_info"
nlines <- length(readLines(NCBI_RPATH))
con  <- file(NCBI_RPATH, open = "r")
gene_list <- vector("list", nlines)

i <- 1
while(length(oneLine <- readLines(con, n = 1)) > 0) {
  if (!startsWith(oneLine, "#")) {
    myVector <- unlist(strsplit(oneLine, "\t"))
    gene_symbol <- as.list(myVector[c(3,5,14)])
    if (grepl("\\|", gene_symbol[[2]])) {
      gene_symbol[2] <- strsplit(gene_symbol[[2]], "\\|")
    }
    if (grepl("\\|", gene_symbol[[3]])) {
      gene_symbol[3] <- strsplit(gene_symbol[[3]], "\\|")
    }
    gene_list[[i]] <- gene_symbol
    i <- i + 1
  }
}
close(con)

# Genes represented in NCBI
all_genes <- unlist(gene_list)
# Identify no length genes that can be converted using NCBI annotation
ncbi_genes <- no_length[no_length %in% all_genes]
# Convert genes using NCBI annotation
fruits <- list(list("apple", "pear"), list("orange", "banana"))

Filter(function(item) ifelse("apple" %in% item[1], TRUE, FALSE), fruits)

#' @param name gene name
#' @return gene symbol
searchReplace <- function(name, gene_list) {
  Filter(function(sublist) {
    ifelse(
    name %in% sublist[[2]] | name %in% sublist[[3]], TRUE, FALSE)},
    gene_list
  )
}

temp_list <- lapply(ncbi_genes, searchReplace, gene_list)
searchReplace("FIGF", gene_list)
ncbi_genes[1]
# Identify genes that are not present in NCBI annotation
# Majority of genes are non-coding RNA or discontinued
# Genes will be removed from count matrix
no_ncbi <- no_length[!(no_length %in% all_genes)]

NCBI_GENES <- "../info/NCBI/refseq-Homo_sapiens.gene_info"
