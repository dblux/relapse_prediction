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

### Compile genes that are represented in NCBI ANNOTATION ###
NCBI_ANNOT <- "../info/NCBI/refseq-Homo_sapiens.gene_info"
nlines <- length(readLines(NCBI_ANNOT))
con  <- file(NCBI_ANNOT, open = "r")
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
ncbi_annot_genes <- no_length[no_length %in% all_genes]
ncbi_annot_genes

# Convert genes using NCBI annotation
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

# Identify genes that are not present in NCBI annotation
# Majority of genes are non-coding RNA or discontinued
# Genes will be removed from count matrix
no_ncbi <- no_length[!(no_length %in% all_genes)]
no_ncbi

NCBI_GENES <- "../info/ref_genome/GRCh38/ncbi_GRCh38_p13-gene_lengths.tsv"
ncbi_genes <- read.table(NCBI_GENES, header = TRUE, row.names = 1,
                         sep = "\t", quote = "", stringsAsFactors = F)
gene_synonyms <- strsplit(ncbi_genes$gene_synonyms, ",")
names(gene_synonyms) <- rownames(ncbi_genes)

# Filter out genes with no synonyms
gene_lookup <- Filter(function(vec) !identical(vec, character(0)),
                      gene_synonyms)

missing_genes <- rownames(raw_rna)[!rownames(raw_rna) %in% rownames(ncbi_genes)]
# List of index for each missing gene
# Return list of index in gene_lookup for each missing gene
index_genes <- lapply(missing_genes, function(gene) {
  which(sapply(gene_lookup, function(list) gene %in% list))
})

names(index_genes) <- missing_genes

# Logical index of genes with no synonyms
discard_genes1 <- names(Filter(function(vec) length(vec) == 0, index_genes))
# Missing genes to be discarded

## CREATING GENE ANNOTATION TABLE
# Save lookup table for missing genes that are recognised synonyms
keep_genes <- names(Filter(function(vec) length(vec) == 1, index_genes))
lookup1 <- sapply(index_genes[keep_genes], function(x) names(x))
lookup1_lines <- paste(names(lookup1), lookup1, sep = "\t")

# How to deal with missing genes with multiple symbols matching
multiple_matches <- Filter(function(vec) length(vec) > 1, index_genes)

# Extract starting string character from all synonyms
# If all are identical then take first one else remove gene

# Convert list of indices to list of substrings
list_substring <- lapply(multiple_matches, function(vec) {
  regmatches(names(vec), regexpr("^[a-z|A-Z]+", names(vec)))
})
# Check if all strings in vector are identical
logi_idx1 <- sapply(list_substring, function(vec) {
  identical(vec, rep(vec[1], length(vec)))
})
# Genes to be discarded
discard_genes2 <- names(list_substring[!logi_idx1])

# Gene synonyms that will be assigned the first gene symbol
identical_symbols <- names(list_substring[logi_idx1])
lookup2 <- sapply(multiple_matches[identical_symbols],
                  function(vec) names(vec[1]))
lookup2_lines <- paste(names(lookup2), lookup2, sep = "\t")

# Genes to be discarded
writeLines(c(discard_genes1, discard_genes2), "dump/discarded_genes.tsv")
# Lookup table for gene synonyms
writeLines(c("gene_synonym\tgene_name", lookup1_lines, lookup2_lines),
           "dump/gene_synonym-lookup.tsv")

## COMPARE TO ENSEMBL
ENSEMBL <- "../info/ref_genome/GRCh38/ensembl_GRCh38_99-gene_lengths.tsv"
ensembl_lengths <- read.table(ENSEMBL, header = T, sep = "\t", row.names = 1)
sum(!(rownames(raw_rna) %in% rownames(ensembl_lengths)))
