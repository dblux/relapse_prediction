# LOADING MTX -------------------------------------------------------------
library(Matrix)

t293_mat <- readMM("data/hie_2019/293t_jurkat/293t/matrix.mtx")
jurkat_mat <- readMM("data/hie_2019/293t_jurkat/jurkat/matrimtx")
mix50_mat <- readMM("data/hie_2019/293t_jurkat/jurkat_293t_50_50/matrix.mtx")
# All matrices have the same rownames
row_genes <- read.table("data/hie_2019/293t_jurkat/293t/genes.tsv",
                        sep = "\t", header = F)
# L2 Normalise
prcomp(t(t293_mat))



# PANCREAS ----------------------------------------------------------------
# celseq1_raw <- read.table("data/scanorama/data/pancreas/pancreas_multi_celseq_expression_matrix.txt.gz",
#                           sep = "\t", header = T, row.names = 1)
# celseq2_raw <- read.table("data/scanorama/data/pancreas/pancreas_multi_celseq2_expression_matrix.txt.gz",
#                           sep = "\t", header = T, row.names = 1)
# indrop_raw <- read.table("data/scanorama/data/pancreas/pancreas_inDrop.txt.gz",
#                          sep = "\t", header = T, row.names = 1)
# smartseq2_raw <- read.table("data/scanorama/data/pancreas/pancreas_multi_smartseq2_expression_matrix.txt.gz",
#                             sep = "\t", header = T, row.names = 1)
# fluidigmc1_raw <- read.table("data/scanorama/data/pancreas/pancreas_multi_fluidigmc1_expression_matrix.txt.gz",
#                              sep = "\t", header = T, row.names = 1)
# pancreas_labels <- read.table("data/scanorama/data/cell_labels/pancreas_cluster.txt",
#                               sep = "\t", header = F)
# head(pancreas_labels)
# colnames(celseq1_raw)
