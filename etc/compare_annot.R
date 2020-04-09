ANNOT_RPATH1 <- "../info/microarray/HG-U133A/annot_entrez-GPL96.tsv"
ANNOT_RPATH2 <- "../info/microarray/HG-U133_Plus_2/annot_entrez-GPL570.tsv"

annot1 <- read.table(ANNOT_RPATH1, sep = "\t", header = T, row.names = 1)
annot2 <- read.table(ANNOT_RPATH2, sep = "\t", header = T, row.names = 1)

common_ps <- intersect(rownames(annot1), rownames(annot2))
head(common_ps)
head(annot1[common_ps,, drop = F])
head(annot2[common_ps,, drop = F])

subset_annot1 <- annot1[common_ps,]
subset_annot2 <- annot2[common_ps,]

temp <- mapply(function(a,b) if(a == b) print(c(a,b)),
       subset_annot1, subset_annot2)

