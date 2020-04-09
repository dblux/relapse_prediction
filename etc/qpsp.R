data(cv)

wm <- apply(RC, 2, gfs)
qpsp_mat <- qpsp(wm, cv)
qpsp_mat

library(NetProt)
library(genefilter)
# Import CORUM df
raw_corum <- read.table("../info/CORUM/entrezId1.txt",
                        sep = "\t", header = T, row.names = 1, stringsAsFactors = F)
# Only human complexes
human_corum <- raw_corum[raw_corum$Organism == "Human",]
list_corum <- strsplit(human_corum[,2], ';')
names(list_corum) <- rownames(human_corum)
