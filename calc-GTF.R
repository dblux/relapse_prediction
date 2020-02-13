RPATH <- "../info/ref_genome/GRCh38/ensembl-Homo_sapiens.GRCh38.99.gtf"

gtf <- read.table(RPATH, sep = "\t")
gtf[1:3,]
