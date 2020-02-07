RPATH1 <- "data/GSE67684/README/batch.info.txt"
RPATH2 <- "data/GSE67684/processed/rna_seq/rna_seq-pid.txt"

batch <- read.table(RPATH1, sep = "\t", header = TRUE)
existing_pid <- sub("_.*", "", batch$ID)

new_pid <- scan(RPATH2, what = "character", nlines = 1, sep = "\t")
new_pid1 <- sub("_.*", "", substring(new_pid[1:234], 4))

new_pid
new_pid[!(new_pid1 %in% existing_pid)]
