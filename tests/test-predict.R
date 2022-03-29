library(testthat)
library(magrittr)

source("../R/calc.R")
source("../R/misc.R")
source("../R/normalise.R")
source("../R/plot.R")
source("../R/predict.R")
source("../R/subset.R")
source("../R/utils.R")


### IMPORT DATA ### 
METADATA_SID <- "../data/GSE67684/processed/metadata/sid-metadata_v2.tsv"
METADATA_PID <- "../data/GSE67684/processed/metadata/pid-metadata_v7.tsv"

metadata_sid <- read.table(METADATA_SID, sep = "\t")
metadata_pid <- read.table(METADATA_PID, sep = "\t", row.names = 1, quote = '"')
metadata_pid$label <- as.factor(metadata_pid$label)

# Removed outliers, patients with timepoints from different batches and batch 5
SUBSET_RPATH <- "../data/GSE67684/processed/subset_yeoh.tsv"
raw_yeoh <- read.table(SUBSET_RPATH, sep = "\t")
# SCALE->REMOVE->FILTER->LOG
scaled_yeoh <- normaliseMeanScaling(raw_yeoh)
selected_yeoh <- removeProbesets(scaled_yeoh)
yeoh <- log2_transform(filterProbesets(selected_yeoh, 0.7, metadata_sid))

yeoh_allps <- log2_transform(scaled_yeoh)
yeoh_unfltr <- log2_transform(selected_yeoh)

### GLOBAL VARIABLES
metadata <- metadata_sid[colnames(yeoh),]
others <- yeoh[, metadata$subtype == "Others"]
others_normal <- yeoh[, metadata$subtype %in% c("Others", "Normal")]

# Define train/test split
sid_mrd_na <- rownames(metadata_pid)[is.na(metadata_pid$d33_mrd)] %>%
  rep(each = 2) %>%
  paste0(c("_D0", "_D8"))

sid_alltrain_local <- rownames(metadata)[
  !(metadata$subtype %in% c("Hypodiploid", "Normal")) &
  !(rownames(metadata) %in% sid_mrd_na)
]
sid_alltrain <- rownames(metadata)[
  !(metadata$subtype %in% c("Hypodiploid", "Hyperdiploid", "Others", "Normal")) &
  !(rownames(metadata) %in% sid_mrd_na)
]
sid_train <- rownames(metadata)[
  metadata$batch_info %in% 1:7 &
  !(metadata$subtype %in% c("Hypodiploid", "Hyperdiploid", "Others", "Normal")) &
  !(rownames(metadata) %in% sid_mrd_na)
]
sid_test <- rownames(metadata)[
  metadata$batch_info %in% 8:10 &
  !(metadata$subtype %in% c("Hypodiploid", "Hyperdiploid", "Others", "Normal")) &
  !(rownames(metadata) %in% sid_mrd_na)
]
sid_remission <- rownames(metadata)[metadata$label == 0]
sid_normal <- paste0("N0", c(1,2,4))

batch_ps <- identify_batch_features(yeoh, metadata, method = 'aov')

X_normal <- yeoh[, sid_normal]
subtype <- "BCR-ABL"
sid_subtype <- sid_alltrain_local[
  metadata_sid[sid_alltrain_local, "subtype"] == subtype
]
X_subtype <- yeoh[, sid_subtype]

print("Testing predict_pipeline...")
test_that("predict_pipeline estimates the correct direction for all features", {
  prediction <- predict_pipeline(
    X_subtype, X_normal, metadata, metadata_pid, batch_ps
  )
  p_sub <- prediction$p_remission_xi %>%
    cbind(label = metadata_pid[rownames(.), "label"]) %>%
    subset(label == 0)
  X_y_sub <- subset(prediction$X_y, label == 0)
  
  expect_equal(
    rank(p_sub$p_erm1_ratio2),
    rank(X_y_sub$erm1_ratio2) 
  )
  expect_equal(
    rank(p_sub$p_l2norm_ratio2),
    rank(X_y_sub$l2norm_ratio2) 
  )
  expect_equal(
    rank(-p_sub$p_angle_d0d8_d0normal),
    rank(X_y_sub$angle_d0d8_d0normal) 
  )
  expect_equal(
    rank(-p_sub$p_log_mrd_d33),
    rank(X_y_sub$log_mrd_d33)
  )
})
