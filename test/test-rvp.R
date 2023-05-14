library(testthat)
library(pryr)

source('../R/rvp.R')


X <- readRDS('../../data_simulator/data/batchqc/sizes/batchqc-2000.rds')
n <- ncol(X)
ncond <- n / 4
batch <- as.factor(rep(1:2, each = ncond * 2))
class <- rep(rep(LETTERS[1:2], each = ncond), 2)
metadata <- data.frame(batch, class)  # assign rownames below

rvp <- RVP(data.frame(t(X)), batch, class)

test_that('rvp gives correct value.', {
  expect_equal(rvp, 0.0442705)
})
