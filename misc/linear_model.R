# # Fitting linear models
# data(mtcars)
# head(mtcars)
# 
# car_size <- as.factor(mtcars$cyl)
# levels(car_size) <- c("small","medium","large")
# 
# model1 <- lm(mpg ~ mod_mat, mtcars)
# model1$coefficients
# abline(model1)
# 
# mod_mat <- model.matrix(~car_size + as.factor(mtcars$gear))
# 
# points(as.factor(car_size), mtcars$mpg)
# abline(model1)

# Simulate single gene data for different configurations
sample_maqc <- read.table("dump/sample_maqc_data.tsv",
                          sep = "\t", header = T, row.names = 1)
batch_info <- as.factor(rep(1:2, each = 10))
class_info <- as.factor(rep(rep(LETTERS[1:2], each = 5), 2))

# Linear modelling
i <- 5
y <- as.numeric(sample_maqc[i,])
fitted_model <- lm(y ~ batch_info + class_info)
# Subtract estimate of batch coefficient
y[11:20] <- y[11:20]- fitted_model$coefficients[2]

# Plot before and after correction
par(mfrow=c(2,1))
plot(as.numeric(sample_maqc[i,]),
     pch = rep(rep(21:22, each = 5), 2),
     cex = 2, bg = rep(2:3, each = 10), ann = FALSE)

plot(y,
     pch = rep(rep(21:22, each = 5), 2),
     cex = 2, bg = rep(2:3, each = 10), ann = FALSE)

# Plot all probesets
par(mfrow=c(4,2))
par(mar=c(1,1,1,1))
idx <- 1

for(i in idx:(idx+7)) {
  y <- as.numeric(sample_maqc[i,])
  plot(1:20, y, pch = rep(rep(21:22, each = 5), 2),
       cex = 2, bg = rep(2:3, each = 10), ann = FALSE)
}

rownames(sample_maqc)
