library(rgl)

### INVESTIGATING DECISION RULES ###
# Set A of 3 features: ERM1, l2norm_d0_d8, d0_normal_proj
# Set B of 3 features: ERM1, l2norm_d0_d8, angle_d0_d8

# Import data -------------------------------------------------------------
features_1 <- read.table("dump/features-cs_quantile.tsv",
                         sep = "\t", header = T, row.names = 1)
# # EDA
# colnames(features_1)
# plot(features_1[,c(1,7,9,11,15)], col = (features_1[,"label"] + 1))

X <- features_1[,c(1,7,9,11)]
print(colnames(X))
y <- features_1$label

# plot3DScatter(X, features_1$label+1)

# Gaussian Naive Bayes ----------------------------------------------------
# Split training set for calculation of parameters
list_X <- split.data.frame(X, y)
# MOM estimate of mean
X_mu <- sapply(list_X, colMeans)
# Biased estimate of standard deviation
X_sigma <- sapply(list_X, apply, 2, sd)

#' # Calculate likelihood ratio for each feature of each sample
#' #' @return vector of likelihood ratios for every feature
#' calcLR <- function(x_vec) {
#'   #' Calculates likelihood ratio for single feature
#'   calcSingleLR <- function(idx) {
#'     x <- x_vec[idx]
#'     dnorm(x, X_mu[idx,2], X_sigma[idx,2])/
#'       dnorm(x, X_mu[idx,1], X_sigma[idx,1])
#'   }
#'   feature_names <- names(x_vec)
#'   sapply(feature_names, calcSingleLR, USE.NAMES = F)
#' }

# Predict likelihood ratios
print(colnames(X_mu)[2]) # Column name
likelihood_ratios <- t(apply(X, 1, calcLR))
X_y <- cbind(likelihood_ratios, y)
X_y
write.table(X_y,
            "dump/global_cs_quantile-lr.tsv",
            quote = F, sep = "\t")

# Plot
par(mfrow=c(5,1))
par(mar=rep(1,4))
for(i in 1:5) {
  plot(likelihood_ratios[,i], col = y+1)  
}

# Experiment with averages of LR
avg_lr <- rowMeans(likelihood_ratios[,1:2])
plot(avg_lr, col = y+1)
par(mfrow=c(1,1))

plot(data.frame(likelihood_ratios[,1:4]), col = y+1)

beta <- c(0.6,0.2,0.1,0.1)
weighted_lr <- sweep(likelihood_ratios[,1:4], 2, beta, "*")
plot(rowSums(weighted_lr), col = y+1)

X_y[,5][X_y[,5] == 0] <- "remission"
X_y[,5][X_y[,5] == 1] <- "relapse"
model.matrix(X_y[,1:4])

# Logistic regression to assign betas
