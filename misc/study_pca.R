source("../functions.R")

# Execute together with set.seed
set.seed(1)
x <- c(rnorm(1000,-5,1), rnorm(1000,6,1)) 
y <- c(rnorm(1000,-5,1), rnorm(1000,6,1))

par(mfrow = c(1,2))
plot(x, y, xlim = c(-20,20), ylim = c(-20,20))

arr <- cbind(x,y)
pca_obj <- prcomp(arr)
# Plot PCA coordinates
plot(pca_obj$x, xlim = c(-10,10), ylim = c(-10,10))

# Total variance before and after PCA remains the same
sum(apply(arr, 2, var))
sum(apply(pca_obj$x, 2, var))

# Eigenvalue represents the variance along each eigenvector
pca_obj$sdev^2 # Eigenvalue
apply(pca_obj$x, 2, var)

# prcomp mean-centers data by default
centered_arr <- scale(arr, scale = F)
# prcomp performs svd on original data
svd_1 <- svd(centered_arr)
# Rotation matrix is given by matrix v
svd_1$v
pca_obj$rotation
# Principal component scores is given by matrix u*d
head(svd_1$u %*% diag(svd_1$d))
head(pca_obj$x)
# Principal component scores can subsequently be given by x*v
head(centered_arr %*% svd_1$v)
# Eigenvalue can be calculated by d^2/(n-1)
pca_obj$sdev^2
svd_1$d^2/1999

# prcomp can be computed more efficiently by:
# Performing svd on t(x)*x
# t(x)*x = v*d^2*t(v)
# Calculate principal component scores using xv
svd_2 <- svd(t(centered_arr) %*% centered_arr)
# v
svd_2$u
# d^2
pca_obj$sdev^2
svd_2$d
# v
svd_2$v
pca_obj$rotation

# Eigenvectors are defined by the rotation matrix
# Eigenvectors all have unit l2norm even if original data is not scaled to have unit variance!
apply(pca_obj$rotation, 2, calc_l2norm)

