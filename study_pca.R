x <- c(rnorm(1000,-5,1), rnorm(1000,6,1)) 
y <- c(rnorm(1000,-5,1), rnorm(1000,6,1))

par(mfrow = c(1,2))
plot(x, y, xlim = c(-20,20), ylim = c(-20,20))

arr <- cbind(x,y)
pca_obj <- prcomp(arr)
plot(pca_obj$x, xlim = c(-10,10), ylim = c(-10,10))
pca_var <- pca_obj$sdev^2

pca_var/sum(pca_var)

old_percentage_var
pca_obj$rotation

M <- matrix(rnorm(50,0,1), 5, 10)
pca_obj <- prcomp(M, center = F, scale. = F)

centered_M <- sweep(M, 2, colMeans(M), "-")
svd_matrices <- svd(M)
M %*% svd_matrices$v

pca_obj$rotation
svd_matrices$v
diag_d <- diag(svd_matrices$d)

pca_obj$x
pca_obj$x[3,2]/
svd_matrices$u %*% svd_matrices$d

M %*% t(svd_matrices$v)
svd_matrices$u

M <- matrix(rnorm(50,0,1), 5, 10)
svd_matrices <- svd(M)
