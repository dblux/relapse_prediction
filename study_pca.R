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
