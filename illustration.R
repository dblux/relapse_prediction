library(RColorBrewer)

A1 <- data.frame(x = rnorm(100, 5), y = rnorm(100,10))
B1 <- data.frame(x = rnorm(100, 15), y = rnorm(100,5))
C1 <- data.frame(x = rnorm(100, 15), y = rnorm(100,15))

B2 <- data.frame(x = rnorm(100, 15), y = rnorm(100,5))
C2 <- data.frame(x = rnorm(100, 15), y = rnorm(100,15))
D2 <- data.frame(x = rnorm(100, 25), y = rnorm(100,10))

B3 <- data.frame(x = rnorm(100, 15), y = rnorm(100,5))
C3 <- data.frame(x = rnorm(100, 15), y = rnorm(100,15))

df_combined <- rbind(A1,B1,C1,B2,C2,D2,B3,C3)
# Shift y-axis of batch 2
df_combined[301:600,2] <- df_combined[301:600,2] + 7
z <- c(rnorm(300, 5, 0.5), rnorm(300, 10, 0.5), rnorm(200, 15, 0.5))
df_combined <- cbind(df_combined, z)

shape_vec <- rep(c(21,22,23), c(300,300,200))
colour_code <- rep(c("lightblue","tomato3","darkolivegreen3",
                     "tomato3","darkolivegreen3","gold",
                     "tomato3","darkolivegreen3"), each = 100)

rgl.open()
rgl.bg(color="white")
rgl.viewpoint(zoom = 0.6)
# rgl.viewpoint(theta = 110, phi = 5, zoom = 0.8)
par3d(windowRect = c(50, 20, 500, 500))
# Plot of MILE dataset
with(df_combined, pch3d(x, z, y, bg = colour_code,
                   pch = shape_vec, cex = 0.3, lwd = 1.5))
axes3d(c('x-', 'y+', 'z+'),
       col = "gray8", labels = F, tick = F)
rgl.postscript("dump/bcm_illustration.pdf", "pdf")

# title3d(xlab = pc_labels[1], ylab = pc_labels[2],
#         zlab = pc_labels[3], col = "black")
# # Plot aspect ratios of axis according to variance
# do.call(aspect3d, ratio_list)

x <- rnorm(1000)
y <- rnorm(1000)
z <- x + y
plot(x,z)
