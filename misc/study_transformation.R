curve(0.5*x^2, add=T)
curve(1/(1+exp(-x)), add = T,  col = 2)
curve(1/(1+exp(-(0.5*x^2))), add = T, col = 3)
curve(1*x, add = T)

curve(x^2, xlim = c(-5,5), ylim = c(0,5))
curve(1/(1+exp(-x)), add = T,  col = 2)
curve(1/(1+exp(-(x^2))), add = T, col = 3)
curve((1/(1+exp(-x)))^2, add = T,  col = 4)
curve(1*x, add = T)

f <- function(x, y) {
  z = ((x^2)+(3*y^2)) * exp(-(x^2)-(y^2))
}

g <- function(x1,x2) x1 + x2
h <- function(x) x^2
hg <- function(x1,x2) h(g(x1,x2))

# Plot a 3D function surface plot
plot3d(g, xlim = c(-3, 3), ylim = c(-3, 3), zlim = c(-3, 3))
plot3d(hg, xlim = c(-3, 3), ylim = c(-3, 3), zlim = c(-3, 3))

# Add axes
axes_limit <- c(-3, 3)
rgl.lines(axes_limit, c(0, 0), c(0, 0), color = "black")
rgl.lines(c(0, 0), axes_limit, c(0, 0), color = "red")
rgl.lines(c(0, 0), c(0, 0), axes_limit, color = "green")