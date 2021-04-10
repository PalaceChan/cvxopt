.libPaths(c("/usr/lib/R/library", "~/rlibs"))
library(data.table)
library(ggplot2)

U <- as.matrix(fread("~/development/cvxopt/hwk7/sphere_fit_data/U.csv"))
plot(U[1,], U[2,], col='blue')

r <- 1.30580895
x <- c(-2.60350921,  6.49561257)

O <- matrix(sapply(seq(0, 2*pi, length.out = 500), function(i) c(r*cos(i), r*sin(i)) + x), nrow=2, byrow = FALSE)
points(O[1,], O[2,])

