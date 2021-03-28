library(data.table)
require(grDevices)

m <- matrix(as.numeric(fread("~/development/cvxopt/hwk5/joint.csv", header=F)), nrow=100, byrow=T)
contour(x, x, m, nlevels = 20, col = terrain.colors(12)[8])
abline(0,-1)
