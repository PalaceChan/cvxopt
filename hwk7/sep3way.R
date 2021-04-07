X <- as.matrix(fread("~/development/cvxopt/hwk7/sep3way_data/X.csv"))
Y <- as.matrix(fread("~/development/cvxopt/hwk7/sep3way_data/Y.csv"))
Z <- as.matrix(fread("~/development/cvxopt/hwk7/sep3way_data/Z.csv"))

a1 <- c(0.28329478, 0.16652294); a2 <- c(-0.20192124, 0.18983291); a3 <- c(-0.10444152, -0.30423805)
b1 <- 0.4585765; b2 <- -0.00016082; b3 <- 0.58926234

plot(X[1,], X[2,], col="red", xlim=c(-7,7), ylim=c(-7,7))
points(Y[1,], Y[2,], col="green")
points(Z[1,], Z[2,], col="blue")
abline(-a1[1]/a1[2], b1/a1[2])
abline(-a2[1]/a2[2], b2/a2[2])
abline(-a3[1]/a3[2], b3/a3[2])


