library(data.table)
library(ggplot2)

n <- 20
M <- 25
K <- 100
D <- -1.7712

y <- fread("~/development/cvxopt/hwk5/y.csv")$V1
X <- t(as.matrix(fread("~/development/cvxopt/hwk5/X.csv")))
c_true <- fread("~/development/cvxopt/hwk5/c_true.csv")$V1

data <- as.data.frame(cbind(X[1:25,], y))
c_ls <- lm("y ~ . - 1", data)$coefficients
r_ls <- sqrt(sum((c_true - c_ls)^2)) / sqrt(sum(c_true^2))
