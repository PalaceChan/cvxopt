.libPaths(c("/usr/lib/R/library", "~/rlibs"))
library(ggplot2)
library(data.table)
## maxIter <- 10
## n <- 3; m <- 5
maxIter <- 1000
n <- 100; m <- 200
set.seed(0)
A <- matrix(rnorm(m*n), nrow=m, ncol=n)
fwrite(A, "~/development/cvxopt/hwk9/grad_newt/A.csv", col.names = FALSE)
ans <- -121.44910007279553

isInDom <- function(A, x) {
    all(A %*% x <= 1) & all(abs(x) <= 1)
}

getObjV <- function(A, x) {
    -sum(log(1 - A %*% x)) - sum(log(1 - x^2))
}

getGrad <- function(A, x) {
    t(A) %*% (1/(1 - A %*% x)) + 2*x / (1-x^2)
}

getHess <- function(A, x) {
    ## sum(1/(1 - A %*% x)^2) * (t(A) %*% A) #+ diag(2*(x^2 + 1)/(1-x^2)^2)
    t(A) %*% diag(as.numeric(1 / (-A %*% x + 1)^2)) %*% A + diag(2*(x^2 + 1)/(1-x^2)^2)
}

getWolfHess <- function(A, x) { #debug hessian with wolfram alpha
    ## Atest <- matrix(c(c(1,2,3),c(4,5,6)), nrow = 2, byrow = TRUE)
    ## {{-16/(1 - 4 x - 5 y - 6 z)^2 - (1 - x - 2 y - 3 z)^(-2), -20/(1 - 4 x - 5 y - 6 z)^2 - 2/(1 - x - 2 y - 3 z)^2, -24/(1 - 4 x - 5 y - 6 z)^2 - 3/(1 - x - 2 y - 3 z)^2}, {-20/(1 - 4 x - 5 y - 6 z)^2 - 2/(1 - x - 2 y - 3 z)^2, -25/(1 - 4 x - 5 y - 6 z)^2 - 4/(1 - x - 2 y - 3 z)^2, -30/(1 - 4 x - 5 y - 6 z)^2 - 6/(1 - x - 2 y - 3 z)^2}, {-24/(1 - 4 x - 5 y - 6 z)^2 - 3/(1 - x - 2 y - 3 z)^2, -30/(1 - 4 x - 5 y - 6 z)^2 - 6/(1 - x - 2 y - 3 z)^2, -36/(1 - 4 x - 5 y - 6 z)^2 - 9/(1 - x - 2 y - 3 z)^2}}

    row1 <- c(-16/(1 - 4*x[1] - 5*x[2] - 6*x[3])^2 - (1 -x[1] - 2*x[2] - 3*x[3])^(-2), -20/(1 - 4*x[1] - 5*x[2] - 6*x[3])^2 - 2/(1 -x[1] - 2*x[2] - 3*x[3])^2, -24/(1 - 4*x[1] - 5*x[2] - 6*x[3])^2 - 3/(1 -x[1] - 2*x[2] - 3*x[3])^2)
    row2 <- c(-20/(1 - 4*x[1] - 5*x[2] - 6*x[3])^2 - 2/(1 -x[1] - 2*x[2] - 3*x[3])^2, -25/(1 - 4*x[1] - 5*x[2] - 6*x[3])^2 - 4/(1 -x[1] - 2*x[2] - 3*x[3])^2, -30/(1 - 4*x[1] - 5*x[2] - 6*x[3])^2 - 6/(1 -x[1] - 2*x[2] - 3*x[3])^2)
    row3 <- c(-24/(1 - 4*x[1] - 5*x[2] - 6*x[3])^2 - 3/(1 -x[1] - 2*x[2] - 3*x[3])^2, -30/(1 - 4*x[1] - 5*x[2] - 6*x[3])^2 - 6/(1 -x[1] - 2*x[2] - 3*x[3])^2, -36/(1 - 4*x[1] - 5*x[2] - 6*x[3])^2 - 9/(1 -x[1] - 2*x[2] - 3*x[3])^2)
    return(-matrix(c(row1, row2, row3), nrow = 3, byrow = TRUE))
}

bakTrak <- function(x, g, dir, alpha, beta) {
    s <- 1
    for (i in seq(0, maxIter)) {
        ## message(sprintf("backtracking iter %s with s=%s", i, s))
        if (i == maxIter -1) stop(sprintf("hit iter %s but s is still %s", i, s))

        if (! isInDom(A, x + s*dir)) {
            s <- beta*s
            ## message(sprintf("x=%s is OOD so shrinking", paste(x+s*dir, collapse=",")))
            next
        }

        fx0 <- getObjV(A, x)
        fx1 <- getObjV(A, x + s*dir)
        thr <- fx0 + alpha * s * (t(g) %*% dir)
        if (fx1 > thr) {
            s <- beta*s
            ## message(sprintf("fx0=%s fx1=%s thr=%s so shrinking", fx0, fx1, thr))
        }
        else {
            ## message(sprintf("fx0=%s fx1=%s thr=%s so done", fx0, fx1, thr))
            break
        }
    }

    s
}

runGradientDescent <- function() {
    x <- rep(0,n)
    eta <- 1e-3
    alpha <- 0.5
    beta <- 0.5
    objIdx <- 1
    objVals <- rep(0, maxIter)

    for (i in seq(0, maxIter)) {    
        g <- getGrad(A, x)
        gNorm <- sqrt(t(g) %*% g)
        if (gNorm < eta) break
        if (i == maxIter -1) stop(sprintf("descent: hit iter %s but gradient is still %s", i, gNorm))    

        s <- bakTrak(x, g, -g, alpha, beta)
        x <- x - s*g

        objIdx <- i
        objVal <- getObjV(A, x)
        objVals[objIdx] <- objVal
        ## message(sprintf("descent: iter %s s=%s objVal=%s", i, s, objVal))
    }

    y <- data.table(i=1:objIdx, v=objVals[1:objIdx])
}

y <- runGradientDescent()
## ggplot(y, aes(x=i, y=v)) + geom_point() + geom_line()
## ggplot(y, aes(x=i, y=abs(ans-v))) + geom_point() + geom_line()

runNewtonDescent <- function() {
    x <- rep(0,n)
    eta <- 1e-3
    alpha <- 0.5
    beta <- 0.5
    objIdx <- 1
    objVals <- rep(0, maxIter)

    for (i in seq(0, maxIter)) {    
        g <- getGrad(A, x)
        H <- getHess(A, x)
        xnt <- solve(-H, g)
        lam2 <- t(g) %*% solve(H) %*% g

        if (lam2/2 < eta) break
        if (i == maxIter -1) stop(sprintf("descent: hit iter %s but lam2 is still %s", i, lam2))    

        s <- bakTrak(x, g, xnt, alpha, beta)
        x <- x + s*xnt

        objIdx <- i
        objVal <- getObjV(A, x)
        objVals[objIdx] <- objVal
        ## message(sprintf("descent: iter %s s=%s objVal=%s", i, s, objVal))
    }

    y <- data.table(i=1:objIdx, v=objVals[1:objIdx])    
}

z <- runNewtonDescent()
## ggplot(z, aes(x=i, y=v)) + geom_point() + geom_line()
## ggplot(z, aes(x=i, y=abs(ans-v))) + geom_point() + geom_line()

# compare methods
ggplot(rbindlist(list(y[,.(i, v, mth="G")], z[,.(i, v, mth="N")], data.table(i=32:544, v=tail(z$v,1), mth="N")))[i < 10], aes(x=i, y=v, col=mth)) + geom_point() + geom_line()

