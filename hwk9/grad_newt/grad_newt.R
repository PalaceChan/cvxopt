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
    sum(1/(1 - A %*% x)^2) * (t(A) %*% A) + diag(2*(x^2 + 1)/(1-x^2)^2)
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
ggplot(y, aes(x=i, y=abs(ans-v))) + geom_point() + geom_line()

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
## ggplot(y, aes(x=i, y=v)) + geom_point() + geom_line()
## ggplot(y, aes(x=i, y=abs(ans-v))) + geom_point() + geom_line()
