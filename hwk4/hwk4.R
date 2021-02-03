## Imports
library(ggplot2)
library(data.table)

## Numerical perturbation analysis
dt <- data.table(x=sort(c(2, seq(2,4,length.out = 100))))
dt[, val := x^2+1][, type := "F0"]
Lfn <- function(lam, x) {(lam + 1)*x^2 - 6 * lam * x + 8 * lam + 1}
dts <- lapply(c(0.1, 1, 2, 3, 4), function(lam) { dt2 <- copy(dt); dt2[, val := Lfn(lam, x)][, type := sprintf("L(%s)", lam)]; dt2 })
dt <- rbindlist(c(list(dt), dts))
ggplot(dt, aes(x=x, y=val, col=type)) + geom_point() + geom_line() + geom_vline(xintercept = 2) + geom_vline(xintercept = 4)
