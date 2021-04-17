.libPaths(c("/usr/lib/R/library", "~/rlibs"))
library(data.table)
library(ggplot2)

v <- fread("~/development/cvxopt/hwk9/rank_aggr/rankings.csv", col.names = c("A","B"))
v <- rbind(data.table(val=v$A, sol="A"), data.table(val=v$B, sol="B"))
ggplot(v, aes(x=val, fill=sol, color=sol)) + geom_histogram(position = "identity", alpha = 0.5)

sum(v[sol == "B"]$val > 0.001) - sum(v[sol == "A"]$val > 0.001)



