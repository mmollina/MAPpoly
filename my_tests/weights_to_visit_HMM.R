rm(list = ls())
v1 <- c(1,1,0,0)
v2 <- c(1,0,1,0)
pl <- length(v1)
a1 <- apply(combn(v1,pl/2), 2, sum)
a2 <- apply(combn(v2,pl/2), 2, sum)
X <- kronecker(a1, t(a2), "+")
y <- min(X):max(X)
I1 <- sapply(y, function(x) apply(X == x, 1, sum))
dimnames(I1) <- list(apply(combn(1:pl,pl/2), 2, paste0, collapse = ""), y)
I1

v1 <- c(1,1,0,0)
v2 <- c(1,0,0,1)
pl <- length(v1)
a1 <- apply(combn(v1,pl/2), 2, sum)
a2 <- apply(combn(v2,pl/2), 2, sum)
X <- kronecker(a1, t(a2), "+")
y <- min(X):max(X)
I2 <- sapply(y, function(x) apply(X == x, 1, sum))
dimnames(I2) <- list(apply(combn(1:pl,pl/2), 2, paste0, collapse = ""), y)
I1;I2

