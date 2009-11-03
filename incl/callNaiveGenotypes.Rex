subplots(3, ncol=1)
par(mar=c(2,4,4,1)+0.1)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# A bimodal distribution
# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
x1 <- rnorm(n=10000, mean=0, sd=1)
x2 <- rnorm(n=10000, mean=4, sd=1)
x <- c(x1,x2)
fit <- findPeaksAndValleys(x)
calls <- callNaiveGenotypes(x, cn=rep(1,length(x)))
xc <- split(x, calls)
print(fit)
print(table(calls))
xx <- c(list(x),xc)
plotDensity(xx, adjust=1.5, lwd=2, col=seq(along=xx), main="c(x1,x2)")
abline(v=fit$x)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# A trimodal distribution
# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
x3 <- rnorm(n=10000, mean=8, sd=1)
x <- c(x1,x2,x3)
fit <- findPeaksAndValleys(x)
calls <- callNaiveGenotypes(x)
xc <- split(x, calls)
print(fit)
print(table(calls))
xx <- c(list(x),xc)
plotDensity(xx, adjust=1.5, lwd=2, col=seq(along=xx), main="c(x1,x2,x3)")
abline(v=fit$x)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# A trimodal distribution with clear separation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
x1b <- rnorm(n=10000, mean=0, sd=0.1)
x2b <- rnorm(n=10000, mean=4, sd=0.1)
x3b <- rnorm(n=10000, mean=8, sd=0.1)
x <- c(x1b,x2b,x3b)
fit <- findPeaksAndValleys(x, tol=0)
calls <- callNaiveGenotypes(x, tol=0)
xc <- split(x, calls)
print(fit)
print(table(calls))
xx <- c(list(x),xc)
plotDensity(xx, adjust=1.5, lwd=2, col=seq(along=xx), main="c(x1b,x2b,x3b)")
abline(v=fit$x)
