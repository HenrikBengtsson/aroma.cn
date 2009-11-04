subplots(3, ncol=1)
par(mar=c(2,4,4,1)+0.1)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# A bimodal distribution
# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
xAA <- rnorm(n=10000, mean=0, sd=0.1)
xBB <- rnorm(n=10000, mean=1, sd=0.1)
x <- c(xAA,xBB)
fit <- findPeaksAndValleys(x)
print(fit)
calls <- callNaiveGenotypes(x, cn=rep(1,length(x)), verbose=-20)
xc <- split(x, calls)
print(table(calls))
xx <- c(list(x),xc)
plotDensity(xx, adjust=1.5, lwd=2, col=seq(along=xx), main="(AA,BB)")
abline(v=fit$x)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# A trimodal distribution with missing values
# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
xAB <- rnorm(n=10000, mean=1/2, sd=0.1)
x <- c(xAA,xAB,xBB)
x[sample(length(x), size=0.05*length(x))] <- NA;
x[sample(length(x), size=0.01*length(x))] <- -Inf;
x[sample(length(x), size=0.01*length(x))] <- +Inf;
fit <- findPeaksAndValleys(x)
print(fit)
calls <- callNaiveGenotypes(x)
xc <- split(x, calls)
print(table(calls))
xx <- c(list(x),xc)
plotDensity(xx, adjust=1.5, lwd=2, col=seq(along=xx), main="(AA,AB,BB)")
abline(v=fit$x)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# A trimodal distribution with clear separation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
xAA <- rnorm(n=10000, mean=0, sd=0.02)
xAB <- rnorm(n=10000, mean=1/2, sd=0.02)
xBB <- rnorm(n=10000, mean=1, sd=0.02)
x <- c(xAA,xAB,xBB)
fit <- findPeaksAndValleys(x)
print(fit)
calls <- callNaiveGenotypes(x)
xc <- split(x, calls)
print(table(calls))
xx <- c(list(x),xc)
plotDensity(xx, adjust=1.5, lwd=2, col=seq(along=xx), main="(AA',AB',BB')")
abline(v=fit$x)
