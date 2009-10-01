# Simulate three shifted tracks of different lengths with same profiles
ns <- c(A=2, B=1, C=0.25)*1000;
xx <- lapply(ns, FUN=function(n) { seq(from=1, to=max(ns), length.out=n) });
zz <- mapply(seq(along=ns), ns, FUN=function(z,n) rep(z,n));

yy <- list(
  A = rnorm(ns["A"], mean=0, sd=0.5),
  B = rnorm(ns["B"], mean=5, sd=0.4),
  C = rnorm(ns["C"], mean=-5, sd=1.1)
);
yy <- lapply(yy, FUN=function(y) {
  n <- length(y);
  y[1:(n/2)] <- y[1:(n/2)] + 2;
  y[1:(n/4)] <- y[1:(n/4)] - 4;
  y;
});

# Shift all tracks toward the first track
yyN <- normalizeDifferencesToAverage(yy, baseline=1);

# The baseline channel is not changed
stopifnot(identical(yy[[1]], yyN[[1]]));

# Get the estimated parameters
fit <- attr(yyN, "fit");

# Plot the tracks
layout(matrix(1:2, ncol=1));
x <- unlist(xx);
col <- unlist(zz);
y <- unlist(yy);
yN <- unlist(yyN);
plot(x, y, col=col, ylim=c(-10,10));
plot(x, yN, col=col, ylim=c(-10,10));
