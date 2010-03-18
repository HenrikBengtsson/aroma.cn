# Load data
pathname <- system.file("exData/TumorBoost,fracB,exampleData.Rbin", package="aroma.cn")
data <- loadObject(pathname)

# Drop loci with missing values
data <- na.omit(data)

attachLocally(data)
pos <- position/1e6

# Call naive genotypes
muN <- callNaiveGenotypes(betaN)

# Genotype classes
isAA <- (muN == 0)
isAB <- (muN == 1/2)
isBB <- (muN == 1)

# Sanity checks
stopifnot(all(muN[isAA] == 0))
stopifnot(all(muN[isAB] == 1/2))
stopifnot(all(muN[isBB] == 1))

# TumorBoost normalization with different flavors
betaTNs <- list()
for (flavor in c("v1", "v2", "v3", "v4")) {
  betaTN <- normalizeTumorBoost(betaT=betaT, betaN=betaN, flavor=flavor)

  # Assert that no non-finite values are introduced
  stopifnot(all(is.finite(betaTN)))

  # Assert that nothing is flipped
  stopifnot(all(betaTN[isAA] < 1/2))
  stopifnot(all(betaTN[isBB] > 1/2))

  betaTNs[[flavor]] <- betaTN
}

# Plot
layout(matrix(1:4, ncol=1))
par(mar=c(2.5,4,0.5,1)+0.1)
ylim <- c(-0.05, 1.05)
col <- rep("#999999", length(muN))
col[muN == 1/2] <- "#000000"
for (flavor in names(betaTNs)) {
  betaTN <- betaTNs[[flavor]]
  ylab <- sprintf("betaTN[%s]", flavor)
  plot(pos, betaTN, col=col, ylim=ylim, ylab=ylab)
}



