setMethodS3("callGenotypes", "RawAlleleBFractions", function(this, ...) {
  beta <- getSignals(this)
  keep <- which(is.finite(beta))
  beta <- beta[keep]
  data <- data.frame(total=2, fracB=beta)
  data <- as.matrix(data)
  data2 <- TotalFracBSnpData(data)
  data3 <- callGenotypes(data2)
  naValue <- NA_real_
  mu <- rep(naValue, times=nbrOfLoci(this))
  mu[keep] <- data3[,2]
  res <- clone(this)
  res$y <- mu
  res <- RawGenotypeCalls(res)
  res
}) # callGenotypes()


setMethodS3("normalizeTumorBoost", "RawAlleleBFractions", function(this, betaN, muN=callGenotypes(betaN), ...) {
  xBetaT <- getSignals(this)
  xBetaN <- getSignals(betaN)
  xMuN <- getSignals(muN)
  xBetaTN <- normalizeTumorBoost(betaT=xBetaT, betaN=xBetaN, muN=xMuN, ...)
  betaTN <- clone(this)
  betaTN$y <- xBetaTN
  betaTN
})
