set.seed(0xBEEF)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Simulate data from the track
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Build up CN profile generation by generation, via aroma.cn::cnr()
pT <- cnr(1,2000, mean=2) +
      cnr(1000,1500) + 
      cnr(1000,1250) +
      cnr(1650,1800) + 
      cnr(200,300) - cnr(650,800)
print(pT)

# Simulate CN signals from this profile
cn <- simulateRawCopyNumbers(pT, n=2000, sd=1/2)
print(cn)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Segment
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
data <- as.data.frame(cn)
fit <- segmentByCBS(data)
plotTracks(fit)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# PruneCNA
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fitList <- pruneCNA(fit)
printf("Number of generations identified: %d\n", nbrOfGenerations(fitList))
plotTracks(fitList)
