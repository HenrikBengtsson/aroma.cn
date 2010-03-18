# Load data
pathname <- system.file("exData/TumorBoost,fracB,exampleData.Rbin", package="aroma.cn")
data <- loadObject(pathname)
attachLocally(data)
pos <- position/1e6
muN <- genotypeN

layout(matrix(1:4, ncol=1))
par(mar=c(2.5,4,0.5,1)+0.1)
ylim <- c(-0.05, 1.05)
col <- rep("#999999", length(muN))
col[muN == 1/2] <- "#000000"

# Allele B fractions for the normal sample
plot(pos, betaN, col=col, ylim=ylim)

# Allele B fractions for the tumor sample
plot(pos, betaT, col=col, ylim=ylim)

# TumorBoost w/ naive genotype calls
betaTN <- normalizeTumorBoost(betaT=betaT, betaN=betaN)
plot(pos, betaTN, col=col, ylim=ylim)

# TumorBoost w/ external multi-sample genotype calls
betaTNx <- normalizeTumorBoost(betaT=betaT, betaN=betaN, muN=muN)
plot(pos, betaTNx, col=col, ylim=ylim)
