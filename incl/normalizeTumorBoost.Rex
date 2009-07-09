# Load data
pathname <- system.file("exData/TumorBoost,fracB,exampleData.Rbin", package="aroma.cn")
data <- loadObject(pathname)
attachLocally(data)
pos <- position/1e6

layout(matrix(1:4, ncol=1))
par(mar=c(2.5,4,0.5,1)+0.1)
ylim <- c(-0.05, 1.05)
col <- (data[,"genotypeN"] == 1/2) + 1L
plot(pos, betaN, col=col, ylim=ylim)
plot(pos, betaT, col=col, ylim=ylim)

betaTN <- normalizeTumorBoost(betaT=betaT, betaN=betaN, muN=genotypeN)
plot(pos, betaTN, col=col, ylim=ylim)

betaTN2 <- normalizeTumorBoost(betaT=betaT, betaN=betaN)
plot(pos, betaTN2, col=col, ylim=ylim)
