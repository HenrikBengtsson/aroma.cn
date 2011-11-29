
if (require("PSCBS")) {

# Load example ASCN data
pathname <- system.file("data-ex/PairedPSCBS,exData,chr01.Rbin", package="PSCBS");
data <- R.utils::loadObject(pathname);

# PSCBS segmentation
fit <- segmentByPairedPSCBS(data, verbose=-10); 

# Normalize
fitN <- normalizeBAFsByRegions(fit, verbose=-10);

set.seed(0xBEEF);
fitBList <- bootstrap(fitN, verbose=-10);

# Extract TCN and DH signals
x <- fitN$data$x;
CT <- fitN$data$CT;
betaTN <- as.vector(fitN$data$betaTN);
cn <- RawCopyNumbers(CT, x=x);
baf <- RawAlleleBFractions(betaTN, x=x);
dh <- extractRawMirroredAlleleBFractions(baf);
dh$y <- 2*dh$y;	
isHet <- (fitN$data$muN == 1/2);
dh$y[!isHet] <- NA;


devNew(width=9, height=5);
grid <- matrix(c(1,1,3,2,2,4), ncol=3, byrow=TRUE);
layout(grid);
par(mar=c(1,3.5,1,0.5)+1);
plot(cn, pch=".", col="gray", ylab="TCN");
dd <- lapply(fitBList, FUN=function(f) { 
  drawLevels(f, col="#ffaaaa", what="tcn", xScale=getXScale(cn));
});
drawLevels(fitN, what="tcn", xScale=getXScale(dh));

plot(dh, pch=".", col="gray", ylim=c(0,1), ylab="DH");
dd <- lapply(fitBList, FUN=function(f) { 
  drawLevels(f, col="#ffaaaa", what="dh", xScale=getXScale(dh));
});
drawLevels(fitN, what="dh", xScale=getXScale(dh));

plotC1C2(fitN, Clim=c(0,3));
dd <- lapply(fitBList, FUN=function(f) { 
  pointsC1C2(f, col="gray");
});
pointsC1C2(fitN);
linesC1C2(fitN);


} # if (require("PSCBS"))
