if (require("psCBS")) {

# Load example ASCN data
pathname <- system.file("data-ex/pcCBS,exData,chr01.Rbin", package="psCBS");
data <- R.utils::loadObject(pathname);
R.oo::attachLocally(data);
x <- position;

# PSCBS segmentation
fit <- segmentByPairedPSCBS(CT, betaT=betaT, betaN=betaN, x=x, verbose=-10); 

# Normalize
fitN <- normalizeBAFsByRegions(fit, verbose=-10);

par(ask=TRUE);

subplots(10, ncol=2, byrow=FALSE);
par(mar=c(1,3.5,1,0.5)+1);
plot(fit, add=TRUE);
plot(fitN, add=TRUE);


subplots(4, ncol=2, byrow=TRUE);
par(mar=c(1,3.5,1,0.5)+1);
plotC1C2(fit, Clim=c(0,3));
linesC1C2(fit);
title(main="(C1,C2)");

plotC1C2(fitN, Clim=c(0,3));
linesC1C2(fitN);
title(main="(C1,C2) - adjusted");

plotC1C2(fit, col="gray", Clim=c(0,3));
linesC1C2(fit, col="gray");
pointsC1C2(fitN);
linesC1C2(fitN);
title(main="(C1,C2) - both");

} # if (require("psCBS"))
