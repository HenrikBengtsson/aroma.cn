if (require("psCBS")) {

# Load example ASCN data
pathname <- system.file("data-ex/pcCBS,exData,chr01.Rbin", package="psCBS");
# pathname <- "wholeGenomeData/GSE20584,ACC,ra,-XY,BPN,-XY,AVG,FLN,-XY/GenomeWideSNP_6/GSE20584,GSM517071vGSM517072,ACC,ra,-XY,BPN,-XY,AVG,FLN,-XY.xdr";
data <- R.utils::loadObject(pathname);
R.oo::attachLocally(data);
x <- position;

# AD HOC: Robustification
CT[CT < 0] <- 0;
CT[CT > 30] <- 30;

# PSCBS segmentation
cc <- 1;
ww <- which(data$chromosome == cc);

fit <- segmentByPairedPSCBS(CT[ww], betaT=betaT[ww], betaN=betaN[ww], x=x[ww], chromosome=cc, verbose=-10);
fit <- postsegmentTCN(fit);

## Normalize
fitN <- normalizeBAFsByRegions(fit, verbose=-10);

# Call allelic balance
fitC <- callAllelicBalanceByBAFs(fitN);
ww <- which(fitC$output$ab.call);
fitC$output[ww, "dh.mean"] <- 0;

# Orthogonalize
fitO <- orthogonalizeC1C2(fitC);

# deShear
fitD <- deShearC1C2(fitC, verbose=verbose);

devSet("tracks");
subplots(10, ncol=2, byrow=FALSE);
par(mar=c(1,3.5,1,0.5)+1);
plot(fit, add=TRUE);
plot(fitD, add=TRUE);

devSet("C1C2");
Clim <- c(0,4);
subplots(4, ncol=2, byrow=TRUE);
par(mar=c(1,3.5,1,0.5)+1);
plotC1C2(fitC, Clim=Clim);
linesC1C2(fitC);
title(main="(C1,C2) - adjusted");

plotC1C2(fitO, Clim=Clim);
linesC1C2(fitO);
title(main="(C1,C2) - orthogonalized");

plotC1C2(fitD, Clim=Clim);
linesC1C2(fitD);
title(main="(C1,C2) - deSheared");

plotC1C2(fitO, col="gray", Clim=Clim);
linesC1C2(fitO, col="gray");
pointsC1C2(fitD);
linesC1C2(fitD);
title(main="(C1,C2) - orthogonalized vs deSheared");

} # if (require("psCBS"))
