
if (Sys.getenv("_R_CHECK_FULL_") != "" && require("PSCBS")) {

# Load example ASCN data
data <- PSCBS::exampleData("paired.chr01");
R.oo::attachLocally(data);

# AD HOC: Robustification
CT[CT < 0] <- 0;
CT[CT > 30] <- 30;

# PSCBS segmentation
cc <- 1;
ww <- which(data$chromosome == cc);

fit <- segmentByPairedPSCBS(CT[ww], betaT=betaT[ww], betaN=betaN[ww],
                            x=x[ww], chromosome=cc, verbose=-10);
fit <- postsegmentTCN(fit);

## Normalize
fitN <- normalizeBAFsByRegions(fit, verbose=-10);

# Call allelic balance
fitC <- callAllelicBalanceByBAFs(fitN);
ww <- which(fitC$output$abCall);
fitC$output[ww, "dhMean"] <- 0;

# deShear
fitD <- deShearC1C2(fitC, verbose=-10);

devSet("tracks");
subplots(2*3, ncol=2, byrow=FALSE);
par(mar=c(1,3.5,1,0.5)+1);
plotTracks(fit, subplots=FALSE);
plotTracks(fitD, subplots=FALSE);

devSet("C1C2");
Clim <- c(0,4);
plotC1C2(fitC, Clim=Clim);
linesC1C2(fitC);
title(main="(C1,C2) - adjusted");

devSet("C1C2,deShear");
Clim <- c(0,4);
plotC1C2(fitD, Clim=Clim);
linesC1C2(fitD);
title(main="(C1,C2) - deSheared");

} # if (require("PSCBS"))
