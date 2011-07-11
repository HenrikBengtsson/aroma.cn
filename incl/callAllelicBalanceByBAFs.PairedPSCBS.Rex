if (require("PSCBS")) {

# Load example ASCN data
pathname <- system.file("data-ex/pcCBS,exData,chr01.Rbin", package="PSCBS");
data <- R.utils::loadObject(pathname);
R.oo::attachLocally(data);
x <- position;

# PSCBS segmentation
fit <- segmentByPairedPSCBS(CT, betaT=betaT, betaN=betaN, x=x, verbose=-10); 

# Normalize
fitN <- normalizeBAFsByRegions(fit, verbose=-10);

fitN2 <- callAllelicBalanceByBAFs(fitN, verbose=-10);
print(fitN2);

} # if (require("PSCBS")) 
