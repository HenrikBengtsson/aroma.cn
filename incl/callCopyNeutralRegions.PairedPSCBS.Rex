
if (require("PSCBS") && Sys.getenv("_R_CHECK_FULL_") != "") {

# Load example ASCN data
pathname <- system.file("data-ex/PairedPSCBS,exData,chr01.Rbin", package="PSCBS");
data <- R.utils::loadObject(pathname);

# PSCBS segmentation
fit <- segmentByPairedPSCBS(data, verbose=-10); 

fit <- bootstrapTCNandDHByRegion(fit, verbose=-10);

fitC <- callCopyNeutralRegions(fit, verbose=-10);
print(fitC);

# Normalize
fitN <- normalizeBAFsByRegions(fitC, verbose=-10);

} # if (require("PSCBS")) 
