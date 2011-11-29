if (require("PSCBS")) {

# Load example ASCN data
pathname <- system.file("data-ex/PairedPSCBS,exData,chr01.Rbin", package="PSCBS");
data <- R.utils::loadObject(pathname);

# PSCBS segmentation
fit <- segmentByPairedPSCBS(data, verbose=-10); 

# Normalize
fitN <- normalizeBAFsByRegions(fit, verbose=-10);

fitN2 <- callCopyNeutralRegions(fitN, verbose=-10);
print(fitN2);

} # if (require("PSCBS")) 
