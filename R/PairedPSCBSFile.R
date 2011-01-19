setConstructorS3("PairedPSCBSFile", function(...) {
  extend(GenericDataFile(...), "PairedPSCBSFile");
})

setMethodS3("loadObject", "PairedPSCBSFile", function(this, ...) {
  pathname <- getPathname(this);
  loadObject(pathname, ...);
})



#############################################################################
# HISTORY:
# 2011-01-18
# o Created.
#############################################################################
