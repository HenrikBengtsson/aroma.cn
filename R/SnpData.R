setConstructorS3("SnpData", function(data=NULL, ...) {
  data <- unclass(data);
  extend(BasicObject(data), "SnpData");
})


setMethodS3("plot", "SnpData", function(x, ...) {
  NextMethod("plot", x, ...);
})

setMethodS3("callGenotypes", "SnpData", function(this, ...) {
  obj <- asTotalFracBSnpData(this);
  res <- callGenotypes(obj, ...);
  asThis(this, res);
})



############################################################################
# HISTORY:
# 2009-03-31
# o Added pairedBoost() for CartesianSnpData.
# 2009-03-30
# o Created.
############################################################################
