# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# DEFUNCT
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# DEPRECATED
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethodS3("nbrOfFiles", "AbstractCurveNormalization", function(this, ...) {
  .Deprecated(msg = "nbrOfFiles(x) is deprecated. Use length(getInputDataSet(x)) instead.")
  ds <- getInputDataSet(this)
  length(ds)
}, protected=TRUE)


setMethodS3("nbrOfFiles", "PairedPscbsModel", function(this, ...) {
  .Deprecated(msg = "nbrOfFiles(x) is deprecated. Use length(getTumorDataSet(x)) instead.")
  dsT <- getTumorDataSet(this)
  length(dsT)
})


setMethodS3("nbrOfFiles", "TumorBoostNormalization", function(this, ...) {
  .Deprecated(msg = "nbrOfFiles(x) is deprecated. Use length(getInputDataSet(x)) instead.")
  ds <- getInputDataSet(this)
  length(ds)
})
