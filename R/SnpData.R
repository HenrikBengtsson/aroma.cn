setConstructorS3("SnpData", function(data=NULL, ...) {
  data <- unclass(data)
  extend(BasicObject(data), "SnpData")
})

setMethodS3("callGenotypes", "SnpData", function(this, ...) {
  obj <- asTotalFracBSnpData(this)
  res <- callGenotypes(obj, ...)
  asThis(this, res)
})
