.setupAromaCn <- function(pkg, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Patches
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # None at the moment.

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Bioconductor package aroma.light
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # require("aroma.light") - install if missing
  ns <- getNamespace("aroma.core")
  .requireBiocPackage <- get(".requireBiocPackage", envir=ns)
  .requireBiocPackage("aroma.light", neededBy=getName(pkg))
} # .setupAromaCn()
