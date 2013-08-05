.setupAromaCn <- function(pkg, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Patches
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # None at the moment.

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Bioconductor package aroma.light
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # require("aroma.light") - install if missing
  aroma.core:::.requireBiocPackage("aroma.light", neededBy=getName(pkg));
} # .setupAromaCn()


############################################################################
# HISTORY:
# 2013-08-04
# o Created from ditto for aroma.affymetrix.
############################################################################
