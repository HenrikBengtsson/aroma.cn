setConstructorS3("PairedPSCBSFileSet", function(...) {
  extend(GenericDataFileSet(...), "PairedPSCBSFileSet")
})

# setMethodS3("byPath", "PairedPSCBSFileSet", function(static, ..., pattern=".*(,PairedPSCBS[.]xdr|[.]RData)$") {
#   # Drop argument 'chipType'
#   args <- list(static=static, ..., pattern=pattern)
#   excl <- which("chipType" == names(args))
#   if (length(excl) > 0) {
#     args <- args[-excl]
#   }
#
#   # Call the "next" method
#   args <- c(list("byPath"), args)
#   do.call("NextMethod", args)
# }, static=TRUE)


#setMethodS3("findByName", "PairedPSCBSFileSet", function(static, ..., chipType=NULL, paths="pscbsData/") {
#  # Argument 'chipType':
#  if (!is.null(chipType)) {
#    chipType <- Arguments$getCharacter(chipType)
#  }
#
#  # Argument 'paths':
#  if (is.null(paths)) {
#    paths <- eval(formals(findByName.PairedPSCBSFileSet)[["paths"]])
#  }
#
#  # Drop argument 'subdirs'
#  args <- list(...)
#  excl <- which("subdirs" == names(args))
#  if (length(excl) > 0) {
#    args <- args[-excl]
#  }
#  args <- c(list(static=static), args, subdirs=chipType, paths=paths)
#
#  # Call the "next" method
#  args <- c(list("findByName"), args)
#  do.call("NextMethod", args)
#}, static=TRUE, protected=TRUE)

setMethodS3("byPath", "PairedPSCBSFileSet", function(static, ..., pattern=".*(,PairedPSCBS[.]xdr|[.]RData)$") {
  NextMethod("byPath")
}, static=TRUE)


setMethodS3("findByName", "PairedPSCBSFileSet", function(static, ..., chipType, paths="pscbsData/") {
  # Argument 'chipType':
  chipType <- Arguments$getCharacter(chipType)

  # Argument 'paths':
  if (!is.null(paths)) paths <- Arguments$getCharacters(paths)

  NextMethod("findByName", subdirs=chipType, paths=paths)
}, static=TRUE, protected=TRUE)


setMethodS3("byName", "PairedPSCBSFileSet", function(static, name, ...) {
 path <- findByName(static, name=name, ...)
 if (is.null(path)) {
   throw("Failed to locate ", class(static)[1L])
 }

  # Call byPath without argument 'chipType'
  args <- list(static, path=path, ...)
  excl <- which("chipType" == names(args))
  if (length(excl) > 0) args <- args[-excl]
  do.call(byPath, args=args)
})


setMethodS3("getPlatform", "PairedPSCBSFileSet", function(this, ...) {
  res <- this$platform
  if (is.null(res)) {
    res <- NA
  }
  res
})


setMethodS3("getDefaultFullName", "PairedPSCBSFileSet", function(this, ...) {
  res <- getPath(this)
  res <- getParent(res)
  res <- basename(res)
  res
}, protected=TRUE)


setMethodS3("getChipType", "PairedPSCBSFileSet", function(this, ...) {
  res <- this$chipType
  if (is.null(res)) {
    path <- getPath(this)
    res <- basename(path)
  }
  res
})
