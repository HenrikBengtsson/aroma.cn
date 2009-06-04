###########################################################################/**
# @RdocClass TotalCnKernelSmoothing
#
# @title "The TotalCnKernelSmoothing class"
#
# \description{
#  @classhierarchy
#
# }
# 
# @synopsis
#
# \arguments{
#  \item{...}{Arguments passed to @see "TotalCnSmoothing".}
#  \item{kernel}{A @character string specifying the type of kernel
#     to be used.}
#  \item{censorH}{}
#  \item{robust}{}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
#*/########################################################################### 
setConstructorS3("TotalCnKernelSmoothing", function( ..., kernel=c("gaussian", "uniform"), censorH=3, robust=FALSE) {
  # Argument 'kernel':
  kernel <- match.arg(kernel);

  # Arguments 'censorH':
  censorH <- Arguments$getDouble(censorH, range=c(0,Inf));

  # Arguments 'robust':
  robust <- Arguments$getLogical(robust);

  extend(TotalCnSmoothing(...), "TotalCnKernelSmoothing",
    .kernel = kernel,
    .censorH = censorH,
    .robust = robust
  );
})


setMethodS3("getParameters", "TotalCnKernelSmoothing", function(this, ...) {
  params <- NextMethod("getParameters", this, ...);
  params$kernel <- this$.kernel;
  params$censorH <- this$.censorH;
  params$robust <- this$.robust;
  params;
}, private=TRUE);


setMethodS3("getAsteriskTags", "TotalCnKernelSmoothing", function(this, collapse=NULL, ...) {
  tags <- NextMethod("getAsteriskTags", this, collapse=NULL, ...);

  # Add class-specific tags

  params <- getParameters(this);
  # "Parameter" 'kernel'
  kernelTag <- params$kernel;
  tags <- c(tags, kernelTag);

  # Parameter 'robust'
  if (params$robust)
    tags <- c(tags, "robust");

  # Collapsed or split?
  if (!is.null(collapse)) {
    tags <- paste(tags, collapse=collapse);
  } 

  tags;
}, protected=TRUE) 


setMethodS3("smoothRawCopyNumbers", "TotalCnKernelSmoothing", function(this, rawCNs, target, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Smoothing one set of copy numbers");
  verbose && print(verbose, rawCNs);

  # Setting up arguments
  params <- getParameters(this);
  args <- c(list(xOut=target$xOut), params, h=params$bandwidth, ...);

  # Keep only known arguments
  knownArguments <- names(formals(colKernelSmoothing.matrix));
  keep <- is.element(names(args), knownArguments);
  args <- args[keep];
  
  args <- c(list(rawCNs), args);

  verbose && cat(verbose, "Calling kernelSmoothing() with arguments:");
  verbose && str(verbose, args);
  args$verbose <- less(verbose, 20);
  smoothCNs <- do.call("kernelSmoothing", args=args);

  verbose && exit(verbose);

  smoothCNs;
}, protected=TRUE)



############################################################################
# HISTORY:
# 2009-02-08
# o Created from TotalCnSmoothing.R.
# 2009-01-26
# o Adopted to the new AromaUnitTotalCnBinarySet.
# o Added Rdoc comments.
# 2008-05-23
# o Created.
############################################################################
