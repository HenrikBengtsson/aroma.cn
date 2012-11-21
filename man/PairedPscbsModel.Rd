%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  PairedPscbsModel.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{PairedPscbsModel}
\docType{class}
\alias{PairedPscbsModel}


 \title{The PairedPscbsModel class}

 \description{
  Package:  aroma.cn \cr
\bold{Class PairedPscbsModel}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.core]{ParametersInterface}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{PairedPscbsModel}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{PairedPscbsModel}\cr
extends \link[aroma.core]{ParametersInterface}\cr



  This class represents the Paired PSCBS method [1], which
  segments matched tumor-normal parental copy-number data into
  piecewise constant segments.
 }

 \usage{PairedPscbsModel(dsT=NULL, dsN=NULL, tags="*", ..., dropTcnOutliers=TRUE, gapMinLength=1e+06, seed=NULL)}

 \arguments{
   \item{dsT, dsN}{The tumor and the normal
     \code{\link[aroma.core]{AromaUnitPscnBinarySet}}.}
   \item{tags}{Tags added to the output data sets.}
   \item{...}{(Optional) Additional arguments passed to
     \code{\link[PSCBS]{segmentByPairedPSCBS}}.}
   \item{dropTcnOutliers}{If \code{\link[base:logical]{TRUE}}, then TCN outliers are dropped using
     \code{\link[PSCBS]{dropSegmentationOutliers}}.}
   \item{gapMinLength}{Genomic regions with no data points that are of
     this length and greater are considered to be "gaps" and are ignored
     in the segmentation.  If +\code{\link[base:is.finite]{Inf}}, no gaps are identified.}
   \item{seed}{An optional \code{\link[base]{integer}} specifying the random seed to be
     used in the segmentation.  Seed needs to be set for exact numerical
     reproducibility.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{fit} \tab  -\cr
 \tab \code{getChipType} \tab  -\cr
 \tab \code{getChromosomes} \tab  -\cr
 \tab \code{getDataSets} \tab  -\cr
 \tab \code{getFullName} \tab  -\cr
 \tab \code{getName} \tab  -\cr
 \tab \code{getNormalDataSet} \tab  -\cr
 \tab \code{getOutputDataSet} \tab  -\cr
 \tab \code{getTags} \tab  -\cr
 \tab \code{getTumorDataSet} \tab  -\cr
 \tab \code{indexOf} \tab  -\cr
 \tab \code{nbrOfFiles} \tab  -\cr
 \tab \code{setTags} \tab  -\cr
}


 \bold{Methods inherited from ParametersInterface}:\cr
getParameters, getParametersAsString

 \bold{Methods inherited from Object}:\cr
asThis, $, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, registerFinalizer, save, print


 }

 \examples{\dontrun{
   dataSet <- "GSE12702"
tags <- "ASCRMAv2"
chipType <- "Mapping250K_Nsp"
ds <- AromaUnitPscnBinarySet$byName(dataSet, tags=tags, chipType=chipType)
print(ds)

# Extract tumors and normals
idxs <- seq(from=1, to=nbrOfFiles(ds), by=2)
dsT <- extract(ds, idxs);
idxs <- seq(from=2, to=nbrOfFiles(ds), by=2)
dsN <- extract(ds, idxs);

# Setup Paired PSCBS model
seg <- PairedPscbsModel(dsT=dsT, dsN=dsN)
print(seg)

# Segment all tumor-normal pairs
fit(seg, verbose=-10)

 }}

 \references{
  [1] ... \cr
 }

 \seealso{
   ...
 }
\keyword{classes}
