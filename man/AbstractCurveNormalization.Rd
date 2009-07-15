%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AbstractCurveNormalization.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{AbstractCurveNormalization}
\docType{class}
\alias{AbstractCurveNormalization}


 \title{The AbstractCurveNormalization class}

 \description{
  Package:  aroma.cn \cr
\bold{Class AbstractCurveNormalization}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\emph{\code{AbstractCurveNormalization}}\cr

\bold{Directly known subclasses:}\cr
\link[aroma.cn]{PrincipalCurveNormalization}, \link[aroma.cn]{XYCurveNormalization}\cr

public abstract static class \bold{AbstractCurveNormalization}\cr
extends \link[R.oo]{Object}\cr


 }

 \usage{AbstractCurveNormalization(dataSet=NULL, targetSet=NULL, subsetToFit=NULL, tags="*", copyTarget=TRUE, ...)}

 \arguments{
  \item{dataSet}{An \code{\link[aroma.core]{AromaUnitTotalCnBinarySet}} of
     "test" samples to be normalized.}
  \item{targetSet}{An \code{\link[aroma.core]{AromaUnitTotalCnBinarySet}} of
     paired target samples.}
  \item{subsetToFit}{The subset of loci to be used to fit the
    normalization functions.
    If \code{\link[base]{NULL}}, loci on chromosomes 1-22 are used, but not on ChrX and ChrY.
  }
  \item{tags}{(Optional) Sets the tags for the output data sets.}
  \item{copyTarget}{If \code{\link[base:logical]{TRUE}}, target arrays are copied to the output
     data set, otherwise not.}
  \item{...}{Not used.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{getFullName} \tab  -\cr
 \tab \code{getInputDataSet} \tab  -\cr
 \tab \code{getName} \tab  -\cr
 \tab \code{getOutputDataSet} \tab  -\cr
 \tab \code{getPath} \tab  -\cr
 \tab \code{getRootPath} \tab  -\cr
 \tab \code{getTags} \tab  -\cr
 \tab \code{getTargetDataSet} \tab  -\cr
 \tab \code{nbrOfFiles} \tab  -\cr
 \tab \code{process} \tab  -\cr
 \tab \code{setTags} \tab  -\cr
}


 \bold{Methods inherited from Object}:\cr
asThis, $, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, registerFinalizer, save


 }

 \section{See also}{
   \code{\link[aroma.light]{fitPrincipalCurve}}
 }

 \author{Henrik Bengtsson}
\keyword{classes}
