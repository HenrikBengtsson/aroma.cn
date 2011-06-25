%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  PrincipalCurveNormalization.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{PrincipalCurveNormalization}
\docType{class}
\alias{PrincipalCurveNormalization}


 \title{The PrincipalCurveNormalization class}

 \description{
  Package:  aroma.cn \cr
\bold{Class PrincipalCurveNormalization}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.cn]{AbstractCurveNormalization}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{PrincipalCurveNormalization}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{PrincipalCurveNormalization}\cr
extends \emph{\link[aroma.cn]{AbstractCurveNormalization}}\cr


 }

 \usage{PrincipalCurveNormalization(..., subset=1/20)}

 \arguments{
  \item{...}{Arguments passed to \code{\link{AbstractCurveNormalization}}.}
  \item{subset}{A \code{\link[base]{double}} in (0,1] specifying the fraction of the
    \code{subsetToFit} to be used for fitting.  Since the fit function
    for this class is rather slow, the default is to use a 1/20:th
    of the default data points.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\emph{No methods defined}.


 \bold{Methods inherited from AbstractCurveNormalization}:\cr
backtransformOne, fitOne, getDataSets, getFullName, getInputDataSet, getName, getOutputDataSet, getPairedDataSet, getPath, getRootPath, getSubsetToFit, getTags, getTargetDataSet, nbrOfFiles, process, setTags

 \bold{Methods inherited from Object}:\cr
asThis, getChecksum, $, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, registerFinalizer, save


 }

 \author{Henrik Bengtsson}
\keyword{classes}
