%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  TotalCnBinnedSmoothing.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{TotalCnBinnedSmoothing}
\docType{class}
\alias{TotalCnBinnedSmoothing}


\title{The TotalCnBinnedSmoothing class}

\description{
 Package:  aroma.cn \cr
\bold{Class TotalCnBinnedSmoothing}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.core]{ParametersInterface}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.core]{AromaTransform}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.cn]{TotalCnSmoothing}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{TotalCnBinnedSmoothing}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{TotalCnBinnedSmoothing}\cr
extends \emph{\link[aroma.cn]{TotalCnSmoothing}}\cr



}

\usage{
TotalCnBinnedSmoothing(..., robust=FALSE)
}

\arguments{
 \item{...}{Arguments passed to \code{\link{TotalCnSmoothing}}.}
 \item{robust}{If \code{\link[base:logical]{TRUE}}, a robust smoother is used, otherwise not.}
}

\details{
 Note that \code{dsS <- TotalCnBinnedSmoothing(ds, targetUgp=ugp)} where
 \code{ugp <- getAromaUgpFile(ds)} returns a data set with an identical
 set of loci as the input data set and identical signals as the
 input ones, \emph{except} for loci with duplicated positions.  If all
 loci have unique positions, the the output is identical to the input.
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\emph{No methods defined}.


 \bold{Methods inherited from TotalCnSmoothing}:\cr
getAsteriskTags, getOutputDataSet0, getOutputFileClass, getOutputFileExtension, getOutputFileSetClass, getOutputFiles, getParameters, getPath, getRootPath, getTargetPositions, getTargetUgpFile, process, smoothRawCopyNumbers

 \bold{Methods inherited from AromaTransform}:\cr
as.character, findFilesTodo, getAsteriskTags, getExpectedOutputFiles, getExpectedOutputFullnames, getFullName, getInputDataSet, getName, getOutputDataSet, getOutputDataSet0, getOutputFiles, getPath, getRootPath, getTags, isDone, process, setTags

 \bold{Methods inherited from ParametersInterface}:\cr
getParameterSets, getParameters, getParametersAsString

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, names, objectSize, print, save, asThis


}

\author{Henrik Bengtsson}
\keyword{classes}
