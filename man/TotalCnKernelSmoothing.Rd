%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  TotalCnKernelSmoothing.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{TotalCnKernelSmoothing}
\docType{class}
\alias{TotalCnKernelSmoothing}


\title{The TotalCnKernelSmoothing class}

\description{
 Package:  aroma.cn \cr
\bold{Class TotalCnKernelSmoothing}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.core]{ParametersInterface}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.core]{AromaTransform}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.cn]{TotalCnSmoothing}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{TotalCnKernelSmoothing}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{TotalCnKernelSmoothing}\cr
extends \emph{\link[aroma.cn]{TotalCnSmoothing}}\cr



}

\usage{
TotalCnKernelSmoothing(..., kernel=c("gaussian", "uniform"), bandwidth=50000, censorH=3,
  robust=FALSE)
}

\arguments{
 \item{...}{Arguments passed to \code{\link{TotalCnSmoothing}}.}
 \item{kernel}{A \code{\link[base]{character}} string specifying the type of kernel
    to be used.}
 \item{bandwidth}{A \code{\link[base]{double}} specifying the bandwidth of the smoothing.}
 \item{censorH}{A positive \code{\link[base]{double}} specifying the bandwidth threshold
    where values outside are ignored (zero weight).}
 \item{robust}{If \code{\link[base:logical]{TRUE}}, a robust smoother is used, otherwise not.}
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
