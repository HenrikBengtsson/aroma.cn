%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  TotalCnSmoothing.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{TotalCnSmoothing}
\docType{class}
\alias{TotalCnSmoothing}


 \title{The abstract TotalCnSmoothing class}

 \description{
  Package:  aroma.cn \cr
\bold{Class TotalCnSmoothing}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.core]{AromaTransform}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\emph{\code{TotalCnSmoothing}}\cr

\bold{Directly known subclasses:}\cr
\link[aroma.cn]{TotalCnKernelSmoothing}\cr

public abstract static class \bold{TotalCnSmoothing}\cr
extends \emph{\link[aroma.core]{AromaTransform}}\cr



 }

 \usage{TotalCnSmoothing(dataSet=NULL, ..., bandwidth=50000, targetUgp=NULL, .reqSetClass="AromaUnitTotalCnBinarySet")}

 \arguments{
  \item{dataSet}{An \code{\link[aroma.core]{AromaUnitTotalCnBinarySet}}.}
  \item{...}{Arguments passed to \code{\link[aroma.core]{AromaTransform}}.}
  \item{targetUgp}{An \code{\link[aroma.core]{AromaUgpFile}} specifying the
    target loci for which smoothed copy-number are generated.}
  \item{bandwidth}{A \code{\link[base]{double}} specify the bandwidth of the smoothing.}
  \item{.reqSetClass}{(internal only)}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{getTargetUgpFile} \tab  -\cr
 \tab \code{process} \tab  -\cr
 \tab \code{smoothRawCopyNumbers} \tab  -\cr
}


 \bold{Methods inherited from AromaTransform}:\cr
getFullName, getInputDataSet, getName, getOutputDataSet, getOutputFiles, getPath, getTags, isDone, process, setTags

 \bold{Methods inherited from Object}:\cr
asThis, $, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
