%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MultiSourceCopyNumberNormalization.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{MultiSourceCopyNumberNormalization}
\docType{class}
\alias{MultiSourceCopyNumberNormalization}


 \title{The MultiSourceCopyNumberNormalization class}

 \description{
  Package:  aroma.cn \cr
\bold{Class MultiSourceCopyNumberNormalization}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{MultiSourceCopyNumberNormalization}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{MultiSourceCopyNumberNormalization}\cr
extends \link[R.oo]{Object}\cr



  The multi-source copy-number normalization (MSCN) method [1] is a
  normalization method that normalizes copy-number estimates measured
  by multiple sites and/or platforms for common samples.  It normalizes the
  estimates toward a common scale such that for any copy-number level
  the mean level of the normalized data are the same.
 }

 \usage{MultiSourceCopyNumberNormalization(dsList=NULL, fitUgp=NULL, subsetToFit=NULL, targetDimension=1, align=c("byChromosome", "none"), tags="*", ...)}

 \arguments{
  \item{dsList}{A \code{\link[base]{list}} of K \code{\link[aroma.core]{AromaUnitTotalCnBinarySet}}:s.}
  \item{fitUgp}{An \code{\link[aroma.core]{AromaUgpFile}} that specifies the
    common set of loci used to normalize the data sets at.}
  \item{subsetToFit}{The subset of loci (as mapped by the \code{fitUgp}
    object) to be used to fit the normalization functions.
    If \code{\link[base]{NULL}}, loci on chromosomes 1-22 are used, but not on ChrX and ChrY.
  }
  \item{targetDimension}{A \code{\link[base]{numeric}} index specifying the data set in
    \code{dsList} to which each platform in standardize towards.
    If \code{\link[base]{NULL}}, the arbitrary scale along the fitted principal curve
    is used.  This always starts at zero and increases.}
  \item{align}{A \code{\link[base]{character}} specifying type of alignment applied, if any.
    If \code{"none"}, no alignment is done.
    If \code{"byChromosome"}, the signals are shifted chromosome
    by chromosome such the corresponding smoothed signals have the same
    median signal across sources.
    For more details, see below.
  }
  \item{tags}{(Optional) Sets the tags for the output data sets.}
  \item{...}{Not used.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[aroma.cn:getAllNames.MultiSourceCopyNumberNormalization]{getAllNames}} \tab Gets the names of all unique samples across all sources.\cr
 \tab \code{getAsteriskTags} \tab  -\cr
 \tab \code{\link[aroma.cn:getInputDataSets.MultiSourceCopyNumberNormalization]{getInputDataSets}} \tab Gets the list of data sets to be normalized.\cr
 \tab \code{getOutputDataSets} \tab  -\cr
 \tab \code{getOutputPaths} \tab  -\cr
 \tab \code{getParameters} \tab  -\cr
 \tab \code{getTags} \tab  -\cr
 \tab \code{nbrOfDataSets} \tab  -\cr
 \tab \code{\link[aroma.cn:process.MultiSourceCopyNumberNormalization]{process}} \tab Normalizes all samples.\cr
}


 \bold{Methods inherited from Object}:\cr
asThis, $, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, registerFinalizer, save


 }

 \details{
   The multi-source normalization method is by nature a single-sample method,
   that is, it normalizes arrays for one sample at the time and independently
   of all other samples/arrays.

   However, the current implementation is such that it first generates
   smoothed data for \emph{all} samples/arrays.  Then, it normalizes the
   sample one by one.
 }

 \section{Different preprocessing methods normalize ChrX \& ChrY differently}{
    Some preprocessing methods estimate copy numbers on sex chromosomes
    differently from the autosomal chromosomes.  The way this is done may
    vary from method to method and we cannot assume anything about what
    approach is.  This is the main reason why the estimation of the
    normalization  function is by default based on signals from autosomal
    chromosomes only; this protects the estimate of the function from
    being biased by specially estimated sex-chromosome signals.
    Note that the normalization function is still applied to all chromosomes.

    This means that if the transformation applied by a particular
    preprocessing method is not the same for the sex chromosomes as the
    autosomal chromosomes, the normalization applied on the sex
    chromosomes is not optimal one.  This is why multi-source
    normalization sometimes fails to bring sex-chromosome signals
    to the same scale across sources.  Unfortunately, there is no
    automatic way to handle this.
    The only way would be to fit a specific normalization function to each
    of the sex chromosomes, but that would require that there exist
    copy-number abberations on those chromosomes, which could be a too
    strong assumption.

    A more conservative approach is to normalize the signals such that
    afterward the median of the smoothed copy-number levels are the same
    across sources for any particular chromosome.
    This is done by setting argument \code{align="byChromosome"}.
 }

 \references{
   [1] H. Bengtsson, A. Ray, P. Spellman & T.P. Speed,
       \emph{A single-sample method for normalizing and combining
         full-resolution copy numbers from multiple platforms,
         labs and analysis methods},
       Bioinformatics 2009. \cr
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
