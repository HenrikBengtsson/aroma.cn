%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  normalizeMirroredBAFsByRegions.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{normalizeMirroredBAFsByRegions.matrix}
\alias{normalizeMirroredBAFsByRegions.matrix}
\alias{normalizeMirroredBAFsByRegions.matrix}

\alias{normalizeMirroredBAFsByRegions}

\title{Normalizes region-level mirrored allele B fractions (mBAFs)}

\description{
 Normalizes region-level mirrored allele B fractions (mBAFs) for heterozygous and homozygous SNPs by rescaling both
 equally much such that the normalized homozygous mBAFs become exactly
 or close to their expected values (which is at zero and one).
}

\usage{\method{normalizeMirroredBAFsByRegions}{matrix}(data, flavor=c("plain", "total"), ..., verbose=FALSE)}

\arguments{
  \item{data}{A \code{\link[base]{numeric}} Kx2 or Kx3 \code{\link[base]{matrix}}, where K is the number of
    segments and the first and the second column contains average
    heterozygous and homozygous mBAF estimates, respectively.
    The third column, which is optional, contains total copy numbers.}
  \item{flavor}{A \code{\link[base]{character}} string specifying how the normalization
    function/scale factors are estimated.}
  \item{...}{Additional arguments passed \code{\link[aroma.light]{fitXYCurve}},
    which is used if \code{flavor="total"}.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
  A \code{\link[base]{numeric}} \code{\link[base]{matrix}} of the same dimensions as argument \code{data}.
}

\details{
  The mBAFs for heterozygous SNPs are also known as the
  Decrease in Heterozygosity signals (DHs).
}

%examples "../incl/normalizeMirroredBAFsByRegions.Rex"

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}


\keyword{methods}
\keyword{internal}
