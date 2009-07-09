%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  callNaiveGenotypes.numeric.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 \name{callNaiveGenotypes.numeric}
\alias{callNaiveGenotypes.numeric}
\alias{callNaiveGenotypes.numeric}

 \alias{callNaiveGenotypes}

 \title{Calls genotypes in a normal sample}

 \description{
   Calls genotypes in a normal sample.
 }

 \usage{\method{callNaiveGenotypes}{numeric}(y, cn=rep(2, length(y)), flavor=c("density"), adjust=1.5, ...)}

 \arguments{
  \item{y}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} of length J containing allele B fractions
    for a normal sample.}
  \item{cn}{An optional \code{\link[base]{numeric}} \code{\link[base]{vector}} of length J specifying the true
    total copy number in \eqn{\{0,1,2,NA\}} at each locus.  This can be
    used to specify which loci are diploid and which are not, e.g.
    autosomal and sex chromosome copy numbers.}
  \item{flavor}{A \code{\link[base]{character}} string specifying the type of algorithm used.}
  \item{adjust}{A postive \code{\link[base]{double}} specifying the amount smoothing for
    the empirical density estimator.}
  \item{...}{Not used.}
 }

 \value{
   Returns a \code{\link[base]{numeric}} \code{\link[base]{vector}} of length J containing the genotype calls
   in allele B fraction space, that is, in [0,1] where 1/2 corresponds
   to a heterozygous call, and 0 and 1 corresponds to homozygous A
   and B, respectively.
   Non called genotypes have value \code{\link[base]{NA}}.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{methods}
