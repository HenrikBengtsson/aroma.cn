%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  testAllelicImbalanceByAlleleBFractions.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{testAllelicImbalanceByAlleleBFractions}
\alias{testAllelicImbalanceByAlleleBFractions.default}
\alias{testAllelicImbalanceByAlleleBFractions}

\title{Tests for allelic imbalance in a genomic region}

\description{
 Tests for allelic imbalance in a genomic region.
}

\usage{\method{testAllelicImbalanceByAlleleBFractions}{default}(betaT, muN, flavor=c("var", "mean"), ..., verbose=FALSE)}

\arguments{
  \item{betaT}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} of J tumor allele B fractions.}
  \item{muN}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} of J normal (germline) genotypes.}
  \item{flavor}{A \code{\link[base]{character}} specifying the type of test to be performed.}
  \item{...}{Not used.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
  A \code{\link[base]{list}} of class "htest".
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

