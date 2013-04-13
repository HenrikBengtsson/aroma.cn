%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  testAllelicBalanceByBAFs.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{testAllelicBalanceByBAFs}
\alias{testAllelicBalanceByBAFs.default}
\alias{testAllelicBalanceByBAFs}

\title{Tests for allelic balance in a genomic region}

\description{
 Tests for allelic balance in a genomic region.
}

\usage{
\method{testAllelicBalanceByBAFs}{default}(betaT, muN, flavor=c("var", "mean"), ..., verbose=FALSE)
}

\arguments{
  \item{betaT}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} of J tumor allele B fractions (BAFs).}
  \item{muN}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} of J normal (germline) genotypes.}
  \item{flavor}{A \code{\link[base]{character}} specifying the type of test to be performed.}
  \item{...}{Not used.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
  A \code{\link[base]{list}} of class "htest".
}

\examples{

if (Sys.getenv("_R_CHECK_FULL_") != "" && require("PSCBS.data") && require("PSCBS")) {

# Load example ASCN data
data <- PSCBS::exampleData("paired.chr01");

# PSCBS segmentation
fit <- segmentByPairedPSCBS(data, verbose=-10);

# Normalize
fitN <- normalizeBAFsByRegions(fit, verbose=-10);

fitN2 <- callAllelicBalanceByBAFs(fitN, verbose=-10);
print(fitN2);

} # if (require("PSCBS"))
}

\author{Henrik Bengtsson, Pierre Neuvial}


\keyword{internal}
