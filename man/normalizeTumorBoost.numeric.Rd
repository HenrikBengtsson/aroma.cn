%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  normalizeTumorBoost.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 \name{normalizeTumorBoost.numeric}
\alias{normalizeTumorBoost.numeric}
\alias{normalizeTumorBoost.numeric}

 \alias{normalizeTumorBoost}

 \title{Normalizes allele B fractions for a tumor given a match normal}

 \description{
  TumorBoost is normalization method that normalizes the allele B fractions
  of a tumor sample given the allele B fractions and genotypes of a
  matched normal.
  The method is a single-sample (single-pair) method.
  It does not require total copy-number estimates.
  The normalization is done such that the total copy number is unchanged
  afterwards.
 }

 \usage{\method{normalizeTumorBoost}{numeric}(betaT, betaN, muN=callNaiveGenotypes(betaN), flavor=c("v4", "v3", "v2", "v1"), ...)}

 \arguments{
  \item{betaT, betaN}{Two \code{\link[base]{numeric}} \code{\link[base]{vector}}s each of length J with
     tumor and normal allele B fractions, respectively.}
  \item{muN}{An optional \code{\link[base]{vector}} of length J containing
     normal genotypes calls in (0,1/2,1,\code{\link[base]{NA}}) for (AA,AB,BB).}
  \item{flavor}{A \code{\link[base]{character}} string specifying the type of
     correction applied.}
  \item{...}{Argument passed to \code{\link{callNaiveGenotypes}}(), if called.}
 }

 \value{
   Returns a \code{\link[base]{numeric}} \code{\link[base]{vector}} of length J containing the normalized
   allele B fractions for the tumor.
   Attribute \code{modelFit} is a \code{\link[base]{list}} containing model fit parameters.
 }

 \details{
   Allele B fractions are defined as the ratio between the allele B signal
   and the sum of both (all) allele signals at the same locus.
   Allele B fractions are typically within [0,1], but may have a slightly
   wider support due to for instance negative noise.
   This is typically also the case for the returned normalized
   allele B fractions.
 }

 \examples{
# Load data
pathname <- system.file("exData/TumorBoost,fracB,exampleData.Rbin", package="aroma.cn")
data <- loadObject(pathname)
attachLocally(data)
pos <- position/1e6

layout(matrix(1:4, ncol=1))
par(mar=c(2.5,4,0.5,1)+0.1)
ylim <- c(-0.05, 1.05)
col <- (data[,"genotypeN"] == 1/2) + 1L
plot(pos, betaN, col=col, ylim=ylim)
plot(pos, betaT, col=col, ylim=ylim)

betaTN <- normalizeTumorBoost(betaT=betaT, betaN=betaN, muN=genotypeN)
plot(pos, betaTN, col=col, ylim=ylim)

betaTN2 <- normalizeTumorBoost(betaT=betaT, betaN=betaN)
plot(pos, betaTN2, col=col, ylim=ylim)
}

 \author{Henrik Bengtsson and Pierre Neuvial}
\keyword{methods}
