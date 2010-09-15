%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  findAtomicAberrations.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{findAtomicAberrations.CopyNumberRegions}
\alias{findAtomicAberrations.CopyNumberRegions}
\alias{CopyNumberRegions.findAtomicAberrations}
\alias{findAtomicAberrations.CopyNumberRegions}
\alias{findAtomicAberrations,CopyNumberRegions-method}

\alias{findAtomicAberrations}

\title{Finds all possible atomic regions}

\description{
  Finds all possible atomic regions of a certain length.
}

\usage{\method{findAtomicAberrations}{CopyNumberRegions}(cnr, data, H=1, alpha=0.02, ..., verbose=FALSE)}

\arguments{
  \item{cnr}{The segments defining the partitioning of the data.}
  \item{data}{The data  used to test for equality.}
  \item{H}{A positive \code{\link[base]{integer}} specifying how many segments each
     atomic abberation should contain.}
  \item{alpha}{A \code{\link[base]{double}} in [0,1] specifying the significance level
     for testing the null-hypothesis that the flanking segments
     are equal.}
  \item{...}{Not used.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
  Returns a \code{\link[base]{data.frame}} with K rows, where K >= 0 is the number
  of atomic aberrations found.
}

\details{
  An \emph{aberration of length H} is defined as any H consecutive segments.
  Each aberrations has two \emph{flanking segments} on each side.
  Regardless of the content of the aberration, it is possible
  to test the null-hypothesis that the two flanking segments are
  equal or not.
  The two flanking regions are said to be \emph{equal}, if the
  null-hypothesis of being equal is \emph{not} rejected.
  If the two flanking regions are called equal, then the contained
  abberation (of length H) is called \emph{atomic}, otherwise not.

  For consistency one may also define atomic aberrations of length H=0.
  Consider that an imaginary aberration of zero length splits a single
  segment into to flanking segments.  Then by construction those two
  segments are equal.  The case where H=0 is still not implemented.
}

\examples{

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Local functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
xAxisBar <- function(x0, x1=NULL, y=par("usr")[3], lwd=4, col=par("col"), length=0.04, xpd=FALSE, ...) {
  if (is.null(x1)) {
    x1 <- x0[2];
    x0 <- x0[1];
  }
  arrows(x0=x0, x1=x1, y0=y, code=3, angle=90, length=length, lwd=lwd, col=col, xpd=xpd, lend=2, ...);
} # xAxisBar()


verbose <- Arguments$getVerbose(-8, timestamp=TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Simulating copy-number data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Build up CN profile generation by generation
pT <- cnr(1,2000, 2) +
     cnr(1000,1500) +
     cnr(1000,1250) +
     cnr(1650,1800) +
     cnr(200,300) - cnr(650,800);
print(pT);

# Simulate data from the track
cn <- simulateRawCopyNumbers(pT, n=2000, sd=1/2);
print(cn);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Find atomic aberrations
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
subplots(4, ncol=1);
par(mar=c(2,4,1,1)+0.1);

# Plot observed and true signals
plot(cn, col="#aaaaaa", ylim=c(0,5));
drawLevels(pT, col="white", lwd=4, lty=1);
drawLevels(pT, col="black", lwd=3, lty=6);
stext(side=3, pos=1, line=-1, "\"H=0\"", cex=1.5);

for (H in 1:3) {
  plot(cn, col="#aaaaaa", ylim=c(0,5));
  drawLevels(t, col="white", lwd=4, lty=1);
  drawLevels(t, col="black", lwd=3, lty=6);
  col <- H+1;
  stext(side=3, pos=1, line=-1, sprintf("H=\%d", H), cex=1.5, col=col);
  par <- par("usr");
  y0 <- par("usr")[3];
  y1 <- par("usr")[4];
  dy <- 0.05*(y1-y0);
 
  fit <- findAtomicAberrations(cnr=t, data=cn, H=H, verbose=verbose);
  df <- fit$res;
 
  for (kk in seq(length=nrow(df))) {
    dfKK <- df[kk,];
    segments <- as.integer(dfKK[,c("firstRegion", "lastRegion")]);
    segments <- segments[1]:segments[2];
    xRange <- as.double(dfKK[,c("start", "stop")]);
    cnrKK <- subset(pT, subset=segments);
    drawLevels(cnrKK, col=col, lwd=3);
    x <- xRange/1e6;
    y <- y0 + 0.2*dy;
    xAxisBar(x0=x[1], x1=x[2], y=y, col=col);
    box();
  }
} # for (H ...)

}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

% \references{
%   [1] \url{http://www.definethat.com/define/7274.htm}
%   [2] \url{http://en.wikipedia.org/wiki/Atomic_(order_theory)}
%   [3] \url{http://en.wikipedia.org/wiki/Atomic_(measure_theory)}
% }

\seealso{
  ...
}


\keyword{internal}
\keyword{methods}
