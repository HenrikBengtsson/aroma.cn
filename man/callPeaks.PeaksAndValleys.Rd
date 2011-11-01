%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  PeaksAndValleys.EXTS.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 \name{callPeaks.PeaksAndValleys}
\alias{callPeaks.PeaksAndValleys}
\alias{PeaksAndValleys.callPeaks}
\alias{callPeaks.PeaksAndValleys}
\alias{callPeaks,PeaksAndValleys-method}


 \title{Calls the peaks in peaks-and-valley estimates}

 \description{
   Calls the peaks in peaks-and-valley estimates to a set of known state.
 }

 \usage{\method{callPeaks}{PeaksAndValleys}(fit, expected=c(-1/2, -1/4, 0, +1/4, +1/2) * pi, flavor=c("decreasing", "all"), verbose=FALSE, ...)}

 \arguments{
  \item{fit}{A KxC \code{\link[base]{data.frame}} of peaks-and-valley estimates.}
  \item{expected}{The expected locations of the peaks to be called.}
  \item{flavor}{A \code{\link[base]{character}} string specifying what flavor of the
    caller to use.}
  \item{verbose}{A \code{\link[base]{logical}} or a \code{\link[R.utils]{Verbose}} object.}
  \item{...}{Not used.}
 }

 \value{
   Returns a Kx(C+2) \code{\link[base]{data.frame}}.
 }

 \section{Flavors}{
  If \code{flavor == "all"}, each peak is called to the state with the
  closest expected value.
  If \code{flavor == "decreasing"}, the strongest peak is called to the
  state with the closest expected value, then the second strongest peak
  is called analogously to one of the remaining states, and so on.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   To get peaks-and-valley estimates, use
   \code{\link[aroma.light]{findPeaksAndValleys}}.
 }
\keyword{internal}
\keyword{methods}
