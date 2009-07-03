%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  TumorBoostNormalization.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{TumorBoostNormalization}
\docType{class}
\alias{TumorBoostNormalization}


 \title{The TumorBoostNormalization class}

 \description{
  Package:  aroma.cn \cr
\bold{Class TumorBoostNormalization}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{TumorBoostNormalization}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{TumorBoostNormalization}\cr
extends \link[R.oo]{Object}\cr



  TumorBoost is normalization method that normalizes the allele B fractions
  of a tumor sample given the allele B fractions and genotype calls for
  a matched normal.
  The method is a single-sample (single-pair) method.  It does not require
  total copy number estimates.
  The normalization is done such that the total copy number is unchanged
  afterwards.
 }

 \usage{TumorBoostNormalization(dsT=NULL, dsN=NULL, gcN=NULL, flavor=c("v1", "v2", "v3"), collapseHomozygotes=FALSE, tags="*", ...)}

 \arguments{
  \item{dsT}{...}
  \item{dsN}{...}
  \item{gcN}{...}
  \item{flavor}{Type of correction applied..}
  \item{collapseHomozygotes}{If \code{\link[base:logical]{TRUE}}, SNPs that are homozygous in the
    matched normal are also called homozygous in the tumor, that is,
    it's allele B fraction is collapsed to either 0 or 1.
    If \code{\link[base:logical]{FALSE}}, the homozygous values are normalized according the
    model.
  }
  \item{tags}{(Optional) Sets the tags for the output data sets.}
  \item{...}{Not used.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{getFullName} \tab  -\cr
 \tab \code{getInputDataSet} \tab  -\cr
 \tab \code{getName} \tab  -\cr
 \tab \code{getNormalDataSet} \tab  -\cr
 \tab \code{getNormalGenotypeCallSet} \tab  -\cr
 \tab \code{getOutputDataSet} \tab  -\cr
 \tab \code{getPath} \tab  -\cr
 \tab \code{getRootPath} \tab  -\cr
 \tab \code{getTags} \tab  -\cr
 \tab \code{nbrOfFiles} \tab  -\cr
 \tab \code{process} \tab  -\cr
 \tab \code{setTags} \tab  -\cr
}


 \bold{Methods inherited from Object}:\cr
asThis, $, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
