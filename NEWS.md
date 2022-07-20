# Version (development version)

## New Features

 * `drawC1C2Density()` gained argument `grid`.
 
## Deprecated & Defunct

 * Deprecated `nbrOfFiles()` methods; please use `length()` on
   corresponding input data sets.
   
 * Removed defunct `callPeaks()` for data.frame:s.

## Miscellaneous

 * CLEANUP: Now `doSegmentByPairedPSCBS()` uses `future_mapply()` of
   **future.apply** instead of deprecated `dsApplyInPairs()` of
   **R.filesets**.
   
 * Utilizing subsetted calculations of **matrixStats** (>= 0.50.0).

 * CLEANUP: Now importing generic `extractC1C2()` from PSCBS to
  avoid creating a new one.


# Version 1.6.1 [2015-10-27]

## Miscellaneous

 * Package now requires R (>= 3.1.1) released July 2014. This allows
   us to use Bioconductor (>= 3.0) (October 2014).
   
 * Bumped package dependencies.
 
 * ROBUSTNESS: Explicitly importing core R functions.

## Deprecated & Defunct

 * Removed `makeSmoothSplinePredict()` defunct since Aug 2013.  Made
   `callPeaks()` for data.frame:s defunct (was deprecated).


# Version 1.6.0 [2015-02-23]

## Miscellaneous

 * ROBUSTNESS: Added the first package tests.
 
 * Bumped package dependencies.
 
 * Package now requires R (>= 3.0.3) and Bioconductor (>= 2.13) which were
   released March 2014 and are in fact old; it's recommended to use
   a more recent version of R.


# Version 1.5.9 [2015-01-06]

## Miscellaneous

 * ROBUSTNESS: Package now does a better job importing objects from
   suggested packages.


# Version 1.5.8 [2014-09-04]

## Bug Fixes

 * It could be that `process()` for `AbstractCurveNormalization` would
   generate an error due to read-only permissions introduced by
   copying the target file without resetting the file permissions.

## Miscellaneous

 * Added a few missing NAMESPACE imports.
   

# Version 1.5.7 [2014-06-14]

## Miscellaneous

 * Package now requires R (>= 3.0.0) and BioC (>= 2.13), which were
   released April 2013 and are in fact old and it's recommended to
   use a more recent version of R.
   
 * Updated package dependencies.


# Version 1.5.6 [2014-03-31]

## New Features

 * Added `doSegmentByPairedPSCBS()` for `AromaUnitPscnBinarySet`.


# Version 1.5.5 [2014-03-09]

## Miscellaneous

 * Updated package dependencies.
 
 * Package requires R (>= 2.15.1) and Bioconductor (>= 2.11.0).


# Version 1.5.4 [2014-02-03]

## Bug Fixes

 * ROBUSTNESS: Now `points()` for `C1C2` passes (modified) argument `x` to
   `NextMethod()` as `object = x`.


# Version 1.5.3 [2014-01-30]

## Miscellaneous

 * Updated package dependencies.


# Version 1.5.2 [2013-12-17]

## Miscellaneous

 * CLEANUP: Package no longer uses `:::` in calls.


# Version 1.5.1 [2013-10-23]

## Miscellaneous

 * CLEANUP: Removed several internal prototype methods no longer
   needed or that have been moved to the (private) **Mikado** package.


# Version 1.5.0 [2013-10-17]

## Miscellaneous

 * Minor tweaks to NAMESPACE.
 
 * Updated package dependencies.
 
 * Package requires R (>= 2.15.0) and Bioconductor (>= 2.10.0).


# Version 1.4.6 [2013-10-07]

## Miscellaneous

 * ROBUSTNESS: Now importing only what needs to be imported and
   formally declaring all S3 methods in NAMESPACE.
   
 * CLEANUP: Dropped obsolete usage of `autoload()`.


# Version 1.4.5 [2013-09-28]

## New Features

 * Now the `**aroma.cn**` Package object is also available when the
   package is only loaded (but not attached).

## Miscellaneous

 * Updated package dependencies.


# Version 1.4.4 [2013-09-26]

## Bug Fixes

 * Forgot to import several functions from **matrixStats**.  These
   went undetected because **aroma.light** (< 1.31.6) attached the
   **matrixStats** in the past.


# Version 1.4.3 [2013-09-20]

## Miscellaneous

 * CLEANUP: Now importing only what is needed from **PSCBS**.
 
 * Updated package dependencies.


# Version 1.4.2 [2013-08-21]

## Miscellaneous

 * More internal updates.


# Version 1.4.1 [2013-08-12]

## Bug Fixes

 * `byPath()`, `byName()`, and `findByPath()` for `PairedPSCBSFileSet`
   was also affected by the bug described in the R-devel thread 'Do
   *not* pass '...' to NextMethod() - it'll do it for you; missing
   documentation, a bug or just me?' on Oct 16, 2012.
   
 * `getPath()` for `PairedPscbsModel` would throw an error on
   `getInputDataSet()` not defined.


# Version 1.4.0 [2013-08-04]

## Deprecated & Defunct

 * Made `makeSmoothSplinePredict()` defunct.
 
## Miscellaneous

 * SPEEDUP: Replaced all `rm()` calls with NULL assignments.  Also
   removed several explicit garbage collector calls.
   
 * CLEANUP: The formal package dependency on Bioconductor package
   **aroma.light** has been relaxed so the package can be installed
   without it.
   
 * CLEANUP: Package now only imports **R.oo**.
 
 * Updated package dependencies.


# Version 1.3.4 [2013-05-20]

## Miscellaneous

 * CRAN POLICY: Now all Rd `\usage{}` lines are at most 90 characters long.
 
 * CRAN POLICY: Now all Rd example lines are at most 100 characters long.


# Version 1.3.3 [2013-04-22]

## Significant Changes

 * `findNeutralCopyNumberState()` is now in **PSCBS**.

## New Features

 * Utilizing new `startupMessage()` of **R.oo**.


# Version 1.3.2 [2013-04-22]

## Miscellaneous

 * Updated package dependencies.
 
 * CLEANUP: No longer using deprecated **PSCBS** methods.
 
 * ROBUSTNESS: `{load,save}Cache()` from **R.cache** are now explicitly
   imported in the namespace.


# Version 1.3.1 [2013-01-17]

## Miscellaneous

 * Updated internal methods for `PairedPSCBS` to recognize when other
   mean-level estimators than the sample mean have been used.


# Version 1.3.0 [2013-01-07]

## Bug Fixes

 * `process()` for `PairedPscbsCaller` used the global `verbose`.

## Miscellaneous

 * Bumped up the package dependencies.


# Version 1.2.20 [2012-12-19]

## Bug Fixes

 * Some `example()` scripts used non-defined values.


# Version 1.2.19 [2012-11-26]

## Miscellaneous

 * Bumped up the package dependencies.


# Version 1.2.18 [2012-11-21]

## New Features

 * Now applicable classes utilize the new `ParametersInterface`.
 
 * DOCUMENTATION: Hiding more internal methods from the help indices.


# Version 1.2.17 [2012-11-13]

## Bug Fixes

 * CLEANUP/FIX: Used `cache:` field modified instead of `cached:`.
   After correction, all `clearCache()` methods could be dropped.


# Version 1.2.16 [2012-11-12]

## Miscellaneous

 * CLEANUP: Now `seq_along(x)` instead of `seq(along = x)` everywhere.
   Similarly, `seq(ds)` where `ds` is `GenericDataFileSet` is now
   replaced by `seq_along(ds)`.  Likewise, `seq_len(x)` replaces
   `seq(length = x)`, and `length(ds)` replaces `nbrOfFiles(ds)`.


# Version 1.2.15 [2012-11-05]

## Miscellaneous

 * CLEANUP: Replaced all `whichVector()` with `which()`, because the
   latter is now the fastest again.


# Version 1.2.14 [2012-11-01]

## Miscellaneous

 * ROBUSTNESS: Now package also imports **PSCBS** to please `R CMD check`.
   The reason was that some of the internal methods call **PSCBS**
   methods, which only happens if **PSCBS** is loaded in the first place
   but `R CMD check` cannot known that.


# Version 1.2.13 [2012-10-29]

## Miscellaneous

 * CLEANUP: Now using `Arguments$get{Read,Writ}ablePath()` instead of
   `filePath(..., expandLinks = "any")`.


# Version 1.2.12 [2012-10-21]

## Miscellaneous

 * ROBUSTNESS: Now using `Arguments$getWritablePath()` everywhere instead
   of `mkdirs()`, because the former will do a better job in creating
   and asserting directories on slow shared file systems, and when it
   fails it gives a more informative error message.


# Version 1.2.11 [2012-10-17]

## Miscellaneous

 * ROBUSTNESS: Now all static `Object` methods that calls "next" methods,
   utilizes `NextMethod()`, which became possible with **R.oo** v1.10.0.


# Version 1.2.10 [2012-10-16]

## Bug Fixes

 * ROBUSTNESS/BUG FIX: No longer passing `...` to `NextMethod()`, cf.
   R-devel thread 'Do *not* pass '...' to NextMethod() - it'll do it
   for you; missing documentation, a bug or just me?' on Oct 16, 2012.


# Version 1.2.9 [2012-10-11]

## New Features

 * Added `getOutputFileClass()` and `getOutputFileExtension()` for
   `TotalCnSmoothing`.


# Version 1.2.8 [2012-09-23]

## Miscellaneous

 * More internal updates.


# Version 1.2.7 [2012-09-20]

## New Features

 * Now `PairedPscbsCaller()` passes `...` to the internal callers,
   which makes it possible to for instance specify the number of
   bootstrap samples done for the AB caller.
   
 * Now `PairedPscbsModel()` excludes the actual gaps from the known
   segments it passes to `segmentByPairedPSCBS()`.


# Version 1.2.6 [2012-09-19]

## New Features

 * Added trial version of `PairedPscbsCaller`.

## Bug Fixes

 * `callPeaks(..., flavor="all")` for `PeaksAndValleys` would return
   an error.

## Miscellaneous

 * Additional internal updates.
 

# Version 1.2.5 [2012-09-16]

## New Features

 * Added `calculateTumorPSCNByGenotype()`.


# Version 1.2.4 [2012-09-15]

## New Features

 * Now `fit()` for `PairedPscbsModel` generates pair names iff tumor and
   normal names don't match, e.g. `GSM517071_vs_GSM517072` (if match
   then just `Patient1`).  It also generated `pair` tags.

## Miscellaneous

 * Bumped up the package dependencies.


# Version 1.2.3 [2012-09-05]

## Miscellaneous

 * Bumped up the package dependencies.


# Version 1.2.2 [2012-08-26]

## Miscellaneous
 
 * DOCUMENTATION: Improved help on `TotalCnBinnedSmoothing`.

 * Bumped up the package dependencies.


# Version 1.2.1 [2012-07-22]

## New Features

 * Added trial version of `PairedPscbsModel`.


# Version 1.2.0 [2012-06-05]

## New Features

 * Adopted `findAtomicAberrations()` for `CBS` from ditto of
   `PairedPSCBS`.
 
 * GENERALIZATION: Now `plotTracks()` for `PruneCNA` supports `CBS`
   segmentation results in additional to `PairedPSCBS` ones.
   
 * GENERALIZATION: Now `pruneCNA()` is implemented for `AbstractCBS`,
   not just `PairedPSCBS` objects.
   
 * Merged updates for `findAtomicAberrations()` for `PairedPSCBS` and
   some additional internal "equality" test functions.

## Miscellaneous
    
 * Updated package dependencies.


# Version 1.1.1 [2012-04-16]

## Miscellaneous
 
 * Updated package dependencies.


# Version 1.1.0 [2012-04-10]

## Miscellaneous
 
 * Updated package dependencies.


# Version 1.0.6 [2012-03-30]

## Miscellaneous
 
 * Added help for `normalizePrincipalCurve()`.


# Version 1.0.5 [2012-03-06]

## Bug Fixes

 * One of the PSCBS examples gave an error.


# Version 1.0.4 [2012-02-27]

## Bug Fixes

 * `drawC1C2Density()` for `PairedPSCBS` would throw an exception if
   there was only one segment, or less than two finite (C1,C2):s.


# Version 1.0.3 [2012-02-24]

## Miscellaneous
 
 * Moved some internal functions to the **PSCBS** package, so package
   dependencies was also updated.


# Version 1.0.2 [2012-02-23]

## Miscellaneous
 
 * Updated package dependencies.
 
 * Additional internal updates.


# Version 1.0.1 [2012-01-16]

## New Features

 * Added `TotalCnBinnedSmoothing()`.


# Version 1.0.0 [2012-01-11]

## Miscellaneous
 
 * CLEANUP: The example code for the internal PairedPSCBS methods now
   only runs if environment variable `_R_CHECK_FULL_` is set. This makes
   the package easier on the CRAN servers.

* ROBUSTNESS: Updated package dependencies.
  

# Version 0.9.5 [2011-12-15]

## New Features

 * ROBUSTNESS: Now `process()` of `TotalCnSmoothing` writes atomically.
 
 * Additional internal updates.


# Version 0.9.4 [2011-11-28]

## Miscellaneous
 
 * Updated package dependencies.
 
 * Additional internal updates.


# Version 0.9.3 [2011-11-12]

## Miscellaneous
 
 * Updated the package dependencies.
 
 * Some internal updates.


# Version 0.9.2 [2011-11-02]

## Bug Fixes

 * Updated the memoiziation keys for some of the PairedPSCBS methods
   so that results prior to **PSCBS** v0.13.3 will not be retrieved.


# Version 0.9.1 [2011-10-31]

## Miscellaneous
 
 * Added Rdoc comments to `callPeaks()` for `PeaksAndValleys`.


# Version 0.9.0 [2011-10-28]

## Miscellaneous
 
 * Added a namespace to the package, which will be more or less
   a requirement starting with R v2.14.0.


# Version 0.8.3 [2011-10-16]

## New Features

 * Now `deShearC1C2()`, `translateC1C2()`, and `transformC1C2()` also
   update C1 and C2 mean levels.
   
 * Now using `getLocusData()` and `getSegments()` internally for all
   `PairedPSCBS` objects wherever applicable.


# Version 0.8.2 [2011-08-07]

## Miscellaneous
 
 * The **aroma.cn** v0.8.1 tarball uploaded to CRAN mistakenly
   contained a NAMESPACE file, which shouldn't have been there.


# Version 0.8.1 [2011-07-27]

## Bug Fixes

 * WORKAROUND: In order for the package to work with the most recent
   version of R devel, which automatically add namespaces to packages
   who do not have one, we explicitly have specify that this package
   should use `cat()` and `getOption()` of **R.utils** (instead of
   **base**).


# Version 0.8.0 [2011-07-10]

## Bug Fixes

 * Forgot to update some examples and test scripts which were still
   referring to the old **psCBS** package (should be **PSCBS**).

 * Updated internal code to work with the new column names in
   **PSCBS** v0.11.0.

## Miscellaneous

 * CLEANUP: Removed some internal functions from the help index.
 

# Version 0.7.3 [2011-06-25]

## Miscellaneous

 * ROBUSTNESS: Updated package dependencies.

## Bug Fixes

## Miscellaneous

 * Now package refers to **PSCBS** package (not old psCBS).


# Version 0.7.2 [2011-04-03]

## Miscellaneous

 * CLEANUP: Utilizing `hpaste()` internally wherever applicable.


# Version 0.7.1 [2011-03-03]

## Miscellaneous

 * Fixed a small code typo that didn't make a difference.


# Version 0.7.0 [2011-01-19]

## New Features

 * Added beta classes `PairedPSCBSFile` and `PairedPSCBSFileSet`.
 
## Miscellaneous

 * Removed a NOTE from `R CMD check`.


# Version 0.6.4 [2010-11-04]

## Miscellaneous

 * ROBUSTNESS: Now all bootstrap methods utilize `resample()`.
 
 * Added more internal utility functions for future usage.


# Version 0.6.3 [2010-10-08]

## Miscellaneous

 * Added more internal utility functions for future usage.


# Version 0.6.2 [2010-09-24]

## Miscellaneous

 * Added more internal utility functions for future usage.


# Version 0.6.1 [2010-09-19]

## Miscellaneous

 * Added more internal utility functions for future usage.


# Version 0.6.0 [2010-09-15]

## Miscellaneous

 * Added internal utility functions for future usage.


# Version 0.5.2 [2010-08-04]

## New Features

 * Added option `preserveScale` to `TumorBoostNormalization` for
   correcting for signal compression in heterozygous SNPs.
   The defaults is to do this correction.


# Version 0.5.1 [2010-07-25]

## Significant Changes

 * `callXXorXY()` no longer calls gender from chr Y, when gender is
   estimated as `XX` from chr X.


# Version 0.5.0 [2010-05-14]

## Miscellaneous

 * Package submitted to CRAN.
 
 * Updated citation information.
 
 * Package now requires **aroma.core** v1.6.0.
 
 * Package pass `R CMD check` on R v2.11.0 and v2.12.0 devel.


# Version 0.4.7 [2010-04-04]

## Significant Changes

 * Moved `normalizeDifferencesToAverage()`, `normalizeTumorBoost()`,
   `callNaiveGenotypes()`, and internal `findPeaksAndValleys()`
   to **aroma.light** v1.5.3.


# Version 0.4.6 [2010-03-18]

## Bug Fixes

 * For flavors "v2" and "v3", `normalizeTumorBoost()` could introduce
   NaN:s if `betaN` was exactly zero or exactly one.


# Version 0.4.5 [2010-01-14]

## New Features

 * Added (for now internal) option to change the degrees of freedom of
   the fitted principal curves in MSCN.
   
 * Added `plotSmoothedPairsOne()` to
   `MultiSourceCopyNumberNormalization`.


# Version 0.4.4 [2010-01-05]

## New Features

 * Added support for transform/untransform functions `h(.)` and `g(.)`
   to `AbstractCurveNormalization`, which allows us to fit say on the
   log scale, e.g. `h(x) = log2(x)`, `g(y) = 2^y`.
   
## Bug Fixes

 * `getOutputDataSet()` of `AbstractCurveNormalization` returned all
   files, not just the ones matching the input set.


# Version 0.4.3 [2010-01-01]

## Miscellaneous

 * ROBUSTNESS: Using new `Arguments$getInstanceOf()` were possible.
 
 * ROBUSTNESS: Now all index arguments are validated correctly
   using the new `max` argument of `Arguments$getIndices()`.  Before
   the case where `"max == 0"` was not handled correctly.


# Version 0.4.2 [2009-12-09]

## Significant Changes

 * Made `flavor = "v4"` of `TumorBoostNormalization` the default, and
   if used then no `"flavor"` tag is added.


# Version 0.4.1 [2009-11-03]

## New Features

 * Now `callXXorXY()` and `callNaiveGenotypes()` handles missing values
   and non-finite values.  They also censor outliers to become
   infinite/extreme values.
   
 * Added `callXXorXY()`.
 
 * Added an `example()` to the Rd help of `callNaiveGenotypes()`.
 
 * Added Rd help to `findPeaksAndValleys()`.
 
 * Now argument `tol` of `findPeaksAndValleys()` can be zero; before it
   had to be at least the smallest possible double.


# Version 0.4.0 [2009-11-01]

## Miscellaneous

 * CLEANUP: Removed suggested dependency on **princurve**, which is
   now indirectly suggested/requested via **aroma.light**.
   
 * More recent dependencies on Bioconductor packages.
 
 * Package passes `R CMD check` on R v2.10.0.


# Version 0.3.8 [2009-10-10]

## New Features

 * Added `normalizeTumorBoost()` for `RawAlleleBFractions`.
 
 * Added `callGenotypes()` for `RawAlleleBFractions`.
 
 * Added `RawGenotypeCalls`.


# Version 0.3.7 [2009-10-02]

## Miscellaneous

 * CLEAN UP: Updated to use `byPath()` instead `fromFiles()`.


# Version 0.3.6 [2009-09-30]

## Significant Changes

 * Renamed argument `alignByChromosome` for the constructor of the
   `MultiSourceCopyNumberNormalization` class to `align` in order to
   allow for more types of aligned.
   
 * The alignment of `MultiSourceCopyNumberNormalization` is now done
   using `normalizeDifferencesToAverage()`, which is robust against
   outliers, waviness, etc.  The previous method which normalized
   toward the same overall median is no longer available.

## New Features

 * Added `normalizeDifferencesToAverage()`.
 
## Bug Fixes

 * `getTags()` of `MultiSourceCopyNumberNormalization` would return
   all asterisk tags as merged, e.g. `c("mscn,align", "tagA",
   "tagB")`.


# Version 0.3.5 [2009-07-15]

## New Features

 * ADDED: `XYCurveNormalization` and `PrincipalCurveNormalization`.
 
## Bug Fixes

 * `TumorBoostNormalization`: the `srcFiles` attribute in file footer
   of the result files contained a duplicated default footer instead
   of the tumor-normal pair.


# Version 0.3.4 [2009-07-08]

## New Features

 * Added low-level `callNaiveGenotype()` and `normalizeTumorBoost()`.


# Version 0.3.3 [2009-07-02]

## New Features

 * Added model flavor `"v4"` which corrects heterozygots according to
   `"v2"` and homozygotes according to `"v1"`.
   
 * Added new model flavor (`"v3"`) of `TumorBoostNormalization` that
   is an extension of last weeks model flavor.


# Version 0.3.2 [2009-06-23]

## New Features

 * Added an optional flavor (`"v2"`) of `TumorBoostNormalization` that
   avoids over correcting (especially at the heterozygotes), but
   adjusting the correction factor.  Use argument `flavor = "v2"`.


# Version 0.3.1 [2009-06-08]

## Significant Changes

 * The constructor of `TumorBoostNormalization` now only takes an
   `AromaUnitGenotypeCallSet` for argument `gcN`.  It no longer takes
   an `AromaUnitFracBCnBinarySet` object, which was only an ad hoc
   solution.


# Version 0.3.0 [2009-05-29]
   
## New Features
   
 * Added argument `alignByChromosomes` to
   `MultiSourceCopyNumberNormalization`.  If TRUE, the signals are
   shifted per chromosome such that the mean of the normalized
   smoothed signals is the same for all sources.  This can for
   instance remove systematic effects on sex chromosomes added by some
   ad hoc preprocessing methods.
  
 * Added a `clearCache()` to `MultiSourceCopyNumberNormalization`.
 
 * ALPHA: Added `TumorBoostNormalization`.
 
 * INTERNAL: Added foundations for TumorBoost, i.e. in memory classes
   such as `TotalAndFracBSnpData`.
 
 * INTERNAL: Added `findPeaksAndValleys()`.

 * ROBUSTNESS: Now all constructors report on unknown arguments.
 
 * ROBUSTNESS: Now `MultiSourceCopyNumberNormalization` first write
   normalized data to a temporary file, which is then renamed. This
   lower the risk for having incomplete data in case of interrupts.
  
 * Now `getOutputDataSets()` of `MultiSourceCopyNumberNormalization`
   only returns output data files with a matching fullname in the
   input set.

## Bug Fixes

 * Added missing argument `verbose` in `getTargetPositions()` of
   `TotalCnSmoothing`.  This caused unwanted verbose output in some
   cases.
  
 * `process()` of `TotalCnSmoothing` would not "recognize" fullname
   translators, that is, the output filenames were always identical to
   the input ones.

## Miscellaneous

 * Package passes `R CMD check` and all redundancy tests.


# Version 0.2.2 [2009-02-23]

## Miscellaneous

 * Minor update in order to work with new `RawGenomicSignals`.


# Version 0.2.1 [2009-02-12]

## Miscellaneous

 * Added redundancy tests to package.
 
 * Further cleanup.  Some functions are now in **aroma.light**.


# Version 0.2.0 [2009-01-26]

## Significant Changes

 * `{fit|backtransform}PrincipalCurve()` were moved to **aroma.light** v1.11.1.

 * Several classes and methods were moved to **aroma.core** v1.0.0.

 * Adopted the package to the new classes of **aroma.core**.


# Version 0.1.7 [2008-10-07]

## New Features

 * ALPHA: Added `backtransformPrincipalCurve()`.


# Version 0.1.6 [2008-08-18]

## New Features

 * Added alpha version of `MultiSourceCopyNumberNormalization`.


# Version 0.1.5 [2008-07-30]

## Miscellaneous

 * Fixed some broken cross links in the Rd help.  Package pass `R CMD
   check` on R v2.7.1 and v2.8.0.


# Version 0.1.4 [2008-06-12]

## New Features

 * Now `extractRawCopyNumbers()` of `AromaTotalCnBinaryFile` adds
   annotation data fields to the returned object, e.g. platform,
   chipType, and the fullname of the source file.


# Version 0.1.3 [2008-05-28]

## New Features

 * ALPHA: Added `normalizePrincipalCurve()` and `fitPrincipalCurve()`.


# Version 0.1.2 [2008-05-22]

## New Features

 * ALPHA: Added `extractRawCopyNumbers()` to `AromaTotalCnBinaryFile`.
 
 * ALPHA: Added `TotalCnSmoothing`.


# Version 0.1.1 [2008-05-18]

## New Features

 * Package now provides platform-independent classes
   `Aroma{Total|FreqB}CnSignal{File|Set}`.   With the more generalized
   **aroma.core** package, it is now possible retrieve the `AromaUgpFile`
   for the above.  This provides the necessary basic methods for
   plotting data along chromosomes.


# Version 0.1.0 [2008-05-09]

 * Created.
