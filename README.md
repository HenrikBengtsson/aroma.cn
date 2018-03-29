# aroma.cn: Copy-Number Analysis of Large Microarray Data Sets


## Installation
R package aroma.cn is available on [CRAN](https://cran.r-project.org/package=aroma.cn) and can be installed in R as:
```r
install.packages('aroma.cn')
```

### Pre-release version

To install the pre-release version that is available in Git branch `develop` on GitHub, use:
```r
source('http://callr.org/install#HenrikBengtsson/aroma.cn@develop')
```
This will install the package from source.  



## Contributions

This Git repository uses the [Git Flow](http://nvie.com/posts/a-successful-git-branching-model/) branching model (the [`git flow`](https://github.com/petervanderdoes/gitflow-avh) extension is useful for this).  The [`develop`](https://github.com/HenrikBengtsson/aroma.cn/tree/develop) branch contains the latest contributions and other code that will appear in the next release, and the [`master`](https://github.com/HenrikBengtsson/aroma.cn) branch contains the code of the latest release, which is exactly what is currently on [CRAN](https://cran.r-project.org/package=aroma.cn).

Contributing to this package is easy.  Just send a [pull request](https://help.github.com/articles/using-pull-requests/).  When you send your PR, make sure `develop` is the destination branch on the [aroma.cn repository](https://github.com/HenrikBengtsson/aroma.cn).  Your PR should pass `R CMD check --as-cran`, which will also be checked by <a href="https://travis-ci.org/HenrikBengtsson/aroma.cn">Travis CI</a> and <a href="https://ci.appveyor.com/project/HenrikBengtsson/aroma-cn">AppVeyor CI</a> when the PR is submitted.


## Software status

| Resource:     | CRAN        | Travis CI       | Appveyor         |
| ------------- | ------------------- | --------------- | ---------------- |
| _Platforms:_  | _Multiple_          | _Linux & macOS_ | _Windows_        |
| R CMD check   | <a href="https://cran.r-project.org/web/checks/check_results_aroma.cn.html"><img border="0" src="http://www.r-pkg.org/badges/version/aroma.cn" alt="CRAN version"></a> | <a href="https://travis-ci.org/HenrikBengtsson/aroma.cn"><img src="https://travis-ci.org/HenrikBengtsson/aroma.cn.svg" alt="Build status"></a>   | <a href="https://ci.appveyor.com/project/HenrikBengtsson/aroma-cn"><img src="https://ci.appveyor.com/api/projects/status/github/HenrikBengtsson/aroma.cn?svg=true" alt="Build status"></a> |
| Test coverage |                     | <a href="https://codecov.io/gh/HenrikBengtsson/aroma.cn"><img src="https://codecov.io/gh/HenrikBengtsson/aroma.cn/branch/develop/graph/badge.svg" alt="Coverage Status"/></a>     |                  |
