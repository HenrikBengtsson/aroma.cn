# CRAN submission aroma.cn 1.7.0

on 2022-07-20

This submission addresses Rd issues detected by R-devel.

Thanks in advance


## Notes not sent to CRAN

### R CMD check validation

The package has been verified using `R CMD check --as-cran` on:

| R version     | GitHub | R-hub    | mac/win-builder |
| ------------- | ------ | -------- | --------------- |
| 4.0.x         | L      |          |                 |
| 4.1.x         | L M W  |          |                 |
| 4.2.x         | L M W  | . M M1 . |  . W            |
| devel         | L M W  | .        |    W            |

*Legend: OS: L = Linux, M = macOS, M1 = macOS M1, W = Windows*


R-hub checks:

```r
res <- rhub::check(platform = c(
#  "debian-clang-devel", "debian-gcc-patched", "linux-x86_64-centos-epel",
  "macos-highsierra-release-cran", "macos-m1-bigsur-release"
#  "windows-x86_64-release"
))
print(res)
```

gives

```
── aroma.cn 1.7.0: WARNING

  Build ID:   aroma.cn_1.7.0.tar.gz-aebdabc8ef73495da51e9e066d3712d9
  Platform:   macOS 10.13.6 High Sierra, R-release, CRAN's setup
  Submitted:  2m 29.7s ago
  Build time: 2m 19.4s

❯ checking whether package ‘aroma.cn’ can be installed ... WARNING
  See below...

0 errors ✔ | 1 warning ✖ | 0 notes ✔

── aroma.cn 1.7.0: OK

  Build ID:   aroma.cn_1.7.0.tar.gz-021f6b9dd67d46f1a6330d017177a007
  Platform:   Apple Silicon (M1), macOS 11.6 Big Sur, R-release
  Submitted:  2m 29.7s ago
  Build time: 1m 21.9s

0 errors ✔ | 0 warnings ✔ | 0 notes ✔
```
