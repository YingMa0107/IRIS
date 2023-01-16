---
layout: page
title: Installation
description: ~
---

`IRIS` is implemented as an Rcpp package, which can be installed from GitHub by:

### Dependencies 
* R version >= 4.2.2.
* Dependent R packages: Rcpp (>= 1.0.9), RcppArmadillo, SingleCellExperiment, SummarizedExperiment, methods, Matrix, MCMCpack, fields, wrMisc, RANN, stats, ggplot2, grDevices, reshape2


#### 1. Install `devtools` if necessary
```r
install.packages('devtools')
```

#### 2. Install `CARD`
```r
devtools::install_github('YingMa0107/IRIS')
```
#### 3. Load package
```r
library(IRIS)
```

This package is supported for Windows 10, MAC and Linux. The package has been tested on the following systems:
- MAC (Catalina 10.15, Monterey 12.4)
- Linux: Ubuntu (16.04.6)
- Windows 10: Home (1903)

#### 4. Some possible issues when installing the package, especially on the MacOS system
(1) Cannot find tools necessary when using R in the MacOS system.
```r
Error: Failed to install 'IRIS' from GitHub:
  Could not find tools necessary to compile a package
Call `pkgbuild::check_build_tools(debug = TRUE)` to diagnose the problem.
``` 
possible solution: in R, type the code ``` options(buildtools.check = function(action) TRUE )```, see the discussion about this error, [link](https://stackoverflow.com/questions/37776377/error-when-installing-an-r-package-from-github-could-not-find-build-tools-neces)

(2) library not found for -lgfortran when using R in the MacOS system.
```r
ld: library not found for -lgfortran
```
It seems the gfortran is not well installed on the MacOS system. Please check this [link](https://thecoatlessprofessor.com/programming/cpp/r-compiler-tools-for-rcpp-on-macos/) for the gfortran installation to see if it helps. 


