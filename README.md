
<!-- README.md is generated from README.Rmd. Please edit that file -->

# scHiCSRS

<!-- badges: start -->

<!-- badges: end -->

The goal of scHiCSRS is to impute single cell HiC data using a self
representation smoothing approach.

## Installation

Rtools4 (https://cran.r-project.org/bin/windows/Rtools/) and Anaconda/Miniconda (https://docs.conda.io/en/latest/miniconda.html) are required for scHiCSRS. It's recommanded to install them first. 

For windows system, one may need to download Microsoft Visual C++ (https://docs.microsoft.com/en-us/cpp/windows/latest-supported-vc-redist?view=msvc-170).

The scHiCSRS package has the following R-package dependencies: SAVER,
keras, tensorflow, Rtsne, ggplot2, ggpubr, and mclust. The dependent
packages will be automatically installed along with scHiCSRS. You can
use the following commands to install scHiCSRS from GitHub.

``` r
# Install and load "devtools" package. 
install.packages("devtools")
library("devtools")

# Install "scHiCSRS" package from github.
install_github("sl-lin/scHiCSRS")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(scHiCSRS)
## data("simudat")
#simudat_res=SRS(simudat, windowsize=2, nbins=61, epochs = 100, estimates.only = FALSE)
```

For more information of functions, please read the vignettes.
