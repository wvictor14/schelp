
<!-- README.md is generated from README.Rmd. Please edit that file -->

# schelp

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/schelp)](https://CRAN.R-project.org/package=schelp)
<!-- badges: end -->

The goal of schelp is to …

## Installation

You can install the development version of schelp from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("wvictor14/schelp")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(schelp)
## basic example code
```

## development: subset data

The package will load up library(scibd) which is a private package on
sonoma’s github containing a small subset \~1% of Smillie2019 gut
scRNAseq data.

This data is loaded into `counts` and `metadata` into the environment
for immediate testing / development.

## development: full data

metadata_cells and counts_lntp10k.rds exists in data/ for development
testing

This is the full smillie 2019 data which is

too large to carry with package

``` r
counts <- readRDS(here::here('data', 'counts_lntp10k.rds'))
metadata <- readRDS(here::here('data', 'metadata_cells.rds'))
```
