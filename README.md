
<!-- README.md is generated from README.Rmd. Please edit that file -->

# turnoverNT

<!-- badges: start -->
<!-- badges: end -->

Contains functions for modeling change in fossil communities with
Hubbell’s neutral theory.

## Installation

You can install the development version of turnoverNT from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("jgsaulsbury/turnoverNT")
```

## Example

Calculating the likelihood of a single observed transition in relative
abundances:

``` r
library(turnoverNT)
xprob(m=c(0.1,0.1,0.8),n=c(0.2,0.3,0.5),tJ=0.1,ss=c(1000,1000))
#> [1] -0.089366
```
