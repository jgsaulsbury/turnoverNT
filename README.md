
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

Calculate the likelihood of a single observed transition in relative
abundances:

``` r
library(turnoverNT)
xprob(m=c(0.1,0.1,0.8),n=c(0.2,0.3,0.5),Jt=10,ss=c(1000,1000))
#> [1] -0.089366
```

Simulate and visualize a timeseries under neutral theory with incomplete
sampling:

``` r
J <- 50000 #50,000 individuals
nsp <- 10 #10 species
tslength <- 5000 #run for 5000 timesteps
every <- 200 #sample every 200 timesteps...
ss <- 1000 #...and sample 1,000 individuals (with replacement) when you do
ages <- seq(0,tslength,every)
X <- simNT(startingabs=rep(J/nsp,nsp),ts=ages,ss=1000)
plot_spindles(X$simulation,X$times)
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />

Find the best-fit J for this simulated timeseries:

``` r
fitJ(occs=X$simulation,ages=X$times,CI=TRUE)
#> $loglik
#> [1] 527.3662
#> 
#> $J
#> [1] 46080.11
#> 
#> $CI
#> [1] 35043.70 61045.01
```
