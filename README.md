
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CBFM (Community-level basis function models)

<!--[![CRAN status](https://www.r-pkg.org/badges/version/CBFM)](https://CRAN.R-project.org/package=CBFM) -->

<!--[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html)-->

`CBFM` is an R package for spatial and spatio-temporal joint species
distributing modeling of multivariate abundance data, based on the idea
of using community-level (species-common) basis functions. It offers an
alternative and connected (but by no means necessarily superior)
approach to the popular spatio-temporal generalized linear latent
variable model (GLLVM) method.

As the name suggests, community-level basis function models (CBFMs) is
built on the idea of fixed-rank kriging (FRK), where spatially- and/or
temporally-indexed basis functions are included to account for
spatio-temporal correlations within and between species. In doing so,
CBFMs bear a lot of similarity to and thus can be set up as a type of
(big) generalized additive model or GAM. This in turns allows CBFMs to
be estimated in a computationally efficient and scalable manner, by
adapting some of the existing machinery available for fitting and
performing inference with GAMs.

The main manuscript introducing CBFMs as a method is available at
[here](https://doi.org/10.1111/2041-210X.14184).

# Installation

Currently, `CBFM` is available and can be installed from github with the
help of `pak` package using:

    pak::pkg_install("fhui28/CBFM")

Alternatively, or if the above does not work, you may download a
(supposedly) stable release of `CBFM` by choosing the latest release on
the right hand side of this Github webpage, and install it manually on
your machine.

# Getting started

**For getting started with `CBFM` though, please have a read of the
manuscript available [here](https://doi.org/10.1111/2041-210X.14184).
After, please install the package then have a look/work through the
examples in the help file for the main workhorse function,
`CBFM::CBFM()`.**

Some broad introductions to joint species distribution modeling,
especially using GLLVMs, please see [So Many Variables: Joint Modeling
in Community Ecology](https://doi.org/10.1016/j.tree.2015.09.007), [How
to make more out of community data? A conceptual framework and its
implementation as models and
software](https://doi.org/10.1111/ele.12757), and the textbook [Joint
Species Distribution Modelling, with Applications in
R](https://doi.org/10.1017/9781108591720).
<!--[Joint dynamic species distribution models: a tool for community ordination and spatio-temporal monitoring](https://doi.org/10.1111/geb.12464) -->

Note there are also some excellent existing packages to fit joint
species distribution models using GLLVMs, including:
[Hmsc](https://cran.r-project.org/package=Hmsc),
[GLLVM](https://cran.r-project.org/package=gllvm),
[sjSDM](https://github.com/TheoreticalEcology/s-jSDM), and
[VAST](https://rdrr.io/github/James-Thorson/VAST/), among others. All of
these have some capacity to handle spatio-temporal multivariate
abundance data using spatio-temporal GLLVMs, or variations thereof.

<!--For general introductions to spatial and/or temporal modeling using basis functions, please check out the excellent [FRK](https://cran.r-project.org/web/packages/FRK/index.html) package for fixed0rank kriging, which heavily inspired this package. Please also see the accompanying software article [FRK: An R Package for Spatial and Spatio-Temporal Prediction with Large Datasets](https://www.jstatsoft.org/article/view/v098i04) and references therein. A more gentle but nevertheless fantastic introduction to basis functions for modeling correlations aimed at ecologists is provided by [The basis function approach for modeling autocorrelation in ecological data](https://esajournals.onlinelibrary.wiley.com/doi/abs/10.1002/ecy.1674). Finally, it would be remiss not to highlight the seminar textbook [Generalized Additive Models: An Introduction with R](https://www.routledge.com/Generalized-Additive-Models-An-Introduction-with-R-Second-Edition/Wood/p/book/9781498728331), and the accompanying [mgcv](https://cran.r-project.org/web/packages/mgcv/index.html) package, which this package both utilizes and takes much inspiration from.-->

## Usage

The below is adapted directly from Example 1a of the help file for the
main workhorse function, `CBFM::CBFM()`, where a CBFM fitted to spatial
multivariate presence-absence data simulated from a spatial latent
variable model. Note this implies the true data generation process
differs from fitted model.

``` r
library(autoFRK)
library(FRK)
library(MASS)
library(mvabund)
library(mvtnorm)
library(ROCR)
library(sp)
library(geoR)
library(tidyverse)
library(CBFM)
```

First, we simulate some spatial multivariate binary (presence-absence)
data

``` r
set.seed(2025)
num_sites <- 500 # Number of (units) sites 
num_spp <- 50 # Number of species
num_X <- 4 # Number of regression slopes

spp_slopes <- matrix(runif(num_spp * num_X, -1, 1), nrow = num_spp)
spp_intercepts <- runif(num_spp, -2, 0)

# Simulate spatial coordinates and environmental covariate components
# We will use this information in later examples as well
xy <- data.frame(x = runif(num_sites, 0, 5), y = runif(num_sites, 0, 5))
X <- mvtnorm::rmvnorm(num_sites, mean = rep(0,4)) 
colnames(X) <- c("temp", "depth", "chla", "O2")
dat <- data.frame(xy, X)
mm <- model.matrix(~ temp + depth + chla + O2 - 1, data = dat) %>% 
scale %>% 
as.matrix

# Simulate latent variable component
# We will use this information in later examples as well
true_lvs <- grf(grid = cbind(xy$x, xy$y), nsim = 2, cov.model = "exponential", 
cov.pars = c(1, 2))$data %>% 
     as.matrix
#> grf: simulation on a set of locations provided by the user
#> grf: process with  1  covariance structure(s)
#> grf: nugget effect is: tausq= 0 
#> grf: covariance model 1 is: exponential(sigmasq=1, phi=2)
#> grf: decomposition algorithm used is:  cholesky 
#> grf: End of simulation procedure. Number of realizations: 2
spp_loadings <- matrix(runif(num_spp * 2, -1, 1), nrow = num_spp) 
set.seed(NULL)

# Simulate spatial multivariate abundance data (presence-absence)
eta <- tcrossprod(cbind(1,mm), cbind(spp_intercepts,spp_slopes)) + 
tcrossprod(true_lvs, spp_loadings)
simy <- matrix(rbinom(num_sites * num_spp, size = 1, 
prob = plogis(eta)), nrow = num_sites)
```

Next, we set up spatial basis functions for CBFM. As discussed in the
help file and in the [article](https://doi.org/10.1111/2041-210X.14184),
there are virtually endless choices for. Below, we will use
multi-resolution thin-plate splines (MRTS) basis functions from the
\[autoFRK()\] package as an (arbitrary) choice.

``` r
num_basisfunctions <- 25 # Number of spatial basis functions to use

train_basisfunctions <- mrts(dat[,c("x","y")], num_basisfunctions) %>% 
as.matrix %>%
{.[,-(1)]} # Remove the first intercept column
```

Fit the CBFM!

``` r
fitcbfm <- CBFM(y = simy, 
                formula = ~ temp + depth + chla + O2, 
                data = dat,
                B_space = train_basisfunctions, 
                family = binomial(), 
                control = list(trace = 1))
#> Compiling TMB C++ file...
#> using C++ compiler: 'g++ (Ubuntu 13.3.0-6ubuntu2~24.04) 13.3.0'
#> Calculating starting values...
#> Commencing model fitting...
#> Updating all coefficients and dispersion/power parameters (this includes running an inner EM algorithm if appropriate).
#> Updating between response correlation (covariance) matrices, G.
#> Updating covariance matrices for basis functions, Sigma, if required.
#> Iteration: 0  Difference in parameter estimates (mean squared error): 0.19559
#> Updating all coefficients and dispersion/power parameters (this includes running an inner EM algorithm if appropriate).
#> Updating between response correlation (covariance) matrices, G.
#> Updating covariance matrices for basis functions, Sigma, if required.
#> Iteration: 1  Difference in parameter estimates (mean squared error): 0.0019
#> Updating all coefficients and dispersion/power parameters (this includes running an inner EM algorithm if appropriate).
#> Updating between response correlation (covariance) matrices, G.
#> Updating covariance matrices for basis functions, Sigma, if required.
#> Iteration: 2  Difference in parameter estimates (mean squared error): 0.00024
#> Updating all coefficients and dispersion/power parameters (this includes running an inner EM algorithm if appropriate).
#> Updating between response correlation (covariance) matrices, G.
#> Updating covariance matrices for basis functions, Sigma, if required.
#> Iteration: 3  Difference in parameter estimates (mean squared error): 7e-05
#> Calculating (components of) the covariance (standard error) matrix...
```

Afterwards, you can a variety of standard modeling functions on the
model fit, such as \[CBFM::summary.CBFM()\], \[CBFM::predict.CBFM()\],
\[CBFM::plot.CBFM()\], among others.

# If you find any bugs and issues…

If you find something that looks like a bug/issue, please let us know
report it, so that we can resolve it and continue to improve this
project (interest and time pending). To report a bug/issue, please make
use of the Github issues and post it up there. As much as possible,
please include in the issue: 1. A description of the bug/issue; 2.
Paste-able code along with some comments that reproduces the problem
e.g., using the [reprex](https://cran.r-project.org/package=reprex)
package. If you also have an idea of how to fix the problem (Francis
tends to make a lot of mistakes in my code, so some may be easy
amendments!), then that is also much appreciated. 3. Required data files
etc…

Thanks heaps!
