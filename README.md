# CBFM (Community-level basis function models)

<!-- badges: start --> [![CRAN status](https://www.r-pkg.org/badges/version/CBFM)](https://CRAN.R-project.org/package=CBFM) [![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental) <!-- badges: end -->

`CBFM` is an R package for spatial and spatio-temporal joint species distributing modeling of multivariate abundance data, based on the idea of using community-level (species-common) basis functions. It offers an alternative and connected (but by no means necessarily superior) approach to the popular spatio-temporal generalized linear latent variable model (GLLVM) method.

As the name suggests, community-level basis function models (CBFMs) is built on the idea of fixed-rank kriging (FRK), where spatially- and/or temporally-indexed basis functions are included to account for spatio-temporal correlations within and between species. In doing so, CBFMs bear a lot of similarity to and thus can be set up as a type of (big) generalized additive model or GAM. This in turns allows CBFMs to be estimated in a computationally efficient and scalable manner, by adapting some of the existing machinery available for fitting and performing inference with GAMs.

The main manuscript introducing CBFMs as a method is available at [here](https://doi.org/10.1111/2041-210X.14184). 


# Installation

````{=html}
<!--Down the road, you can install the package on CRAN using:
```
install.packages("CBFM")
```
-->
````

Currently, `CBFM` is available and can be installed from github with the help of `devtools` package using:

```         
devtools::install_github("fhui28/CBFM")
```

Alternatively, or if the above does not work, you may download a (supposedly) stable release of `CBFM` by choosing the latest release on the right hand side of this Github webpage, and install it manually on your machine.

# Getting started

For getting started with `CBFM` though, please have a read of the manuscript available  [here](https://doi.org/10.1111/2041-210X.14184). We also recommend installing it and reading the help file for the main CBFM function. 

Some broad introductions to joint species distribution modeling, especially using GLLVMs, please see [So Many Variables: Joint Modeling in Community Ecology](https://doi.org/10.1016/j.tree.2015.09.007), [How to make more out of community data? A conceptual framework and its implementation as models and software](https://doi.org/10.1111/zele.12757), 
<!--[Joint dynamic species distribution models: a tool for community ordination and spatio-temporal monitoring](https://doi.org/10.1111/geb.12464) -->
and the textbook [Joint Species Distribution Modelling, with Applications in R](https://doi.org/10.1017/9781108591720).

Note there are also some excellent existing \texttt{R} packages to fit joint species distribution models using GLLVMs, including: [Hmsc](https://cran.r-project.org/web/packages/Hmsc/index.html), [GLLVM](https://cran.r-project.org/web/packages/gllvm/index.html), [sjSDM](https://github.com/TheoreticalEcology/s-jSDM), and [VAST](https://rdrr.io/github/James-Thorson/VAST/), among others. All of these have some capacity to handle spatio-temporal multivariate abundance data using spatio-temporal GLLVMs, or variations thereof.

<!--For general introductions to spatial and/or temporal modeling using basis functions, please check out the excellent [FRK](https://cran.r-project.org/web/packages/FRK/index.html) package for fixed0rank kriging, which heavily inspired this package. Please also see the accompanying software article [FRK: An R Package for Spatial and Spatio-Temporal Prediction with Large Datasets](https://www.jstatsoft.org/article/view/v098i04) and references therein. A more gentle but nevertheless fantastic introduction to basis functions for modeling correlations aimed at ecologists is provided by [The basis function approach for modeling autocorrelation in ecological data](https://esajournals.onlinelibrary.wiley.com/doi/abs/10.1002/ecy.1674). Finally, it would be remiss not to highlight the seminar textbook [Generalized Additive Models: An Introduction with R](https://www.routledge.com/Generalized-Additive-Models-An-Introduction-with-R-Second-Edition/Wood/p/book/9781498728331), and the accompanying [mgcv](https://cran.r-project.org/web/packages/mgcv/index.html) package, which this package both utilizes and takes much inspiration from.-->

# If you find any bugs and issues...

If you find something that looks like a bug/issue, please let us know report it, so that we can resolve it and continue to improve this project (interest and time pending). To report a bug/issue, please make use of the Github issues and post it up there. As much as possible, please include in the issue: 1. A description of the bug/issue; 2. Paste-able code along with some comments that reproduces the problem e.g., using the [reprex](https://cran.r-project.org/web/packages/reprex/index.html) package. If you also have an idea of how to fix the problem (Francis tends to make a lot of mistakes in my code, so some may be easy amendments!), then that is also much appreciated. 3. Required data files etc...

Thanks heaps!

Finally, please note Github issues is a place to post any bugs/issues/unusual behavior in relation to the package. **Please try to refrain from posting any feature requests of general statistical modeling questions**. They will very likely be ignored and/or deleted without prior consent. You are probably better off emailing one of the package maintainers instead. Much appreciated!
