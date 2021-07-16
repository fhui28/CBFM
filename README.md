# CBFM (Community-level basis function models)

`CBFM` is an R package for spatio-temporal joint species distributing modeling of multivariate abundance data, based on the use of community-level basis functions. As an alternative (but by no means necessarily superior) approach to the increasingly popular framework generalized linear latent variable model (GLLVM), and as the name suggests, community-level basis function models (CBFMs) build on the idea of fixed-rank kriging (FRK), where spatial and/or temporally indexed basis functions are included to account for spatio-temporal correlations both within and between species. In doing, CBFMs bear a lot of similarity to and thus can be set up as a type of generalized additive model (GAM).


# Installation

Down the road, you can install the package on CRAN using:
```
install.packages("CBFM")
```
Currently, the developmental version of `CBFM` is available and can be installed from github with the help of `devtools` package using:
```
devtools::install_github("fhui28/CBFM")
```

# Getting started

We are in the process of producing the first manuscript introducing and explaining CBFMs to community ecologists. For getting started with `CBFM` though, we recommend installing it and reading the help file for the main CBFM function. For broad introductions to joint species distribution modeling, especially using GLLVMs, we recommend reading [So Many Variables: Joint Modeling in Community Ecology](https://www.sciencedirect.com/science/article/pii/S0169534715002402?casa_token=me_1KcIBbeMAAAAA:9-EUCdI5o5e1g5pSk5biiKO9zKj-wdwxtc4yNcgtRhFrPSQeXLl9a3n1DE_Furrnigb5i5PzbrM), [How to make more out of community data? A conceptual framework and its implementation as models and software](https://onlinelibrary.wiley.com/doi/full/10.1111/ele.12757), [Joint dynamic species distribution models: a tool for community ordination and spatio-temporal monitoring](https://onlinelibrary.wiley.com/doi/abs/10.1111/geb.12464) and [Understanding co-occurrence by modelling species simultaneously with a Joint Species Distribution Model (JSDM)](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12180). 

There are some existing and good \texttt{R} packages to fit joint species distribution models using the GLLVMs, including [Hmsc](https://cran.r-project.org/web/packages/Hmsc/index.html), [GLLVM](https://cran.r-project.org/web/packages/gllvm/index.html),  [boral](https://cran.r-project.org/web/packages/boral/index.html), and [VAST](https://rdrr.io/github/James-Thorson/VAST/).

For general introductions to spatial and/or temporal modeling using basis functions, please check out the excellent [FRK](https://cran.r-project.org/web/packages/FRK/index.html) package for fixed rank kriging, which heavily inspired this package. Please also see the accompanying software article [FRK: An R Package for Spatial and Spatio-Temporal Prediction with Large Datasets](https://www.jstatsoft.org/article/view/v098i04) and references therein. A more gentle but nevertheless fantastic introduction to basis functions for modeling correlations aimed at ecologists is provided  [The basis function approach for modeling autocorrelation in ecological data](https://esajournals.onlinelibrary.wiley.com/doi/abs/10.1002/ecy.1674)



