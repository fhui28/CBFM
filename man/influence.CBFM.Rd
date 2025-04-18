% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/influence.CBFM.R
\name{influence.CBFM}
\alias{influence.CBFM}
\alias{influence}
\title{Regression diagnostics from a CBFM fit}
\usage{
\method{influence}{CBFM}(object, ncores = NULL, ...)
}
\arguments{
\item{object}{An object of class \code{CBFM}.}

\item{ncores}{To speed up calculation of the influence measures, parallelization can be performed, in which case this argument can be used to supply the number of cores to use in the parallelization. Defaults to \code{detectCores()-1}.}

\item{...}{Not used.}
}
\value{
A list containing two elements:
\describe{
\item{hat: }{A matrix of estimated hat values i.e., diagonal elements of the influence/hat matrix. The dimensions of this matrix should be the same as \code{object$y}.}
\item{cooks: }{A matrix of estimated and approximate Cook's distances. The dimensions of this matrix should be the same as \code{object$y}. }
}
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#experimental}{\figure{lifecycle-experimental.svg}{options: alt='[Experimental]'}}}{\strong{[Experimental]}}

Takes a fitted \code{CBFM} object and calculates hat values i.e., diagonal elements of the influence/hat matrix and an approximate Cook's distance.
}
\details{
This is a relatively basic function that calculates the so-called "hat values" and an approximate Cook's distance. Classically, hat values are used to identify potentially high-leverage observations, and Cook's distance is a well-known regression-deletion diagnostic measuring the overall influence of a point on the fit of a model; we refer the reader to \code{\link[stats:influence.measures]{stats::influence.measures()}} for more information and reference pertaining to its use in standard linear and other regression models. Because both concepts of both can be (sort of) carried over to a generalized additive model or GAM e.g., see \code{\link[mgcv:influence.gam]{mgcv::influence.gam()}}, then analogous diagnostics can be produced for a CBFM.

All measures are obtained on a \emph{per-species basis}. This makes sense since the influence/leverage of a unit e.g., space-time coordinate, will differ for different species and their relationships to the measured covariates.

** We leave it up to the user to decide how to use the regression diagnostics, if at all.** In a GAM let alone a CBFM for spatio-temporal multivariate abundance data these diagnostics may only be approximate, and so providing rules-of-thumb related to their usage is challenging. Besides, we echo the sentiment provided in Chapter 4.4 of Fox et al., (2019) that cutoffs and rules-of-thumb should not be given too much weight, with more attention placed on graphical displays and assessing \emph{relative influence} of observations.
}
\examples{
\dontrun{
library(autoFRK)
library(FRK)
library(MASS)
library(mvabund)
library(mvtnorm)
library(ROCR)
library(sp)
library(geoR)
library(tidyverse)

##------------------------------
## **Example 1: Fitting a CBFM to spatial multivariate presence-absence data**
## simulated from a spatial latent variable model
## Please note the data generation process (thus) differs from CBFM.
##------------------------------
set.seed(2021)
num_sites <- 500 # 500 (units) sites 
num_spp <- 50 # Number of species
num_X <- 4 # Number of regression slopes

spp_slopes <- matrix(runif(num_spp * num_X, -1, 1), nrow = num_spp)
spp_intercepts <- runif(num_spp, -2, 0)

# Simulate spatial coordinates and environmental covariate components
xy <- data.frame(x = runif(num_sites, 0, 5), y = runif(num_sites, 0, 5))
X <- mvtnorm::rmvnorm(num_sites, mean = rep(0,4))
colnames(X) <- c("temp", "depth", "chla", "O2")
dat <- data.frame(xy, X)
mm <- model.matrix(~ temp + depth + chla + O2 - 1, data = dat) \%>\% 
scale \%>\% 
as.matrix

# Simulate latent variable component
true_lvs <- grf(grid = cbind(xy$x, xy$y), nsim = 2, cov.model = "exponential",
cov.pars = c(1, 2))$data \%>\%
     as.matrix
spp_loadings <- matrix(runif(num_spp * 2, -1, 1), nrow = num_spp)
set.seed(NULL)

# Simulate spatial multivariate abundance data (presence-absence)
eta <- tcrossprod(cbind(1,mm), cbind(spp_intercepts,spp_slopes)) + 
tcrossprod(true_lvs, spp_loadings)
simy <- matrix(rbinom(num_sites * num_spp, size = 1, 
prob = binomial()$linkinv(eta)), nrow = num_sites)
rm(X, mm, spp_loadings, true_lvs, xy, eta)


# Set up spatial basis functions for CBFM -- Most users will start here! 
num_basisfunctions <- 25 # Number of spatial basis functions to use
basisfunctions <- mrts(dat[,c("x","y")], num_basisfunctions) \%>\% 
as.matrix \%>\%
{.[,-(1)]} # Remove the first intercept column

# Fit CBFM 
useformula <- ~ temp + depth + chla + O2
fitcbfm <- CBFM(y = simy, formula = useformula, data = dat, 
B_space = basisfunctions, family = binomial(), control = list(trace = 1))

influence(fitcbfm)
}

}
\references{
Fox, J. (2019). Regression diagnostics: An introduction. Sage publications.
}
\seealso{
\code{\link[=CBFM]{CBFM()}} for fitting CBFMs, \code{\link[=edf.CBFM]{edf.CBFM()}} for extracting estimated degrees of freedom associated with a CBFM fit, \code{\link[=residuals.CBFM]{residuals.CBFM()}} for calculating different types of residuals from a CBFM fit, and \code{\link[=plot.CBFM]{plot.CBFM()}} for basic residual diagnostics.
}
\author{
Francis K.C. Hui \href{mailto:fhui28@gmail.com}{fhui28@gmail.com}, Chris Haak
}
