% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate.CBFM.R
\name{simulate.CBFM}
\alias{simulate.CBFM}
\alias{simulate.CBFM_hurdle}
\alias{simulate}
\title{Simulate data from (hurdle) CBFM fit}
\usage{
\method{simulate}{CBFM}(
  object,
  nsim = 1,
  seed = NULL,
  max_resp = Inf,
  conditional = TRUE,
  ...
)

\method{simulate}{CBFM_hurdle}(
  object,
  nsim = 1,
  seed = NULL,
  max_resp = Inf,
  conditional = TRUE,
  ...
)

simulate(object, ...)
}
\arguments{
\item{object}{An object of class \code{CBFM} or \code{CBFM_hurdle}.}

\item{nsim}{A positive integer specifying the number of simulated datasets. Defaults to 1.}

\item{seed}{An integer to set seed number. Defaults to a random seed number.}

\item{max_resp}{A upper bound to limit the maximum value of responses obtained. This is useful if the user wants, say, all counts to not exceed a particular value. In such case, the function will attempt to simulate counts that do not \code{max_resp}. Note it only \emph{attempts} this: it will give up after 10 unsuccessful attempts and then return whatever is simulated on the 10-th attempt.}

\item{conditional}{If \code{conditional = TRUE}, the data are simulated conditional on the estimated species-specific regression coefficients associated with the basis functions. Otherwise if \code{conditional = FALSE} then new species-specific regression coefficients are generated from the estimated values of the \eqn{\Sigma}'s and \eqn{G}'s, and their corresponding random effects distribution. Please see the details section in \code{\link[=CBFM]{CBFM()}} for more details. Note with CBFM being set up much a like generalized additive model or GAM, then simulating conditionally is generally what most users will required, if they need to simulate.}

\item{...}{not used.}
}
\value{
A three dimensional array of dimension \eqn{N} by \eqn{m} by \code{nsim} is returned, where the number of simulated spatio-temporal multivariate abundance data sets is given by the last index.
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}

Simulate new spatio-temporal multivariate abundance data based on a fitted \code{CBFM} or \code{CBFM_hurdle} object.
}
\details{
Simulates spatio-temporal multivariate abundance data from a fitted community-level basis function model (CBFM). For the purposes of the package, the CBFM is characterized by the following mean regression model: for observational unit \eqn{i=1,\ldots,N} and species \eqn{j=1,\ldots,m}, we have

\deqn{g(\mu_{ij}) = \eta_{ij} = x_i^\top\beta_j + b_i^\top a_j,}

where \eqn{g(.)} is a known link function, \eqn{x_i} denotes a vector of predictors for unit \eqn{i} i.e., the \eqn{i}-th row from the created model matrix, \eqn{\beta_j} denotes the corresponding regression coefficients for species \eqn{j}, \eqn{b_i} denotes a vector of spatial, temporal, and/or spatio-temporal basis functions for unit \eqn{i} , and \eqn{a_j} denotes the corresponding regression coefficients for species \eqn{j}. In the function, \eqn{x_i} is created based on the \code{formula} and \code{data} arguments, \eqn{\beta_j} is supplied as part of the \code{betas} argument, and \eqn{b_i} is formed from the \code{B_space}, \code{B_time} and \code{B_spacetime} arguments. Finally, \eqn{a_j} is either supplied directly as part of \code{basis_effects_mat} argument, or generated based on the \code{Sigma} and \code{G} arguments. For hurdle CBFM models, there are two component models corresponding to the modeling the probability of presence and modeling the distribution of the data conditional on the species being present; please see \code{\link[=makeahurdle]{makeahurdle()}} for more details.

Now, as an example suppose we have a CBFM which involves spatial and temporal (but no spatio-temporal) basis functions. Then \eqn{b_i = (b_{i,space}, b_{i,time})} is formed from the \eqn{i}-th rows of \code{B_space} and \code{B_time} (which were inputs into the fitted CBFM), while \eqn{a_j = (a_{j,space}, a_{j,time})} comes from the \eqn{j}-th row of the fitted CBFM i.e., from \code{object$basis_effects_mat}. If \code{object$basis_effects_mat} is not supplied i.e., \code{conditional = FALSE}, then it is instead obtain by simulating

\deqn{(a_{1,space}, \ldots, a_{m,space}) \sim N(0, kronecker(G_{space}, \Sigma_{space})),}

where \eqn{G_{space}} and \eqn{\Sigma_{space}} are supplied from the corresponding estimates from the fitted CBFM, and \eqn{kronecker(\cdot)} is the Kronecker product operator. Similarly, we have \eqn{(a_{1,time}, \ldots, a_{m,time}) \sim N(0, kronecker(G_{time}, \Sigma_{time}))}.

By plugging estimated (or simulated) values of the parameters from the fitted CBFM into the mean model given above, responses \eqn{y_{ij}} are then subsequently simulated from the assumed distribution given this mean value, along with any estimated values of the dispersion/power/zero-inflation probability parameters from the fitted CBFM as appropriate. For hurdle CBFMs, the simulation process is analogous to this, except it consists of two steps: 1) simulating presence-absence responses first, and then simulating count responses for the presences only.
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

# Conditional on estimated basis function coefficients
simulate(fitcbfm, nseed = 5)

# Generate new basis function coefficients as part of the simulation process
simulate(fitcbfm, nseed = 1, conditional = FALSE)


# See also the examples in the help file for the makeahurdle function.
}


}
\seealso{
\code{\link[=create_CBFM_life]{create_CBFM_life()}} for manually simulating spatio-temporal multivariate abundance data.
}
\author{
Francis K.C. Hui \href{mailto:fhui28@gmail.com}{fhui28@gmail.com}, Chris Haak
}
