#' @title Extract the number of observational units from a (hurdle) CBFM fit
#' 
#' @description 
#' `r lifecycle::badge("stable")`
#' 
#' Extract the number of observational units from a fitted \code{CBFM} or \code{CBFM_hurdle} object. This is principally intended to be used in computing information criterion such as in the function [AIC.CBFM()].
#'
#' @param object An object of class \code{CBFM} or \code{CBFM_hurdle}.
#' @param ... Not used in this case.
#'
#' @details 
#' With CBFM being set up much a like generalized additive model or GAM, then this function returns the number of observational units \eqn{N} e.g., the number of rows in \code{object$y}, as opposed to the total number of elements in the response matrix \code{obejct$y}.
#' 
#' 
#' @return A positive integer.
#'
#'
#' @author Francis K.C. Hui <fhui28@gmail.com>, Chris Haak
#'
#' @seealso [CBFM()] for fitting CBFMs, [AIC.CBFM()] and [AICc.CBFM()] for calculating various information criteria from a CBFM fit, and [logLik.CBFM()] for extracting the log-likelihood value from a CBFM fit.
#'
#' @examples
#' \dontrun{
#' library(autoFRK)
#' library(FRK)
#' library(MASS)
#' library(mvabund)
#' library(mvtnorm)
#' library(ROCR)
#' library(sp)
#' library(geoR)
#' library(tidyverse)
#' 
#' ##------------------------------
#' ## **Example 1: Fitting a CBFM to spatial multivariate presence-absence data** 
#' ## simulated from a spatial latent variable model
#' ## Please note the data generation process (thus) differs from CBFM.
#' ##------------------------------
#' set.seed(2021)
#' num_sites <- 500 # 500 (units) sites 
#' num_spp <- 50 # Number of species
#' num_X <- 4 # Number of regression slopes
#' 
#' spp_slopes <- matrix(runif(num_spp * num_X, -1, 1), nrow = num_spp)
#' spp_intercepts <- runif(num_spp, -2, 0)
#' 
#' # Simulate spatial coordinates and environmental covariate components
#' xy <- data.frame(x = runif(num_sites, 0, 5), y = runif(num_sites, 0, 5))
#' X <- mvtnorm::rmvnorm(num_sites, mean = rep(0,4))
#' colnames(X) <- c("temp", "depth", "chla", "O2")
#' dat <- data.frame(xy, X)
#' mm <- model.matrix(~ temp + depth + chla + O2 - 1, data = dat) %>% 
#' scale %>% 
#' as.matrix
#' 
#' # Simulate latent variable component
#' true_lvs <- grf(grid = cbind(xy$x, xy$y), nsim = 2, cov.model = "exponential",
#' cov.pars = c(1, 2))$data %>%
#'      as.matrix
#' spp_loadings <- matrix(runif(num_spp * 2, -1, 1), nrow = num_spp)
#' set.seed(NULL)
#' 
#' # Simulate spatial multivariate abundance data (presence-absence)
#' eta <- tcrossprod(cbind(1,mm), cbind(spp_intercepts,spp_slopes)) + 
#' tcrossprod(true_lvs, spp_loadings)
#' simy <- matrix(rbinom(num_sites * num_spp, size = 1, 
#' prob = binomial()$linkinv(eta)), nrow = num_sites)
#' rm(X, mm, spp_loadings, true_lvs, xy, eta)
#' 
#' 
#' # Set up spatial basis functions for CBFM -- Most users will start here! 
#' num_basisfunctions <- 25 # Number of spatial basis functions to use
#' basisfunctions <- mrts(dat[,c("x","y")], num_basisfunctions) %>% 
#' as.matrix %>%
#' {.[,-(1)]} # Remove the first intercept column
#' 
#' # Fit CBFM 
#' useformula <- ~ temp + depth + chla + O2
#' fitcbfm <- CBFM(y = simy, formula = useformula, data = dat, 
#' B_space = basisfunctions, family = binomial(), control = list(trace = 1))
#' 
#' nobs(fitcbfm)
#'}
#'
#' @export
#' @md

nobs.CBFM <- function(object, ...) {
    if(!inherits(object, "CBFM")) 
        stop("`object' is not of class \"CBFM\"")
     
     n <- nrow(object$y)
     return(n)
     }


#' @rdname nobs.CBFM
#' @export
nobs.CBFM_hurdle <- function(object, ...) {
    if(!inherits(object, "CBFM_hurdle")) 
        stop("`object' is not of class \"CBFM_hurdle\"")
     
     n <- nrow(object$pa_fit$y)
     return(n)
     }


#' @rdname nobs.CBFM
#' @export
nobs <- function(object, ...) {
    UseMethod("nobs")
    }
