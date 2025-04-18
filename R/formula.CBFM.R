#' @title Extract the formula used from a CBFM fit
#' 
#' @description 
#' `r lifecycle::badge("stable")`
#' 
#' Extract the \code{formula}, and \code{ziformula} if appropriate, argument from a fitted \code{CBFM} object. 
#'
#' @param x An object of class \code{CBFM}.
#' @param ... Not used in this case.
#'
#' @details 
#' This function works in a similar manner to functions such as [stats::formula.lm()] and [mgcv::formula.gam()]. Recall in the main CBFM fitting function, the argument \code{formula} is a symbolic description of the model matrix of covariates to be created, while for zero-inflated distributions \code{ziformula} is a symbolic description of the model matrix to be created for the zero-inflation component.
#' 
#' Formulas based on generalized additive models or GAMs are permitted (at least, for the smoothing terms we have tried so far!); please see [mgcv::formula.gam()] for more details. 
#' 
#' @return A object of class \code{formula}. Note there will be nothing on the left hand side of the "~". 
#'
#' @author Francis K.C. Hui <fhui28@gmail.com>, Chris Haak
#'
#' @seealso [CBFM()] for fitting CBFMs.
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
#' formula(fitcbfm)
#'}
#'
#' @export
#' @md
#' 
formula.CBFM <- function(x, ...) {
    if(!inherits(x, "CBFM")) 
        stop("`x' is not of class \"CBFM\"")
     
     out <- x$formula
     
     if(x$family$family[1] %in% c("zipoisson", "zinegative.binomial"))
          out <- list(formula = x$formula, ziformula = x$ziformula)
     
     return(out)
     }

