#' @title Extract fitted values from a (hurdle) CBFM fit
#' 
#' @description 
#' `r lifecycle::badge("stable")`
#' 
#' Extracts fitted mean values from a \code{CBFM} or \code{CBFM_hurdle} fit.
#' 
#' @param object An object of class \code{CBFM} or \code{CBFM_hurdle}.
#' @param ... Not used.
#' 
#' @details 
#' To clarify, the returned fitted values are on the response scale i.e., a matrix of the estimated means \eqn{\hat{\mu}_{ij}} after model fitting. Note that for zero-inflated distributions, while the mean of the non-zero-inflated component is modeled in CBFM, the fitted values are the *actual expected mean values* i.e., it returns estimated values of \eqn{(1-\pi_j)*\mu_{ij}} where \eqn{\pi_j} is the species-specific probability of zero inflation and \eqn{\mu_{ij}} is the mean of the non-zero-inflated component. 
#' 
#' Similarly, for zero-truncated count distributions, while the mean of the base count distribution is modeled in CBFM, the fitted values are the *actual expected mean values* i.e., it returns estimated values of \eqn{\mu_{ij}/(1-p(0,\mu_{ij}))} where \eqn{\mu_{ij}} is the mean of the base count distribution component and \eqn{p(0,\mu_{ij})} generically denotes the probability of observing a zero count for the base count distribution (and it returns \code{NA} values for elements corresponding to zero counts in \code{object$y}.
#' 
#' For hurdle CBFMs, the function returns the estimated values of \eqn{\pi_{ij}*\mu_{ij}} where \eqn{\pi_{ij}} is the probability of observing a presence based on the presence-absence component of the model, and \eqn{\mu_{ij}} is the mean of the zero-truncated component of the model.  
#' 
#' @return A matrix of fitted values.
#' 
#' @author Francis K.C. Hui <fhui28@gmail.com>, Chris Haak
#' 
#' @seealso [CBFM()] for fitting CBFMs, [model.matrix.CBFM()] for extracting the model matrix from a CBFM fit, [predict.CBFM()] for constructing predictions from a CBFM fit, and [residuals.CBFM()] for calculating various types of residuals.
#' 
#' @examples
#' \donttest{
#' library(autoFRK)
#' library(FRK)
#' library(MASS)
#' library(mvabund)
#' library(mvtnorm)
#' library(ROCR)
#' library(sp)
#' library(RandomFields)
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
#' X <- rmvnorm(num_sites, mean = rep(0,4)) 
#' colnames(X) <- c("temp", "depth", "chla", "O2")
#' dat <- data.frame(xy, X)
#' mm <- model.matrix(~ temp + depth + chla + O2 - 1, data = dat) %>% 
#' scale %>% 
#' as.matrix
#' 
#' # Simulate latent variable component
#' true_lvs <- RFsimulate(model = RMexp(var=1, scale=2), 
#' x = xy$x, y = xy$y, n = 2)@data %>% 
#' as.matrix
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
#' fitcbfm <- CBFM(y = simy, formula_X = useformula, data = dat, 
#' B_space = basisfunctions, family = binomial(), control = list(trace = 1))
#' 
#' fitted(fitcbfm)
#' 
#' 
#' # See also the examples in the help file for the makeahurdle function.
#' }
#' 
#' @aliases fitted.CBFM fitted.CBFM_hurdle
#' @export
#' @md

fitted.CBFM <- function(object, ...) {
    if(!inherits(object, "CBFM")) 
        stop("`object' is not of class \"CBFM\"")
    
    return(object$fitted)
     }


#' @rdname fitted.CBFM
#' @export
fitted.CBFM_hurdle <- function(object, ...) {
     if(!inherits(object, "CBFM_hurdle")) 
        stop("`object' is not of class \"CBFM_hurdle\"")

     if(object$count_fit$family$family[1] %in% c("ztpoisson", "ztnegative.binomial"))
          out <- object$pa_fit$fitted * predict.CBFM(object$count_fit, type = "response")
     
     return(out)
     }

