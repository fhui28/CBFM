#' @title Akaike's An Information Criterion (AIC) for a CBFM fit
#' 
#' @description 
#' `r lifecycle::badge("experimental")`
#' 
#' Calculates Akaike's "An Information Criterion" (AIC, Akaike, 1974) from a fitted \code{CBFM} object. This can also be generalized to other information criterion by modifying the model complexity penalty \code{k}.
#'
#' @param object An object of class "CBFM".
#' @param k The model complexity penalty to use in the calculation of the information criterion. Defaults to \code{k = 2}, which is the classical AIC.
#' @param ... Not used in this case.
#'
#' @details 
#' While the default returns the much celebrated AIC for a fitted CBFM, using a default complexity penalty of 2, the user is free to modify the model complexity penalty according to whatever information criterion they wish to calculate. Another common choice is the Bayesian Information criterion (BIC, Schwarz, 1978), where \code{k = log(nobs(object))}. 
#' 
#' The generic form of the information criterion this function uses is \eqn{-2 \times \ell + k \times df}, where \eqn{\ell} is the maximized log-likelihood value of the \code{CBFM} object at convergence and \eqn{df} is the (estimated) degrees of freedom; please [logLik.CBFM()] for more details.
#' 
#' 
#' @return A numeric value of the calculated information criterion.
#'
#' @author Francis K.C. Hui <fhui28@gmail.com>, Chris Haak
#' 
#' @references 
#' Akaike, H. (1974). A new look at the statistical model identification. IEEE Transactions on Automatic Control, 19, 716-723.
#' 
#' Schwarz, G. (1978). Estimating the dimension of a model. The Annals of Atatistics, 6, 461-464.
#' 
#' 
#' @seealso [AICc.CBFM()] for calculating the corrected AIC specifically, [logLik.CBFM()] for extracting the log-likelihood value from a CBFM fit, and [nobs.CBFM()] for extracting the number of observational units from a CBFM fit.
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
#' ## Example 1: Fitting a CBFM to spatial multivariate presence-absence data 
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
#' # Set up spatial basis functions for CBFM -- Most practitioners will start here! 
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
#' 
#' AIC(fitcbfm) # AIC
#' 
#' AIC(fitcbfm, k = log(nobs(fitcbfm))) # BIC
#'}
#'
#' @export
#' @md
#' 
AIC.CBFM <- function(object, k = 2, ...) {
     if(!inherits(object, "CBFM")) 
          stop("`object' is not of class \"CBFM\"")

     get_dfs <- attributes(logLik.CBFM(object))$df
     
     IC <- -2 * object$logLik + k * get_dfs
     return(IC)
     }


# #' @export
#AIC <- function(object, ...) {
#     UseMethod(generic = "AIC")
#     }
