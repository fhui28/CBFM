#' @title Extract parameter (estimator) covariance matrix from CBFM fit
#' 
#' @description 
#' `r lifecycle::badge("stable")`
#'
#' Extracts the so-called Bayesian posterior covariance matrix of the parameters estimators from a \code{CBFM} fit.
#' 
#' @param object An object of class \code{CBFM}.
#' @param ... Not used.
#' 
#' @details 
#' 
#' This basically extracts \code{object$covar_components}, assuming \code{stderrors = TRUE} in the original CBFM fit. Please note that to save memory, rather than returning the full Bayesian posterior covariance matrix, only relevant sub-blocks that are needed for standard errors and inference/prediction are returned. 
#' 
#' 
#' @return A list containing with the following components: 
#' \item{topleft: }{A matrix corresponding to the top-left block of the full Bayesian posterior covariance matrix. The top-left block specifically relates to the regression coefficients associated with the measured predictors i.e., the covariance matrix associated with \code{object$betas}, and the species-specific zero-inflated probabilities on the logit scale if the response distribution involved one.}

#' \item{topright: }{matrix of the top-right block of the full Bayesian posterior covariance matrix. The top-right block specifically relates to the cross-covariance of the regression coefficients associated with the measured predictors (plus the species-specific zero-inflated probabilities on the logit scale) and the basis functions i.e., the cross-covariance matrix between \code{object$betas} and \code{object$basis_effects_mat}.}

#' \item{bottomright: }{A matrix containing components of the bottom-right block of the full Bayesian posterior covariance matrix. The bottom-left block specifically relates to the regression coefficients associated with the basis functions i.e., the covariance matrix associated with \code{object$basis_effects_mat}.} 
#' 
#' 
#' @author Francis K.C. Hui <fhui28@gmail.com>, Chris Haak
#' 
#' @seealso [CBFM()] for fitting CBFMs, [coef.CBFM()] to obtain standard errors and confidence interval limits in a more user-friendly form.
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
#' fitcbfm <- CBFM(y = simy, formula = useformula, data = dat, 
#' B_space = basisfunctions, family = binomial(), control = list(trace = 1))
#' 
#' vcov(fitcbfm)
#' }
#' 
#' @export
#' @md
#' 
vcov.CBFM <- function(object, ...) {
    if(!inherits(object, "CBFM")) 
        stop("`object' is not of class \"CBFM\"")
    
    if(object$stderrors == FALSE)
        stop("Standard errors were not asked for and thus not calculated as part of the the CBFM object. Sorry!")

    return(object$covar_components)
    }

# #' @export
#vcov <- function(object, ...) {
#    UseMethod(generic = "vcov")
#    }     
