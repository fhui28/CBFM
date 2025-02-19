#' @title Extract model coefficients from a CBFM fit
#' 
#' @description 
#' `r lifecycle::badge("stable")`
#' 
#' Extracts the estimated regression coefficients associated with the covariates from a fitted \code{CBFM} object.
#'
#' @param object An object of class \code{CBFM}.
#' @param ... Not used in this case.
#'
#' @details 
#' For the purposes of the package, the CBFM is characterized by the following mean regression model: for observational unit \eqn{i=1,\ldots,N} and species \eqn{j=1,\ldots,m}, we have
#' 
#' \deqn{g(\mu_{ij}) = \eta_{ij} = x_i^\top\beta_j + b_i^\top a_j,}
#'
#' where \eqn{g(.)} is a known link function, \eqn{x_i} denotes a vector of predictors for unit \eqn{i} i.e., the \eqn{i}-th row from the created model matrix, \eqn{\beta_j} denotes the corresponding regression coefficients for species \eqn{j}, \eqn{b_i} denotes a vector of spatial, temporal, and/or spatio-temporal basis functions for unit \eqn{i} , and \eqn{a_j} denotes the corresponding regression coefficients for species \eqn{j}. 
#' 
#' This function will extract the estimated coefficients \eqn{\hat{\beta}_j}'s from the fitted CBFM, noting that this may included the estimated smoothing coefficients if any smoothers were included. For zero-inflated distributions, it will also return the estimated coefficients associated with modeling the probability of zero-inflation, noting that this may included the estimated smoothing coefficients if any smoothers were included.
#' 
#' This function does \emph{not} return the estimated regression coefficients associated with the basis functions i.e., the \eqn{\hat{a}_j}'s. These can be obtained from \code{object$basis_effects_mat}.
#'
#' @return A matrix of estimated species-specific regression coefficients corresponding to the model matrix created, where the number of rows is equal to the number of species. For zero-inflated distributions, it returns a list containing both the matrix of estimated species-specific regression coefficients corresponding to the model matrix created, and a vector of estimated species-specific probabilities of zero-inflation on the logit scale.
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
#' prob = plogis(eta)), nrow = num_sites)
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
#' coef(fitcbfm)
#'}
#'
#' @export
#' @md

coef.CBFM <- function(object, ...) {
    if(!inherits(object, "CBFM")) 
        stop("`object' is not of class \"CBFM\"")

     out <- object$betas
     if(object$family$family[1] %in% c("zipoisson", "zinegative.binomial"))
          out <- list(betas = object$betas, zibetas = object$zibetas)
     return(out)
     }
