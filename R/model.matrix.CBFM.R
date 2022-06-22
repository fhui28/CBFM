#' @title Extract model matrices (otherwise known as design matrices) from a CBFM fit
#' 
#' @description 
#' `r lifecycle::badge("stable")`
#' 
#' Obtains the model matrix from a fitted \code{CBFM} object. This is especially useful when the fitted CBFM includes smoothing terms, say, in which case, the function will return the precise model matrix used. 
#' 
#' @param object An object of class \code{CBFM}.
#' @param zi For zero-inflated distributions, set this to \code{TRUE} if the model matrix associated with modeling the probabilities of zero-inflation is desired.
#' @param ... Not used.
#' 
#' @details 
#' Similar to how [mgcv::model.matrix.gam()] works, it calls [mgcv::predict.gam()] with no \code{newdata} argument and \code{type = "lpmatrix"} in order to obtain the model matrix. Note this is the model matrix associated with the covariates i.e., based on arguments \code{object$formula}, and **not** the basis functions. 
#' 
#' For zero-inflated distributions, it can also obtain the model matrix associated with modeling the probabilities of zero-inflation.
#' 
#' @return A model matrix.
#' 
#' @author Francis K.C. Hui <fhui28@gmail.com>, Chris Haak
#' 
#' @seealso [CBFM()] for fitting CBFMs, [fitted.CBFM()] for extracting fitted values from a CBFM fit, [residuals.CBFM()] for calculating various types of residuals, and [predict.CBFM()] for constructing predictions from a CBFM fit.
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
#' model.matrix(fitcbfm)
#' }
#' 
#' @export
#' @importFrom mgcv gam predict.gam
#' @md

model.matrix.CBFM <- function(object, zi = FALSE, ...) {
     if(!inherits(object, "CBFM")) 
          stop("`object' is not of class \"CBFM\"")
        
     tmp_formula <- as.formula(paste("response", paste(as.character(object$formula),collapse="") ) )
     nullfit <- gam(tmp_formula, data = data.frame(response = runif(nrow(object$y)), object$data), fit = TRUE, control = list(maxit = 1))
     MM <- predict.gam(nullfit, type = "lpmatrix", ...)
        
     if(zi) {
          if(!(object$family$family[1] %in% c("zipoisson","zinegative.binomial")))
               stop("A model model associated with the probability of zero-inflation can be only be obtained for zero-inflated CBFMs.")
          
          tmp_formula <- as.formula(paste("response", paste(as.character(object$ziformula),collapse="") ) )
          nullfit <- gam(tmp_formula, data = data.frame(response = runif(nrow(object$y)), object$data), fit = TRUE, control = list(maxit = 1))
          MM <- predict.gam(nullfit, type = "lpmatrix", ...)
          }
        return(MM)
        }
     
