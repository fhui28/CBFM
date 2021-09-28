#' @title Log-likelihood of a CBFM fit
#' 
#' @description 
#' `r lifecycle::badge("experimental")`
#' 
#' Extracts the log-likelihood from a fitted \code{CBFM} object.
#'
#' @param object An object of class "CBFM".
#' @param use_edf If \code{TRUE}, then the estimated degrees of freedom for the species-specific coefficients related to the spatial and/temporal basis functions is used instead. Defaults to \code{FALSE}, in which case species-specific coefficients related to the basis functions are regarded as fixed effects.
#' @param ... Not used in this case.
#'
#' @details 
#' In this package, CBFM are fitted using maximized penalized quasi-likelihood or PQL estimation. In turn, this function returns the maximized log-likelihood value (*excluding* the quadratic penalty term in the PQL) of the \code{CBFM} object at convergence. 
#' 
#' By defautl, the degrees of freedom calculated as part of this function is based on: 1) using the estimated of effective degrees of freedom for the component of the model related to \code{object$formula_X} plus any nuisance parameters (see [mgcv::logLik.gam()] for details about estimated degrees of freedom when smoothing terms are involved); 2) treating the species-specific regression coefficients related to the spatial and/or temporal basis functions as fixed effects. Overall, this means the calculation of degrees of freedom only involves the number of regression coefficients in the model (see Hui et al., 2017, Hui, 2021, for justificatons for these in the context of mixed models). Alternatively, if \code{use_edf = TRUE}, then point 2 is modified to instead use the estimated degrees of freedom instead. This is done by making a call to [edf.CBFM()] and we refer to the corresponding help file for more information.
#' 
#' The community-level covariance matrices are not involved in the calculation of the degrees of freedom. 
#' 
#' @return Returns an object of class "logLik".  This is a number with at least one attribute, "df" (degrees of freedom), giving the number of (estimated) parameters in the CBFM model; please see see [stats::logLik()].
#'
#' @author Francis K.C. Hui <fhui28@gmail.com>, Chris Haak
#'
#' @references 
#' Hui, F. K. C., Mueller, S., and Welsh, A. H. (2017). Joint selection in mixed models using regularized PQL. Journal of the American Statistical Association, 112, 1323-1333.
#' 
#' Hui, F. K. C.(2021). On the use of a penalized quasilikelihood information criterion for generalized linear mixed models. Biometrika, 108, 353-365.
#' 
#' @seealso [CBFM()] for fitting CBFMs, and [AIC.CBFM()] and [AICc.CBFM()] for calculation various information criteria from a CBFM fit.
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
#' logLik(fitcbfm)
#' 
#' logLik(fitcbfm, use_edf = TRUE) # Degrees of freedom for this are much smaller
#'}
#'
#' @export
#' @md

logLik.CBFM <- function(object, use_edf = FALSE, ...) {
    if(!inherits(object, "CBFM")) 
        stop("`object' is not of class \"CBFM\"")

     logL <- object$logLik
     num_params <- min(sum(object$edf), nrow(object$betas)*ncol(object$betas)) + nrow(object$basis_effects_mat)*ncol(object$basis_effects_mat)
     if(use_edf) {
         getedf <- edf.CBFM(object)
         getedf <- sum(getedf[(nrow(getedf) - sum(object$which_B_used) + 1):nrow(getedf),])
         num_params <- min(sum(object$edf), nrow(object$betas)*ncol(object$betas)) + min(getedf, nrow(object$basis_effects_mat)*ncol(object$basis_effects_mat))
        }
     if(object$family$family %in% c("Beta","gaussian","Gamma","negative.binomial","tweedie","ztnegative.binomial"))                       
          num_params <- num_params + length(object$dispparam)
     if(object$family$family %in% c("tweedie"))                        
          num_params <- num_params + length(object$powerparam)
     if(object$family$family %in% c("zipoisson", "zinegative.binomial"))                        
          num_params <- num_params + length(object$zeroinfl_prob)
     
     
     attributes(logL)$df <- num_params
     attributes(logL)$nobs <- dim(object$y)[1]
     class(logL) <- "logLik"
     return(logL)
     }
