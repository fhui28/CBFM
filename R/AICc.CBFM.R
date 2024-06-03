#' @title Corrected An Information Criterion (AICc) for a (hurdle) CBFM fit
#' 
#' @description 
#' `r lifecycle::badge("experimental")`
#' 
#' Calculates the corrected AIC (AICc, Cavanaugh, 1997; Burnham and Anderson, 2002) from a fitted \code{CBFM} or \code{CBFM_hurdle} object. 
#'
#' @param object An object of class \code{CBFM} or \code{CBFM_hurdle}.
#' @param use_edf If \code{TRUE}, then the estimated degrees of freedom for the species-specific coefficients related to the spatial and/temporal basis functions is used instead. Defaults to \code{FALSE}, in which case species-specific coefficients related to the basis functions are regarded as fixed effects.
#' @param ... Not used in this case.
#'
#' @details 
#' The corrected Akaike Information Criterion (AICc) is a correction of AIC for small sample sizes. For a CBFM, and adapting the idea from [gllvm::AICc()], it takes the form
#' \deqn{-2 \times \ell + k \times df + \frac{2 \times df \times (df+1)}{Nm - df - 1},}
#' where \eqn{\ell} is the maximized log-likelihood value (*excluding* the quadratic penalty term in the PQL) of the \code{CBFM} object at convergence \eqn{df} is the (estimated) degrees of freedom (please see [logLik.CBFM()] for more details, especially regarding the use of the argument \code{use_edf}), \eqn{N} is the number of observational units, and \eqn{m} is the number of species. In light of the use of the maximized log-likelihood value, one may consider that this function constructs something more akin to a conditional AICc; see the discussion in [mgcv::logLik.gam()] and references therein.
#' 
#' Basically, AICc is essentially AIC with an extra penalty term for the number of parameters in the model. Note that as the number of observational units increase, this extra penalty term converges to zero, and thus AICc converges to AIC. 
#' 
#' As an alternative to using information criteria, CBFM also has available built-in approaches for smoothing term (but not parametric term) selection via shrinkage smoothers or null space penalization; please see the \code{select} argument in the [CBFM()] help file as well as [mgcv::gam.selection()] for more information.  
#' 
#' 
#' @return A numeric value of the calculated corrected AIC.
#'
#'
#' @author Francis K.C. Hui <fhui28@gmail.com>, Chris Haak
#' 
#' @references 
#' Burnham, K. P., and Anderson, D. R. (2002). Model selection and multimodel inference: A practical information-theoretic approach (2nd ed.). Springer-Verlag.
#' 
#' Cavanaugh, J. E. (1997). Unifying the derivations for the Akaike and corrected Akaike information criteria. Statistics and Probability Letters, 33, 201-208.
#' 
#' 
#' @seealso [AIC.CBFM()] for calculating AIC and other information criteria and [logLik.CBFM()] for extracting the log-likelihood value from a CBFM fit.
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
#' X <- rmvnorm(num_sites, mean = rep(0,4)) 
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
#' AICc(fitcbfm)
#' 
#' 
#' # See also the examples in the help file for the makeahurdle function.
#'}
#'
#' @aliases AICc AICc.CBFM AICc.CBFM_hurdle 
#' @export
#' @export AICc.CBFM
#' @md

AICc.CBFM <- function(object, use_edf = FALSE, ...) {
     if(!inherits(object, "CBFM")) 
          stop("`object' is not of class \"CBFM\"")

     N <- nrow(object$y)
     m <- ncol(object$y)
     get_dfs <- attributes(logLik.CBFM(object, use_edf = use_edf))$df
     
     IC <- -2 * object$logLik + 2 * get_dfs + 2 * get_dfs * (get_dfs + 1) / (N*m - get_dfs - 1)
     return(IC)
     }

#' @rdname AICc.CBFM
#' @export
AICc.CBFM_hurdle <- function(object, use_edf = FALSE, ...) {
        if(!inherits(object, "CBFM_hurdle")) 
                stop("`object' is not of class \"CBFM_hurdle\"")

        N <- nrow(object$pa_fit$y)
        N2 <- sum(object$count_fit$y > 0)
        m <- ncol(object$pa_fit$y)
        IC <- logLik.CBFM_hurdle(object, use_edf = use_edf)
        get_dfs <- attributes(IC)$df
     
     IC <- -2 * as.vector(IC) + 2 * get_dfs + 2 * get_dfs * (get_dfs + 1) / (N*m + N2 - get_dfs - 1)
     return(IC)
     }


#' @rdname AICc.CBFM
#' @export AICc 
AICc <- function(object, ...) {
     UseMethod("AICc")
     }
