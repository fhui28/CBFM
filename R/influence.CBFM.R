#' @title Regression diagnostics from a CBFM fit
#' 
#' @description 
#' `r lifecycle::badge("stable")`
#'
#' Takes a fitted \code{CBFM} object and calculates both the hat values i.e., diagonal elements of the influence/hat matrix, and an approximate Cook's distance.  
#'
#' @param object An object of class "CBFM".
#' @param ncores To speed up calculation of the influence measures, parallelization can be performed, in which case this argument can be used to supply the number of cores to use in the parallelization. Defaults to \code{detectCores()-1}.
#'
#' @details 
#' This is a pretty simple function that calculates the so-called "hat values" and an approximate Cook's distance from a fitted CBFM. Classically, the former is used to identify potentially high-leverage observations, while the latter is a well-known regression-deletion diagnostic that measures the overall influence of a point on the fit of a model; we refer the reader to [stats::influence.measures()] for more information and reference pertaining to its use in standard linear and other regression models. Because both concepts of both can be (sort of) carried over to a generalized additive model or GAM e.g., see [mgcv::influence.gam()], then analogous diagnostics can be produced for a CBFM. 
#' 
#' Note hat values and Cook's distance are obtained on a *per-species basis*. This makes sense since the influence/leverage of a unit e.g., space-time coordinate, will differ for different species and their relationships to the measured covariates. 
#' 
#' ** We leave it up to the practitioner to decide how to use the regression diagnostics, if at all.** In a GAM let alone a CBFM for spatio-temporal multivariate abundance data these diagnostics may only be approximate, and so providing rules-of-thumb related to their usage is challenging. Besides, we echo the sentiment provided in Chapter 4.4 of Fox et al., (2019) that cutoffs and rules-of-thumb should not be given too much weight, with more attention placed on graphical displays and assessing *relative infuence* of observations.   
#'  
#'  
#' @return A list containing two elements:
#' \describe{
#' \item{hat: }{A matrix of estimated hat values i.e., diagonal elements of the influence/hat matrix. The dimensions of this matrix should be the same as \code{object$y}.}
#' \item{cooks: }{A matrix of estimated and approximate Cook's distances. The dimensions of this matrix should be the same as \code{object$y}. }
#' }
#'
#' @author Francis K.C. Hui <fhui28@gmail.com>, Chris Haak
#'
#' @references 
#'Fox, J. (2019). Regression diagnostics: An introduction. Sage publications.
#'   
#' @seealso [CBFM()] for fitting CBFMs, [residuals.CBFM()] for calculating different types of residuals from a CBFM fit, and [plot.CBFM()] for basic residual diagnostics.
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
#' influence(fitcbfm)
#' }
#'
#' @export
#' @import foreach  
#' @import Matrix
#' @importFrom doParallel registerDoParallel
#' @md

# Interesting note that in simulations, most of the time, the results from this are not too far from (but very slightly more conservative than) just ad-hoc applying summary.gam to the last iteration of the PQL estimation algorithm in CBFM! 
influence_CBFM <- function(object, ncores = NULL) {
  if(!inherits(object, "CBFM")) 
    stop("`object' is not of class \"CBFM\"")
  
  if(is.null(ncores))
    registerDoParallel(cores = detectCores()-1)
  if(!is.null(ncores))
    registerDoParallel(cores = ncores)

  if(!object$stderrors)
    stop("Standard errors must have been produced from `object' for influence measures to work. Sorry!")

  # (X^T W X + S)^{-1}== Bayesian posterior covariance matrix
  bigV <- cbind(object$covar_components$topleft, object$covar_components$topright)
  bigV <- rbind(bigV, cbind(t(object$covar_components$topright), object$covar_components$bottomright))
  
  # Square root of Weights
  num_spp <- nrow(object$betas)
  num_units <- nrow(object$B)
  if(is.null(object$dispparam))
      object$dispparam <- rep(1, num_spp)
  if(is.null(object$powerparam))
      object$powerparam <- rep(0, num_spp)
  if(is.null(object$zeroinfl_prob_intercept))
      object$zeroinfl_prob_intercept <- rep(0, num_spp)
  
  weights_mat <- .neghessfamily(family = object$family, eta = object$linear_predictor, y = object$y, 
                                phi = matrix(object$dispparam, num_units, num_spp, byrow = TRUE), 
                                powerparam = matrix(object$powerparam, num_units, num_spp, byrow = TRUE),
                                zeroinfl_prob_intercept = matrix(object$zeroinfl_prob_intercept, num_units, num_spp, byrow = TRUE), 
                                trial_size = object$trial_size, domore = TRUE)
  if(!(object$family$family[1] %in% c("zipoisson","zinegative.binomial")))
    weights_mat <- matrix(weights_mat$out, nrow = num_units, ncol = num_spp) # Overwrite weights_mat since only one quantity needed
  if(object$family$family[1] %in% c("zipoisson","zinegative.binomial"))
    weights_mat_betabeta <- matrix(weights_mat$out, nrow = num_units, ncol = num_spp)

  # [W^{1/2}X, W^{1/2}B]
  X <- model.matrix.CBFM(object)
  WsqrtXB <- function(j) {                
    if(!(object$family$family[1] %in% c("zipoisson","zinegative.binomial"))) {
      WsqrtX <- X*sqrt(weights_mat[,j])
      WsqrtB <- object$B * sqrt(weights_mat[,j])
      return(list(WsqrtX = WsqrtX, WsqrtB = WsqrtB))
      }
    
    if(object$family$family[1] %in% c("zipoisson","zinegative.binomial")) {
      Xi <- bdiag(matrix(1, num_units, 1), X)
      bigW <- cbind(Diagonal(x = sqrt(weights_mat$out_zeroinflzeroinfl[,j])),  Diagonal(x = sqrt(weights_mat$out_zeroinflbetas[,j])))
      bigW <- rbind(bigW, cbind(Diagonal(x = sqrt(weights_mat$out_zeroinflbetas[,j])),  Diagonal(x = sqrt(weights_mat_betabeta[,j]))))
      WsqrtX <- Xi %*% sqrt(bigW)
      WsqrtB <- object$B * sqrt(weights_mat_betabeta[,j])
      return(list(WsqrtX = WsqrtX, WsqrtB = WsqrtB))
      }
    }
  
  getall_WsqrtXB <- foreach(j = 1:num_spp) %dopar% WsqrtXB(j = j)
  gc()
  rm(X)
  
  bigsqrtWXB <- cbind(bdiag(lapply(getall_WsqrtXB, function(x) x$WsqrtX)), bdiag(lapply(getall_WsqrtXB, function(x) x$WsqrtB))) 
  
  
  # hat values -- This is a big bottleneck!
  gethatvals_j <- function(j) {
    sel_units <- (j*num_units - num_units + 1):(j*num_units)
    hatvals_j <- rowSums((bigsqrtWXB[sel_units,] %*% bigV) * bigsqrtWXB[sel_units,]) 
    
    return(hatvals_j)
    }
  #hatvals <- rowSums((bigsqrtWXB %*% bigV) * bigsqrtWXB)
  #hatvals <- matrix(hatvals, nrow = num_units, ncol = num_spp)
  hatvals <- foreach(j = 1:num_spp, .combine = cbind) %dopar% gethatvals_j(j = j)
  rm(bigsqrtWXB, bigV, weights_mat)    

  
  # Cook's distance
  res <- residuals.CBFM(object, type = "pearson")
  final_dfs <- colSums(object$edf) + nrow(object$basis_effects_mat)
  cookD <- (res / (1 - hatvals))^2 * hatvals / matrix(final_dfs, nrow = num_units, ncol = num_spp, byrow = TRUE)
  
  rownames(hatvals) <- rownames(cookD) <- rownames(object$y)
  colnames(hatvals) <- colnames(cookD) <- colnames(object$y)
  
  return(list(hat = hatvals, cooks = cookD))
  }

