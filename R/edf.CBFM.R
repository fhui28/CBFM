#' @title Extract estimated/effective degrees of freedom associated with a CBFM fit 
#' 
#' @description 
#' `r lifecycle::badge("experimental")`
#'
#' Takes a fitted \code{CBFM} object and calculates the estimated or effective degrees of freedom associated with each smoothing terms included in the model, as well as for each of the spatial and/or temporal basis functions included in the model.
#'
#' @param object An object of class \code{CBFM}.
#' @param ncores To speed up calculation of the estimated degrees of freedom, parallelization can be performed, in which case this argument can be used to supply the number of cores to use in the parallelization. Defaults to \code{detectCores()-1}.
#' @param ... Not used.
#'
#' @details 
#' For the estimated or effective of freedom (EDF) associated with any smoothing terms included in the model as part of \code{object$formula_X}, the function directly returns what is available from \code{object$pen.edf}. Note that as stated in the [CBFM()], these values are pulled straight from the GAM part of the estimation algorithm and consequently may only be *very* approximate. 
#' 
#' For the EDF associated with the spatial and/or temporal basis function coefficients, similar to [influence.CBFM()]  for each species up to three EDFs are given depending on which of \code{B_space/B_time/B_spacetime} are included in the model. Note that because of the way the CBFM is set up, there is usually a considerable amount of penalization taking place for regression coefficients corresponding to the spatial and/or temporal basis functions, and so one should expect these value to usually be *much* smaller than the corresponding number of basis functions included in the model. 
#'  
#' @return A matrix of species-specific EDFs, with the the number of columns equal to the number of columns in \code{object$y}, while the number of rows depends on the number of smoothing terms included in \code{object$formula_X} and which of \code{B_space/B_time/B_spacetime} were included in the model.
#' 
#'
#' @author Francis K.C. Hui <fhui28@gmail.com>, Chris Haak
#'
#' @seealso [CBFM()] for fitting CBFMs, [influence.CBFM()] for computing some influence measures including the estimates degrees of freedom for the species-specific regression coefficients corresponding to spatial and/or temporal basis functions included.
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
#' useformula <- ~ temp + s(depth) + s(chla) + s(O2) # Temperature is not "smoothed"
#' fitcbfm <- CBFM(y = simy, formula_X = useformula, data = dat, 
#' B_space = basisfunctions, family = binomial(), control = list(trace = 1))
#' 
#' edf(fitcbfm)
#' }
#'
#' @aliases edf edf.CBFM
#' @export
#' @export edf.CBFM
#' @importFrom foreach foreach %dopar%
#' @import Matrix 
#' @importFrom doParallel registerDoParallel
#' @md

edf.CBFM <- function(object, ncores = NULL, ...) {
  if(!inherits(object, "CBFM")) 
    stop("`object' is not of class \"CBFM\"")
  
  if(!object$stderrors)
    stop("Standard errors must have been produced from `object' for influence measures to work.")
  
  if(is.null(ncores))
    registerDoParallel(cores = detectCores()-1)
  if(!is.null(ncores))
    registerDoParallel(cores = ncores)

  if(is.null(object$pen_edf[[1]]))
    edf_part1 <- NULL
  if(!is.null(object$pen_edf[[1]]))
    edf_part1 <- simplify2array(object$pen_edf) # Estimated degrees of freedom for any smoothing terms in the model
  edf_part2 <- .edf_basisfunctionpart(object = object, ncores = ncores) # Estimated degrees of freedom for basis function component
  edf_part2 <- edf_part2[!apply(edf_part2, 1, function(x) all(is.na(x))),, drop = FALSE]
  final_edf <- rbind(edf_part1, edf_part2)
  return(final_edf)
  }


.edf_basisfunctionpart <- function(object, ncores = NULL) {
  if(!inherits(object, "CBFM")) 
    stop("`object' is not of class \"CBFM\"")
  if(!object$stderrors)
    stop("Standard errors must have been produced from `object' for influence measures to work.")
  
  if(is.null(ncores))
    registerDoParallel(cores = detectCores()-1)
  if(!is.null(ncores))
    registerDoParallel(cores = ncores)

  ##------------------
  ## Calculate building block quantities
  ##------------------
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
  
  weights_mat <- .neghessfamily(family = object$family, eta = object$linear_predictors, y = object$y, 
                                phi = matrix(object$dispparam, num_units, num_spp, byrow = TRUE), 
                                powerparam = matrix(object$powerparam, num_units, num_spp, byrow = TRUE),
                                zeroinfl_prob_intercept = matrix(object$zeroinfl_prob_intercept, num_units, num_spp, byrow = TRUE), 
                                trial_size = object$trial_size, domore = TRUE)
  if(object$family$family[1] %in% c("ztpoisson", "ztnegative.binomial"))
    weights_mat$out[is.na(weights_mat$out)] <- 0
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
  
  
  ##------------------
  ## Now the EDF
  ##------------------
  # Estimated degrees of freedom. Actually EDF is calculated for all coefficients, but here we only make available those for basis functions 
  # In fact, in additional simulations we found that the EDFs calculated for the terms in formula_X are typically very very close to those from edf and the last update using GAM in the PQL estimation algorithm. 
  edfs <- diag(bigV %*% crossprod(bigsqrtWXB))
  names(edfs) <- colnames(bigV)
  if(object$num_B == 0)
    edfs_out <- NULL
  if(object$num_B > 0) {
    edfs_out <- matrix(NA, nrow = num_spp, ncol = 3)
    rownames(edfs_out) <- colnames(object$y)
    colnames(edfs_out) <- c("B_space", "B_time", "B_spacetime")
    
    for(j in 1:num_spp) { 
      sub_edfs <- edfs[grep(paste0(colnames(object$y)[j],"$"), names(edfs))]
      sub_edfs <- sub_edfs[-(1:(nrow(object$covar_components$topleft)/num_spp))] 
      
      if(object$which_B_used[1]) {
        edfs_out[j,1] <- sum(sub_edfs[1:object$num_B_space])
        sub_edfs <- sub_edfs[-(1:object$num_B_space)] 
        }
      if(object$which_B_used[2]) {
        edfs_out[j,2] <- sum(sub_edfs[1:object$num_B_time])
        sub_edfs <- sub_edfs[-(1:object$num_B_time)] 
        }
      if(object$which_B_used[3]) {
        edfs_out[j,3] <- sum(sub_edfs[1:object$num_B_spacetime])
        sub_edfs <- sub_edfs[-(1:object$num_B_spacetime)] 
        }
      }
    }
  
  return(t(edfs_out))
  }


#' @method edf CBFM
#' @export edf 
edf <- function(object, ...) {
     UseMethod("edf")
     }
