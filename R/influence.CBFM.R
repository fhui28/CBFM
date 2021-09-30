#' @title Regression diagnostics from a CBFM fit
#' 
#' @description 
#' `r lifecycle::badge("experimental")`
#'
#' Takes a fitted \code{CBFM} object and calculates hat values i.e., diagonal elements of the influence/hat matrix, an approximate Cook's distance, and estimated or effective degrees of freedom for the species-specific regression coefficients corresponding to spatial and/or temporal basis functions included in the model.
#'
#' @param object An object of class \code{CBFM}.
#' @param ncores To speed up calculation of the influence measures, parallelization can be performed, in which case this argument can be used to supply the number of cores to use in the parallelization. Defaults to \code{detectCores()-1}.
#' @param ... Not used.
#'
#' @details 
#' This is a relatively basic function that calculates the so-called "hat values", an approximate Cook's distance, and some estimated degrees of freedom from a fitted CBFM. Classically, hat values are used to identify potentially high-leverage observations, and Cook's distance is a well-known regression-deletion diagnostic that measures the overall influence of a point on the fit of a model; we refer the reader to [stats::influence.measures()] for more information and reference pertaining to its use in standard linear and other regression models. Because both concepts of both can be (sort of) carried over to a generalized additive model or GAM e.g., see [mgcv::influence.gam()], then analogous diagnostics can be produced for a CBFM. 
#' 
#' Regarding the estimated degrees of freedom (EDF) for the basis function coefficients, for each species up to three EDFs are given depending on which of \code{B_space/B_time/B_spacetime} are included in the model. These degrees of freedom values are analogous to what are available in \code{object$edf/object$edf1}, which are the EDFs for each model parameter in \code{formula_X}; see [CBFM()] for more information. Note however that because of the way the CBFM is set up, there is usually a considerable amount of penalization taking place for regression coefficients corresponding to the spatial and/or temporal basis functions, and so one should expect these value to usually be *much* smaller than the corresponding number of basis functions included in the model. On their own, the EDFs are not really of much use at the moment, especially since they are at the moment **not** used in calculating overall degrees of freedom e.g., [logLik.CBFM()] or in information criteria e.g., [AIC.CBFM()] and [AICc.CBFM()]. This might change at some point down the road though... 
#' 
#' All measures are obtained on a *per-species basis*. This makes sense since the influence/leverage of a unit e.g., space-time coordinate, will differ for different species and their relationships to the measured covariates. 
#' 
#' ** We leave it up to the practitioner to decide how to use the regression diagnostics, if at all.** In a GAM let alone a CBFM for spatio-temporal multivariate abundance data these diagnostics may only be approximate, and so providing rules-of-thumb related to their usage is challenging. Besides, we echo the sentiment provided in Chapter 4.4 of Fox et al., (2019) that cutoffs and rules-of-thumb should not be given too much weight, with more attention placed on graphical displays and assessing *relative influence* of observations.   
#'  
#'  
#' @return A list containing two elements:
#' \describe{
#' \item{hat: }{A matrix of estimated hat values i.e., diagonal elements of the influence/hat matrix. The dimensions of this matrix should be the same as \code{object$y}.}

#' \item{cooks: }{A matrix of estimated and approximate Cook's distances. The dimensions of this matrix should be the same as \code{object$y}. }

#' \item{B_edf: }{A matrix of estimated degrees of freedom corresponding to the spatial and/or temporal basis functions included in the model. The number of columns of this matrix should be the same as \code{object$y}, while it always has three rows. Any elements in the matrix equal to \code{NA} imply that that set "type" of basis function was not included in the model.}
#' }
#'
#' @author Francis K.C. Hui <fhui28@gmail.com>, Chris Haak
#'
#' @references 
#'Fox, J. (2019). Regression diagnostics: An introduction. Sage publications.
#'   
#' @seealso [CBFM()] for fitting CBFMs, [edf.CBFM()] for extracting estimated degrees of freedom associated with a CBFM fit, [residuals.CBFM()] for calculating different types of residuals from a CBFM fit, and [plot.CBFM()] for basic residual diagnostics.
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
#' influence(fitcbfm)
#' }
#'
#' @aliases influence influence.CBFM
#' @export 
#' @export influence.CBFM
#' @importFrom foreach foreach %dopar%
#' @import Matrix 
#' @importFrom doParallel registerDoParallel
#' @md

influence.CBFM <- function(object, ncores = NULL, ...) {
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
  
   weights_mat <- .neghessfamily(family = object$family, eta = object$linear_predictor, y = object$y, 
                                 phi = matrix(object$dispparam, num_units, num_spp, byrow = TRUE), 
                                 powerparam = matrix(object$powerparam, num_units, num_spp, byrow = TRUE),
                                 zeroinfl_prob_intercept = matrix(object$zeroinfl_prob_intercept, num_units, num_spp, byrow = TRUE), 
                                 trial_size = object$trial_size, domore = TRUE)
   if(object$family$family[1] %in% c("ztpoisson"))
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
   ## Now the measures
   ##------------------
   # hat values 
   gethatvals_j <- function(j) {
      sel_units <- (j*num_units - num_units + 1):(j*num_units)
      hatvals_j <- rowSums((bigsqrtWXB[sel_units,] %*% bigV) * bigsqrtWXB[sel_units,]) 
    
      return(hatvals_j)
      }
    #hatvals <- rowSums((bigsqrtWXB %*% bigV) * bigsqrtWXB)
   #hatvals <- matrix(hatvals, nrow = num_units, ncol = num_spp)
   hatvals <- foreach(j = 1:num_spp, .combine = cbind) %dopar% gethatvals_j(j = j)
   if(object$family$family[1] %in% c("ztpoisson"))
      hatvals[which(object$y == 0)] <- NA
   
  
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
  
  # Experimentation. But most of the time frequentist either has smaller standard errors, so not used
  #freq_cov_stderrs <- sqrt(rowSums((bigV %*% crossprod(bigsqrtWXB)) * bigV))
  #freq_cov_stderrs <- matrix(freq_cov_stderrs[1:(num_spp*ncol(model.matrix(object)))], ncol = num_spp)
  #bayes_cov_stderrs <- matrix(sqrt(diag(object$covar_components$topleft)), ncol = num_spp)
  #matplot(bayes_cov_stderrs, freq_cov_stderrs)
  #abline(0,1)
  
   # Cook's distance
   rm(bigsqrtWXB, bigV, weights_mat, edfs)    
   res <- residuals.CBFM(object, type = "pearson")
   final_dfs <- colSums(object$edf) + nrow(object$basis_effects_mat)
   cookD <- (res / (1 - hatvals))^2 * hatvals / matrix(final_dfs, nrow = num_units, ncol = num_spp, byrow = TRUE)
  
   rownames(hatvals) <- rownames(cookD) <- rownames(object$y)
   colnames(hatvals) <- colnames(cookD) <- colnames(object$y)
  

   return(list(hat = hatvals, cooks = cookD, B_edf = t(edfs_out)))
   }


#' @method influence CBFM
#' @export influence
influence <- function(object, ...) {
     UseMethod("influence")
     }
