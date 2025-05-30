#' @title Regression diagnostics from a CBFM fit
#' 
#' @description 
#' `r lifecycle::badge("experimental")`
#'
#' Takes a fitted \code{CBFM} object and calculates hat values i.e., diagonal elements of the influence/hat matrix and an approximate Cook's distance.
#'
#' @param object An object of class \code{CBFM}.
#' @param ncores To speed up calculation of the influence measures, parallelization can be performed, in which case this argument can be used to supply the number of cores to use in the parallelization. Defaults to \code{detectCores()-1}.
#' @param ... Not used.
#'
#' @details 
#' This is a relatively basic function that calculates the so-called "hat values" and an approximate Cook's distance. Classically, hat values are used to identify potentially high-leverage observations, and Cook's distance is a well-known regression-deletion diagnostic measuring the overall influence of a point on the fit of a model; we refer the reader to [stats::influence.measures()] for more information and reference pertaining to its use in standard linear and other regression models. Because both concepts of both can be (sort of) carried over to a generalized additive model or GAM e.g., see [mgcv::influence.gam()], then analogous diagnostics can be produced for a CBFM.
#' 
#' All measures are obtained on a *per-species basis*. This makes sense since the influence/leverage of a unit e.g., space-time coordinate, will differ for different species and their relationships to the measured covariates. 
#' 
#' ** We leave it up to the user to decide how to use the regression diagnostics, if at all.** In a GAM let alone a CBFM for spatio-temporal multivariate abundance data these diagnostics may only be approximate, and so providing rules-of-thumb related to their usage is challenging. Besides, we echo the sentiment provided in Chapter 4.4 of Fox et al., (2019) that cutoffs and rules-of-thumb should not be given too much weight, with more attention placed on graphical displays and assessing *relative influence* of observations.   
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
#' @seealso [CBFM()] for fitting CBFMs, [edf.CBFM()] for extracting estimated degrees of freedom associated with a CBFM fit, [residuals.CBFM()] for calculating different types of residuals from a CBFM fit, and [plot.CBFM()] for basic residual diagnostics.
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
  
     # Weights W
     num_spp <- nrow(object$betas)
     num_units <- nrow(object$B)
     if(is.null(object$dispparam))
          object$dispparam <- rep(1, num_spp)
     if(is.null(object$powerparam))
          object$powerparam <- rep(0, num_spp)
     if(object$family$family[1] %in% c("zipoisson","zinegative.binomial")) {                        
          zieta <- as.vector(tcrossprod(model.matrix.CBFM(object, zi = TRUE), object$zibetas))
          }
   
     weights_mat <- .neghessfamily(family = object$family, eta = object$linear_predictors, y = object$y, 
                                   phi = matrix(object$dispparam, num_units, num_spp, byrow = TRUE), 
                                   powerparam = matrix(object$powerparam, num_units, num_spp, byrow = TRUE),
                                   zieta = zieta, trial_size = object$trial_size, domore = TRUE)
     if(object$family$family[1] %in% c("ztpoisson"))
          weights_mat$out[is.na(weights_mat$out)] <- 0
     
     if(!(object$family$family[1] %in% c("zipoisson","zinegative.binomial"))) {
          weights_mat <- matrix(weights_mat$out, nrow = num_units, ncol = num_spp) # Overwrite weights_mat since only one quantity needed
          weights_mat[is.na(object$y)] <- 0
          }
     if(object$family$family[1] %in% c("zipoisson","zinegative.binomial")) {
          weights_mat_betabeta <- matrix(weights_mat$out, nrow = num_units, ncol = num_spp)
          weights_mat_betabeta[is.na(object$y)] <- 0
          weights_mat$out_zeroinflzeroinfl[is.na(object$y)] <- 0
          weights_mat$out_zeroinflbetas[is.na(object$y)] <- 0
          }
   
   
     # [W^{1/2}X, W^{1/2}B]
     X <- model.matrix.CBFM(object)
     if(!(object$family$family[1] %in% c("zipoisson","zinegative.binomial"))) {
          WsqrtXB <- function(j) {                
               return(list(WsqrtX = X*sqrt(weights_mat[,j]), WsqrtB = object$B * sqrt(weights_mat[,j])))
               }
          getall_WsqrtXB <- foreach(j = 1:num_spp) %dopar% WsqrtXB(j = j)
          bigsqrtWXB <- cbind(bdiag(lapply(getall_WsqrtXB, function(x) x$WsqrtX)), bdiag(lapply(getall_WsqrtXB, function(x) x$WsqrtB))) 
          }
   
     if(object$family$family[1] %in% c("zipoisson","zinegative.binomial")) {
          ziX <- model.matrix.CBFM(object, zi = TRUE)
          WsqrtXB <- function(j) {                
               #Xi <- bdiag(model.matrix.CBFM(object, zi = TRUE), X)
               #bigW <- cbind(Diagonal(x = weights_mat$out_zeroinflzeroinfl[,j]),  Diagonal(x = weights_mat$out_zeroinflbetas[,j]))
               #bigW <- rbind(bigW, cbind(Diagonal(x = weights_mat$out_zeroinflbetas[,j]),  Diagonal(x = weights_mat_betabeta[,j])))
               # XTWX <- crossprod(Xi, bigW) %*% Xi
               # BTWX <- crossprod(object$B, cbind(Diagonal(x = weights_mat$out_zeroinflbetas[,j]), Diagonal(x = weights_mat_betabeta[,j]))) %*% Xi
               # BTWB <- crossprod(object$B*sqrt(weights_mat_betabeta[,j]))
               # return(list(XTWX = XTWX, BTWX = BTWX, BTWB = BTWB))

               ## For ZI distributions, the weight matrix is not diagonal due to the non-orthgonality of the coefficients related to the count and zero-inflation components. For hat matrix calculations though, focus and use only the diagonal elements    
               return(list(WsqrtX = cbind(ziX*sqrt(weights_mat$out_zeroinflzeroinfl[,j]), X*sqrt(weights_mat_betabeta[,j])), 
                           WsqrtB = object$B*sqrt(weights_mat_betabeta[,j])))
               }
          getall_WsqrtXB <- foreach(j = 1:num_spp) %dopar% WsqrtXB(j = j)
          bigsqrtWXB <- cbind(bdiag(lapply(getall_WsqrtXB, function(x) x$WsqrtX)), bdiag(lapply(getall_WsqrtXB, function(x) x$WsqrtB))) 
          rm(ziX)
          }
     gc()
     rm(X, getall_WsqrtXB)
   
  
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
    hatvals <- foreach(j = 1:num_spp, .combine = cbind) %dopar% gethatvals_j(j = j) # Occassionally you might some get hat-values to be exactly zero due to weights doing funky things, but on the whole they should be rare
   
  
  # Experimentation. But most of the time frequentist either has smaller standard errors, so not used
  #freq_cov_stderrs <- sqrt(rowSums((bigV %*% crossprod(bigsqrtWXB)) * bigV))
  #freq_cov_stderrs <- matrix(freq_cov_stderrs[1:(num_spp*ncol(model.matrix(object)))], ncol = num_spp)
  #bayes_cov_stderrs <- matrix(sqrt(diag(object$covar_components$topleft)), ncol = num_spp)
  #matplot(bayes_cov_stderrs, freq_cov_stderrs)
  #abline(0,1)
  
    # Cook's distance -- The formula for this is taken from some hidden code in mgcv, and credit goes to Simon Wood for this [dont' know if this code is still hidden in later versions of mgcv]
    rm(bigV, weights_mat)    
    res <- residuals.CBFM(object, type = "pearson")
    get_B_edf <- edf.CBFM(object)
    final_dfs <- colSums(object$edf) + get_B_edf[((nrow(get_B_edf)-sum(object$which_B_used)+1):nrow(get_B_edf)),,drop=FALSE]
    if(object$family$family[1] %in% c("zipoisson","zinegative.binomial"))
         final_dfs <- final_dfs + colSums(object$ziedf)
     cookD <- (res / (1 - hatvals))^2 * hatvals / matrix(final_dfs, nrow = num_units, ncol = num_spp, byrow = TRUE)
     rm(get_B_edf, final_dfs)
   
     rownames(hatvals) <- rownames(cookD) <- rownames(object$y)
     colnames(hatvals) <- colnames(cookD) <- colnames(object$y)
     hatvals[is.na(object$y)] <- NA
     cookD[is.na(object$y)] <- NA
     if(object$family$family[1] %in% c("ztpoisson", "ztnegative.binomial")) {
          hatvals[which(object$y == 0)] <- NA
          cookD[which(object$y == 0)] <- NA
          }
  

   return(list(hat = hatvals, cooks = cookD))
   }


#' @method influence CBFM
#' @export influence
influence <- function(object, ...) {
     UseMethod("influence")
     }
