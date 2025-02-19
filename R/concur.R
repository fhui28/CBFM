#' @title Concurvity measures for a CBFM fit
#' 
#' @description 
#' `r lifecycle::badge("experimental")`
#' 
#' Calculates measures of concurvity from a fitted \code{CBFM} object. The measures are adapted from those found in [mgcv::concurvity()], and expanded to work on individual parametric terms as well on the spatial and/or temporal basis functions included in the CBFM.  
#' 
#' @param object An object of class \code{CBFM}.
#' @param zi For zero-inflated distributions, set this to \code{TRUE} if measures of concurvity associated with modeling the probabilities of zero-inflation are desired.
#' @param ... Not used.
#' 
#' @details 
#' As explained in [mgcv::concurvity()], the concurvity can be viewed as a generalization of multicollinearity sometimes/often seen in parametric regression models, where smooth terms in a model could be approximated by one or more of the other smooth terms in the model. Like collinearity, concurvity **may** lead to issues in terms of interpretation of response-covariate relationship, and potentially make estimates unstable with potentially inflated standard errors. It is important to emphasize the "**may**": depending on the scientific question/s of interest and the type/s of inference the user is interested in, concurvity like collinearity may or may not be an issue.    
#' 
#' Concurvity is something that is perhaps particularly worth keeping in mind in the context of CBFMs: many of the measured covariates included are likely to be spatially and/or temporally indexed, in which case their inclusion (whether as a smooth or parametric) maybe exhibit concurvity with the "smooths" resulting from the spatial and/or temporal basis functions also included to account for residual spatio-temporal correlations both between- and within-species. Again however, we stress the **perhaps**, as extreme concurvity may be something to think/worry about but otherwise it may not be something to stress over depending on the goal of the CBFM. Concurvity is also closely related (I suspect?!) to the issue of spatial confounding (Hodges et al., 2010, Hanks et al., 2015; Lany et al., 2020) 
#' 
#' For each species, this function presents two measures of concurvity, bounded between zero and one (zero indicating no problem and one indicating a total lack of identifiability in the CBFM), which are based on those found in [mgcv::concurvity()]. Both are based on something of the form \eqn{\|g\|^2 / \|f\|^2}, where \eqn{\|\cdot\|^2} is the squared norm, \eqn{f} is the smooth term of focus, and \eqn{g} is the "component" of that smooth which can be represented by other remaining terms in the model. The two measures specifically are 
#' 1. \code{Observed}, which is based on the corresponding linear predictors and their Euclidean norm. Basically, it is a measure of how well the linear predictor for a specific covariate can be represented by a linear predictor formed from the other terms in the model. According to [mgcv::concurvity()], the measure can be a bit over-optimistic about the potential for a problem in some cases. 
#' 2. \code{Estimate}, which is based on the the covariate only and its representation in terms of the other terms in the model. Basically, it is a measure of how well a specific covariate can be represented by the other terms in the model (in terms of breaking things down to basis functions). According to [mgcv::concurvity()], it does not suffer from the pessimism or potential for over-optimism than the first measure.
#' 
#' Note concurvity to measured on a per-species basis, although this does not matter for the second measure due to the covariates community-level nature of the basis functions.
#' 
#' Finally, note this function only provides an indicator of potential collinearity/concurvity; it does not provide a remedy to the problem if it arises. 
#' 
#' @return A data frame providing the species-specific values of the two concurvity measures for each term in the CBFM.
#' 
#' @details # Warning
#' 
#' This function currently does not work well with factor variables. However it should be kept in mind that concurvity and collinearity with factor variables is generally more challenging to diagnose and deal with. 
#' 
#' @author Francis K.C. Hui <fhui28@gmail.com>, Chris Haak
#' 
#' @references
#' Hanks, E. M., Schliep, E. M., Hooten, M. B., and Hoeting, J. A. (2015). Restricted spatial regression in practice: geostatistical models, confounding, and robustness under model misspecification. Environmetrics, 26, 243-254.
#' 
#' Hodges, J. S., and Reich, B. J. (2010). Adding spatially-correlated errors can mess up the fixed effect you love. The American Statistician, 64, 325-334.
#' 
#' Lany, N. K., Zarnetske, P. L., Finley, A. O., and McCullough, D. G. (2020). Complementary strengths of spatially‐explicit and multi‐species distribution models. Ecography, 43, 456-466.
#' 
#' #' @seealso [CBFM()] for fitting CBFMs, [plot.CBFM()] for basic residual diagnostics from a CBFM fit, [residuals.CBFM()] for calculating various types of residuals.
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
#' concur(fitcbfm)
#' 
#' 
#' ##------------------------------
#' ## **Example 1b: Repeat Example 1a but illustrate the use of smoothing terms in CBFM**
#' ## Since the true model only involves parametric terms, then we do not expect its performance
#' ## to be as good as assuming the right form for the mean model.
#' ## It is purely for illustration purposes. 
#' ## Please note this will take a while to run...get a cuppa and stretch your legs! 
#' ##------------------------------
#' # Set up spatial basis functions for CBFM -- Most users will start here! 
#' num_basisfunctions <- 25 # Number of spatial basis functions to use
#' basisfunctions <- mrts(dat[,c("x","y")], num_basisfunctions) %>% 
#' as.matrix %>%
#' {.[,-(1)]} # Remove the first intercept column
#' 
#' # Fit CBFM 
#' tic <- proc.time()
#' useformula <- ~ s(temp) + s(depth) + s(chla) + s(O2)
#' fitcbfm_gam <- CBFM(y = simy, formula = useformula, 
#' data = dat, B_space = basisfunctions, family = binomial(), control = list(trace = 1))
#' toc <- proc.time()
#' toc-tic
#' 
#' concur(fitcbfm_gam)
#' }
#' 
#' @export
#' @importFrom mgcv gam 
#' @md

concur <- function(object, zi = FALSE, ...) {
     if(!inherits(object, "CBFM")) 
        stop("`object' is not of class \"CBFM\"")
     
     num_spp <- ncol(object$y)
     tmp_formula <- as.formula(paste("response", paste(as.character(object$formula), collapse = " ") ) )
     nullfit <- gam(tmp_formula, data = data.frame(response = runif(nrow(object$y)), object$data), knots = object$knots, fit = TRUE, control = list(maxit = 1))
     
     if(zi) {
          if(!(object$family$family[1] %in% c("zipoisson","zinegative.binomial"))) 
               stop("Measures of concurvity associated with the probability of zero-inflation can be only be obtained for zero-inflated CBFMs.")
          if(object$family$family[1] %in% c("zipoisson","zinegative.binomial")) {
               tmp_formula <- as.formula(paste("response", paste(as.character(object$ziformula), collapse = " ") ) )
               nullfit <- gam(tmp_formula, data = data.frame(response = runif(nrow(object$y)), object$data), knots = object$ziknots, fit = TRUE, control = list(maxit = 1))
               }
          }
          
     
     # Start with parametric terms, but exclude intercept. Then include smoothing terms. Then basis functions.
     if(length(nullfit$smooth) == 0) {
          if(length(nullfit$coefficients) == 1)  
               start <- stop <- labs <- NULL
          if(length(nullfit$coefficients) > 1) {
               start <- stop <- 2:length(nullfit$coefficients)
               labs <- names(nullfit$coefficients)[-1]
               }
          }
     if(length(nullfit$smooth) > 0) {
          find_last_parametric_term <- nullfit$smooth[[1]]$first.para-1
          if(find_last_parametric_term == 1) 
                start <- stop <- labs <- NULL
          if(find_last_parametric_term > 1) {
                start <- stop <- 2:find_last_parametric_term
                labs <- names(nullfit$coefficients)[2:find_last_parametric_term]
                }
          for(i in 1:length(nullfit$smooth)) { 
               start <- c(start, nullfit$smooth[[i]]$first.para)
               stop <- c(stop, nullfit$smooth[[i]]$last.para)
               labs <- c(labs, nullfit$smooth[[i]]$label)
               }
          }
     
     num_X_terms <- length(start)
     
     if(!zi) {
          if(object$num_B_space > 0) {
               start <- c(start, stop[length(stop)]+1)
               stop <- c(stop, stop[length(stop)] + object$num_B_space)
               labs <- c(labs, "B_space")
               }
          if(object$num_B_time > 0) {
               start <- c(start, stop[length(stop)]+1)
               stop <- c(stop, stop[length(stop)] + object$num_B_time)
               labs <- c(labs, "B_time")
               }
          if(object$num_B_spacetime > 0) {
               start <- c(start, stop[length(stop)]+1)
               stop <- c(stop, stop[length(stop)] + object$num_B_spacetime)
               labs <- c(labs, "B_spacetime")
               }
          }
     
     num_terms <- length(start)
     full_concur <- array(0, dim = c(num_spp, 2, num_terms), dimnames = list(species = colnames(object$y), measure = c("Observed", "Estimate"), covariates = labs))
     rm(labs, tmp_formula, nullfit)
     # So start and stop capture all things implied by formula. num_terms counts this along with the any basis coefficients in model

     
     X <- cbind(model.matrix.CBFM(object), object$B)
     if(zi) {
          X <- model.matrix.CBFM(object, zi = TRUE)
          } 
     X <- base::qr.R(base::qr(X)) # Using a QR decomposition speeds up remaining computation for calculating measures?
     
     # Doing concurvity measures for terms in formula
     for(k0 in 1:(num_terms-sum(object$which_B_used))) {
          Xj <- X[, start[k0]:stop[k0], drop = FALSE]
          Xi <- X[, -(start[k0]:stop[k0]), drop = FALSE]
          sel_num_cols <- ncol(Xi) 
          R <- base::qr.R(base::qr(cbind(Xi,Xj), tol = 0))[, -(1:sel_num_cols), drop = FALSE] 
          Rt <- base::qr.R(base::qr(R, tol = 0)) 
          
          for(k1 in 1:num_spp) { 
               # Observed 
               # Consider a QR decompostion [everything else, covariate] = Q R = Q [R_everythingelse, R_covariate]. Thus covariate = Q R_covariate = Q_everythingelse*R_everythingelse + Q_covariate*R covariate. The first part lies entirely in the space of one or more other covariates in the model. The remainder part that is completely within the covariate's own space. 
               # Numerator is ||Q_everything*R_everythingelse %*% beta ||^2 
               # Denominator is ||X%*%beta||^2 for that covariate
               sel_betas <- object$betas[k1, start[k0]:stop[k0]]
               if(zi)
                    sel_betas <- object$zibetas[k1, start[k0]:stop[k0]]
               full_concur[k1,1,k0] <- sum((R[1:sel_num_cols,,drop = FALSE] %*% sel_betas)^2) / sum((Rt %*% sel_betas)^2)
               sel_betas <- object$betas[k1, start[k0]:stop[k0]]
               
               # Less pessimistic estimate
               # Numerator is ||Q_everything*R_everythingelse||^2 
               # Denominator is ||X||^2 for that covariate
               full_concur[k1,2,k0] <- sum(R[1:sel_num_cols,]^2) / sum(R^2)
               }
          }
     
     
     # Doing concurvity measures for basis function terms
     if(!zi) {
          for(k0 in (num_X_terms+1):num_terms) {
               Xj <- X[, start[k0]:stop[k0], drop = FALSE]
               Xi <- X[, -(start[k0]:stop[k0]), drop = FALSE]
               sel_num_cols <- ncol(Xi) 
               
               R <- base::qr.R(base::qr(cbind(Xi,Xj), tol = 0))[, -(1:sel_num_cols), drop = FALSE] 
               Rt <- base::qr.R(base::qr(R, tol = 0)) 
               
               for(k1 in 1:num_spp) { 
                    # Observed 
                    sel_betas <- cbind(object$betas, object$basis_effects_mat)[k1, start[k0]:stop[k0]]
                    full_concur[k1,1,k0] <- sum((R[1:sel_num_cols,,drop = FALSE] %*% sel_betas)^2) / sum((Rt %*% sel_betas)^2)
     
                    # Less pessimistic estimate
                    full_concur[k1,2,k0] <- sum(R[1:sel_num_cols,]^2) / sum(R^2)
                    }
               }
          }
     
     full_concur <- as.data.frame.table(full_concur)
     colnames(full_concur) <- c("Species", "Measure", "Term", "Value")
     return(full_concur)
     }

     
# @method concur CBFM
# @export concur 
# concur <- function(object, ...) {
#      UseMethod("concur")
#      }
