#' @title Construct predictions for a CBFM fit
#' 
#' @description 
#' `r lifecycle::badge("experimental")`
#' 
#' Takes a fitted \code{CBFM} object and produces predictions given (potentially) a new set of observational units with their corresponding covariate and basis function functions. Predictions can be accompanied by standard errors, based on the Bayesian estimated covariance matrix of the parameter estimates.
#' 
#' @param object An object of class "CBFM".
#' @param newdata A data frame containing the values of the covariates at which predictions are to be calculated. If this is not provided, then predictions corresponding to the original data are returned. If \code{newdata} is provided then it should contain all the variables needed for prediction, that is, it can construct a model matrix from this as \code{object$formula_X}.
#' @param manualX A manually supplied model matrix at which predictions are to be calculated. This can be used if for some reason the user wants to supply a very custom model matrix for calculating predictions. Note supply of this overrides any supplied \code{newdata} argument. The number of columns in \code{manualX} should equal to \code{ncol(object$betas)}.
#' @param new_B_space A matrix of new spatial basis functions at which predictions are to be calculated. If this is not provided, then predictions corresponding to the original \code{B_space} argument are returned. Please note this should only be supplied if \code{B_space} was supplied in the original CBFM fit.  
#' @param new_B_time A matrix of new temporal basis functions at which predictions are to be calculated. If this is not provided, then predictions corresponding to the original \code{B_time} argument are returned. Please note this should only be supplied if \code{B_time} was supplied in the original CBFM fit.  
#' @param new_B_spacetime A matrix of new spatio-temporal basis functions at which predictions are to be calculated. If this is not provided, then predictions corresponding to the original \code{B_spacetime} argument are returned. Please note this should only be supplied if \code{B_spacetime} was supplied in the original CBFM fit.  
#' @param type The type of prediction required. The default \code{type = "link"} is on the scale of the linear predictors. Alternatively, \code{type = "response"} returns predictions on the scale of the response variable. For example, the predictions for a binomial CBFM are the predicted probabilities.
#' @param se_fit When this is set to \code{TRUE} (not default), then standard error estimates are returned for each predicted value.
#' @param coverage The coverage probability of the uncertainty intervals for prediction. Defaults to 0.95, which corresponds to 95% uncertainty intervals.
#' @param ncores To speed up calculation of the standard error estimates, parallelization can be performed, in which case this argument can be used to supply the number of cores to use in the parallelization. Defaults to \code{detectCores()-1}.
#' @param ... Not used.
#' 
#' @details 
#' The standard errors produced by \code{predict.CBFM} are based on the Bayesian posterior covariance matrix of the estimated parameters from the fitted \code{CBFM} object, and associated uncertainty intervals are obtained on the associated large sample normality result of the estimated parameters. Both of these are similar to [mgcv::predict.gam()].
#' 
#' @return If \code{se_fit = TRUE}, then a list with the following components (if applicable) is returned:
#' \item{fit: }{A matrix of predicted values.}
#' \item{stderr: }{A matrix of estimated standard errors associated with the predictions.}
#' \item{lower: }{A matrix of the lower bound of the uncertainty intervals for the predictions.}
#' \item{upper: }{A matrix of the upper bound of the uncertainty intervals for the predictions.}
#' Otherwise if \code{se_fit = FALSE}, then a matrix of predicted values is returned. 
#' 
#' @author Francis K.C. Hui <fhui28@gmail.com>, Chris Haak
#' 
#' @seealso [CBFM()] for fitting CBFMs, [fitted.values.CBFM()] for obtaining fitted values from a CBFM fit, and [residuals.CBFM()] for calculating various types of residuals.
#' 
#' @examples
#' \donttest{
#' # Please see examples in the help file for the main CBFM function 
#' }
#' 
#' @export
#' 
#' @import foreach 
#' @import Matrix
#' @importFrom doParallel registerDoParallel
#' @importFrom mgcv gam model.matrix.gam predict.gam
#' @importFrom parallel detectCores
#' @importFrom stats qnorm
#' @md

predict.CBFM <- function(object, newdata = NULL, manualX = NULL, new_B_space = NULL, new_B_time = NULL, new_B_spacetime = NULL, type = c("link", "response"), se_fit = FALSE, coverage = 0.95, ncores = NULL, ...) {
     if(is.null(ncores))
          registerDoParallel(cores = detectCores()-1)
     if(!is.null(ncores))
          registerDoParallel(cores = ncores)
     
     
     if(!is.null(manualX)) {
          new_X <- manualX
          warning("manualX has been supplied. This overrides the creation of a model matrix based on object$formula_X and/or newdata.")
          }
     if(is.null(manualX)) {
          tmp_formula <- as.formula(paste("response", paste(as.character(object$formula_X),collapse="") ) )
          nullfit <- gam(tmp_formula, data = data.frame(response = rnorm(nrow(object$data)), object$data), fit = TRUE, control = list(maxit = 1))
          if(is.null(newdata))
               new_X <- predict.gam(nullfit, type = "lpmatrix")
          if(!is.null(newdata))
               new_X <- predict.gam(nullfit, newdata = data.frame(newdata), type = "lpmatrix")
          rm(tmp_formula, nullfit)
          }
     if(ncol(new_X) != ncol(object$betas))
          stop("Number of columns in new_X should match the number of columns in new_X.")

     newB <- NULL
     if(!is.null(new_B_space)) {
          if(object$num_B_space != ncol(new_B_space))
               stop("The number of columns of new_B_space does not object$num_B_space.")
          newB <- cbind(newB, new_B_space)
          }
     if(!is.null(new_B_time)) {
          if(object$num_B_time != ncol(new_B_time))
               stop("The number of columns of new_B_time does not object$num_B_time.")
          newB <- cbind(newB, new_B_time)
          }
     if(!is.null(new_B_spacetime)) {
          if(object$num_B_spacetime != ncol(new_B_spacetime))
               stop("The number of columns of new_B_spacetime does not object$num_B_spacetime.")
          newB <- cbind(newB, new_B_spacetime)
          }
     if(is.null(newB))
          newB <- object$B
     if(!is.null(newB)) {
          if(nrow(newB) != nrow(new_X))
               stop("newB does not contain the same number of rows as the model matrix created based on object and/or newdata.")
          if(ncol(newB) != ncol(object$basis_effects_mat))
               stop("newB does not contain the same number of columns as the number of columns in object$basis_effects_mat.")
          }
          
     type <- match.arg(type, choices = c("link", "response"))
     num_spp <- nrow(object$betas)
     num_X <- ncol(new_X)
     num_basisfns <- ncol(newB)
     
     
     if(se_fit == TRUE & object$stderrors == FALSE)
          stop("Standard errors can not be calculated since the covariance matrix estimate was not detected to be available in object.")

     ptpred <- tcrossprod(new_X, object$betas) + tcrossprod(newB, object$basis_effects_mat)
     ptpred <- as.matrix(ptpred)
     colnames(ptpred) <- rownames(object$betas)          

     if(!se_fit) {
          if(type == "response") {
               if(object$family$family != "ztnegative.binomial")
                    ptpred <- object$family$linkinv(ptpred)
               # if(object$family$family == "ztnegative.binomial")
               #      ptpred <- object$family$linkinv(eta = ptpred, phi = matrix(object$dispparam, nrow = nrow(y), ncol = ncol(y), byrow = TRUE))
               }
          return(ptpred)
          }
     
     if(se_fit) {
          ci_alpha <- qnorm((1-coverage)/2, lower.tail = FALSE)
          
          stderr1_fun <- function(j) {
               out1 <- object$covar_components$topleft[(num_X*j - num_X + 1):(num_X*j), (num_X*j - num_X + 1):(num_X*j),drop=FALSE]
               return(diag(new_X %*% tcrossprod(out1, new_X)))
               }
          stderr2_fun <- function(j) {
               out1 <- object$covar_components$topright[[j]]
               return(diag(new_X %*% tcrossprod(out1, newB)))
               }
          stderr3_fun <- function(j) {
               out1 <- object$covar_components$bottomright[[j]]
               return(diag(newB %*% tcrossprod(out1, newB)))
               }

          stderr <- foreach(j = 1:num_spp, .combine = "cbind") %dopar% stderr1_fun(j = j)
          stderr <- stderr + 2 * foreach(j = 1:num_spp, .combine = "cbind") %dopar% stderr2_fun(j = j)
          stderr <- stderr + foreach(j = 1:num_spp, .combine = "cbind") %dopar% stderr3_fun(j = j)          
          stderr <- as.matrix(stderr)
          colnames(stderr) <- rownames(object$betas)    
          
          alllower <- ptpred - ci_alpha * sqrt(stderr) 
          allupper <- ptpred + ci_alpha * sqrt(stderr)
          if(type == "response") {
               if(object$family$family != "ztnegative.binomial") {
                    ptpred <- object$family$linkinv(ptpred)
                    alllower <- object$family$linkinv(alllower)
                    allupper <- object$family$linkinv(allupper)
                    }
               # if(object$family$family == "ztnegative.binomial") {
               #      ptpred <- object$family$linkinv(eta = ptpred, phi = matrix(object$dispparam, nrow = nrow(y), ncol = ncol(y), byrow = TRUE))
               #      alllower <- object$family$linkinv(eta = alllower, phi = matrix(object$dispparam, nrow = nrow(y), ncol = ncol(y), byrow = TRUE))
               #      allupper <- object$family$linkinv(eta = allupper, phi = matrix(object$dispparam, nrow = nrow(y), ncol = ncol(y), byrow = TRUE))
               #      }
               }
          
          return(list(fit = ptpred, stderr = sqrt(stderr), lower = alllower, upper = allupper))
          }
     }
     
