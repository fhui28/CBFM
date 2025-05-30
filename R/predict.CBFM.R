#' @title Construct predictions for a (hurdle) CBFM fit
#' 
#' @description
#' `r lifecycle::badge("experimental")`
#' 
#' Takes a fitted \code{CBFM} or \code{CBFM_hurdle} object and produces predictions given (potentially) a new set of observational units with their corresponding covariate and basis function functions. Predictions can be accompanied by standard errors, based on the Bayesian covariance matrix of the parameter estimates. As another option, the function can return the model matrix of the covariates constructed (potentially) using at the new set of observational units; in [mgcv::predict.gam()] this is also known as the linear predictor matrix.  
#' 
#' @param object An object of class \code{CBFM} or \code{CBFM_hurdle}.
#' @param newdata A data frame containing the values of the covariates at which predictions are to be calculated. If this is not provided, then predictions corresponding to the original data are returned. If \code{newdata} is provided then it should contain all the variables needed for prediction, that is, it can construct a model matrix from this as \code{object$formula}.
#' @param manualX A manually supplied model matrix at which predictions are to be calculated. This can be used if for some reason the user wants to supply a very custom model matrix for calculating predictions. Note supply of this overrides any supplied \code{newdata} argument. The number of columns in \code{manualX} should equal to \code{ncol(object$betas)}.
#' @param manualziX For zero-inflated distributions, a manually supplied model matrix associated with modeling the probability of zero-inflation, at which predictions are to be calculated. This can be used if for some reason the user wants to supply a very custom model matrix for calculating predictions. Note supply of this overrides any supplied \code{newdata} argument. The number of columns in \code{manualziX} should equal to \code{ncol(object$zibetas)}.
#' @param new_B_space A matrix of new spatial basis functions at which predictions are to be calculated. If this is not provided, then predictions corresponding to the original \code{B_space} argument are returned. Please note this should only be supplied if \code{B_space} was supplied in the original CBFM fit.  
#' @param new_B_time A matrix of new temporal basis functions at which predictions are to be calculated. If this is not provided, then predictions corresponding to the original \code{B_time} argument are returned. Please note this should only be supplied if \code{B_time} was supplied in the original CBFM fit.  
#' @param new_B_spacetime A matrix of new spatio-temporal basis functions at which predictions are to be calculated. If this is not provided, then predictions corresponding to the original \code{B_spacetime} argument are returned. Please note this should only be supplied if \code{B_spacetime} was supplied in the original CBFM fit.  
#' @param newdata_pa For hurdle CBFM models, similar to \code{newdata} above but applied to the presence-absence component model.
#' @param manualX_pa For hurdle CBFM models, similar to \code{manual_X} above but applied to the presence-absence component model. 
#' @param new_B_space_pa For hurdle CBFM models, similar to \code{new_B_space} above but applied to the presence-absence component model. 
#' @param new_B_time_pa For hurdle CBFM models, similar to \code{new_B_time} above but applied to the presence-absence component model. 
#' @param new_B_spacetime_pa For hurdle CBFM models, similar to \code{new_B_spacetime} above but applied to the presence-absence component model.
#' @param newdata_count For hurdle CBFM models, similar to \code{newdata} above but applied to the zero-truncated count component model.
#' @param manualX_count For hurdle CBFM models, similar to \code{manual_X} above but applied to the zero-truncated count component model. 
#' @param new_B_space_count For hurdle CBFM models, similar to \code{new_B_space} above but applied to the zero-truncated count component model. 
#' @param new_B_time_count For hurdle CBFM models, similar to \code{new_B_time} above but applied to the zero-truncated count component model. 
#' @param new_B_spacetime_count For hurdle CBFM models, similar to \code{new_B_spacetime} above but applied to the zero-truncated count component model.
#' @param type The type of prediction required: 
#' The default \code{type = "link"} is on the scale of the linear predictors. Alternatively, \code{type = "response"} returns predictions on the scale of the response variable. For example, the predicted response for a binomial CBFM are the predicted probabilities.
#' Note for zero-inflated distributions, \code{type = "link"} returns the predicted linear predictor of the count component, consistent with what \code{object$linear_predictors} returns. But \code{type = "response"} returns the *actual predicted mean values* of the distribution, consistent with what \code{object$fitted} returns. If the linear predictor associated with modeling the probability of zero-inflation is desired, then use \code{type = "zilink"}.
#' Similarly, for zero-truncated distributions, \code{type = "link"} returns the predicted linear predictor of the base count distribution, consistent with what \code{object$linear_predictors} returns. But \code{type = "response"} returns the *actual predicted mean values* of the zero-truncated distribution, consistent with what \code{object$fitted} returns. 
#' Another option is given by \code{type  = "lpmatrix"}, which returns the model matrix of the covariates constructing (potentially) using \code{newdata} and based on \code{object$formula}. In [mgcv::predict.gam()], this is also known as the linear predictor matrix. Note under this option, arguments such as \code{new_B_space/new_B_time/new_B_spacetime/se_fit/coverage} are irrelevant and hence ignored. 
#' Another is given by \code{type  = "zilpmatrix"}, which returns the model matrix of the covariates constructing (potentially) using \code{newdata} and based on \code{object$ziformula} i.e., the linear predictor matrix associated with modeling the probability of zero inflation.
#' Finally, note for hurdle CBFMS, this argument is not available, as the only option which is really of use here (and which you can not obtain based on applying the \code{predict} to the components of the hurdle CBFM itself) is the predicted mean values i.e., \code{type = "response"}.
#' @param se_fit When this is set to \code{TRUE} (not default), then standard error estimates are returned for each predicted value.
#' @param coverage The coverage probability of the uncertainty intervals for prediction. Defaults to 0.95, which corresponds to 95% uncertainty intervals.
#' @param ncores To speed up calculation of the standard error estimates, parallelization can be performed, in which case this argument can be used to supply the number of cores to use in the parallelization. Defaults to \code{detectCores()-1}.
#' @param num_sims If simulation is required for constructing uncertainty intervals, then this specifies the number of Monte-Carlo examples to simulate.
#' @param return_internals Not used. Please ignore this argument!
#' @param ... Not used.
#' 
#' @details 
#' The standard errors produced by \code{predict.CBFM} are based on the Bayesian posterior covariance matrix of the estimated parameters from the fitted \code{CBFM} object, and associated uncertainty intervals are obtained based on the associated large sample normality result of i.e., basically a Gaussian approximation to the posterior distribution of the parameters. This construction is similar to [mgcv::predict.gam()].
#' 
#' The functions tries to avoid using simulation (also known as single-fit bootstrapping, Fletcher and Jowett, 2022) for constructing the uncertainty intervals for the predictions (to minimize computational burden) However, in some cases it will fallback to doing so when it can not find a simple way to construct these intervals  (blame Francis!) e.g., for the zero-inflated Poisson and negative-binomial distribution when \code{type = "response"}. 
#' 
#' 
#' @return If \code{se_fit = TRUE}, then a list with the following components (if applicable) is returned:
#' \item{fit: }{A matrix of predicted values.}

#' \item{stderr: }{A matrix of standard errors associated with the uncertainty intervals for the predictions.}

#' \item{lower: }{A matrix of the lower bound of the uncertainty intervals for the predictions.}

#' \item{upper: }{A matrix of the upper bound of the uncertainty intervals for the predictions.}
#' Otherwise if \code{se_fit = FALSE}, then a matrix of predicted values is returned. 
#' 
#' Other if \code{type = "lpmatrix"}, then a model matrix with the number of rows equal to either the number of rows if \code{object$y} (if \code{newdata} is not supplied) or the number of rows in \code{newdata} when it is supplied.
#' 
#' 
#' @author Francis K.C. Hui <fhui28@gmail.com>, Chris Haak
#' 
#' @references
#' Fletcher, D., and Jowett, T. (2022). Single-fit bootstrapping: A simple alternative to the delta method. Methods in Ecology and Evolution, 13, 1358-1367.
#' 
#' @seealso [CBFM()] for fitting CBFMs, [fitted.CBFM()] for obtaining fitted values from a CBFM fit, and [residuals.CBFM()] for calculating various types of residuals.
#' 
#' @examples
#' \dontrun{
#' # Please see examples in the help file for the CBFM and makeahurdle functions 
#' }
#' 
#' @export
#' 
#' @importFrom foreach foreach %dopar%
#' @import Matrix
#' @importFrom abind abind
#' @importFrom doParallel registerDoParallel
#' @importFrom mgcv gam model.matrix.gam predict.gam
#' @importFrom parallel detectCores
#' @importFrom stats qnorm plogis
#' @md

predict.CBFM <- function(object, newdata = NULL, manualX = NULL, manualziX = NULL, new_B_space = NULL, new_B_time = NULL, new_B_spacetime = NULL, 
                         type = "link", se_fit = FALSE, coverage = 0.95, ncores = NULL, num_sims = 500, return_internals = FALSE, ...) {
     if(!inherits(object, "CBFM")) 
          stop("`object' is not of class \"CBFM\"")

     if(is.null(ncores))
          registerDoParallel(cores = detectCores()-1)
     if(!is.null(ncores))
          registerDoParallel(cores = ncores)
        
     ## Construct X and ziX, and offsets
     if(!is.null(manualX)) {
          new_X <- manualX
          warning("manualX has been supplied. This overrides the creation of a model matrix based on object$formula and/or newdata.")
          }
     if(is.null(manualX)) {
          tmp_formula <- as.formula(paste("response", paste(as.character(object$formula),collapse = " ") ) )
          nullfit <- gam(tmp_formula, data = data.frame(response = runif(nrow(object$y)), object$data), knots = object$knots, fit = TRUE, control = list(maxit = 1))
          if(is.null(newdata)) {
               new_X <- predict.gam(nullfit, type = "lpmatrix")
               }
          if(!is.null(newdata)) {
               new_X <- predict.gam(nullfit, newdata = data.frame(newdata), type = "lpmatrix")
               }
          
          formula_offset <- attributes(new_X)$model.offset
          if(all(formula_offset == 0))
               formula_offset <- numeric(nrow(new_X))
          
          rm(tmp_formula, nullfit)
          }
     if(ncol(new_X) != ncol(object$betas))
          stop("Something went wrong! The number of columns in constructed model matrix should match the number of columns in object$betas.")

     if(!(object$family$family[1] %in% c("zipoisson","zinegative.binomial")) & !is.null(manualziX))
          stop("manualziX should only be supplied for zero-inflated CBFMs.")
     if(object$family$family[1] %in% c("zipoisson","zinegative.binomial")) {
          if(!is.null(manualziX)) {
               new_ziX <- manualziX
               warning("manualziX has been supplied. This overrides the creation of a model matrix based on object$ziformula and/or newdata.")
               }
          if(is.null(manualziX)) {
               tmp_formula <- as.formula(paste("response", paste(as.character(object$ziformula),collapse = " ") ) )
               nullfit <- gam(tmp_formula, data = data.frame(response = runif(nrow(object$y)), object$data), knots = object$ziknots, fit = TRUE, control = list(maxit = 1))
               if(is.null(newdata))
                    new_ziX <- predict.gam(nullfit, type = "lpmatrix")
               if(!is.null(newdata))
                    new_ziX <- predict.gam(nullfit, newdata = data.frame(newdata), type = "lpmatrix")
               
               ziformula_offset <- attributes(new_ziX)$model.offset
               if(all(ziformula_offset == 0))
                    ziformula_offset <- numeric(nrow(new_ziX))
               
               rm(tmp_formula, nullfit)
               }
          if(ncol(new_ziX) != ncol(object$zibetas))
               stop("Something went wrong! The number of columns in constructed model matrix should match the number of columns in object$zibetas.")
     }
     
          
     ## Construct B
     new_B <- NULL
     if(!is.null(new_B_space)) {
          if(object$num_B_space != ncol(new_B_space))
               stop("The number of columns of new_B_space does not object$num_B_space.")
          new_B <- cbind(new_B, new_B_space)
          }
     if(!is.null(new_B_time)) {
          if(object$num_B_time != ncol(new_B_time))
               stop("The number of columns of new_B_time does not object$num_B_time.")
          new_B <- cbind(new_B, new_B_time)
          }
     if(!is.null(new_B_spacetime)) {
          if(object$num_B_spacetime != ncol(new_B_spacetime))
               stop("The number of columns of new_B_spacetime does not object$num_B_spacetime.")
          new_B <- cbind(new_B, new_B_spacetime)
          }
     if(is.null(new_B))
          new_B <- object$B
     if(nrow(new_B) != nrow(new_X))
          stop("new_B does not contain the same number of rows as the model matrix created based on object and/or newdata.")
     if(ncol(new_B) != ncol(object$basis_effects_mat))
          stop("new_B does not contain the same number of columns as the number of columns in object$basis_effects_mat.")

     ## Final checks
     type <- match.arg(type, choices = c("link", "zilink", "response", "lpmatrix", "zilpmatrix"))
     if(!(object$family$family[1] %in% c("zipoisson","zinegative.binomial")) & type %in% c("zilink", "zilpmatrix"))
          stop("zilink/zilpmatrix can be only be obtained for zero-inflated CBFMs.")
     
     num_spp <- nrow(object$betas)
     num_X <- ncol(new_X)
     num_basisfns <- ncol(new_B)
     
     ## Start producing outputs
     if(type == "lpmatrix")
          return(new_X)
     if(type == "zilpmatrix")
          return(new_ziX)
     
     
     if(se_fit == TRUE & object$stderrors == FALSE)
          stop("Standard errors can not be calculated since the covariance matrix estimate was not detected to be available in object.") 

     ptpred <- tcrossprod(new_X, object$betas) + tcrossprod(new_B, object$basis_effects_mat)
     ptpred <- ptpred + formula_offset
     ptpred <- as.matrix(ptpred)
     colnames(ptpred) <- rownames(object$betas)          
     if(object$family$family[1] %in% c("zipoisson", "zinegative.binomial")) {
          ziptpred <- tcrossprod(new_ziX, object$zibetas) 
          ziptpred <- ziptpred + ziformula_offset
          ziptpred <- as.matrix(ziptpred)
          colnames(ziptpred) <- rownames(object$zibetas)          
          }
     
     if(!se_fit) {
          if(type == "zilink")
               return(ziptpred)
          if(type == "response") {
               if(!(object$family$family[1] %in% c("zipoisson", "zinegative.binomial","ztpoisson","ztnegative.binomial"))) 
                    ptpred <- object$family$linkinv(ptpred)
               if(object$family$family[1] %in% c("zipoisson", "zinegative.binomial"))
                    ptpred <- object$family$linkinv(ptpred) * (1-plogis(ziptpred))
               if(object$family$family[1] %in% c("ztpoisson")) {
                    ptpred <- exp(ptpred) / (1-exp(-exp(ptpred)))
                    ptpred[ptpred < 1 | ptpred == Inf] <- 1 # Predictions less than 1 will be almost certainly due to the linear predictor being so close to zero that you get underflow issues
                    }
               if(object$family$family[1] %in% c("ztnegative.binomial")) {
                    ptpred <- exp(ptpred) / (1-dnbinom(0, mu = exp(ptpred), size = matrix(1/object$dispparam, nrow(new_X), num_spp, byrow = TRUE)))                             
                    ptpred[ptpred < 1 | ptpred == Inf] <- 1 # Predictions less than 1 will be almost certainly due to the linear predictor being so close to zero that you get underflow issues
                    }
          }
     return(ptpred)
     }
     
     
     if(se_fit) {
          # This is needed *only* for predict.CBFM_hurdle downstream
          if(return_internals) {
               mu_vec <- as.vector(t(cbind(object$betas, object$basis_effects_mat)))
               bigcholcovar <- as.matrix(rbind(cbind(object$covar_components$topleft, object$covar_components$topright),
                                               cbind(t(object$covar_components$topright), object$covar_components$bottomright)))
               bigcholcovar <- t(chol(bigcholcovar))

               inner_linpredfn <- function(j) {
                    parameters_sim <- matrix(mu_vec + as.vector(bigcholcovar %*% rnorm(length(mu_vec))), nrow = num_spp, byrow = TRUE)
                    betas_sim <- parameters_sim[,1:num_X, drop=FALSE]
                    basiseff_sim <- parameters_sim[,-(1:num_X), drop=FALSE]
              
                    linpred <- tcrossprod(new_X, betas_sim) + tcrossprod(new_B, basiseff_sim)
                    linpred <- linpred + formula_offset
                    return(as.matrix(linpred))
                    }
                    
               alllinpreds <- foreach(j = 1:num_sims) %dopar% inner_linpredfn(j = j)
               rm(bigcholcovar, mu_vec)
               return(abind(alllinpreds, along = 3))
               }
                
                
           ci_alpha <- qnorm((1-coverage)/2, lower.tail = FALSE)
           need_sim <- FALSE
           if(object$family$family[1] %in% c("zipoisson","zinegative.binomial")) 
                need_sim <- TRUE
           if(object$family$family[1] %in% c("ztpoisson","ztnegative.binomial") & type == "response") 
                   need_sim <- TRUE
           
           if(!need_sim) {
                stderr1_fun <- function(j) {
                    out1 <- object$covar_components$topleft[(num_X*j - num_X + 1):(num_X*j), (num_X*j - num_X + 1):(num_X*j),drop=FALSE]
                    #return(diag(new_X %*% tcrossprod(out1, new_X)))
                    return(rowSums((new_X %*% out1) * new_X))
                    }
                stderr2_fun <- function(j, num_X, num_basisfns) {
                    out1 <- .extractcovarblocks_topright(j = j, Q = object$covar_components$topright, num_X = num_X, num_basisfns = num_basisfns)
                    #return(diag(new_X %*% tcrossprod(out1, new_B)))
                    return(rowSums((new_X %*% out1) * new_B))
                    }
                stderr3_fun <- function(j, num_basisfns) {
                    out1 <- .extractcovarblocks_bottomright(j = j, Q = object$covar_components$bottomright, num_basisfns = num_basisfns)
                    #return(diag(new_B %*% tcrossprod(out1, new_B)))
                    return(rowSums((new_B %*% out1) * new_B))
                    }

               stderr <- foreach(j = 1:num_spp, .combine = "cbind") %dopar% stderr1_fun(j = j)
               stderr <- stderr + 2 * foreach(j = 1:num_spp, .combine = "cbind") %dopar% stderr2_fun(j = j, num_X = num_X, num_basisfns = num_basisfns)
               stderr <- stderr + foreach(j = 1:num_spp, .combine = "cbind") %dopar% stderr3_fun(j = j, num_basisfns = num_basisfns)          
               stderr <- as.matrix(stderr)
               colnames(stderr) <- rownames(object$betas)    
          
               alllower <- ptpred - ci_alpha * sqrt(stderr) 
               allupper <- ptpred + ci_alpha * sqrt(stderr)
               if(type == "response") {
                    ptpred <- object$family$linkinv(ptpred)
                    alllower <- object$family$linkinv(alllower)
                    allupper <- object$family$linkinv(allupper)
                    }
               }

          if(need_sim) {
               message("Simulation required for prediction uncertainty intervals. This could take a while...grab a cuppa while you're waiting uwu")
               mu_vec <- as.vector(t(cbind(object$zibetas, object$betas, object$basis_effects_mat)))
               bigcholcovar <- as.matrix(rbind(cbind(object$covar_components$topleft, object$covar_components$topright),
                                               cbind(t(object$covar_components$topright), object$covar_components$bottomright)))
               check_eigen <- eigen(bigcholcovar, only.values = TRUE)
               if(any(check_eigen$values <= 0)) {
                    bigcholcovar <- Matrix::nearPD(bigcholcovar)
                    bigcholcovar <- bigcholcovar$mat
                    }
               rm(check_eigen)
               bigcholcovar <- t(chol(bigcholcovar))
                        
               if(object$family$family[1] %in% c("zipoisson","zinegative.binomial")) {
                    innerzip_predfn <- function(j) {
                         parameters_sim <- matrix(mu_vec + as.vector(bigcholcovar %*% rnorm(length(mu_vec))), nrow = num_spp, byrow = TRUE)
                         zibetas_sim <- parameters_sim[,1:ncol(object$zibetas)]
                         betas_sim <- parameters_sim[,ncol(object$zibetas)+(1:num_X), drop=FALSE]
                         basiseff_sim <- parameters_sim[,(ncol(object$zibetas)+num_X+1):ncol(parameters_sim), drop=FALSE]
                         
                         ptpred <- tcrossprod(new_X, betas_sim) + tcrossprod(new_B, basiseff_sim)
                         if(!is.null(formula_offset))0
                              ptpred <- ptpred + formula_offset
                         ziptpred <- tcrossprod(new_ziX, zibetas_sim) 
                         ziptpred <- ziptpred + ziformula_offset
                         if(type == "link")
                              return(as.matrix(ptpred))
                         if(type == "zilink")
                              return(as.matrix(ziptpred))
                         if(type == "response") {
                              ptpred <- as.matrix(object$family$linkinv(ptpred) * (1-plogis(ziptpred)))
                              return(ptpred)
                              }
                         }
                    
                         allpreds <- foreach(j = 1:num_sims) %dopar% innerzip_predfn(j = j)
                         }
               if(object$family$family[1] %in% c("ztpoisson") & type == "response") {
                    inner_predfn <- function(j) {
                         parameters_sim <- matrix(mu_vec + as.vector(bigcholcovar %*% rnorm(length(mu_vec))), nrow = num_spp, byrow = TRUE)
                         betas_sim <- parameters_sim[,1:num_X, drop=FALSE]
                         basiseff_sim <- parameters_sim[,-(1:num_X), drop=FALSE]
                        
                         ptpred <- tcrossprod(new_X, betas_sim) + tcrossprod(new_B, basiseff_sim)
                         ptpred <- ptpred + formula_offset
                         ptpred <- exp(ptpred) / (1-exp(-exp(ptpred)))
                         return(as.matrix(ptpred))
                         }
                    
                    allpreds <- foreach(j = 1:num_sims) %dopar% inner_predfn(j = j)
                    }
               if(object$family$family[1] %in% c("ztnegative.binomial") & type == "response") {
                    inner_predfn <- function(j) {
                         parameters_sim <- matrix(mu_vec + as.vector(bigcholcovar %*% rnorm(length(mu_vec))), nrow = num_spp, byrow = TRUE)
                         betas_sim <- parameters_sim[,1:num_X, drop=FALSE]
                         basiseff_sim <- parameters_sim[,-(1:num_X), drop=FALSE]
                        
                         ptpred <- tcrossprod(new_X, betas_sim) + tcrossprod(new_B, basiseff_sim)
                         ptpred <- ptpred + formula_offset
                         ptpred <- exp(ptpred) / (1-dnbinom(0, mu = exp(ptpred), size = matrix(1/object$dispparam, nrow(new_X), num_spp, byrow = TRUE)))
                         return(as.matrix(ptpred))
                         }
                    
                    allpreds <- foreach(j = 1:num_sims) %dopar% inner_predfn(j = j)
                    }
                        
               rm(bigcholcovar, mu_vec)
               allpreds <- abind(allpreds, along = 3)
               ptpred <- apply(allpreds, c(1,2), mean)
               stderr <- apply(allpreds, c(1,2), var)
               alllower <- apply(allpreds, c(1,2), quantile, prob = (1-coverage)/2)
               allupper <- apply(allpreds, c(1,2), quantile, prob = coverage + (1-coverage)/2) 
               rm(allpreds)
               gc()
               }

          return(list(fit = ptpred, stderr = sqrt(stderr), lower = alllower, upper = allupper))
          }
     }



#' @rdname predict.CBFM
#' @export
predict.CBFM_hurdle <- function(object, 
                         newdata_pa = NULL, manualX_pa = NULL, new_B_space_pa = NULL, new_B_time_pa = NULL, new_B_spacetime_pa = NULL, 
                         newdata_count = NULL, manualX_count = NULL, new_B_space_count = NULL, new_B_time_count = NULL, new_B_spacetime_count = NULL, 
                         se_fit = FALSE, coverage = 0.95, ncores = NULL, num_sims = 500, ...) {
        if(!inherits(object, "CBFM_hurdle")) 
                stop("`object' is not of class \"CBFM_hurdle\"")
        
        if(!se_fit) {
                preds_pa <- predict.CBFM(object$pa_fit, newdata = newdata_pa, manualX = manualX_pa,  
                                         new_B_space = new_B_space_pa, new_B_time = new_B_time_pa, new_B_spacetime = new_B_spacetime_pa,
                                         type = "response", se_fit = FALSE, coverage = coverage, ncores = ncores)
                preds_counts <- predict.CBFM(object$count_fit, newdata = newdata_count, manualX = manualX_count,  
                                         new_B_space = new_B_space_count, new_B_time = new_B_time_count, new_B_spacetime = new_B_spacetime_count,
                                         type = "response", se_fit = FALSE, coverage = coverage, ncores = ncores)
                
                return(preds_pa * preds_counts)
                }

        if(se_fit) {
                message("Simulation required for prediction uncertainty intervals. This could take a while...grab a cuppa while you're waiting uwu")
                preds_pa <- predict.CBFM(object$pa_fit, newdata = newdata_pa, manualX = manualX_pa,  
                                         new_B_space = new_B_space_pa, new_B_time = new_B_time_pa, new_B_spacetime = new_B_spacetime_pa,
                                         type = "response", se_fit = TRUE, coverage = coverage, ncores = ncores, num_sims = num_sims, return_internals = TRUE)
                preds_pa <- object$pa_fit$family$linkinv(preds_pa)
                preds_counts <- predict.CBFM(object$count_fit, newdata = newdata_count, manualX = manualX_count,  
                                         new_B_space = new_B_space_count, new_B_time = new_B_time_count, new_B_spacetime = new_B_spacetime_count,
                                         type = "response", se_fit = TRUE, coverage = coverage, ncores = ncores, num_sims = num_sims, return_internals = TRUE)
                if(object$count_fit$family$family[1] == "ztpoisson")
                        preds_counts <- exp(preds_counts) / (1-exp(-exp(preds_counts)))
                if(object$count_fit$family$family[1] == "ztnegative.binomial")
                        preds_counts <- exp(preds_counts) / (1-dnbinom(0, mu = exp(preds_counts), size = matrix(1/object$count_fit$dispparam, nrow(preds_counts), ncol(preds_counts), byrow = TRUE)))
                final_preds <- preds_pa * preds_counts 
                rm(preds_pa, preds_counts)
                
                ptpred <- apply(final_preds, c(1,2), mean)
                stderr <- apply(final_preds, c(1,2), var)
                alllower <- apply(final_preds, c(1,2), quantile, prob = (1-coverage)/2)
                allupper <- apply(final_preds, c(1,2), quantile, prob = coverage + (1-coverage)/2) 
                rm(final_preds)
                gc()
                return(list(fit = ptpred, stderr = sqrt(stderr), lower = alllower, upper = allupper))
                }
        }
