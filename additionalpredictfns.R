predict_PA_CBFM <- function(object, newdata = NULL, manualX = NULL, new_B_space = NULL, new_B_time = NULL, new_B_spacetime = NULL, 
                         type = "link", ncores = NULL, parameter_uncertainty = FALSE, num_sims = 100, ...) {
     
     se_fit <- FALSE # Switching off standard errors for this function
     
     if(!inherits(object, "CBFM")) 
          stop("`object' is not of class \"CBFM\"")
     if(object$family$family != "binomial")
          stop("This prediction function is specifically designed for binary responses!")

     if(is.null(ncores))
          registerDoParallel(cores = detectCores()-1)
     if(!is.null(ncores))
          registerDoParallel(cores = ncores)
        
     if(!is.null(manualX)) {
          new_X <- manualX
          warning("manualX has been supplied. This overrides the creation of a model matrix based on object$formula_X and/or newdata.")
          }
     if(is.null(manualX)) {
          tmp_formula <- as.formula(paste("response", paste(as.character(object$formula),collapse="") ) )
          nullfit <- gam(tmp_formula, data = data.frame(response = runif(nrow(object$y)), object$data), fit = TRUE, control = list(maxit = 1))
          if(is.null(newdata))
               new_X <- predict.gam(nullfit, type = "lpmatrix")
           if(!is.null(newdata))
                new_X <- predict.gam(nullfit, newdata = data.frame(newdata), type = "lpmatrix")
           rm(tmp_formula, nullfit)
           }
     if(ncol(new_X) != ncol(object$betas))
          stop("Number of columns in new_X should match the number of columns in object$betas.")

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
     if(!is.null(new_B)) {
           if(nrow(new_B) != nrow(new_X))
               stop("new_B does not contain the same number of rows as the model matrix created based on object and/or newdata.")
           if(ncol(new_B) != ncol(object$basis_effects_mat))
               stop("new_B does not contain the same number of columns as the number of columns in object$basis_effects_mat.")
           }
          
     type <- match.arg(type, choices = c("link", "response", "lpmatrix", "y"))
     num_spp <- nrow(object$betas)
     num_X <- ncol(new_X)
     num_basisfns <- ncol(new_B)
     
     if(type == "lpmatrix")
           return(new_X)
        
     
     ptpred <- tcrossprod(new_X, object$betas) + tcrossprod(new_B, object$basis_effects_mat)
     ptpred <- as.matrix(ptpred)
     colnames(ptpred) <- rownames(object$betas)          
          
     if(!parameter_uncertainty) {
          if(type == "response") {
               ptpred <- object$family$linkinv(ptpred)
               }
          if(type == "y") {
               ptpred <- object$family$linkinv(ptpred)
               ptpred <- replicate(num_sims, matrix(rbinom(length(ptpred), size = object$trial_size, prob = ptpred), ncol = num_spp))
               dimnames(ptpred) <- list(units = rownames(new_X), spp = rownames(object$betas), replications = 1:num_sims)
               }
          
          return(ptpred)
          }
     
        
     if(parameter_uncertainty) {
          mu_vec <- as.vector(t(cbind(object$betas, object$basis_effects_mat)))
          bigcholcovar <- as.matrix(rbind(cbind(object$covar_components$topleft, object$covar_components$topright),
                                          cbind(t(object$covar_components$topright), object$covar_components$bottomright)))
          check_eigen <- eigen(bigcholcovar, only.values = TRUE)
          if(any(check_eigen$values <= 0)) {
               bigcholcovar <- Matrix::nearPD(bigcholcovar)
               bigcholcovar <- bigcholcovar$mat 
               }
          rm(check_eigen)
          bigcholcovar <- t(chol(bigcholcovar))

          innerzip_predfn <- function(j) {
               parameters_sim <- matrix(mu_vec + as.vector(bigcholcovar %*% rnorm(length(mu_vec))), nrow = num_spp, byrow = TRUE)
               betas_sim <- parameters_sim[,1:num_X, drop=FALSE]
               basiseff_sim <- parameters_sim[,-(1:num_X), drop=FALSE]
                        
               ptpred <- tcrossprod(new_X, betas_sim) + tcrossprod(new_B, basiseff_sim)
               ptpred <- as.matrix(ptpred)
               return(ptpred)
               }
                    
          allpreds <- foreach(j = 1:num_sims) %dopar% innerzip_predfn(j = j)
          rm(bigcholcovar, mu_vec)
          allpreds <- abind(allpreds, along = 3)
          
          if(type == "response") {
               allpreds <- object$family$linkinv(allpreds)
               }
          if(type == "y") {
               allpreds <- object$family$linkinv(allpreds)
               ptpreds <- allpreds
               for(k0 in 1:num_sims) {
                    ptpreds[,,k0] <- matrix(rbinom(length(allpreds[,,k0]), size = object$trial_size, prob = allpreds[,,k0]), ncol = num_spp)
                    }
               allpreds <- ptpreds
               }
          
          dimnames(allpreds) <- list(units = rownames(new_X), spp = rownames(object$betas), replications = 1:num_sims)
          gc()
          return(allpreds)
          }
        }

predict_ztcount_CBFM <- function(object, newdata = NULL, manualX = NULL, new_B_space = NULL, new_B_time = NULL, new_B_spacetime = NULL, B_dampen = NULL,
                            type = "link", ncores = NULL, parameter_uncertainty = FALSE, num_sims = 100, ...) {
     
     se_fit <- FALSE # Switching off standard errors for this function
     
     if(!inherits(object, "CBFM")) 
          stop("`object' is not of class \"CBFM\"")
     if(!(object$family$family %in% c("ztpoisson","ztnegative.binomial")))
          stop("This prediction function is specifically designed for zero-truncated count responses!")
     
     if(is.null(ncores))
          registerDoParallel(cores = detectCores()-1)
     if(!is.null(ncores))
          registerDoParallel(cores = ncores)
     
     if(!is.null(manualX)) {
          new_X <- manualX
          warning("manualX has been supplied. This overrides the creation of a model matrix based on object$formula_X and/or newdata.")
     }
     if(is.null(manualX)) {
          tmp_formula <- as.formula(paste("response", paste(as.character(object$formula),collapse="") ) )
          nullfit <- gam(tmp_formula, data = data.frame(response = runif(nrow(object$y)), object$data), fit = TRUE, control = list(maxit = 1))
          if(is.null(newdata))
               new_X <- predict.gam(nullfit, type = "lpmatrix")
          if(!is.null(newdata))
               new_X <- predict.gam(nullfit, newdata = data.frame(newdata), type = "lpmatrix")
          rm(tmp_formula, nullfit)
     }
     if(ncol(new_X) != ncol(object$betas))
          stop("Number of columns in new_X should match the number of columns in object$betas.")
     
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
     if(!is.null(new_B)) {
          if(nrow(new_B) != nrow(new_X))
               stop("new_B does not contain the same number of rows as the model matrix created based on object and/or newdata.")
          if(ncol(new_B) != ncol(object$basis_effects_mat))
               stop("new_B does not contain the same number of columns as the number of columns in object$basis_effects_mat.")
          }
     
      if(!is.null(B_dampen)) {
          if(nrow(B_dampen) != nrow(B_dampen))
               stop("B_dampen should be a matrix with the same number of rows as the corresponding basis function matrix.")
          }


     type <- match.arg(type, choices = c("link", "response", "lpmatrix", "y"))
     num_spp <- nrow(object$betas)
     num_X <- ncol(new_X)
     num_basisfns <- ncol(new_B)
     
     if(type == "lpmatrix")
          return(new_X)
     
     ptpred <- tcrossprod(new_X, object$betas)
     if(is.null(B_dampen))
          ptpred <- ptpred  + tcrossprod(new_B, object$basis_effects_mat)
     if(!is.null(B_dampen)) {
          for(j in 1:nrow(object$betas))
               ptpred[,j] <- ptpred[,j] + as.vector((new_B * as.vector(unlist(B_dampen[,j]))) %*% object$basis_effects_mat[j,])
          }
     ptpred <- as.matrix(ptpred)
     colnames(ptpred) <- rownames(object$betas)          
     
     ztpR <- gamlss.tr::trun.r(par = 0, family = gamlss.dist::PO()$family[1], type = "left") 
     ztnbR <- gamlss.tr::trun.r(par = 0, family = gamlss.dist::NBI()$family[1], type = "left") 

          
     if(!parameter_uncertainty) {
          if(type == "response") {
               if(object$family$family[1] %in% c("ztpoisson")) {
                    ptpred <- exp(ptpred) / (1-exp(-exp(ptpred)))
                    ptpred[ptpred < 1 | ptpred == Inf] <- 1 # Predictions less than 1 will be almost certainly due to the linear predictor being so close to zero that you get underflow issues
                    }
               if(object$family$family[1] %in% c("ztnegative.binomial")) {
                    ptpred <- exp(ptpred) / (1-dnbinom(0, mu = exp(ptpred), size = matrix(1/object$dispparam, nrow = nrow(ptpred), ncol = num_spp, byrow = TRUE)))                             
                    ptpred[ptpred < 1 | ptpred == Inf] <- 1 # Predictions less than 1 will be almost certainly due to the linear predictor being so close to zero that you get underflow issues
                    }
               }
          if(type == "y") {
               ptpred <- exp(ptpred)
               if(object$family$family[1] == "ztpoisson")
                    ptpred <- replicate(num_sims, matrix(ztpR(length(ptpred), mu = ptpred), ncol = num_spp))
               if(object$family$family[1] == "ztnegative.binomial")
                    ptpred <- replicate(num_sims, matrix(ztnbR(length(ptpred), mu = ptpred, 
                                                               sigma = matrix(object$dispparam, nrow = nrow(ptpred), ncol = num_spp, byrow = TRUE)), 
                                                         ncol = num_spp))
               
               dimnames(ptpred) <- list(units = rownames(new_X), spp = rownames(object$betas), replications = 1:num_sims)
               }
          
          return(ptpred)
          }
     
     
     if(parameter_uncertainty) {
          mu_vec <- as.vector(t(cbind(object$betas, object$basis_effects_mat)))
          bigcholcovar <- as.matrix(rbind(cbind(object$covar_components$topleft, object$covar_components$topright),
                                          cbind(t(object$covar_components$topright), object$covar_components$bottomright)))
          check_eigen <- eigen(bigcholcovar, only.values = TRUE)
          if(any(check_eigen$values <= 0)) {
               bigcholcovar <- Matrix::nearPD(bigcholcovar)
               bigcholcovar <- bigcholcovar$mat
               }
          rm(check_eigen)
          bigcholcovar <- t(chol(bigcholcovar))
          
          innerzip_predfn <- function(j) {
               parameters_sim <- matrix(mu_vec + as.vector(bigcholcovar %*% rnorm(length(mu_vec))), nrow = num_spp, byrow = TRUE)
               betas_sim <- parameters_sim[,1:num_X, drop=FALSE]
               basiseff_sim <- parameters_sim[,-(1:num_X), drop=FALSE]
               
               ptpred <- tcrossprod(new_X, betas_sim)
               if(is.null(B_dampen))
                    ptpred <- ptpred  + tcrossprod(new_B, basiseff_sim)
               if(!is.null(B_dampen)) {
                    for(j in 1:nrow(object$betas))
                         ptpred[,j] <- ptpred[,j] + (new_B * as.vector(unlist(B_dampen[,j]))) %*% basiseff_sim[j,]
                    }

               ptpred <- as.matrix(ptpred)
               return(ptpred)
               }
          
          allpreds <- foreach(j = 1:num_sims) %dopar% innerzip_predfn(j = j)
          rm(bigcholcovar, mu_vec)
          allpreds <- abind(allpreds, along = 3)
          
          if(type == "response") {
               if(object$family$family[1] %in% c("ztpoisson")) {
                    allpreds <- exp(allpreds) / (1-exp(-exp(allpreds)))
                    allpreds[allpreds < 1 | allpreds == Inf] <- 1 # Predictions less than 1 will be almost certainly due to the linear predictor being so close to zero that you get underflow issues
                    }
               if(object$family$family[1] %in% c("ztnegative.binomial")) {
                    for(k0 in 1:num_sims) {
                         allpreds[,,k0] <- exp(allpreds[,,k0]) / (1-dnbinom(0, mu = exp(allpreds[,,k0]), 
                                                                            size = matrix(1/object$dispparam, nrow = nrow(allpreds[,,k0]), ncol = num_spp, byrow = TRUE)))                         
                         }
                    allpreds[allpreds < 1 | allpreds == Inf] <- 1 # Predictions less than 1 will be almost certainly due to the linear predictor being so close to zero that you get underflow issues
                    }
               }
          if(type == "y") {
               ptpreds <- allpreds
               allpreds <- exp(allpreds)
               for(k0 in 1:num_sims) {
                    if(object$family$family[1] == "ztpoisson")
                         ptpreds[,,k0] <- matrix(ztpR(length(allpreds[,,k0]), mu = allpreds[,,k0]+1e-12), ncol = num_spp)
                    if(object$family$family[1] == "ztnegative.binomial")
                         ptpreds[,,k0] <- matrix(ztnbR(length(allpreds[,,k0]), mu = allpreds[,,k0]+1e-12, 
                                                       sigma = matrix(object$dispparam, nrow = nrow(allpreds[,,k0]), ncol = num_spp, byrow = TRUE)), 
                                                 ncol = num_spp)
                    }
               allpreds <- ptpreds
               }
          
          dimnames(allpreds) <- list(units = rownames(new_X), spp = rownames(object$betas), replications = 1:num_sims)
          gc()
          return(allpreds)
          }
     }



predict_hurdle_CBFM <- function(object, prediction_form = "mean", linex_a = 0.001, quantile_threshold = 0.9,
                                  newdata_pa = NULL, manualX_pa = NULL, new_B_space_pa = NULL, new_B_time_pa = NULL, new_B_spacetime_pa = NULL, 
                                  newdata_count = NULL, manualX_count = NULL, new_B_space_count = NULL, new_B_time_count = NULL, new_B_spacetime_count = NULL, 
                                  se_fit = FALSE, ncores = NULL, parameter_uncertainty = FALSE, num_sims = 100, ...) {
     
     ##-------------------
     #' # Checks and balances
     ##-------------------
     if(!inherits(object, "CBFM_hurdle"))  
          stop("`object' is not of class \"CBFM_hurdle\"")
     if(object$count_fit$family[1] != "ztnegative.binomial")  
          stop("The current function is designed only of zero-truncated negative binomial hurdle CBFMs.")
     
     prediction_form <- match.arg(prediction_form, choices = c("mean", "median", "linex", "threshold"))
     if(linex_a == 0) 
          stop("For linex prediction, the tuning parameter a must not equal to zero.")
     if(quantile_threshold <= 0 | quantile_threshold >= 1) 
          stop("For threshold prediction, the quantile threshold must be between zero and one.")
     se_fit <- FALSE # Switching off standard errors for this function
     
     
     ##-------------------
     #' # PA component prediction
     ##-------------------
     preds_pa <- predict_PA_CBFM(object$pa_fit, newdata = newdata_pa, manualX = manualX_pa,  
                                 new_B_space = new_B_space_pa, new_B_time = new_B_time_pa, new_B_spacetime = new_B_spacetime_pa,
                                 type = "response", coverage = coverage, parameter_uncertainty = parameter_uncertainty, num_sims = num_sims, ncores = ncores)

     all_dispparam <- matrix(object$count_fit$dispparam, nrow = nrow(preds_pa), ncol = ncol(preds_pa), byrow = TRUE)
     
     ##-------------------
     #' # If parameter uncertainty is not required
     ##-------------------
     if(length(dim(preds_pa)) == 2) {
          if(prediction_form %in% c("mean", "threshold")) {
               preds_counts <- predict_ztcount_CBFM(object$count_fit, newdata = newdata_count, manualX = manualX_count,  
                                            new_B_space = new_B_space_count, new_B_time = new_B_time_count, new_B_spacetime = new_B_spacetime_count,
                                            type = "response", se_fit = FALSE, coverage = coverage, parameter_uncertainty = parameter_uncertainty, 
                                            num_sims = num_sims, ncores = ncores)
               
               final_prediction <- preds_pa * preds_counts
               }
     
          if(prediction_form == "median") {
               preds_counts <- predict_ztcount_CBFM(object$count_fit, newdata = newdata_count, manualX = manualX_count,  
                                            new_B_space = new_B_space_count, new_B_time = new_B_time_count, new_B_spacetime = new_B_spacetime_count,
                                            type = "link", se_fit = FALSE, coverage = coverage, parameter_uncertainty = parameter_uncertainty, 
                                            num_sims = num_sims, ncores = ncores)
               
               final_prediction <- gamlss.dist::qZANBI(0.5, 
                                                       mu = exp(preds_counts)+1e-12,  
                                                       sigma = all_dispparam,
                                                       nu = 1 - preds_pa)
               final_prediction <- matrix(final_prediction, nrow = nrow(preds_pa), ncol = ncol(preds_pa))
               }
          
          if(prediction_form == "linex") {
               preds_counts <- predict_ztcount_CBFM(object$count_fit, newdata = newdata_count, manualX = manualX_count,  
                                            new_B_space = new_B_space_count, new_B_time = new_B_time_count, new_B_spacetime = new_B_spacetime_count,
                                            type = "link", se_fit = FALSE, coverage = coverage, parameter_uncertainty = parameter_uncertainty, 
                                            num_sims = num_sims, ncores = ncores)
               
               final_prediction <- .linexpred_hurdlenb(a = linex_a, 
                                                       p = preds_pa, 
                                                       mu = exp(preds_counts)+1e-12,
                                                       phi = all_dispparam)
               }
          
          if(prediction_form == "threshold") {
               preds_counts <- predict_ztcount_CBFM(object$count_fit, newdata = newdata_count, manualX = manualX_count,  
                                            new_B_space = new_B_space_count, new_B_time = new_B_time_count, new_B_spacetime = new_B_spacetime_count,
                                            type = "link", se_fit = FALSE, coverage = coverage, parameter_uncertainty = parameter_uncertainty, 
                                            num_sims = num_sims, ncores = ncores)
               
               final_prediction <- pmin(final_prediction,
                                        gamlss.dist::qZANBI(quantile_threshold, 
                                                            mu = exp(preds_counts)+1e-12,  
                                                            sigma = all_dispparam,
                                                            nu = 1-preds_pa)) 
               final_prediction <- matrix(final_prediction, nrow = nrow(preds_pa), ncol = ncol(preds_pa))
               }
          
          
          return(final_prediction)
          }
          
     ##-------------------
     #' # If parameter uncertainty is needed required
     ##-------------------
     if(length(dim(preds_pa)) == 3) {
          if(prediction_form %in% c("mean", "threshold")) {
               preds_counts <- predict_ztcount_CBFM(object$count_fit, newdata = newdata_count, manualX = manualX_count,  
                                            new_B_space = new_B_space_count, new_B_time = new_B_time_count, new_B_spacetime = new_B_spacetime_count,
                                            type = "response", se_fit = FALSE, coverage = coverage, parameter_uncertainty = parameter_uncertainty, 
                                            num_sims = num_sims, ncores = ncores)
               
               final_prediction <- array(NA, dim = dim(preds_pa))
               for(k0 in 1:dim(preds_pa)[3])
                    final_prediction[,,k0] <- preds_pa[,,k0] * preds_counts[,,k0]
          }
          
          if(prediction_form == "median") {
               preds_counts <- predict_ztcount_CBFM(object$count_fit, newdata = newdata_count, manualX = manualX_count,  
                                            new_B_space = new_B_space_count, new_B_time = new_B_time_count, new_B_spacetime = new_B_spacetime_count,
                                            type = "link", se_fit = FALSE, coverage = coverage, parameter_uncertainty = parameter_uncertainty, 
                                            num_sims = num_sims, ncores = ncores)
       
               final_prediction <- array(NA, dim = dim(preds_pa))
               for(k0 in 1:dim(preds_pa)[3]) {
                    final_prediction[,,k0] <- matrix(gamlss.dist::qZANBI(0.5, 
                                                       mu = exp(preds_counts[,,k0])+1e-12,  
                                                       sigma = all_dispparam,
                                                       nu = 1 - preds_pa[,,k0]), nrow = nrow(preds_pa), ncol = ncol(preds_pa))
                    }
          }
          
          if(prediction_form == "linex") {
               preds_counts <- predict_ztcount_CBFM(object$count_fit, newdata = newdata_count, manualX = manualX_count,  
                                            new_B_space = new_B_space_count, new_B_time = new_B_time_count, new_B_spacetime = new_B_spacetime_count,
                                            type = "link", se_fit = FALSE, coverage = coverage, parameter_uncertainty = parameter_uncertainty, 
                                            num_sims = num_sims, ncores = ncores)
               
               final_prediction <- array(NA, dim = dim(preds_pa))
               for(k0 in 1:dim(preds_pa)[3]) {
                    final_prediction[,,k0] <- .linexpred_hurdlenb(a = linex_a, 
                                                            p = preds_pa[,,k0], 
                                                            mu = exp(preds_counts[,,k0])+1e-12,
                                                            phi = all_dispparam)
                    }
               }
          
          if(prediction_form == "threshold") {
               preds_counts <- predict_ztcount_CBFM(object$count_fit, newdata = newdata_count, manualX = manualX_count,  
                                            new_B_space = new_B_space_count, new_B_time = new_B_time_count, new_B_spacetime = new_B_spacetime_count,
                                            type = "link", se_fit = FALSE, coverage = coverage, parameter_uncertainty = parameter_uncertainty, 
                                            num_sims = num_sims, ncores = ncores)
               
               final_prediction <- array(NA, dim = dim(preds_pa))
               for(k0 in 1:dim(preds_pa)[3]) {
                    final_prediction[,,k0] <- matrix(pmin(final_prediction[,,k0],
                                                   gamlss.dist::qZANBI(quantile_threshold, 
                                                                       mu = exp(preds_counts[,,k0])+1e-12,  
                                                                       sigma = all_dispparam,
                                                                       nu = 1-preds_pa[,,k0])),
                                                   nrow = nrow(preds_pa), ncol = ncol(preds_pa))
                    }
               }
          
          
          return(final_prediction)
          }
     
     }



.linexpred_hurdlenb <- function(a = 0.01, p, phi, mu) { 
     out <- (1-p) + p * (1/(1+phi*mu))^(1/phi) * (1 - (1/(1+phi*mu))^(1/phi))^(-1) * ((1 - mu*exp(-a)/(mu+1/phi))^(-1/phi) - 1) # E(exp(-aY)) where Y ~ HurdleNB(p, mu, phi), and p is probability of presence; from https://etd.ohiolink.edu/apexprod/rws_etd/send_file/send?accession=osu1543573678017356&disposition=inline
     return(-1/a * log(out))
     }


#' @title Construct species-specific dampening weights for the basis functions.
#' @abtract After distances between the prediction spatial locations and the confidence polygon are calculated, the dampening weights are constructed as \code{exp(- dampening_scale * distance / mean(distance) )}. The argument can either be a scalar or a vector of length equal to \code{ncol(object$y)} i.e., the number of species.
#'
#' This is an **experimental and ad-hoc** function/method for altering the predictions, and should not be used without first consulting with authors i.e., Chris Haak!

#' @param object A CBFM object. For the current set up, it more or less has to be a zero-truncated CBFM fit e.g., forming the count component of a hurdle CBFM.
#' @param new_spatial_locations A matrix of spatial locations on which the dampening weights are to be construct. Defaults to \code{NULL}, in which case it is set to \code{object$data[,coordinate_columns]}.
#' @param method Which method by which to construct the confidence polygons. Currently accepted method are \code{data_ellipse}, which uses the [car::dataEllipse()] function, or \code{kde}, which uses the [adehabitatHR::kernelUD()] function. The latter is preferred as it tends to be more stable.
#' @param level The percent level by which to construct the confidence polygons.
#' @param dampening_scale The dampening factor. Essentially, after distances between the prediction spatial locations and the confidence polygon are calculated, the dampening weights are constructed as \code{exp(- dampening_scale * distance / mean(distance) )}. The argument can either be a scalar or a vector of length equal to \code{ncol(object$y)} i.e., the number of species.
#' @param coordinate_columns The columns of \code{object$data} which contain the spatial coordinates. These must be coordinates for which Euclidean distances are applicable.
#' @param ncores The number of cores to use. Defaults to \code{NULL}, in which it uses \code{detectCores()-2}.

construct_dampening_weights <- function(object, new_spatial_locations = NULL, method = "kde", level = 0.95, dampening_scale,
                                        coordinate_columns = c("UTM_X", "UTM_Y"), ncores = NULL) {
    if(!inherits(object, "CBFM"))
        stop("`object' is not of class \"CBFM\"")

    method <- match.arg(method, choices = c("data_ellipse", "kde"))

    if(is.null(ncores))
        ncores <- detectCores()-2
    registerDoParallel(ncores)

    if(length(dampening_scale) == 1)
        dampening_scale <- rep(dampening_scale, ncol(object$y))
    if(!(length(dampening_scale) %in% c(1, ncol(object$y))))
        stop("dampening_scale must have length equal to 1 or the number of species i.e., ncol(object$y).")

    if(is.null(new_spatial_locations))
        new_spatial_locations <- object$data[,coordinate_columns]
    if(!is.null(new_spatial_locations))
        colnames(new_spatial_locations) <- coordinate_columns



    # Make sf object from all points
    mypoints <- sf::st_as_sf(new_spatial_locations, coords = coordinate_columns)


    create_poly_fn <- function(sel_spp) {
        #' Use dataellipse from car package to generate ellipse enclosing 95% of the data, and turn this into an sf object
        if(method == "data_ellipse") {
            myelips <- car::dataEllipse(
                x = as.vector(unlist(object$data[coordinate_columns[1]]))[(which(object$y[,sel_spp] > 0))],
                y = as.vector(unlist(object$data[coordinate_columns[2]]))[(which(object$y[,sel_spp] > 0))],
                levels = level,
                draw = FALSE)
            mypoly <- sf::st_polygon(x = list(myelips))
        }

        #' ## Alternative approach based on Kernel density estimation
        if(method == "kde") {
            #' Use kernelUD to estimate a probability density distribution of the data
            #kernel_ref <- kernelUD(nzpoints_sp, h = "href", grid=100, kern="bivnorm")
            kernel_ref <- adehabitatHR::kernelUD(SpatialPoints(coords = cbind(
                as.vector(unlist(object$data[coordinate_columns[1]]))[(which(object$y[,sel_spp] > 0))],
                as.vector(unlist(object$data[coordinate_columns[2]]))[(which(object$y[,sel_spp] > 0))])),
                h = "href",
                grid = 100,
                kern = "epa")
            #kernel_ref <- kernelUD(nzpoints_sp, h = "LSCV", grid=100)

            #' Get vertices
            mykernelpoly <- adehabitatHR::getverticeshr(kernel_ref, percent = 100*level)
            #plot(mykernelpoly)
            #plot(nzpoints_sp, add = TRUE)

            #' Convert poly to sf
            mypoly <- sf::st_as_sf(mykernelpoly)
        }

        #' Find the points that fall inside polygon/kernel, then take inverse (this assigns NAs to pts outside the polygon)
        outpoints <- mypoints[is.na(as.numeric(sf::st_intersects(mypoints, mypoly, sparse = TRUE))),]
        #outpoints<- mypoints[lengths(st_intersects(mypoints, mypoly))==0,]
        #plot(outpoints)
        #plot(mypoly, add = TRUE)

        #' Get distance for all points, and assigns 0 to any points inside the poly
        mypoints$distance <- sf::st_distance(mypoints, mypoly)

        return(mypoints$distance)
    }

    all_distances <- foreach(k = 1:ncol(object$y), .combine = "cbind") %dopar% create_poly_fn(sel_spp = k)
    all_distances <- Matrix::Matrix(all_distances, sparse = TRUE)


    #' Now calculate weights
    average_nonzero_distance <- apply(all_distances, 2, function(x) mean(x[x>0]))
    BFweights <- exp(-all_distances / matrix(average_nonzero_distance/dampening_scale, nrow = nrow(all_distances), ncol = ncol(all_distances), byrow = TRUE))
    #summary(mypoints$BFweight)
    #plot(mypoints["BFweight"])

    return(BFweights)
    }

