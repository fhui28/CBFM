predict_PA_CBFM <- function(object, newdata = NULL, manualX = NULL, new_B_space = NULL, new_B_time = NULL, new_B_spacetime = NULL, 
                         type = "link", coverage = 0.95, ncores = NULL, parameter_uncertainty = FALSE, num_sims = 10, ...) {
     
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

predict_ztcount_CBFM <- function(object, newdata = NULL, manualX = NULL, new_B_space = NULL, new_B_time = NULL, new_B_spacetime = NULL, 
                            type = "link", coverage = 0.95, ncores = NULL, parameter_uncertainty = FALSE, num_sims = 10, ...) {
     
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
     
     type <- match.arg(type, choices = c("link", "response", "lpmatrix", "y"))
     num_spp <- nrow(object$betas)
     num_X <- ncol(new_X)
     num_basisfns <- ncol(new_B)
     
     if(type == "lpmatrix")
          return(new_X)
     
     
     ptpred <- tcrossprod(new_X, object$betas) + tcrossprod(new_B, object$basis_effects_mat)
     ptpred <- as.matrix(ptpred)
     colnames(ptpred) <- rownames(object$betas)          
     
     ztpR <- gamlss.tr::trun.r(par = 0, family = PO()$family[1], type = "left") 
     ztnbR <- gamlss.tr::trun.r(par = 0, family = NBI()$family[1], type = "left") 

          
     if(!parameter_uncertainty) {
          if(type == "response") {
               ptpred <- object$family$linkinv(ptpred)
               }
          if(type == "y") {
               ptpred <- object$family$linkinv(ptpred)
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
                    if(object$family$family[1] == "ztpoisson")
                         ptpreds[,,k0] <- matrix(ztpR(length(allpreds[,,k0]), mu = allpreds[,,k0]), ncol = num_spp)
                    if(object$family$family[1] == "ztnegative.binomial")
                         ptpreds[,,k0] <- matrix(ztnbR(length(allpreds[,,k0]), mu = allpreds[,,k0], 
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


## To get a hurdle prediction, do something like
#' pa_predictions <- predict_PA_CBFM(object = fit$pa_fit, xxx)
#' ztcount_predictions <- predict_ztcount_CBFM(object = fit$count_fit, xxx)
#' hurdle_predictions <- pa_predictions * ztcount_predictions