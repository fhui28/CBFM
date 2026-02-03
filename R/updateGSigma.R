## Functions for updating G and Sigma
## Hidden and not exported

#' @noRd

.update_G_fn <- function(Ginv, basis_effects_mat, Sigmainv, B, X, ziX = NULL, zioffset = NULL, 
                         y_vec, linpred_vec, dispparam, powerparam, zibetas,
                         trial_size, family, G_control, use_rank_element, return_correlation) {
     
     num_spp <- nrow(basis_effects_mat)
     num_basisfns <- ncol(Sigmainv)
     trial_size <- as.vector(trial_size)
     use_structure <- G_control$structure[use_rank_element]
     #use_min_sp <- G_control$min_sp[use_rank_element]
     
     #if(num_spp == 1) 
     #        return(solve(Ginv))
     
     G_control$method <- match.arg(G_control$method, choices = c("simple","REML","ML")) 
     use_structure <- match.arg(use_structure, choices = c("unstructured","identity","homogeneous")) 
     G_control$inv_method <- "chol2inv" #match.arg(G_control$inv_method, choices = c("chol2inv","schulz"))
     A_Sigmain_AT <- basis_effects_mat %*% tcrossprod(Sigmainv, basis_effects_mat)

     if(family$family == "binomial") {
          linpred_vec[linpred_vec > 5] <- 5
          linpred_vec[linpred_vec < -5] <- -5
          }

     if(G_control$method == "simple") {
          if(use_structure == "unstructured") 
               new_G <- A_Sigmain_AT / num_basisfns
          if(use_structure == "homogeneous") #' Pass the estimation of lambda in G = lambda * R into estimation of Sigma
               new_G <- A_Sigmain_AT / num_basisfns
          if(use_structure == "identity")  #' Pass the estimation of lambda in G = lambda * R into estimation of Sigma
               new_G <- Matrix::Diagonal(n = num_spp)
          }
     
     if(G_control$method %in% c("REML","ML")) {
          zieta <- NULL
          if(family$family[1] %in% c("zipoisson","zinegative.binomial")) {                        
               zieta <- tcrossprod(ziX, zibetas)
               if(!is.null(zioffset))
                    zieta <- zieta + zioffset
               zieta <- as.vector(zieta)
               }
          
          ## Note weights are on a per-species basis i.e., site runs faster than species
          num_units <- nrow(B)
          weights_mat <- .neghessfamily(family = family,
                                        eta = linpred_vec,
                                        y = y_vec,
                                        phi = rep(dispparam, each = num_units),
                                        powerparam = rep(powerparam, each = num_units),
                                        zieta = zieta,
                                        trial_size = trial_size)

          ## Set up REML and ML
          weights_mat[is.na(y_vec)] <- 0
          weights_mat <- matrix(weights_mat, nrow = nrow(B), ncol = num_spp, byrow = FALSE)
          if(G_control$method == "REML") {
               inner_fn <- function(j) {
                    w_j <- weights_mat[,j]
                    sqrt_w <- sqrt(w_j)

                    XTX_inv <- chol2inv(chol(crossprod(X*sqrt_w) + Diagonal(x = 1e-8, n = ncol(X))))
                    BTWX <- crossprod(B, X*w_j)
                    return(crossprod(B*sqrt_w) - BTWX %*% tcrossprod(XTX_inv, BTWX))
                    }
               }
          if(G_control$method == "ML") {
               inner_fn <- function(j) {
                    sqrt_w <- sqrt(weights_mat[,j])
                    return(crossprod(B * sqrt_w))
                    }
               }
          
          if(use_structure %in% c("unstructured", "homogeneous")) {
               BtKB <- foreach(j = 1:num_spp) %dopar% inner_fn(j = j) 
               BtKB <- Matrix::bdiag(BtKB)
               gc()
      
               vecPsi <- do.call(rbind, lapply(1:num_basisfns, function(k) kronecker(Matrix::Diagonal(n = num_spp), Sigmainv[,k]))) 
               Q2 <- kronecker(Matrix::Diagonal(n = num_spp), vecPsi)
               needonly_cols <- unlist(sapply(1:num_spp, function(j) return((j-1)*num_spp + j:num_spp)))
               Q2 <- Q2[, needonly_cols, drop = FALSE]
               rm(vecPsi, weights_mat)
               gc()
      
               counter <- 0
               diff <- 10
               cw_G <- chol2inv(chol(Ginv))
      
               while(diff > G_control$tol & counter < G_control$maxit) {
                    cw_Ginv_Sigmainv <- Matrix::Matrix(kronecker(chol2inv(chol(cw_G)), Sigmainv), sparse = TRUE)
           
                    if(G_control$inv_method == "chol2inv")
                         tic <- proc.time()
                         Q1 <- as.vector(chol2inv(chol(Matrix::forceSymmetric(BtKB + cw_Ginv_Sigmainv)))) ## THIS IS THE BOTTLENECK
           #                if(G_control$inv_method == "schulz") {
           #                     mat <- forceSymmetric(BtKB + cw_Ginv_Sigmainv)
           #                     Q1 <- as.vector(.schulz_inversion_cmpfn(mat = mat))
           #                     rm(mat)
           #                     }
                         toc <- proc.time()
                         
                    new_G <- matrix(0, nrow = num_spp, ncol = num_spp)
                    if(num_spp > 1)
                         new_G[lower.tri(new_G, diag = TRUE)] <- crossprod(Q2, Q1)
           #                if(num_spp == 1)
           #                     new_G <- matrix(as.vector(crossprod(Q2, Q1)), 1, 1)
                    new_G <- new_G + t(new_G) - diag(x = diag(new_G), nrow = num_spp)
                    new_G <- (new_G + A_Sigmain_AT)/num_basisfns
                    new_G <- Matrix::forceSymmetric(new_G) 
           
                    diff <- 0.5 * mean(as.vector((new_G - cw_G)^2))
                    if(G_control$trace > 0)
                         message("Inner iteration: ", counter, "\t Difference: ", round(diff,5))
                    cw_G <- new_G
                    counter <- counter + 1
                    }
               }
          
          #' Pass the estimation of lambda in G = lambda * R into estimation of Sigma
          if(use_structure == "identity") {
               new_G <- Matrix::Diagonal(n = num_spp)
               }     
          }
          
     if(return_correlation)
          new_G <- stats::cov2cor(as.matrix(new_G))
     
     return(as.matrix(new_G))
     }
     
     
#' @noRd

.update_LoadingG_fn <- function(G, G_control, use_rank_element, correlation) {
     num_spp <- nrow(G)
     num_rank <- G_control$rank[use_rank_element]
     use_structure <- G_control$structure[use_rank_element]
     if(num_rank != "full")
             num_rank <- as.numeric(num_rank)
     
     if(correlation) {
          if(any(G_control$nugget_profile < 0) || any(G_control$nugget_profile > 1))
               stop("All values in nugget_profile should be between 0 and 1.")
          }
     
     if(use_structure == "identity") {
          out <- list(Loading = NULL, nugget = NULL, cov = G)
          out$invcov <- chol2inv(chol(out$cov))
          return(out)
          }
     
     if(num_spp <= 2) {
        out <- list(Loading = NULL, nugget = NULL, cov = G)
        out$invcov <- chol2inv(chol(out$cov))
        return(out)
        }
     
     if(num_rank == "full") {
        out <- list(Loading = NULL, nugget = NULL, cov = G)
        out$invcov <- chol2inv(chol(out$cov))
        return(out)
        }     
     
     if(correlation) {
          min_err <- Inf
          fit_eigen <- eigen(G) 
          for(k0 in 1:length(G_control$nugget_profile)) {     
               Z <- fit_eigen$vectors[,1:num_rank,drop=FALSE] %*% diag(x = fit_eigen$values[1:num_rank]^0.5, nrow = num_rank)
               start_Loading <- Z/matrix(sqrt(rowSums(Z^2)), nrow = num_spp, ncol = num_rank, byrow = FALSE)
               start_Loading <- start_Loading * sqrt(1 - G_control$nugget_profile[k0])
               rm(Z)
               
               counter <- 0
               diff <- 10
               err <- Inf
               new_Loading <- start_Loading

               while(diff > G_control$tol & counter < G_control$maxit) {
                    Rtilde <- G - diag(x = G_control$nugget_profile[k0], nrow = num_spp)
               
                    for(j in 1:num_spp) {
                         Bj <- rowSums(matrix(apply(new_Loading[-j,,drop=FALSE], 1, tcrossprod), nrow = num_rank^2))               
                         Bj <- matrix(Bj, nrow = num_rank)
                         max_evalue_Bj <- eigen(Bj, only.values = TRUE)$values[1]
                         
                         Ej <- max_evalue_Bj * as.vector(new_Loading[j,,drop=FALSE]) - Bj %*% as.vector(new_Loading[j,,drop=FALSE]) + colSums(Rtilde[j,-j] * new_Loading[-j,,drop=FALSE])
                         if(!all(Ej == 0))
                              new_Loading[j,] <- as.vector(Ej / sqrt(sum(Ej^2)) * sqrt(1 - G_control$nugget_profile[k0]))
                         }
               
                    err <- c(err, 0.5 * sum((G - tcrossprod(new_Loading))^2))
                    diff <- err[length(err)-1]/err[length(err)] - 1
                    if(G_control$trace > 0)
                         message("Inner iteration: ", counter, "\t Difference: ", round(diff,5))
                    counter <- counter + 1
                    }
                    
               if(err[length(err)] < min_err) {
                    best_Loading <- new_Loading
                    best_nugget <- G_control$nugget_profile[k0]
                    min_err <- err[length(err)]
                    }
               }
          
          out <- list(Loading = best_Loading,
                      nugget = best_nugget,
                      cov = tcrossprod(best_Loading) + diag(x = best_nugget, nrow = num_spp))
          out$invcov <- chol2inv(chol(out$cov))
          return(out)
          }
     
     if(!correlation) {
          min_err <- Inf
          counter <- 0
          diff <- 10
          err <- Inf
          cw_nugget <- min(diag(G))*0.1
          while(diff > G_control$tol & counter < G_control$maxit) {
               Gtilde <- G - diag(x = cw_nugget, nrow = num_spp)
               
               do_svd <- svd(Gtilde)
               new_approx <- do_svd$u[, 1:num_rank,drop=FALSE] %*% tcrossprod(diag(x = do_svd$d[1:num_rank], nrow = num_rank), do_svd$v[, 1:num_rank,drop=FALSE])
               new_nugget <- mean(diag(G - new_approx))
                         
               err <- c(err, 0.5 * mean((G - new_approx - diag(x = new_nugget, nrow = num_spp))^2))
               diff <- err[length(err)-1]/err[length(err)] - 1
               if(G_control$trace)
                    message("Inner iteration: ", counter, "\t Difference: ", round(diff,5))
               cw_nugget <- new_nugget
               counter <- counter + 1
               }
               rm(do_svd)
          
          do_svd <- svd(new_approx)
          new_Loading <- do_svd$u[, 1:num_rank,drop=FALSE] %*% diag(x = sqrt(do_svd$d[1:num_rank]), nrow = num_rank)
          rm(do_svd)
          
          out <- list(Loading = new_Loading, nugget = new_nugget, cov = tcrossprod(new_Loading) + diag(x = new_nugget, nrow = num_spp))
          out$invcov <- chol2inv(chol(out$cov))
          return(out)
          }
     }
     

#' @noRd

.update_Sigma_fn <- function(Sigmainv, lambdas = NULL, basis_effects_mat, Ginv,
                             B, X, ziX = NULL, zioffset = NULL,
                             y_vec, linpred_vec, dispparam, powerparam, zibetas,
                             trial_size, family, Sigma_control, estimate_lambda, which_B) {

     num_spp <- nrow(basis_effects_mat)
     num_basisfns <- ncol(Sigmainv)
     trial_size <- as.vector(trial_size)

     Sigma_control$method <- match.arg(Sigma_control$method, choices = c("simple","REML","ML")) 
     Sigma_control$inv_method <- "chol2inv" #match.arg(Sigma_control$inv_method, choices = c("chol2inv","schulz")) 
     
     AT_Ginv_A <- crossprod(basis_effects_mat, Ginv) %*% basis_effects_mat
     
     if(family$family == "binomial") {
          linpred_vec[linpred_vec > 5] <- 5
          linpred_vec[linpred_vec < -5] <- -5
          }

          
     if(Sigma_control$method == "simple") { 
          if(estimate_lambda)
               stop("The simple and fast method for estimating Sigma can not be used when the structure of G is set to either \"identity\" or \"homogeneous\".")
          
          new_Sigma <- AT_Ginv_A / num_spp
          return(as.matrix(new_Sigma))
          }

     if(Sigma_control$method %in% c("REML","ML")) { 
          zieta <- NULL
          if(family$family[1] %in% c("zipoisson","zinegative.binomial")) {                        
               zieta <- tcrossprod(ziX, zibetas)
               if(!is.null(zioffset))
                    zieta <- zieta + zioffset
               zieta <- as.vector(zieta)
               }
          ## Note weights are on a per-species basis i.e., site runs faster than species
          num_units <- nrow(B)
          weights_mat <- .neghessfamily(family = family,
                                        eta = linpred_vec,
                                        y = y_vec,
                                        phi = rep(dispparam, each = num_units),
                                        powerparam = rep(powerparam, each = num_units),
                                        zieta = zieta, trial_size = trial_size)

          ## Set up to REML and ML
          weights_mat[is.na(y_vec)] <- 0
          weights_mat <- matrix(weights_mat, nrow = nrow(B), ncol = num_spp, byrow = FALSE)

          if(Sigma_control$method == "REML") {
               inner_fn <- function(j) {
                    w_j <- weights_mat[,j]
                    sqrt_w <- sqrt(w_j)

                    XTX_inv <- chol2inv(chol(crossprod(X*sqrt_w) + Diagonal(x = 1e-8, n = ncol(X))))
                    BTWX <- crossprod(B, X*w_j)
                    return(crossprod(B*sqrt_w) - BTWX %*% tcrossprod(XTX_inv, BTWX))
                    }
               }
          if(Sigma_control$method == "ML") {
               inner_fn <- function(j) {
                    return(crossprod(B*sqrt(weights_mat[,j])))
                    }
               }
                
          if(estimate_lambda == 0) {
               BtKB <- foreach(j = 1:num_spp) %dopar% inner_fn(j = j) 
               BtKB <- Matrix::bdiag(BtKB)
               gc()               

               Q2 <- do.call(rbind, lapply(1:num_spp, function(k) kronecker(Matrix::Diagonal(n = num_basisfns), as(kronecker(Ginv[,k], diag(nrow = num_basisfns)), "sparseMatrix"))) )
               #Q2 <- do.call(rbind, lapply(1:num_spp, function(k) kronecker(Matrix::Diagonal(n = num_basisfns), kronecker(Ginv[,k], Diagonal(n = num_basisfns)))) )
               needonly_cols <- unlist(sapply(1:num_basisfns, function(j) (j-1)*num_basisfns + j:num_basisfns))
               Q2 <- Q2[, needonly_cols, drop = FALSE]
               rm(weights_mat)
               gc()
     
               counter <- 0
               diff <- 10
               cw_Sigma <- chol2inv(chol(Sigmainv))
               while(diff > Sigma_control$tol & counter < Sigma_control$maxit) {
                       cw_Ginv_Sigmainv <- Matrix(kronecker(Ginv, chol2inv(chol(cw_Sigma))), sparse = TRUE)    
     
                       if(Sigma_control$inv_method == "chol2inv")
                            Q1 <- as.vector(chol2inv(chol(Matrix::forceSymmetric(BtKB + cw_Ginv_Sigmainv)))) ## THIS IS THE BOTTLENECK
                              
     #                if(Sigma_control$inv_method == "schulz") {
     #                     mat <- forceSymmetric(BtKB + cw_Ginv_Sigmainv)
     #                     Q1 <- as.vector(.schulz_inversion_cmpfn(mat = mat))
     #                     rm(mat)
     #                     }
                    
                       new_Sigma <- matrix(0, nrow = num_basisfns, ncol = num_basisfns)
                       new_Sigma[lower.tri(new_Sigma, diag = TRUE)] <- crossprod(Q2, Q1)
                       new_Sigma <- new_Sigma + t(new_Sigma) - diag(x = diag(new_Sigma), nrow = num_basisfns)
                       new_Sigma <- (new_Sigma + AT_Ginv_A)/num_spp
                       new_Sigma <- Matrix::forceSymmetric(new_Sigma) 
                    
                       diff <- 0.5 * mean(as.vector((new_Sigma - cw_Sigma)^2))
                       if(Sigma_control$trace > 0)
                            message("Inner iteration: ", counter, "\t Difference: ", round(diff,5))
                       cw_Sigma <- new_Sigma
                       counter <- counter + 1
                       }
               return(as.matrix(new_Sigma))
               }
          
          if(estimate_lambda == 1) {
               BtKB <- foreach(j = 1:num_spp) %dopar% inner_fn(j = j) 
               BtKB <- Matrix::bdiag(BtKB)
               gc()
               
               #' Parametrizing in terms of rho = notLog(1/lambda), where lambda is the variance component. 
               #' 1. When there is a single lambda, then the covariance matrix of the B's is parametrized as (lambda x G) \otimes Sigma = G \otimes (lambda x Sigma)
               #' 2. When there are multiple lambdas, then the covariance matrix of the B's generalizes the form on the RHS of the above as = G \otimes (\sum_r (lambda_r x Sigma_r)^{-1})^{-1} = G \otimes (\sum_r Sigma_r^{-1}/lambda_r)^{-1}. This form is chosen to be consistent with how GAMs in mgcv set up their tensor product penalties.  
               
               if(which_B == 1) {
                    if(any(Sigma_control$upper_space_lambdas > 1e8))
                         Sigma_control$upper_space_lambdas[which(Sigma_control$upper_space_lambdas > 1e8)] <- 1e8
                    if(any(Sigma_control$lower_space_lambdas < 1e-8))
                         Sigma_control$lower_space_lambdas[which(Sigma_control$lower_space_lambdas < 1e-8)] <- 1e-8
                    
                    if(is.matrix(Sigma_control[["custom_space"]])) {
                        Sigmainv_space <- .pinv(Sigma_control[["custom_space"]])
                        trace_quantity <- sum(diag(Sigmainv_space %*% AT_Ginv_A))
                        KronGinvSigmainv <- kronecker(Ginv, Sigmainv_space)

                        fn <- function(x) {
                            expx <- mgcv::notExp(x)
                            out <- 0.5*num_spp*num_basisfns*log(expx) - 0.5*expx*trace_quantity
                            out <- out - 0.5*determinant(BtKB + expx*KronGinvSigmainv)$mod
                            return(as.vector(-out))
                            }
                        
                        update_lambda <- stats::optimise(f = fn, 
                                                         interval = c(mgcv::notLog(1/Sigma_control$upper_space_lambdas), mgcv::notLog(1/Sigma_control$lower_space_lambdas)), 
                                                         maximum = FALSE)
                        new_lambda <- 1/mgcv::notExp(update_lambda$m)
                        # update_lambda <- optim(par = mgcv::notLog(1/lambdas), 
                        #                        fn = fn, 
                        #                        control = list(trace = 0, maxit = 1000), 
                        #                        method = "BFGS")
                        # new_lambda <- 1/mgcv::notExp(update_lambda$par)
                        }
                    if(is.list(Sigma_control[["custom_space"]])) {
                         fn_lambda <- function(x) { 
                              expx <- mgcv::notExp(x)
                              Sigmainv_space <- Reduce("+", lapply(1:length(lambdas), function(j2) .pinv(Sigma_control[["custom_space"]][[j2]]) * expx[j2])) 
                              trace_quantity <- sum(diag(Sigmainv_space %*% AT_Ginv_A))
                              e <- eigen(Sigmainv_space, only.values = TRUE, symmetric = TRUE)$values
                              
                              out <- 0.5*num_spp*sum(log(e[e > .Machine$double.eps])) - 0.5*trace_quantity
                              out <- out - 0.5*determinant(BtKB + kronecker(Ginv, Sigmainv_space))$mod
                              return(as.vector(-out))
                              }
                         update_lambda <- optim(par = mgcv::notLog(1/lambdas), 
                                                fn = fn_lambda, 
                                                control = list(trace = 0, maxit = 1000), 
                                                method = "L-BFGS-B",
                                                lower = mgcv::notLog(1/Sigma_control$upper_space_lambdas), 
                                                upper = mgcv::notLog(1/Sigma_control$lower_space_lambdas))
                         #update_lambda <- nlminb(start = mgcv::notLog(1/lambdas), objective = fn_lambda, control = list(trace = 0, maxit = 1000))
                         new_lambda <- 1/mgcv::notExp(update_lambda$par)
                         }
                    }
               if(which_B == 2) {
                    if(any(Sigma_control$upper_time_lambdas > 1e8))
                         Sigma_control$upper_time_lambdas[which(Sigma_control$upper_time_lambdas > 1e8)] <- 1e8
                    if(any(Sigma_control$lower_time_lambdas < 1e-8))
                         Sigma_control$lower_time_lambdas[which(Sigma_control$lower_time_lambdas < 1e-8)] <- 1e-8
                    
                    if(is.matrix(Sigma_control[["custom_time"]])) {
                         Sigmainv_time <- .pinv(Sigma_control[["custom_time"]])
                         trace_quantity <- sum(diag(Sigmainv_time %*% AT_Ginv_A))
                         KronGinvSigmainv <- kronecker(Ginv, Sigmainv_time)
                         fn <- function(x) { 
                              expx <- mgcv::notExp(x)
                              out <- 0.5*num_spp*num_basisfns*log(expx) - 0.5*expx*trace_quantity
                              out <- out - 0.5*determinant(BtKB + expx*KronGinvSigmainv)$mod
                              return(as.vector(-out))
                              }
                    
                         update_lambda <- stats::optimise(f = fn, 
                                                          interval = c(mgcv::notLog(1/Sigma_control$upper_time_lambdas), mgcv::notLog(1/Sigma_control$lower_time_lambdas)), 
                                                          maximum = FALSE)
                         new_lambda <- 1/mgcv::notExp(update_lambda$m)
                         # update_lambda <- optim(par = mgcv::notLog(1/lambdas), 
                         #                        fn = fn, 
                         #                        control = list(trace = 0, maxit = 1000), 
                         #                        method = "BFGS")
                         # new_lambda <- 1/mgcv::notExp(update_lambda$par)
                         }
                    if(is.list(Sigma_control[["custom_time"]])) {
                         fn_lambda <- function(x) { 
                              expx <- mgcv::notExp(x)
                              Sigmainv_time <- Reduce("+", lapply(1:length(lambdas), function(j2) .pinv(Sigma_control[["custom_time"]][[j2]]) * expx[j2])) 
                              trace_quantity <- sum(diag(Sigmainv_time %*% AT_Ginv_A))
                              e <- eigen(Sigmainv_time, only.values = TRUE, symmetric = TRUE)$values
                              
                              out <- 0.5*num_spp*sum(log(e[e > .Machine$double.eps])) - 0.5*trace_quantity
                              out <- out - 0.5*determinant(BtKB + kronecker(Ginv, Sigmainv_time))$mod
                              return(as.vector(-out))
                              }
                         update_lambda <- optim(par = mgcv::notLog(1/lambdas), 
                                                fn = fn_lambda, 
                                                control = list(trace = 0, maxit = 1000), 
                                                method = "L-BFGS-B",
                                                lower = mgcv::notLog(1/Sigma_control$upper_time_lambdas), 
                                                upper = mgcv::notLog(1/Sigma_control$lower_time_lambdas))
                         new_lambda <- 1/mgcv::notExp(update_lambda$par)
                         }
                    }
               if(which_B == 3) {
                    if(any(Sigma_control$upper_spacetime_lambdas > 1e8))
                         Sigma_control$upper_spacetime_lambdas[which(Sigma_control$upper_spacetime_lambdas > 1e8)] <- 1e8
                    if(any(Sigma_control$lower_spacetime_lambdas < 1e-8))
                         Sigma_control$lower_spacetime_lambdas[which(Sigma_control$lower_spacetime_lambdas < 1e-8)] <- 1e-8
                    
                    if(is.matrix(Sigma_control[["custom_spacetime"]])) {
                         Sigmainv_spacetime <- .pinv(Sigma_control[["custom_spacetime"]])
                         trace_quantity <- sum(diag(Sigmainv_spacetime %*% AT_Ginv_A))
                         KronGinvSigmainv <- kronecker(Ginv, Sigmainv_spacetime)
                         fn <- function(x) { 
                              expx <- mgcv::notExp(x)
                              out <- 0.5*num_spp*num_basisfns*log(expx) - 0.5*expx*trace_quantity
                              out <- out - 0.5*determinant(BtKB + expx*KronGinvSigmainv)$mod
                              return(as.vector(-out))
                              }
                         
                         update_lambda <- stats::optimise(f = fn, 
                                                          interval = c(mgcv::notLog(1/Sigma_control$upper_spacetime_lambdas), mgcv::notLog(1/Sigma_control$lower_spacetime_lambdas)), 
                                                          maximum = FALSE)
                         new_lambda <- 1/mgcv::notExp(update_lambda$m)
                         # update_lambda <- optim(par = mgcv::notLog(1/lambdas), 
                         #                        fn = fn, 
                         #                        control = list(trace = 0, maxit = 1000), 
                         #                        method = "BFGS")
                         # new_lambda <- 1/mgcv::notExp(update_lambda$par)
                         }
                    if(is.list(Sigma_control[["custom_spacetime"]])) {
                         fn_lambda <- function(x) { 
                              expx <- mgcv::notExp(x)
                              Sigmainv_spacetime <- Reduce("+", lapply(1:length(lambdas), function(j2) .pinv(Sigma_control[["custom_spacetime"]][[j2]]) * expx[j2])) 
                              trace_quantity <- sum(diag(Sigmainv_spacetime %*% AT_Ginv_A))
                              e <- eigen(Sigmainv_spacetime, only.values = TRUE, symmetric = TRUE)$values
                              
                              out <- 0.5*num_spp*sum(log(e[e > .Machine$double.eps])) - 0.5*trace_quantity
                              out <- out - 0.5*determinant(BtKB + kronecker(Ginv, Sigmainv_spacetime))$mod
                              return(as.vector(-out))
                              }
                         
                         update_lambda <- optim(par = mgcv::notLog(1/lambdas), 
                                                fn = fn_lambda, 
                                                control = list(trace = 0, maxit = 1000),
                                                method = "L-BFGS-B",
                                                lower = mgcv::notLog(1/Sigma_control$upper_spacetime_lambdas), 
                                                upper = mgcv::notLog(1/Sigma_control$lower_spacetime_lambdas))
                         new_lambda <- 1/mgcv::notExp(update_lambda$par)
                         }
                    }
               
               
               return(as.vector(new_lambda))
               }
          } 
     }
	

#' @noRd
.update_LoadingSigma_fn <- function(Sigma, Sigma_control, use_rank_element, estimate_lambda_not_Sigma, which_B) {
     if(!estimate_lambda_not_Sigma) {
          num_basisfns <- nrow(Sigma)
          num_rank <- Sigma_control$rank[use_rank_element]
     
          if(num_rank == "full") {
             out <- list(Loading = NULL, nugget = NULL, cov = Sigma)
             out$invcov <- chol2inv(chol(out$cov))
             return(out)
             }
          
          min_err <- Inf
          counter <- 0
          diff <- 10
          err <- Inf
          cw_nugget <- min(diag(Sigma))*0.1
          num_rank <- as.numeric(num_rank)
          while(diff > Sigma_control$tol & counter < Sigma_control$maxit) {
               Sigmatilde <- Sigma - diag(x = cw_nugget, nrow = num_basisfns)
               
               do_svd <- svd(Sigmatilde)
               new_approx <- do_svd$u[, 1:num_rank,drop=FALSE] %*% tcrossprod(diag(x = do_svd$d[1:num_rank], nrow = num_rank), do_svd$v[, 1:num_rank,drop=FALSE])
               new_nugget <- mean(diag(Sigma - new_approx))
                         
               err <- c(err, 0.5 * mean((Sigma - new_approx)^2))
               diff <- err[length(err)-1]/err[length(err)] - 1
               if(Sigma_control$trace)
                    message("Inner iteration: ", counter, "\t Difference: ", round(diff,5))
               cw_nugget <- new_nugget
               counter <- counter + 1
               }
               
          rm(do_svd)
          do_svd <- svd(new_approx)
          new_Loading <- do_svd$u[, 1:num_rank,drop=FALSE] %*% diag(x = sqrt(do_svd$d[1:num_rank]), nrow = num_rank)
          rm(do_svd)
     
          out <- list(Loading = new_Loading, nugget = new_nugget, cov = tcrossprod(new_Loading) + diag(x = new_nugget, nrow = num_basisfns))
          out$invcov <- chol2inv(chol(out$cov))
          }
     
     #' After lambda has been estimated in the updating Sigma part of the algorithm, Sigma itself is adjusted based on this i.e., kronecker(lambda x R, Sigma) = kronecker(R, lambda x Sigma)
     if(estimate_lambda_not_Sigma) {
          if(which_B == 1) {
               if(is.matrix(Sigma_control[["custom_space"]]))
                    out <- list(Loading = NULL, nugget = NULL, lambdas = Sigma,
                                cov = Sigma_control[["custom_space"]] * Sigma, invcov = .pinv(Sigma_control[["custom_space"]]) / Sigma)
               if(is.list(Sigma_control[["custom_space"]]))
               out <- list(Loading = NULL, nugget = NULL, lambdas = Sigma,
                           invcov = Reduce("+", lapply(1:length(Sigma), function(x) .pinv(Sigma_control[["custom_space"]][[x]]) / Sigma[x])) )
               out$cov <- .pinv(out$invcov)
               }
          if(which_B == 2) {
               if(is.matrix(Sigma_control[["custom_time"]]))
                    out <- list(Loading = NULL, nugget = NULL, lambdas = Sigma,
                                cov = Sigma_control[["custom_time"]] * Sigma, invcov = .pinv(Sigma_control[["custom_time"]]) / Sigma)
               if(is.list(Sigma_control[["custom_time"]]))
                    out <- list(Loading = NULL, nugget = NULL, lambdas = Sigma,
                                invcov = Reduce("+", lapply(1:length(Sigma), function(x) .pinv(Sigma_control[["custom_time"]][[x]]) / Sigma[x])) )
               out$cov <- .pinv(out$invcov)
               }
          if(which_B == 3) {
               if(is.matrix(Sigma_control[["custom_spacetime"]]))
                    out <- list(Loading = NULL, nugget = NULL, lambdas = Sigma,
                                cov = Sigma_control[["custom_spacetime"]] * Sigma, invcov = .pinv(Sigma_control[["custom_spacetime"]]) / Sigma)
               if(is.list(Sigma_control[["custom_spacetime"]]))
                    out <- list(Loading = NULL, nugget = NULL, lambdas = Sigma,
                                invcov = Reduce("+", lapply(1:length(Sigma), function(x) .pinv(Sigma_control[["custom_spacetime"]][[x]]) / Sigma[x])) )
                    out$cov <- .pinv(out$invcov)
                    }
          }
     
     return(out)
     }
