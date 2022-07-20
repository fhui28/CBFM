update_G_fn <- function(Ginv, basis_effects_mat, Sigmainv, B, X, ziX = NULL, y_vec, linpred_vec, dispparam, powerparam, zibetas,  
                        trial_size, family, G_control, return_correlation = TRUE) {
        
        num_spp <- nrow(basis_effects_mat)
        num_basisfns <- ncol(Sigmainv)
        trial_size <- as.vector(trial_size)

        # if(num_spp == 1) 
        #         return(solve(Ginv))
     
        G_control$method <- match.arg(G_control$method, choices = c("simple","REML","ML")) 
        G_control$inv_method <- "chol2inv" #match.arg(G_control$inv_method, choices = c("chol2inv","schulz"))
        A_Sigmain_AT <- basis_effects_mat %*% tcrossprod(Sigmainv, basis_effects_mat)

        if(family$family == "binomial") {
                linpred_vec[linpred_vec > 5] <- 5
                linpred_vec[linpred_vec < -5] <- -5
                }
     
        if(G_control$method == "simple") {
                new_G <- A_Sigmain_AT / num_basisfns
                }
     
        if(G_control$method %in% c("REML","ML")) {
             zieta <- NULL
             if(family$family[1] %in% c("zipoisson","zinegative.binomial")) {                        
                  zieta <- as.vector(tcrossprod(ziX, zibetas))
                  }
             ## Note weights are on a per-species basis i.e., site runs faster than species
             weights_mat <- .neghessfamily(family = family, eta = linpred_vec, y = y_vec, phi = rep(dispparam, each = nrow(B)), 
                                           powerparam = rep(powerparam, each = nrow(B)), zieta = zieta, trial_size = trial_size)

             ## Set up REML and ML
             weights_mat[is.na(y_vec)] <- 0
             weights_mat <- matrix(weights_mat, nrow = nrow(B), ncol = num_spp, byrow = FALSE)
             if(G_control$method == "REML") {
                  inner_fn <- function(j) {
                       XTX_inv <- chol2inv(chol(crossprod(X*sqrt(weights_mat[,j]))))
                       BTWX <- crossprod(B, X*weights_mat[,j])               
                       return(crossprod(B*sqrt(weights_mat[,j])) - BTWX %*% tcrossprod(XTX_inv, BTWX))
                    }
                  }
             if(G_control$method == "ML") {
                  inner_fn <- function(j) {
                       return(crossprod(B*sqrt(weights_mat[,j])))
                    }
                  }
                
             BWB_minus_BWX_XWXinv_XWB <- foreach(j = 1:num_spp) %dopar% inner_fn(j = j) 
             BWB_minus_BWX_XWXinv_XWB <- Matrix::bdiag(BWB_minus_BWX_XWXinv_XWB)
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
                       Q1 <- as.vector(chol2inv(chol(Matrix::forceSymmetric(BWB_minus_BWX_XWXinv_XWB + cw_Ginv_Sigmainv)))) ## THIS IS THE BOTTLENECK

#                if(G_control$inv_method == "schulz") {
#                     mat <- forceSymmetric(BWB_minus_BWX_XWXinv_XWB + cw_Ginv_Sigmainv)
#                     Q1 <- as.vector(.schulz_inversion_cmpfn(mat = mat))
#                     rm(mat)
#                     }

                  new_G <- matrix(0, nrow = num_spp, ncol = num_spp)
                  new_G[lower.tri(new_G, diag = TRUE)] <- crossprod(Q2, Q1)

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

        if(return_correlation)
             new_G <- stats::cov2cor(as.matrix(new_G))
          
        return(as.matrix(new_G))
        }

## SEPARARE TESTING CODE TO ESTIMATE THE SPECIAL CASE OF A DIAGONAL G. DO DIRECT NUMERICAL OPTIMIZATION SINCE IT IS FASTER
# objfn <- function(x) {
#      M <- matrix(0, num_spp, num_spp)
#      diag(M) <- exp(x)
#      out <- 0.5*num_basisfns*determinant(M, logarithm = TRUE)$mod - 0.5*sum(diag(M %*% basis_effects_mat %*% tcrossprod(Sigmainv, basis_effects_mat)))
#      out <- out - 0.5*determinant(BWB_minus_BWX_XWXinv_XWB + kronecker(M, Sigmainv), logarithm = TRUE)$mod
#      -out
#      }
# 
# check_optim <- optim(par = numeric(num_spp), fn = objfn, control = list(trace = 1, maxit = 1000), method = "BFGS")
# 
# new_Ginv <- matrix(0, num_spp, num_spp)
# diag(new_Ginv) <- exp(check_optim$par)
# new_G <- chol2inv(chol(new_Ginv))
     
     
update_LoadingG_fn <- function(G, G_control, use_rank_element, correlation = TRUE) {
     num_spp <- nrow(G)
     num_rank <- G_control$rank[use_rank_element]
     if(num_rank != "full")
             num_rank <- as.numeric(num_rank)
     
     if(correlation) {
          if(any(G_control$nugget_profile < 0) || any(G_control$nugget_profile > 1))
               stop("All values in nugget_profile should be between 0 and 1.")
          }
     if(num_spp <= 2) {
        out <- list(Loading = NULL, nugget = NULL, cov = G)
        out$covinv <- chol2inv(chol(out$cov))
        return(out)
        }
     
     if(num_rank == "full") {
        out <- list(Loading = NULL, nugget = NULL, cov = G)
        out$covinv <- chol2inv(chol(out$cov))
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
          
          out <- list(Loading = best_Loading, nugget = best_nugget, cov = tcrossprod(best_Loading) + diag(x = best_nugget, nrow = num_spp))
          out$covinv <- chol2inv(chol(out$cov))
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
                         
               err <- c(err, 0.5 * mean((G - new_approx)^2))
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
          out$covinv <- chol2inv(chol(out$cov))
          return(out)
          }
     }
     

update_Sigma_fn <- function(Sigmainv, basis_effects_mat, Ginv, B, X, ziX = NULL, y_vec, linpred_vec, dispparam, powerparam, zibetas, 
                            trial_size, family, Sigma_control) {
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
                new_Sigma <- AT_Ginv_A / num_spp
                }

        if(Sigma_control$method %in% c("REML","ML")) {
             zieta <- NULL
             if(family$family[1] %in% c("zipoisson","zinegative.binomial")) {                        
                  zieta <- as.vector(tcrossprod(ziX, zibetas))
                  }
             ## Note weights are on a per-species basis i.e., site runs faster than species
                weights_mat <- .neghessfamily(family = family, eta = linpred_vec, y = y_vec, phi = rep(dispparam, each = nrow(B)), 
                                              powerparam = rep(powerparam, each = nrow(B)),  zieta = zieta, trial_size = trial_size)

                ## Set up to REML and ML
                weights_mat[is.na(y_vec)] <- 0
                weights_mat <- matrix(weights_mat, nrow = nrow(B), ncol = num_spp, byrow = FALSE)
                if(Sigma_control$method == "REML") {
                     inner_fn <- function(j) {
                          XTX_inv <- chol2inv(chol(crossprod(X*sqrt(weights_mat[,j]))))
                          BTWX <- crossprod(B, X*weights_mat[,j])               
                          return(crossprod(B*sqrt(weights_mat[,j])) - BTWX %*% tcrossprod(XTX_inv, BTWX))
                         }
                    }
                if(Sigma_control$method == "ML") {
                     inner_fn <- function(j) {
                          return(crossprod(B*sqrt(weights_mat[,j])))
                          }
                    }
                
               
                BWB_minus_BWX_XWXinv_XWB <- foreach(j = 1:num_spp) %dopar% inner_fn(j = j) 
                BWB_minus_BWX_XWXinv_XWB <- Matrix::bdiag(BWB_minus_BWX_XWXinv_XWB)
                gc()               

                Q2 <- do.call(rbind, lapply(1:num_spp, function(k) kronecker(Matrix::Diagonal(n = num_basisfns), kronecker(Ginv[,k], Matrix::Diagonal(n = num_basisfns)))) )
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
                                Q1 <- as.vector(chol2inv(chol(Matrix::forceSymmetric(BWB_minus_BWX_XWXinv_XWB + cw_Ginv_Sigmainv)))) ## THIS IS THE BOTTLENECK
                         
#                if(Sigma_control$inv_method == "schulz") {
#                     mat <- forceSymmetric(BWB_minus_BWX_XWXinv_XWB + cw_Ginv_Sigmainv)
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
                }
        
        return(as.matrix(new_Sigma))
        }
	

update_LoadingSigma_fn <- function(Sigma, Sigma_control, use_rank_element) {
     num_basisfns <- nrow(Sigma)
     num_rank <- Sigma_control$rank[use_rank_element]
     
     if(num_rank == "full") {
        out <- list(Loading = NULL, nugget = NULL, cov = Sigma)
        out$covinv <- chol2inv(chol(out$cov))
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
     out$covinv <- chol2inv(chol(out$cov))
     return(out)
     }
     