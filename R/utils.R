## Calculate -0.5*tr(G^{-1} * A * Sigma^{-1} * t(A)), which is the quadratic term in the PQL function of the CBFM related to the basis functions coefficients
.calc_pqlquadraticterm_basiseffects <- function(basis_effects_mat, Ginv, Sigmainv) {
     out <- Ginv %*% basis_effects_mat %*% tcrossprod(Sigmainv, basis_effects_mat)
     return(-0.5*sum(diag(out)))
     }
          
     
## Two functions used to extract relevant principle submatrices of the Bayesian posterior covariance matrix from the CBFM fit. 
## This was initially done to save memory and because only certain components of the matrices are needed later on, but has since been abandoned since some predictions need to simulate from the whole kitchen sink!
.extractcovarblocks_topright <- function(j, Q, num_X, num_basisfns) { 
     return(Q[(num_X*j - num_X + 1):(num_X*j), (num_basisfns*j - num_basisfns + 1):(num_basisfns*j)]) 
     }

.extractcovarblocks_bottomright <- function(j, Q, num_basisfns) { 
     return(Q[(num_basisfns*j - num_basisfns + 1):(num_basisfns*j), (num_basisfns*j - num_basisfns + 1):(num_basisfns*j)]) 
     }

               
## E-step functions for zero-inflated and zero-truncated distributions -- calculate the posterior probability of being in the zero-inflation component/posterior probability of observing a zero
.estep_fn <- function(family, cwfit, y, X, B) {
   num_units <- nrow(y)
   num_spp <- ncol(y)
   out <- Matrix::Matrix(0, nrow = num_units, ncol = num_spp, sparse = TRUE)
   if(family$family %in% c("zipoisson","zinegative.binomial")) {
                fitvals <- exp(tcrossprod(X, cwfit$betas) + tcrossprod(B, cwfit$basis_effects_mat))
                zeroinfl_prob <- plogis(cwfit$zeroinfl_prob_intercept)
                
                for(j in 1:num_spp) {
                        sel_zerospp <- which(y[,j] == 0)
                        if(family$family[1] == "zipoisson")
                                out[sel_zerospp,j] <- zeroinfl_prob[j] / (zeroinfl_prob[j] + (1-zeroinfl_prob[j])*dpois(0, lambda = fitvals[sel_zerospp,j]))
                        if(family$family[1] == "zinegative.binomial")
                                out[sel_zerospp,j] <- zeroinfl_prob[j] / (zeroinfl_prob[j] + (1-zeroinfl_prob[j])*dnbinom(0, mu = fitvals[sel_zerospp,j], size = 1/cwfit$dispparam[j]))
                        }
                }
        
   # if(family$family %in% c("ztnegative.binomial")) {
   #    fitvals <- exp(tcrossprod(X, cwfit$betas) + tcrossprod(B, cwfit$basis_effects_mat))
   #    for(j in 1:num_spp) {
   #       out[,j] <- dnbinom(0, mu = fitvals[,j], size = 1/cwfit$dispparam[j]) / (1-dnbinom(0, mu = fitvals[,j], size = 1/cwfit$dispparam[j]))
   #       }
   #    }
        
   return(out)
   }



## Get the full S matrix from GAMs. Relies on the fact that gam always move the parametric terms first
.get_bigS <- function(fit_gam, num_X) {
   # if(class(fit_gam)[1] == "gamlss")
   #    fit_gam <- getSmo(fit_gam)
     
   bigS <- Matrix::Matrix(0, num_X, num_X, sparse = TRUE)
   num_smooth_terms <- length(fit_gam$smooth)
   if(num_smooth_terms == 0)
      return(bigS)
          
   num_Smatrices_per_smooth <- sapply(fit_gam$smooth, function(x) length(x$S)) # The sum of this should equal length(fit_gam$sp)
   sp_index <- split(1:length(fit_gam$sp), rep(1:num_smooth_terms, num_Smatrices_per_smooth)) # Because te and ti smooths have multiple S and smoothing parameters, then this tells you how many S/sp's are within each smooth term  
   rm(num_Smatrices_per_smooth)
   num_smooth_cols <- sum(sapply(fit_gam$smooth, function(x) x$df)) # According to ?smooth.construct, this is the degrees of freedom associated with this term when unpenalized and unconstrained
   num_parametric_cols <- num_X - num_smooth_cols

   subS <- lapply(1:num_smooth_terms, function(j) {
      out <- fit_gam$sp[sp_index[[j]][1]] * fit_gam$smooth[[j]]$S[[1]]
      if(length(sp_index[[j]]) > 1) { # To deal with smooths that have multiple S matrices and smoothing parameters
         for(l0 in 2:length(sp_index[[j]]))
            out <- out + fit_gam$sp[sp_index[[j]][l0]] * fit_gam$smooth[[j]]$S[[l0]]
         }
      return(out)
      })
   subS <- Matrix::bdiag(subS)
   bigS[-(1:num_parametric_cols), -(1:num_parametric_cols)] <- subS
          
   return(bigS)
   }
   
   
## This function is used specifically when an additive form of the CBFM is used, for use in the construction of the Bayesian posterior covariance matrix for standard errors.
## For example, given G_space, Sigma_space, G_time, Sigma_time, it forms for covariance matrix for vector of the random slopes (a_{space,1},a_{time,1}, a_{space,2},a_{time,2},...a_{space,m},a_{time,m}). 
.kkproduct <- function(G1, G2, G3 = NULL, Sigma1, Sigma2, Sigma3 = NULL, inverse = TRUE) {
     G1 <- as.matrix(G1)     
     G2 <- as.matrix(G2)               
     Sigma1 <- as.matrix(Sigma1)     
     Sigma2 <- as.matrix(Sigma2)     
     num_spp <- nrow(G1)
     num_basisfns_Sigma1 <- nrow(Sigma1)
     num_basisfns_Sigma2 <- nrow(Sigma2)
          
     if(is.null(G3)) {     
          GSigma1 <- kronecker(G1, Sigma1)
          GSigma2 <- kronecker(G2, Sigma2)
          num_basisfns <- num_basisfns_Sigma1+num_basisfns_Sigma2
          
          out <- matrix(0, nrow = num_spp*num_basisfns, ncol = num_spp*num_basisfns)
          sel_spaceseq <- (0:(num_spp-1))*num_basisfns
          sel_spaceseq <- as.vector(sapply(sel_spaceseq, function(x) x + 1:num_basisfns_Sigma1))
          out[sel_spaceseq, sel_spaceseq] <- GSigma1
          out[-sel_spaceseq, -sel_spaceseq] <- GSigma2
          }
          
     if(!is.null(G3)) {     
          G3 <- as.matrix(G3)     
          Sigma3 <- as.matrix(Sigma3)     
          num_basisfns_Sigma3 <- nrow(Sigma3)
          num_basisfns <- num_basisfns_Sigma1+num_basisfns_Sigma2+num_basisfns_Sigma3
          
          out <- matrix(0, nrow = num_spp*num_basisfns, ncol = num_spp*num_basisfns) # Do not use sparse matrices here as it is much slower making a sparse-matrix non-sparse!
          for(j in 1:num_spp) { for(k in 1:j) {
               sel_rows <- num_basisfns*j - num_basisfns + 1:num_basisfns
               sel_cols <- num_basisfns*k - num_basisfns + 1:num_basisfns
               
               out[sel_rows, sel_cols] <- as.matrix(bdiag(G1[j,k]*Sigma1, G2[j,k]*Sigma2, G3[j,k]*Sigma3))
               } }                    
          
          out <- 0.5*(out + t(out))
          }
     
     out <- Matrix::Matrix(out, sparse = TRUE)
     if(!inverse)
          return(out)
     if(inverse)
          return(.cholthenpinv(out))
     }


## A local pseudo-inverse function -- straight from summary.gam in mgcv package. Full credit goes to Simon Wood for this!
.pinv <- function(V, M, rank.tol = 1e-6) {
    D <- eigen(V, symmetric = TRUE)
    M1 <- length(D$values[D$values > rank.tol * D$values[1]])
    if(M > M1)
        M<-M1 # avoid problems with zero eigen-values
    if(M+1 <= length(D$values))
        D$values[(M+1):length(D$values)]<-1
    D$values<- 1/D$values
    if(M+1 <= length(D$values))
        D$values[(M+1):length(D$values)]<-0
    res <- D$vectors %*% diag(x = D$values, nrow = length(D$values)) %*% t(D$vectors)

    return(res)
    }


# Tries chol2inv then, which will fail for singular matrices. If it does fail then go to use .pinv. 
# One could use .pinv directly, but I suspect chol2inv is more scalable (when it works) and so is preferred as a first option?
.cholthenpinv <- function(V) {
   out <- try(chol2inv(chol(V)), silent = TRUE)
   if(inherits(out, "try-error"))
      out <- .pinv(V = V, M = ncol(V))
   
   return(out)
   }

# schulz_inversion_fn <- function(mat, max_iter = 100) {
#      diff <- 10
#      counter <- 0
#      
#      cwV <-  Diagonal(n = nrow(mat)) / norm(mat, "F") ## Initial guess
# 
#      while(diff > 1e-3 & counter < max_iter) {
#           ## Hotelling-Bodewig algorithm 
#           #newV <- cwV %*% (Diagonal(x = 2, n = num_spp*num_basisfns) - mat %*% cwV)
#           ## Li et al., 2011
#           I_minus_AV <- as.matrix(Diagonal(n = num_spp*num_basisfns) - mat %*% cwV)
#           newV <- cwV %*% (Diagonal(n = num_spp*num_basisfns) + I_minus_AV + I_minus_AV %^% 2)
#           
#           diff <- norm(newV - cwV, "F")
#           message("Schulz Inversion: Iteration ", counter, "\t Diff:", round(diff,6))
#           cwV <- newV
#           counter <- counter + 1
#           }
#      return(newV)
#      }
# 
# .schulz_inversion_cmpfn <- compiler::cmpfun(schulz_inversion_fn)
# rm(schulz_inversion_fn)
# 
