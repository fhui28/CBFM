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

               
## E-step functions for zero-inflated distributions -- calculate the posterior probability of being in the zero-inflation component
## Hidden and not exported
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
        
        return(out)
        }



## Get the full S matrix from GAMs. Relies on the fact that gam always move the parametric terms first
.get_bigS <- function(fit_gam, num_X) {
   # if(class(fit_gam)[1] == "gamlss")
   #    fit_gam <- getSmo(fit_gam)
     
   bigS <- Matrix::Matrix(0, num_X, num_X, sparse = TRUE)
   num_smooth_terms <- length(fit_gam$sp)
   if(num_smooth_terms == 0)
      return(bigS)
          
   num_smooth_cols <- sum(sapply(fit_gam$smooth, function(x) ncol(x$S[[1]])))
   num_parametric_cols <- num_X - num_smooth_cols

   subS <- lapply(1:num_smooth_terms, function(j) {
      fit_gam$sp[j] * fit_gam$smooth[[j]]$S[[1]]
      })
   subS <- Matrix::bdiag(subS)
   bigS[-(1:num_parametric_cols), -(1:num_parametric_cols)] <- subS
          
   return(bigS)
   }
   
   
## This function is used specifically when an additive form of the CBFM is used, for use in the construction of the Bayesian posterior covariance matrix for standard errors.
## For example, given G_space, Sigma_space, G_time, Sigma_time, it forms for covariance matrix for vecttor of the random slopes (a_{space,1},a_{time,1}, a_{space,2},a_{time,2},...a_{space,m},a_{time,m}). 
.kkproduct <- function(G1, G2, G3 = NULL, Sigma1, Sigma2, Sigma3 = NULL, inverse = TRUE) {
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
          num_basisfns_Sigma3 <- nrow(Sigma3)
          num_basisfns <- num_basisfns_Sigma1+num_basisfns_Sigma2+num_basisfns_Sigma3
          
          out <- matrix(0, nrow = num_spp*num_basisfns, ncol = num_spp*num_basisfns)
          for(j in 1:num_spp) { for(k in 1:j) {
               sel_rows <- num_basisfns*j - num_basisfns + 1:num_basisfns
               sel_cols <- num_basisfns*k - num_basisfns + 1:num_basisfns
               
               out[sel_rows, sel_cols] <- bdiag(G1[j,k]*Sigma1, G2[j,k]*Sigma2, G3[j,k]*Sigma3)
               } }                    
          
          out <- 0.5*(out + t(out))
          }
     
     out <- Matrix::Matrix(out, sparse = TRUE)
     if(!inverse)
          return(out)
     if(inverse)
          return(chol2inv(chol(out)))
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
