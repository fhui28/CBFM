## Calculate -0.5*tr(G^{-1} * A * Sigma^{-1} * t(A)), which is the quadratic term in the PQL function of the CBFM related to the basis functions coefficients
.calc_pqlquadraticterm_basiseffects <- function(basis_effects_mat, Ginv, Sigmainv) {
     out <- Ginv %*% basis_effects_mat %*% tcrossprod(Sigmainv, basis_effects_mat)
     return(-0.5*sum(diag(out)))
     }
          

## Function to trick mgcv and subsequently gratia so the right standard errors are obtained, along with everything else, when applying gratia::parametric_effects
.calc_parametric_effects <- function(j, object) {
     tmp_formula <- as.formula(paste("response", paste(as.character(object$formula),collapse = " ") ) )
     nulldat <- data.frame(response = object$y[,j], object$data)
     nullfit <- mgcv::gam(tmp_formula, data = nulldat, knots = object$knots, fit = TRUE, control = list(maxit = 1))
     num_X <- ncol(model.matrix(nullfit))
     
     if(object$family$family[1] %in% c("zipoisson","zinegative.binomial")) {
          tmp_ziformula <- as.formula(paste("response", paste(as.character(object$ziformula),collapse = " ") ) )
          zinullfit <- mgcv::gam(tmp_ziformula, data = nulldat, knots = object$ziknots, fit = TRUE, control = list(maxit = 1))
          num_ziX <- ncol(model.matrix(zinullfit))
          }
     
     if(length(attr(nullfit$pterms, "term.labels")) == 0) 
          out <- NULL
     if(length(attr(nullfit$pterms, "term.labels")) > 0) {
          nullfit$coefficients <- object$betas[j,]
          sel_rowcols <- grep(paste0(colnames(object$y)[j],"$"), rownames(object$covar_components$topleft))
          if(object$family$family[1] %in% c("zipoisson","zinegative.binomial")) {
               sel_rowcols <- sel_rowcols[-(1:num_ziX)]
               }
          nullfit$Vp <- as.matrix(object$covar_components$topleft[sel_rowcols, sel_rowcols, drop = FALSE])
          out <- suppressMessages(parametric_effects(object = nullfit))
          out$species <- colnames(object$y)[j]
          }
     
     # Parametric effects for zero-inflation is appropriate
     ziout <- NULL
     if(object$family$family[1] %in% c("zipoisson","zinegative.binomial")) {
          if(length(attr(zinullfit$pterms, "term.labels")) == 0) 
               ziout <- NULL
          if(length(attr(zinullfit$pterms, "term.labels")) > 0) {
               zinullfit$coefficients <- object$zibetas[j,]
               sel_rowcols <- grep(paste0(colnames(object$y)[j],"$"), rownames(object$covar_components$topleft))
               sel_rowcols <- sel_rowcols[(1:num_ziX)]
               nullfit$Vp <- as.matrix(object$covar_components$topleft[sel_rowcols, sel_rowcols, drop = FALSE])
               ziout <- suppressMessages(parametric_effects(object = zinullfit))
               ziout$species <- colnames(object$y)[j]
               }
          }
     
     return(list(out = out, ziout = ziout))
     }


## Function to trick mgcv and subsequently gratia so the right standard errors are obtained, along with everything else, when applying gratia::smooth_estimates
.calc_smooth_estimates <- function(j, object) {
     tmp_formula <- as.formula(paste("response", paste(as.character(object$formula),collapse = " ") ) )
     nulldat <- data.frame(response = object$y[,j], object$data)
     nullfit <- mgcv::gam(tmp_formula, data = nulldat, knots = object$knots, fit = TRUE, control = list(maxit = 1))
     num_X <- ncol(model.matrix(nullfit))
     
     if(object$family$family[1] %in% c("zipoisson","zinegative.binomial")) {
          tmp_ziformula <- as.formula(paste("response", paste(as.character(object$ziformula),collapse = " ") ) )
          zinullfit <- mgcv::gam(tmp_ziformula, data = nulldat, knots = object$ziknots, fit = TRUE, control = list(maxit = 1))
          num_ziX <- ncol(model.matrix(zinullfit))
          }
     
     if(length(nullfit$smooth) == 0) 
          out <- NULL
     if(length(nullfit$smooth) > 0) {
          nullfit$coefficients <- object$betas[j,]
          sel_rowcols <- grep(paste0(colnames(object$y)[j],"$"), rownames(object$covar_components$topleft))
          if(object$family$family[1] %in% c("zipoisson","zinegative.binomial")) {
               sel_rowcols <- sel_rowcols[-(1:num_ziX)]
               }
          nullfit$Vp <- nullfit$Ve <- nullfit$Vc <- as.matrix(object$covar_components$topleft[sel_rowcols, sel_rowcols,drop=FALSE])
          out <- suppressMessages(smooth_estimates(object = nullfit))
          out$species <- colnames(object$y)[j]
          }
     
     # Smooth estimates for zero-inflation is appropriate
     ziout <- NULL
     if(object$family$family[1] %in% c("zipoisson","zinegative.binomial")) {
          if(length(zinullfit$smooth) == 0) 
               ziout <- NULL
          if(length(zinullfit$smooth) > 0) {
               zinullfit$coefficients <- object$zibetas[j,]
               sel_rowcols <- grep(paste0(colnames(object$y)[j],"$"), rownames(object$covar_components$topleft))
               sel_rowcols <- sel_rowcols[(1:num_ziX)]
               zinullfit$Vp <- zinullfit$Ve <- zinullfit$Vc <- as.matrix(object$covar_components$topleft[sel_rowcols, sel_rowcols,drop=FALSE])
               ziout <- suppressMessages(smooth_estimates(object = zinullfit))
               ziout$species <- colnames(object$y)[j]
               }
          }
     
     return(list(out = out, ziout = ziout))
     }


## Two functions used to extract relevant principle submatrices of the Bayesian posterior covariance matrix from the CBFM fit. 
## This was initially done to save memory and because only certain components of the matrices are needed later on, but has since been abandoned since some predictions need to simulate from the whole kitchen sink!
.extractcovarblocks_topright <- function(j, Q, num_X, num_basisfns) { 
     return(Q[(num_X*j - num_X + 1):(num_X*j), (num_basisfns*j - num_basisfns + 1):(num_basisfns*j)]) 
     }

.extractcovarblocks_bottomright <- function(j, Q, num_basisfns) { 
     return(Q[(num_basisfns*j - num_basisfns + 1):(num_basisfns*j), (num_basisfns*j - num_basisfns + 1):(num_basisfns*j)]) 
     }

               
## E-step functions for zero-inflated distributions -- calculate the posterior probability of being in the zero-inflation component/posterior probability of observing a zero
.estep_fn <- function(family,
                      cwfit,
                      y,
                      X,
                      B,
                      offset,
                      formula_offset,
                      ziX = NULL,
                      zioffset = NULL) {

   num_units <- nrow(y)
   num_spp <- ncol(y)
   out <- Matrix::Matrix(0, nrow = num_units, ncol = num_spp, sparse = TRUE)
   
   if(family$family %in% c("zipoisson","zinegative.binomial")) {
                fitvals <- exp(tcrossprod(X, cwfit$betas) + tcrossprod(B, cwfit$basis_effects_mat) + offset + formula_offset)
                zieta <- tcrossprod(ziX, cwfit$zibetas)
                if(!is.null(zioffset))
                    zieta <- zieta + zioffset
                zeroinfl_prob <- plogis(zieta)
                rm(zieta)
                
                for(j in 1:num_spp) {
                        sel_zerospp <- which(y[,j] == 0)
                        if(family$family[1] == "zipoisson")
                                out[sel_zerospp,j] <- zeroinfl_prob[sel_zerospp,j] / (zeroinfl_prob[sel_zerospp,j] + (1-zeroinfl_prob[sel_zerospp,j])*dpois(0, lambda = fitvals[sel_zerospp,j]))
                        if(family$family[1] == "zinegative.binomial")
                                out[sel_zerospp,j] <- zeroinfl_prob[sel_zerospp,j] / (zeroinfl_prob[sel_zerospp,j] + (1-zeroinfl_prob[sel_zerospp,j])*dnbinom(0, mu = fitvals[sel_zerospp,j], size = 1/cwfit$dispparam[j]))
                        }
                
                }
        
   return(out)
   }


## Get the full S matrix from GAMs. Relies on the fact gam always move the parametric terms first
.get_bigS <- function(fit_gam, num_X) {
   # if(class(fit_gam)[1] == "gamlss")
   #    fit_gam <- getSmo(fit_gam)
     
   bigS <- Matrix::Matrix(0, num_X, num_X, sparse = TRUE)
   num_smooth_terms <- length(fit_gam$smooth)
   if(num_smooth_terms == 0)
      return(bigS)
          
   num_Smatrices_per_smooth <- lapply(fit_gam$smooth, function(x) length(x$S)) # The sum of this should equal length(fit_gam$sp)
   sp_index <- split(1:length(fit_gam$sp), rep(1:num_smooth_terms, num_Smatrices_per_smooth)) # Because fs, te, and ti smooths have multiple S and smoothing parameters, then this tells you how many and indexes the S/sp's within each smooth term. This is very similar to extracting first.sp and last.sp from each smooth
   rm(num_Smatrices_per_smooth)
   #num_smooth_cols <- sum(sapply(fit_gam$smooth, function(x) x$df)) # According to ?smooth.construct, this is the degrees of freedom associated with this term when unpenalized and unconstrained
   num_smooth_cols <- sum(sapply(fit_gam$smooth, function(x) x$last.para - x$first.para + 1)) 
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
   
  
## Function for calculating some starting lambda values, which are then used into .update_Sigma_fn() when multiple lambda's need to be estimated.
## Based on and acknowledgments go to the initial.sp function in the mgcv package.
## Not currently not actually used!
# .get_initial_lambdas <- function(BtKB, custom_spacetime) {
#      starting_lambdainv <- array(0, length(custom_spacetime))
#      ldxx <- diag(BtKB)
#      ldss <- ldxx * 0
#      pen <- rep(FALSE, length(ldxx))
#      
#      S <- lapply(custom_spacetime, .pinv)
#      
#      for(j in 1:length(custom_spacetime)) {
#           rsS <- rowMeans(abs(S[[j]]))
#           csS <- colMeans(abs(S[[j]]))
#           dS <- diag(abs(S[[j]]))
#           thresh <- .Machine$double.eps^0.8 * max(abs(S[[j]]))
#           ind <- rsS > thresh & csS > thresh & dS > thresh
#           ss <- diag(S[[j]])[ind]
#           xx <- ldxx
#           xx <- xx[ind]
#           pen <- pen | ind
#           sizeXX <- mean(xx)
#           sizeS <- mean(ss)
#           starting_lambdainv[j] <- sizeXX / sizeS
#           ldss <- ldss + starting_lambdainv[j] * diag(S[[j]])
#           }
#           
#      ind <- ldss > 0 & pen & ldxx > 0
#      ldxx <- ldxx[ind]
#      ldss <- ldss[ind]
#      while(mean(ldxx/(ldxx + ldss)) > 0.4) {
#           starting_lambdainv <- starting_lambdainv * 10
#           ldss <- ldss * 10
#           }
#      while (mean(ldxx/(ldxx + ldss)) < 0.4) {
#           starting_lambdainv <- starting_lambdainv/10
#           ldss <- ldss/10
#           }
#      
#      return(1/starting_lambdainv)
#      }


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
          return(.pinv(out))
     }


## A local pseudo-inverse function -- straight from summary.gam in mgcv package. Full credit goes to Simon Wood for this!
.pinv <- function(V, M, rank.tol = 1e-6) {
     if(missing(M))
          M <- ncol(V)
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



#' ##------------------------------
#' ## **Example 1h: Repeat Example 1a but illustrate applications to zero-truncated count data** 
#' ## This one could take a while...grab a cuppa while you wait!
#' ##------------------------------
#' library(gamlss)
#' library(gamlss.tr)
#' gen.trun(0, family = "NBI")
#'
#' set.seed(2021)
#' num_sites <- 1000 # 500 (units) sites for training set + 500 sites for testing.
#' num_spp <- 50 # Number of species
#' num_X <- 4 # Number of regression slopes
#' 
#' spp_slopes <- matrix(runif(num_spp * num_X, -1, 1), nrow = num_spp)
#' spp_intercepts <- runif(num_spp, -1, 1)
#' 
#' # Simulate spatial coordinates and environmental covariate components
#' # We will use this information in later examples as well
#' xy <- data.frame(x = runif(num_sites, 0, 5), y = runif(num_sites, 0, 5))
#' X <- mvtnorm::rmvnorm(num_sites, mean = rep(0,4)) 
#' colnames(X) <- c("temp", "depth", "chla", "O2")
#' dat <- data.frame(xy, X)
#' mm <- model.matrix(~ temp + depth + chla + O2 - 1, data = dat) %>% 
#' scale %>% 
#' as.matrix
#' 
#' # Simulate latent variable component
#' # We will use this information in later examples as well
#' true_lvs <- grf(grid = cbind(xy$x, xy$y), nsim = 2, cov.model = "exponential", 
#' cov.pars = c(1, 2))$data %>% 
#'      as.matrix
#' spp_loadings <- matrix(runif(num_spp * 2, -1, 1), nrow = num_spp) 
#' set.seed(NULL)
#' 
#' # Simulate spatial multivariate abundance data
#' # Note the double loop is needed as rNBItr behaves oddly when using vectorized arguments
#' eta <- tcrossprod(cbind(1,mm), cbind(spp_intercepts,spp_slopes)) + 
#' tcrossprod(true_lvs, spp_loadings)
#' spp_dispersion <- runif(num_spp)
#' simy <- eta*0 
#' for(i in 1:num_sites) {
#'      for(j in 1:num_spp) {
#'           simy[i,j] <- rNBItr(1, mu = exp(eta[i,j]), sigma = spp_dispersion[j])
#'      } }
#' 
#' # Form training and test sets
#' dat_train <- dat[1:500,]
#' dat_test <- dat[501:1000,]
#' simy_train <- simy[1:500,]
#' simy_test <- simy[501:1000,]
#' 
#' 
#' # Fit stacked zero-truncated regression models as a baseline
#' # Note gamlss may sometimes fail...
#' fitstacked <- NULL 
#' for(j in 1:num_spp) {
#' fitstacked[[j]] <- gamlss(resp ~ temp + depth + chla + O2, 
#' data = data.frame(resp = simy_train[,j], dat_train), family = NBItr)
#' }
#' 
#' 
#' # Set up spatial basis functions for CBFM -- Most users will start here! 
#' # This is the same set up as examples above
#' num_basisfunctions <- 25 # Number of spatial basis functions to use
#' # Training set basis functions
#' train_basisfunctions <- mrts(dat_train[,c("x","y")], num_basisfunctions) %>% 
#' as.matrix %>%
#' {.[,-(1)]} # Remove the first intercept column
#' # Testing set basis functions
#' test_basisfunctions <- mrts(dat_train[,c("x","y")], num_basisfunctions) %>% 
#' predict(newx = dat_test[,c("x","y")]) %>% 
#' as.matrix %>%
#' {.[,-c(1)]} 
#' 
#' # Fit zero-truncated negative binomial CBFM
#' tic <- proc.time()
#' useformula <- ~ temp + depth + chla + O2
#' fitcbfm <- CBFM(y = simy_train, formula = useformula, data = dat_train, 
#' B_space = train_basisfunctions, family = ztnb2(), control = list(trace = 1))
#' toc <- proc.time()
#' toc - tic
#' 
#' summary(fitcbfm) %>% 
#' str
#' 
#' 
#' # Calculate predictions onto test dataset
#' # Note extra step needed for stacked models from gamlss to get the actual fitted values 
#' predictions_stacked <- sapply(1:num_spp, function(j) predict(fitstacked[[j]], 
#' newdata = dat_test, type = "response"))
#' stacked_dispparam <- sapply(1:num_spp, function(j) exp(fitstacked[[j]]$sigma.coefficient))
#' stacked_dispparam <- matrix(stacked_dispparam, nrow(dat_test), num_spp, byrow = TRUE)
#' predictions_stacked <- predictions_stacked / 
#' (1-dnbinom(0, mu = predictions_stacked, size = 1/stacked_dispparam))
#' predictions_cbfm <- predict(fitcbfm, newdata = dat_test, type = "response", 
#' new_B_space = test_basisfunctions)
#' 
#' # Evaluation predictions
#' # Pseudo R-squared across species
#' pseudoR2 <- data.frame(
#' stacked = sapply(1:num_spp, function(j) { 
#' out <- cor(predictions_stacked[,j], simy_test[,j], method = "spearman")
#' out^2 * sign(out)     
#' }),
#' cbfm = sapply(1:num_spp, function(j) { 
#' out <- cor(predictions_cbfm[,j], simy_test[,j], method = "spearman")
#' out^2 * sign(out)     
#' })
#' )
#' 
#' boxplot(pseudoR2, main = "Pseudo-R2", names = c("Stacked GLM", "CBFM"))
#' 
#' ggplot(pseudoR2, aes(x = stacked, y = cbfm)) +
#' geom_point() +
#' geom_abline(intercept = 0, slope = 1, linetype = 2) +
#' labs(x = "Stacked SDM", y = "CBFM", main = "Pseudo-R2") +
#' theme_bw()
#' 
#' # Predictive deviance across species (lower is better)
#' # Need to define density of zero-inflated NB distribution first (or get it from a package)
#' preddeviance <- data.frame(
#' stacked = sapply(1:num_spp, function(j) { 
#' -2*sum(dNBItr(simy_test[,j], mu = predictions_stacked[,j], 
#' sigma = exp(fitstacked[[j]]$sigma.coefficient), log = TRUE))
#' }),
#' cbfm = sapply(1:num_spp, function(j) { 
#' -2*sum(dNBItr(simy_test[,j], mu = predictions_cbfm[,j], 
#' sigma = fitcbfm$dispparam[j], log = TRUE))
#' })
#' )
#' 
#' boxplot(preddeviance, main = "Deviance", names = c("Stacked GLM", "CBFM"))
#' 
#' ggplot(preddeviance, aes(x = stacked, y = cbfm)) +
#' geom_point() +
#' geom_abline(intercept = 0, slope = 1, linetype = 2) +
#' labs(x = "Stacked SDM", y = "CBFM", main = "Deviance") +
#' theme_bw()
#' 
#' 
#' ## Please see the makeahurdle function for uses of the above for hurdle count models
#' 
