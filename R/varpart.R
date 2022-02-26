#' @title Variance partitioning for a CBFM fit
#' 
#' @description 
#' `r lifecycle::badge("stable")`
#' 
#' For each species, partition the variance of the linear predictor from a fitted \code{CBFM} object into components associated with (groups of) the covariates, and the basis function included.
#' 
#' @param object An object of class \code{CBFM}.
#' @param groupX A vector of group indicator variables, which allows the variance partitioning to be done for groups of covariates (including the intercept) i.e., calculating how much of the total variation does a certain subset of the covariates explain. The length of \code{groupX} must be equal to \code{ncol(object$betas)}. Defaults to \code{NULL}, in whih case all the covariates are treated as a single group.
#' 
#' @details 
#' Variance partitioning can be performed for a communiy-level basis function model (CBFM) fit in a similar manner to how it is performed for latent variables models e.g., by [boral::calc.varpart()]. For the purposes of the package, the CBFM is characterized by the following mean regression model: for observational unit \eqn{i=1,\ldots,N} and species \eqn{j=1,\ldots,m}, we have
#' 
#' \deqn{g(\mu_{ij}) = \eta_{ij} = x_i^\top\beta_j + b_i^\top a_j,}
#'
#' where \eqn{g(.)} is a known link function, \eqn{x_i^\top\beta_j} is the component of the linear predictor due to the explained covariates (however they end up being included in the CBFM fit), and \eqn{b_i^\top a_j} is the component due to the inclusion of basis functions, noting that that may include spatial, temporal, and/or spatio-temporal basis functions. Note that for zero-inflated distributions, this \eqn{\mu_{ij}} is the mean of the non-zero-inflated component. 
#'     
#' For each species, variation partitioning is performed by calculating the variance due to each component in \eqn{\eta_{ij}}, and then rescaling them to ensure that they sum to one. The general details of this type of variation partitioning is given in Ovaskainen et al., (2017) and Bjork et al. (2018) among others for latent variables models. Note that the variance due the basis functions is itself broken up into components depending on what basis functions are included in the model. For example, if a fitted \code{CBFM} object included basis functions for \code{B_space} and \code{B_space} (but not \code{B_spacetime}), then the variance is partitioned into separate components for the spatial and temporal basis functions. 
#'
#' If \code{groupX} is supplied, the variance due to the included covariates is done based on subsets of covariates (including the intercept) as identified by \code{groupX}, and then rescaled correspondingly. This is useful if one was to, for example, quantify the proportion of variation in each species which is explained by each specific environmental covariate and there are multiple terms associated with each covariate e.g., due to polynomial or smoothing terms.  
#'
#' Finally, note that if there are missing values in the response matrix, then the variance partitioning calculation will ignore all values of \eqn{\eta_{ij}} corresponding to this.
#'
#' @return A list with the following components (if applicable):
#' \describe{
#' \item{varpart_X}{Vector containing the proportion of variance for each species explained by measure predictors, and grouped according to \code{goupX} if supplied.}

#' \item{varpart_B_space}{Vector containing the proportion of variance for each species explained by spatial basis functions.}

#' \item{varpart_B_time}{Vector containing the proportion of variance for each species explained by temporal basis functions.}

#' \item{varpart_B_spacetime}{Vector containing the proportion of variance for each species explained by spatio-temporal basis functions.}
#' }
#' 
#' @section Warning:
#' There is some controversy over what quantities proportion of variance explained exactly mean in the context of latent variable models, and this also carries over to the case of CBFMs. Please interpret the numbers with a grain of salt e.g., there may not be a consistent let alone sensible relationship between proportion of variance and estimated parameters in the CBFM fit and whether they are statistically significantly or not.
#'
#'
#' @author Francis K.C. Hui <fhui28@gmail.com>, Chris Haak
#' 
#' @references 
#' Bjork, J. R., Hui, F. K. C., O'Hara, R. B., and Montoya, J. M. (2018). Uncovering the drivers of host-associated microbiota with joint species distribution modelling. Molecular Ecology, 27, 2714-2724.
#' 
#' Ovaskainen, O., Tikhonov, G., Norberg, A., Guillaume Blanchet, F., Duan, L., Dunson, D., and Abrego, N. (2017). How to make more out of community data? A conceptual framework and its implementation as models and software. Ecology Letters, 20, 561-576.
#' 
#' @seealso [CBFM()] for fitting CBFMs.
#' 
#' @examples
#' \donttest{
#' library(autoFRK)
#' library(FRK)
#' library(MASS)
#' library(mvabund)
#' library(mvtnorm)
#' library(ROCR)
#' library(sp)
#' library(RandomFields)
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
#' X <- rmvnorm(num_sites, mean = rep(0,4)) 
#' colnames(X) <- c("temp", "depth", "chla", "O2")
#' dat <- data.frame(xy, X)
#' mm <- model.matrix(~ temp + depth + chla + O2 - 1, data = dat) %>% 
#' scale %>% 
#' as.matrix
#' 
#' # Simulate latent variable component
#' true_lvs <- RFsimulate(model = RMexp(var=1, scale=2), 
#' x = xy$x, y = xy$y, n = 2)@data %>% 
#' as.matrix
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
#' fitcbfm <- CBFM(y = simy, formula_X = useformula, data = dat, 
#' B_space = basisfunctions, family = binomial(), control = list(trace = 1))
#' 
#' varpart(fitcbfm)
#' 
#' # Intercept, temp and depth considered as one group; chla and O2 in the other group
#' varpart(fitcbfm, groupX = c(1,1,1,2,2))  
#' }
#' 
#' @export
#' 
#' @importFrom mgcv gam
#' @importFrom stats cor cov var
#' @md

varpart <- function(object, groupX = NULL) {
    if(!inherits(object, "CBFM")) 
        stop("`object' is not of class \"CBFM\"")

    num_spp <- nrow(object$betas)
     
     #tmp_formula <- as.formula(paste("response", paste(as.character(object$formula_X),collapse="") ) )
     #nullfit <- gam(tmp_formula, data = data.frame(response = runif(nrow(object$y)), object$data), fit = TRUE, control = list(maxit = 1))
     X <- model.matrix.CBFM(object)
     #rm(tmp_formula, nullfit)
     rownames(X) <- rownames(object$data)
    
     if(!is.null(groupX)) { 
          if(length(groupX) != ncol(object$betas)) 
               stop("If groupX is supplied, then it must be a vector with the length equal to ncol(object$betas).") 
        }

    #if(!is.null(object$traits)) {
    #    trait.linpred <- array(0, dim = c(nrow(fit.mcmc), object$n, num_spp))
    #    all_cors_spp <- matrix(0, nrow(fit.mcmc), num_spp)
    #    }

    X_var <- tcrossprod(X, object$betas)
    X_var[is.na(object$y)] <- NA
    X_var <- apply(tcrossprod(X, object$betas), 2, var, na.rm = TRUE)
    if(!is.null(groupX)) {
          X_var <- matrix(0, nrow = length(unique(groupX)), ncol = num_spp)        
          for(k0 in 1:length(unique(groupX))) {
              cw_X_var <- tcrossprod(X[,which(groupX==k0),drop=FALSE], object$betas[,which(groupX==k0),drop=FALSE])
              cw_X_var[is.na(object$y)] <- NA
              X_var[k0,] <- apply(cw_X_var, 2, var, na.rm = TRUE)
               }
          }
          
     B_space_var <- B_time_var <- B_spacetime_var <- rep(0, num_spp)
     if(object$which_B_used[1]) {
        cw_B_var <- tcrossprod(object$B[,1:object$num_B_space,drop=FALSE], object$basis_effects_mat[,1:object$num_B_space,drop=FALSE])
        cw_B_var[is.na(object$y)] <- NA
        B_space_var <- apply(cw_B_var, 2, var, na.rm = TRUE)
        }
     if(object$which_B_used[2]) {
        cw_B_var <- tcrossprod(object$B[,object$num_B_space + (1:object$num_B_time),drop=FALSE], object$basis_effects_mat[,object$num_B_space + (1:object$num_B_time),drop=FALSE])
        cw_B_var[is.na(object$y)] <- NA
         B_time_var <- apply(cw_B_var, 2, var, na.rm = TRUE)
        }
     if(object$which_B_used[3]) {
        cw_B_var <- tcrossprod(object$B[,object$num_B_space + object$num_B_time + (1:object$num_B_spacetime),drop=FALSE], object$basis_effects_mat[,object$num_B_space + object$num_B_time + (1:object$num_B_spacetime),drop=FALSE])
        cw_B_var[is.na(object$y)] <- NA
        B_spacetime_var <- apply(cw_B_var, 2, var, na.rm = TRUE)
        }
     #if(!is.null(object$traits)) {
     #     cw.traits.coefs <- cbind(fit.mcmc[k, grep("traits.int",colnames(fit.mcmc))], matrix(fit.mcmc[k, grep("traits.coefs",colnames(fit.mcmc))], nrow = ncol(object$X)+1))
     #     rownames(cw.traits.coefs) <- c("beta0", colnames(object$X))
     #     trait.X.coefs <- tcrossprod(cbind(1,object$traits), cw.traits.coefs) ## beta = intercept + trait %*% trait.coefs
     #     cw.trait.linpred <- tcrossprod(cbind(1,object$X), trait.X.coefs)
     #     all_cors_spp[k,] <-  sapply(1:num_spp, function(i) cor(cw.X.linpred[,i], cw.trait.linpred[,i])^2)
     #     }
    
     B_total_var <- B_space_var + B_time_var + B_spacetime_var
     if(is.null(groupX)) {
          var_X <- X_var / (X_var + B_total_var) 
          var_Bspace <- B_space_var / (X_var + B_total_var)
          var_Btime <- B_time_var / (X_var + B_total_var)
          var_Bspacetime <- B_spacetime_var / (X_var + B_total_var)
          names(var_X) <- names(var_Bspace) <- names(var_Btime) <- names(var_Bspacetime) <- rownames(object$betas)
          }

     if(!is.null(groupX)) {
          var_X <- matrix(0, nrow = length(unique(groupX)), ncol = num_spp) 
          for(k0 in 1:length(unique(groupX))) 
               var_X[k0,] <- X_var[k0,] / (colSums(X_var) + B_total_var)
          rownames(var_X) <- paste0("Group", unique(groupX))
          
          var_Bspace <- B_space_var / (colSums(X_var) + B_total_var)
          var_Btime <- B_time_var / (colSums(X_var) + B_total_var)
          var_Bspacetime <- B_spacetime_var / (colSums(X_var) + B_total_var)
          colnames(var_X) <- names(var_Bspace) <- names(var_Btime) <- names(var_Bspacetime) <- rownames(object$betas)
          }

          
    out <- list(varpart_X = var_X)
     if(object$which_B_used[1])
          out$varpart_B_space <- var_Bspace
     if(object$which_B_used[2])
          out$varpart_B_time <- var_Btime
     if(object$which_B_used[3])
          out$varpart_B_spacetime <- var_Bspacetime
    #if(!is.null(object$traits)) {
    #    out$R2.traits <- colMeans(all_cors_spp)
    #    names(all_cors_spp) <- colnames(object$y) 
    #    }
    return(out)
    }


