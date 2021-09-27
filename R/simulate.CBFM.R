#' @title Simulate data from CBFM fit
#' 
#' @description 
#' `r lifecycle::badge("experimental")`
#'
#' Simulate new spatio-temporal multivariate abundance data based on a fitted \code{CBFM} object.
#'
#' @param object An object of class "CBFM".
#' @param nsim A positive integer specifying the number of simulated datasets. Defaults to 1.
#' @param seed An integer to set seed number. Defaults to a random seed number.
#' @param max_resp A upper bound to limit the maximum value of responses obtained. This is useful if the user wants, say, all counts to not exceed a particular value. In such case, the function will attempt to simulate counts that do not \code{max_resp}. Note it only \emph{attempts} this: it will give up after 10 unsuccessful attempts and then return whatever is simulated on the 10-th attempt.
#' @param conditional If \code{conditional = TRUE}, the data are simulated conditional on the estimated species-specific regression coefficients associated with the basis functions. Otherwise if \code{conditional = FALSE} then new species-specific regression coefficients are generated from the estimated values of the \eqn{\Sigma}'s and \eqn{G}'s, and their corresponding random effects distribution. Please see the details section in [CBFM()] for more details.
#' With CBFM being set up much a like generalized additive model or GAM, then simulating conditionally is generally what most practitioners will required...if they need to simulate.
#' @param ... not used.
#'
#' @details 
#' Simulates spatio-temporal multivariate abundance data from a fitted community-level basis function model (CBFM). For the purposes of the package, the CBFM is characterized by the following mean regression model: for observational unit \eqn{i=1,\ldots,N} and species \eqn{j=1,\ldots,m}, we have
#' 
#' \deqn{g(\mu_{ij}) = \eta_{ij} = x_i^\top\beta_j + b_i^\top a_j,}
#'
#' where \eqn{g(.)} is a known link function, \eqn{x_i} denotes a vector of predictors for unit i i.e., the i-th row from the created model matrix, \eqn{\beta_j} denotes the corresponding regression coefficients for species j, \eqn{b_i} denotes a vector of spatial, temporal, and/or spatio-temporal basis functions for unit i , and \eqn{a_j} denotes the corresponding regression coefficients for species j. In the function, \eqn{x_i} is created based on the \code{formula_X} and \code{data} arguments, \eqn{\beta_j} is supplied as part of the code{betas} argument, and \eqn{b_i} is formed from the \code{B_space}, \code{B_time} and \code{B_spacetime} arguments. Finally, \eqn{a_j} is either supplied directly as part of \code{basis_effects_mat} argument, or generated based on the \code{Sigma} and \code{G} arguments. 
#' 
#' As an example, suppose we have a CBFM which involves spatial and temporal (but no spatio-temporal) basis functions. Then \eqn{b_i = (b_{i,space}, b_{i,time})} is formed from the i-th rows of \code{B_space} and \code{B_time}, while \eqn{a_j = (a_{j,space}, a_{j,time})} comes from the j-th row \code{basis_effects_mat}. If \code{basis_effects_mat} is not supplied, then it is instead obtain by simulating 
#' 
#' \deqn{(a_{1,space}, \ldots, a_{m,space}) \sim N(0, kronecker(G_{space}, \Sigma_{space})),} 
#' 
#' where \eqn{G_{space}} and \eqn{\Sigma_{space}} are supplied from \code{G$space} and \code{Sigma$space} respectively, and \eqn{kronecker(\cdot)} is the Kronecker product operator. Similarly, we have \eqn{(a_{1,time}, \ldots, a_{m,time}) \sim N(0, kronecker(G_{time}, \Sigma_{time}))}. 
#' 
#' Based on the mean model given above, responses \eqn{y_{ij}} are then simulated from the assumed distribution, using the additional dispersion/power/zero inflation probability parameters as appropriate. In this function, all parameter are replaced by their estimated values from the fitted CBFM, as appropriate.
#' 
#' @return A three dimensional array of dimension \eqn{N} by \eqn{m} by \code{nsim} is returned, where the number of simulated spatio-temporal multivariate abundance data sets is given by the last index.
#' 
#' @author Francis K.C. Hui <fhui28@gmail.com>, Chris Haak
#' 
#' @seealso [create_CBFM_life()] for manually simulating spatio-temporal multivariate abundance data.
#'
#' @examples
#'  \donttest{
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
#' # Set up spatial basis functions for CBFM -- Most practitioners will start here! 
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
#' # Conditional on estimated basis function coefficients
#' simulate(fitcbfm, nseed = 5)
#' 
#' # Generate new basis function coefficients as part of the simulation process
#' simulate(fitcbfm, nseed = 1, conditional = FALSE)
#' }
#'
#'
#' @export
#' @importFrom stats plogis
#' @md
#' 
simulate.CBFM <- function (object, nsim = 1, seed = NULL, max_resp = Inf, conditional = TRUE, ...) {
     ## Code chunk from simulate.lm to sort out the seed thing
     if(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) 
          runif(1)
     if(is.null(seed)) 
          RNGstate <- get(".Random.seed", envir = .GlobalEnv)
     else {
          R.seed <- get(".Random.seed", envir = .GlobalEnv)
          set.seed(seed)
          RNGstate <- structure(seed, kind = as.list(RNGkind()))
          on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
          }
  
    if(!inherits(object, "CBFM")) 
        stop("`object' is not of class \"CBFM\"")
        

     ##-----------------------
     ## Set up true model
     ##-----------------------
     true_model <- list(
          family = object$family, 
          formula_X = object$formula_X,
          data = object$data, 
          B_space = NULL, 
          B_time = NULL, 
          B_spacetime = NULL, 
          offset = object$offset,
          betas = object$betas,
          basis_effects_mat = NULL,
          Sigma = NULL,
          G = NULL,
          trial_size = object$trial_size,
          dispparam = object$dispparam, 
          powerparam = object$powerparam, 
          max_resp = max_resp
          )
     if(!is.null(object$zeroinfl_prob_intercept))
       true_model$zeroinfl_prob <- plogis(object$zeroinfl_prob_intercept)
     if(object$which_B_used[1] == 1)
          true_model$B_space <- object$B[,1:object$num_B_space,drop=FALSE]
     if(object$which_B_used[2] == 1)
          true_model$B_time <- object$B[,object$num_B_space + (1:object$num_B_time),drop=FALSE]
     if(object$which_B_used[3] == 1)
          true_model$B_spacetime <- object$B[,object$num_B_space + object$num_B_time + (1:object$num_B_spacetime),drop=FALSE]
        
     if(conditional == TRUE) {
          true_model$basis_effects_mat <- object$basis_effects_mat
          }
     
     if(conditional == FALSE) {
          true_model$Sigma <- list(space = object$Sigma_space, time = object$Sigma_time, spacetime = object$Sigma_spacetime)
          true_model$G <- list(space = object$G_space, time = object$G_time, spacetime = object$G_spacetime)
          }
     
     ##-------------------
     ## Simulate response matrices
     ##-------------------
     out <- replicate(nsim, create_CBFM_life(
               family = true_model$family, 
               formula_X = true_model$formula_X, 
               data = true_model$data, 
               B_space = true_model$B_space, 
               B_time = true_model$B_time, 
               B_spacetime = true_model$B_spacetime, 
               offset = true_model$offset,
               betas = true_model$betas,
               basis_effects_mat = true_model$basis_effects_mat,
               Sigma = true_model$Sigma,
               G = true_model$G,
               trial_size = true_model$trial_size,
               dispparam = true_model$dispparam, 
               powerparam = true_model$dispparam, 
               zeroinfl_prob = true_model$zeroinfl_prob, 
               max_resp = max_resp,
               only_y = TRUE)
          )

     return(out)
     }

#' @rdname simulate.CBFM
#' @export 
simulate <- function(object, ...) {
    UseMethod(generic = "simulate")
    }
