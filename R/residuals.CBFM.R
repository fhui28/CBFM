#' @title Extract residuals from a CBFM fit
#' 
#' @description 
#' `r lifecycle::badge("experimental")`
#' 
#' Calculate various types of residuals from a fitted \code{CBFM} object, including probability integral transform (PIT) residuals and Dunn-Smyth residuals.
#' 
#' @param object An object of class "CBFM".
#' @param type The type of residuals which should be returned. Currently the options available are: "repose" (default), "pearson", "PIT", "dunnsmyth", and "partial". Can be abbreviated.
#' @param seed This can be used set the seed when constructing the PIT and Dunn-Smyth residuals, which for discrete responses involve some degree of jittering.  
#' @param ... Not used.
#' 
#' @details 
#' Suppose that the fitted values from a CBFM fit are denoted by \eqn{\hat{\mu}_{ij}} for observational unit \eqn{i = 1,\ldots,N} and species \eqn{j=1,\ldots,m}. Then:
#' 
#' For \code{type = "response"}, this returns the raw residuals \eqn{y_{ij} - \hat{\mu}_{ij}}. Note for binomial responses what is returned is \eqn{y_{ij}/N_{trial,ij} - \hat{\mu}_{ij}} where \eqn{N_{trial,ij}} is the corresponding trial size and \eqn{\hat{\mu}_{ij}} is the fitted probability of "success". 
#' For zero-inflated distributions, this returns  \eqn{y_{ij} - (1-\hat{\pi}_j)\hat{\mu}_{ij}} where \eqn{\hat{\pi}_j} is the estimated species-specific probability of zero-inflation and \eqn{\hat{\mu}_{ij}} is the estimated mean of the non-zero-inflated component.
#' 
#' For \code{type = "pearson"}, this returns the Pearson residuals, which are calculated by standardizing the raw residuals by the square root of their corresponding variance. Note the variance incorporates any scaled parameters if present, and so sometimes this type of residual is referred to as scaled Pearson residuals.
#' 
#' For \code{type = "PIT"}, this returns the probability integral transform residuals that are also used in [DHARMa::simulateResiduals()] and [mpcmp::rPIT()], among other packages. If the (estimated) model is correct, then these residuals should behave as random variables from a standard uniform distribution (Dunn and Smyth, 1996). Note there is a level of jitting used in producing the PIT residuals.
#' 
#' For \code{type = "dunnsmyth"}, this returns the Dunn-Smyth residuals that are also used in [boral::ds.residuals()] and[gllvm::residuals.gllvm()], among other packages. If the (estimated) model is correct, then these residuals should behave as random variables from a standard normal distribution (Dunn and Smyth, 1996). Note there is a level of jitting used in producing the Dunn-Smyth residuals.
#' 
#' @return A matrix of residuals.
#' 
#' @author Francis K.C. Hui <fhui28@gmail.com>, Chris Haak
#' 
#' @references 
#' Dunn, P. K., and Smyth, G. K. (1996). Randomized quantile residuals. Journal of Computational and Graphical Statistics, 5, 236-244.
#' 
#' @seealso [CBFM()] for fitting CBFMs and [fitted.values.CBFM()] for calculating fitted values from a CBFM fit.
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
#' residuals(fitcbfm)
#' 
#' residuals(fitcbfm, type = "dunnsmyth")
#' }
#' 
#' @export
#' 
#' @importFrom gamlss.tr trun.p
#' @importFrom stats runif qnorm pbeta pbinom pgamma plogis pnorm ppois pnbinom
#' @importFrom tweedie ptweedie 
#' @md

residuals.CBFM <- function(object, type = "response", seed = NULL, ...) {
        type <- match.arg(type, choices = c("response", "pearson", "dunnsmyth", "PIT"))
        num_units <- nrow(object$y)
        num_spp <- ncol(object$y)

     
        out <- object$y - object$fitted
        if(object$family$family[1] == "binomial")
                out <- object$y/object$trial_size - object$fitted
     
        if(type == "response")
                out <- out
     
        if(type == "pearson") {
                if(object$family$family[1] %in% c("binomial")) 
                        out <- out / sqrt(object$family$variance(object$fitted)/object$trial_size)     
                if(object$family$family[1] %in% c("gaussian")) 
                        out <- out / sqrt(matrix(object$dispparam, nrow = num_units, ncol = num_spp, byrow = TRUE))
                if(object$family$family[1] %in% c("Gamma")) 
                        out <- out / sqrt(matrix(object$dispparam, nrow = num_units, ncol = num_spp, byrow = TRUE) * object$family$variance(object$fitted))
                if(object$family$family[1] %in% c("negative.binomial", "ztnegative.binomial","Beta")) 
                        out <- out / sqrt(object$family$variance(object$fitted, phi = matrix(object$dispparam, nrow = num_units, ncol = num_spp, byrow = TRUE)))
                if(object$family$family[1] %in% c("poisson", "ztpoisson")) 
                        out <- out / sqrt(object$family$variance(object$fitted))     
                if(object$family$family[1] %in% c("tweedie")) 
                        out <- out / sqrt(object$family$variance(object$fitted, phi = matrix(object$dispparam, nrow = num_units, ncol = num_spp, byrow = TRUE),
                                                                 power = matrix(object$powerparam, nrow = num_units, ncol = num_spp, byrow = TRUE)))
                if(object$family$family[1] %in% c("zipoisson")) 
                        out <- out / sqrt(object$family$actual_variance(object$fitted, zeroinfl_prob = matrix(plogis(object$zeroinfl_prob_intercept), nrow = num_units, ncol = num_spp, byrow = TRUE)))
                if(object$family$family[1] %in% c("zinegative.binomial")) 
                        out <- out / sqrt(object$family$actual_variance(object$fitted, zeroinfl_prob = matrix(plogis(object$zeroinfl_prob_intercept), nrow = num_units, ncol = num_spp, byrow = TRUE), phi = matrix(object$dispparam, nrow = num_units, ncol = num_spp, byrow = TRUE))
                            )
                if(object$family$family[1] %in% c("ztpoisson")) 
                        out <- out / sqrt(object$family$actual_variance(mu = object$fitted, lambda = exp(object$linear_predictor)))
          }
          
        if(type %in% c("PIT","dunnsmyth")) {
                set.seed(seed)
                if(object$family$family[1] %in% c("beta")) 
                        out <- pbeta(object$y, shape1 = matrix(object$dispparam, nrow = num_units, ncol = num_spp, byrow = TRUE) * object$fitted, 
                                     shape2 = matrix(object$dispparam, nrow = num_units, ncol = num_spp, byrow = TRUE) * (1-object$fitted))
                if(object$family$family[1] %in% c("binomial")) 
                        out <- runif(length(object$y), 
                                     min = pbinom(object$y-1, size = object$trial_size, prob = object$fitted), 
                                     max = pbinom(object$y, size = object$trial_size, prob = object$fitted)
                                     )
                if(object$family$family[1] %in% c("Gamma"))
                        out <- pgamma(object$y, scale = matrix(object$dispparam, nrow = num_units, ncol = num_spp, byrow = TRUE) * object$fitted, 
                                      shape = matrix(1/object$dispparam, nrow = num_units, ncol = num_spp, byrow = TRUE))
                if(object$family$family[1] %in% c("gaussian")) 
                        out <- pnorm(object$y, mean = object$fitted, sd = sqrt(matrix(object$dispparam, nrow = num_units, ncol = num_spp, byrow = TRUE)))
                if(object$family$family[1] %in% c("poisson"))
                        out <- runif(length(object$y), min = ppois(object$y-1, lambda = object$fitted), max = ppois(object$y, lambda = object$fitted))
                if(object$family$family[1] %in% c("negative.binomial")) 
                        out <- runif(length(object$y), 
                                     min = pnbinom(object$y-1, mu = object$fitted, size = matrix(1/object$dispparam, nrow = num_units, ncol = num_spp, byrow = TRUE)),
                                     max = pnbinom(object$y, mu = object$fitted, size = matrix(1/object$dispparam, nrow = num_units, ncol = num_spp, byrow = TRUE))
                                     )
                if(object$family$family[1] %in% c("tweedie")) {
                        a <- b <- ptweedie(object$y, mu = object$fitted, phi = matrix(object$dispparam, nrow = num_units, ncol = num_spp, byrow = TRUE),
                                           power = matrix(object$powerparam, nrow = num_units, ncol = num_spp, byrow = TRUE))
                        a[object$y == 0] <- 0
                        out <- runif(length(object$y), min = a, max = b)
                        }
                if(object$family$family[1] %in% c("zipoisson")) {
                        out <- runif(length(object$y), 
                                     min = .pzipois(object$y-1, lambda = exp(object$linear_predictor), zeroinfl_prob = matrix(plogis(object$zeroinfl_prob_intercept), nrow = num_units, ncol = num_spp, byrow = TRUE)),
                                     max = .pzipois(object$y, lambda = exp(object$linear_predictor), zeroinfl_prob = matrix(plogis(object$zeroinfl_prob_intercept), nrow = num_units, ncol = num_spp, byrow = TRUE))
                                     )
                        }
                if(object$family$family[1] %in% c("zinegative.binomial")) {
                        out <- runif(length(object$y), 
                                     min = .pzinegativebinomial(object$y-1, lambda = exp(object$linear_predictor), zeroinfl_prob = matrix(plogis(object$zeroinfl_prob_intercept), nrow = num_units, ncol = num_spp, byrow = TRUE), phi = matrix(object$dispparam, nrow = num_units, ncol = num_spp, byrow = TRUE)),
                                     max = .pzinegativebinomial(object$y, lambda = exp(object$linear_predictor), zeroinfl_prob = matrix(plogis(object$zeroinfl_prob_intercept), nrow = num_units, ncol = num_spp, byrow = TRUE), phi = matrix(object$dispparam, nrow = num_units, ncol = num_spp, byrow = TRUE))
                                     )
               }
          if(object$family$family[1] %in% c("ztpoisson")) {
                  tmp_trun <- trun.p(par = 0, family = "PO", type = "left")
                  out <- matrix(NA, nrow = num_units, ncol = num_spp)
                  a <- tmp_trun((object$y-1)[which(object$y>0)], mu = exp(object$linear_predictor[which(object$y>0)]))
                  b <- tmp_trun(object$y[which(object$y>0)], mu = exp(object$linear_predictor[which(object$y>0)]))
                  out[which(object$y>0)] <- runif(length(a), min = a, max = b)
                  }
#           if(object$family$family[1] %in% c("ztnegative.binomial"))
#                out <- runif(length(object$y), min = pztnbinom(y-1, mu = object$fitted, size = matrix(1/object$dispparam, nrow = num_units, ncol = num_spp, byrow = TRUE)), 
#                     max = pztnbinom(object$y, mu = object$fitted, size = matrix(1/object$dispparam, nrow = num_units, ncol = num_spp, byrow = TRUE)))
          
          if(type == "dunnsmyth")
               out <- qnorm(out)
          
          out <- matrix(out, nrow = num_units, ncol = num_spp)
          rownames(out) <- rownames(object$y)
          colnames(out) <- colnames(object$y)
          }
          
     
     set.seed(NULL)
     return(out)
     }

