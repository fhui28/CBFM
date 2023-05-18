#' @title Use the gratia package to construct smooth estimates and parametric effect tables from a CBFM fit
#' 
#' @description 
#' `r lifecycle::badge("stable")`
#' 
#' Calculates parametric effects and smooth estimates for all species from a \code{CBFM} fit, based on applying the corresponding functions [gratia::parametric_effects()] and [gratia::smooth_estimates()].
#' 
#' @param object An object of class \code{CBFM}.
#' @param ncores To speed up calculation of the smooth estimates and parametric effects, parallelization can be performed, in which case this argument can be used to supply the number of cores to use in the parallelization. Defaults to \code{detectCores()-1}.
#' @param ... Not used.
#' 
#' @details 
#' This function is more or less direct wrapper for [gratia::parametric_effects()] and [gratia::smooth_estimates()], applied to every species; please see the help files for those two functions for more details. Please note the returned outputs \code{all_parametric_estimates} and \code{all_smooth_effects} can be *very* large data frames!
#' 
#' Acknowledgments to Gavin Simpson for building and assisting with the implementation of the \code{gratia} package!
#' 
#' 
#' @return A list with the following components (as appropriate):
#' \item{all_parametric_effects: }{If \code{formula} included any parametric terms excluding the intercept, then a long format data frame is returned containing each estimated parametric effect for each species, which is then primarily used for visualizing the estimated parametric model terms. The data frame is effectively constructed by applying [gratia::parametric_effects()] for each species. If no smoothing terms are included in \code{formula}, then this will equal to \code{NULL}.}

#' \item{allzi_parametric_effects: }{If \code{ziformula} included any parametric terms excluding the intercept, then a long format data frame is returned containing each estimated parametric effect for each species in relation to the probability of zero-inflation, which is then primarily used for visualizing the estimated parametric model terms. The data frame is effectively constructed by applying [gratia::parametric_effects()] for each species. If no smoothing terms are included in \code{ziformula}, then this will equal to \code{NULL}.}

#' \item{all_smooth_estimates: }{If \code{formula} included any smoothing terms excluding the intercept, then a long format data frame is returned containing each estimated smoothed effect (evaluated on a grid of evenly spaced values over the range of each corresponding covariate) for each species, which is then primarily used for visualizing the smooth model terms. The data frame is effectively constructed by applying [gratia::smooth_estimates()] for each species. If no smoothing terms are included in \code{formula}, then this will equal to \code{NULL}.}

#' \item{allzi_smooth_estimates: }{If \code{ziformula} included any smoothing terms excluding the intercept, then a long format data frame is returned containing each estimated smoothed effect (evaluated on a grid of evenly spaced values over the range of each corresponding covariate) for each species in relation to the probability of zero-inflation, which is then primarily used for visualizing the smooth model terms. The data frame is effectively constructed by applying [gratia::smooth_estimates()] for each species. If no smoothing terms are included in \code{ziformula}, then this will equal to \code{NULL}.}

#' 
#' @author Francis K.C. Hui <fhui28@gmail.com>, Chris Haak
#' 
#' @seealso [CBFM()] for fitting CBFMs, [model.matrix.CBFM()] for extracting the model matrix from a CBFM fit, [predict.CBFM()] for constructing predictions from a CBFM fit, and [residuals.CBFM()] for calculating various types of residuals.
#' 
#' @examples
#' \dontrun{
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
#' fitcbfm <- CBFM(y = simy, formula = useformula, data = dat, 
#' B_space = basisfunctions, family = binomial(), control = list(trace = 1))
#' 
#' # Example of plotting parametric model terms
#' geteffects <- gratia_effects(fitcbfm)
#' geteffects$all_parametric_effects$species <- geteffects$all_parametric_effects$species %>%
#' fct_inorder
#' ggplot(geteffects$all_parametric_effects, aes(x = value, y = partial, color = species)) +
#' geom_line() +
#' geom_rug(aes(x = value), sides = "b", show.legend = FALSE, color = "black") +
#' facet_wrap(. ~ term, nrow = 2) +
#' labs(x = "Covariate", y = "Effect") +
#' theme_bw() +
#' theme(legend.position = "bottom")
#' 
#'  
#' ##------------------------------
#' ## **Example 1b: Repeat Example 1a but illustrate the use of smoothing terms in CBFM**
#' ## Since the true model only involves parametric terms, then we do not expect its performance
#' ## to be as good as assuming the right form for the mean model.
#' ## It is purely for illustration purposes. 
#' ## Please note this will take a while to run...get a cup of tea and stretch your legs! 
#' ##------------------------------
#' 
#' # Fit CBFM 
#' tic <- proc.time()
#' useformula <- ~ temp + s(depth) + chla + s(O2)
#' fitcbfm_gam <- CBFM(y = simy, formula = useformula, 
#' data = dat, B_space = basisfunctions, family = binomial(), control = list(trace = 1))
#' toc <- proc.time()
#' toc-tic
#' 
#' summary(fitcbfm_gam) %>% 
#' str
#' 
#' # Example of plotting smooth model terms
#' geteffects <- gratia_effects(fitcbfm_gam)
#' geteffects$all_smooth_estimates$species <- geteffects$all_smooth_estimates$species %>%
#' fct_inorder
#' ggplot() +
#' geom_line(data = geteffects$all_smooth_estimates %>% subset(smooth == "s(depth)"), 
#' aes(x = depth, y = est, color = species), show.legend = FALSE) +
#' geom_rug(aes(x = depth), data = dat, sides = "b", color = "black") +
#' labs(x = "depth", y = "Effect") +
#' theme_bw()
#' 
#' 
#' }
#' 
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar% 
#' @importFrom gratia smooth_estimates parametric_effects
#' 
#' @export
#' @md

gratia_effects <- function(object, ncores = NULL, ...) {
     if(!inherits(object, "CBFM")) 
          stop("`object' is not of class \"CBFM\"")
     if(!object$stderrors) 
          stop("Parametric estimates and smooth effects can not be constructed if the CBFM object has stderrors = FALSE.")
     
     if(is.null(ncores))
          registerDoParallel(cores = detectCores()-1)
     if(!is.null(ncores))
          registerDoParallel(cores = ncores)
     
     
     num_spp <- nrow(object$betas)
     out <- list()
     
     all_parametric_effects <- foreach(j = 1:num_spp) %dopar% .calc_parametric_effects(j = j, object = object)
     out$all_parametric_effects <- do.call(rbind, lapply(all_parametric_effects, function(x) x$out))
     out$allzi_parametric_effects <- do.call(rbind, lapply(all_parametric_effects, function(x) x$ziout))
     rm(all_parametric_effects)
     
     all_smooth_estimates <- foreach(j = 1:num_spp) %dopar% .calc_smooth_estimates(j = j, object = object)
     out$all_smooth_estimates <- do.call(rbind, lapply(all_smooth_estimates, function(x) x$out))
     out$allzi_smooth_estimates <- do.call(rbind, lapply(all_smooth_estimates, function(x) x$ziout))
     rm(all_smooth_estimates)
     
     return(out)
     }


