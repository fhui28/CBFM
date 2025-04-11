#' @title Construct spatio-temporal ordinations from a CBFM fit
#' 
#' @description 
#' `r lifecycle::badge("experimental")`
#' 
#' Calculates observational unit-specific scores and corresponding species-specific loadings from a fitted \code{CBFM} object, potentially at new values of the basis functions. which in turn can be used as a means of spatio-temporal ordination; see Hui et al., (2023) for more details.
#' 
#' @param object An object of class \code{CBFM}.
#' @param num_comp The number of ordination axes to construct. 
#' @param new_B_space A matrix of new spatial basis functions at which ordinations are to be constructed. If this is not provided, then ordinations based on the original set of basis functions are returned. If at least one of \code{new_B_space/new_B_time/new_B_spacetime} are provided, then it is \emph{assumed} ordinations at new basis function values are desired. 
#' @param new_B_time A matrix of new temporal basis functions at which ordinations are to be constructed. If this is not provided, then ordinations based on the original set of basis functions are returned. If at least one of \code{new_B_space/new_B_time/new_B_spacetime} are provided, then it is \emph{assumed} ordinations at new basis function values are desired. 
#' @param new_B_spacetime A matrix of new spatio-temporal basis functions at which ordinations are to be constructed. If this is not provided, then ordinations based on the original set of basis functions are returned. If at least one of \code{new_B_space/new_B_time/new_B_spacetime} are provided, then it is \emph{assumed} ordinations at new basis function values are desired. 
#' @param alpha A numeric scalar between 0 and 1 that controls the relative scaling of the score and corresponding loadings, when constructing a biplot. Defaults to 0.5, which gives equal weighting to above.
#' @param ... Not used.
#' 
#' @details 
#' For a CBFM, one way of performing spatio-temporal ordination is as follows. For observational unit \eqn{i=1,\ldots,N} and species \eqn{j=1,\ldots,m}, recall for the purposes of this package, the CBFM is defined by the mean model:
#' 
#' \deqn{g(\mu_{ij}) = \eta_{ij} = x_i^\top\beta_j + b_i^\top a_j,}
#'
#' where \eqn{g(.)} is a known link function, \eqn{x_i} denotes a vector of predictors for unit \eqn{i} i.e., the \eqn{i}-th row from the created model matrix, \eqn{\beta_j} denotes the corresponding regression coefficients for species \eqn{j}, \eqn{b_i} denotes a vector of spatial and/or temporal basis functions for unit \eqn{i} , and \eqn{a_j} denotes the corresponding regression coefficients for species \eqn{j}. 
#'
#' Consider now the estimates values of \eqn{b_i^\top a_j} and construct an \eqn{N \times m} matrix of estimated linear predictors with this. We can then construct an ordination from this matrix, albeit in a bit of an ad-hoc manner, by performing a singular value decomposition to obtain matrices of left and and right singular vectors, along with a diagonal matrix of singular values in descending order. The left singular vectors correspond to observational unit-specific scores (similar to latent variables from an latent variable model; see for example Hui et al., 2015, Warton et al., 2015, Niku et al., 2019, van der Veen, 2021) while the right singular vectors correspond to the corresponding species-specific loadings (similar to latent variables from an latent variable model). By taking only the first \code{num_comp} of the left and right singular vectors and scaling according to the singular values, we can then use these as a means of spatial and/or temporal ordination (Hui et al., 2023)
#'
#' The scores and loadings and thus ordinations can also be constructed at a new set of observational units. This is achieved by supplying the values of \code{new_B_space/new_B_time/new_B_spacetime} as appropriate, similar to how [predict.CBFM()] functions. Ordinations at new observational units are often desired for spatial temporal ordinations e.g., if the ordination is to be constructed and visualized on a grid of spatial locations, and potentially over time.
#'
#' Please note the current implementation does *not* produce plots automatically: it only returns the scores and loadings, which the user then has to manually use to construct plots. Future versions of the package may change this! 
#'
#' @return A list with the following items:
#' \item{scores: }{A matrix of observational unit-specific scores.}

#' \item{loadings: }{A matrix of species-specific loadings.}
#' 
#' @details #Warning
#' For zero-truncated count CBFMs, observational unit-specific scores are produced for all observational units, although perhaps these should be treated cautiously if the response matrix itself contains zero counts.
#' 
#' @author Francis K.C. Hui <fhui28@gmail.com>, Chris Haak
#' 
#' @references
#' Hui, F. K. C., Taskinen, S., Pledger, S., Foster, S. D., and Warton, D. I. (2015). Model‐based approaches to unconstrained ordination. Methods in Ecology and Evolution, 6, 399-411.
#'
#' Hui, F. K. C., Warton, D. I., Foster, S. D., & Haak, C. R. (2023). Spatiotemporal joint species distribution modelling: A basis function approach. Methods in Ecology and Evolution, 14(8), 2150-2164.
#' 
#' Niku, J., Hui, F. K. C., Taskinen, S., and Warton, D. I. (2019). gllvm: Fast analysis of multivariate abundance data with generalized linear latent variable models in R. Methods in Ecology and Evolution, 10, 2173-2182.
#'
#' van der Veen, B., Hui, F. K. C., Hovstad, K. A., Solbu, E. B., & O'Hara, R. B. (2021). Model‐based ordination for species with unequal niche widths. Methods in Ecology and Evolution. In press.
#' 
#' Warton, D. I., Blanchet, F. G., O'Hara, R. B., Ovaskainen, O., Taskinen, S., Walker, S. C., and Hui, F. K. C. (2015). So many variables: joint modeling in community ecology. Trends in Ecology and Evolution, 30, 766-779.
#'
#'
#' @seealso [CBFM()] for fitting CBFMs, [corB()] for calculating residual between-species (cross-)correlations due to the basis functions, [predict.CBFM()] for constructing predictions from a CBFM fit.
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
#' library(geoR)
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
#' X <- mvtnorm::rmvnorm(num_sites, mean = rep(0,4))
#' colnames(X) <- c("temp", "depth", "chla", "O2")
#' dat <- data.frame(xy, X)
#' mm <- model.matrix(~ temp + depth + chla + O2 - 1, data = dat) %>% 
#' scale %>% 
#' as.matrix
#' 
#' # Simulate latent variable component
#' true_lvs <- grf(grid = cbind(xy$x, xy$y), nsim = 2, cov.model = "exponential",
#' cov.pars = c(1, 2))$data %>%
#'      as.matrix
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
#' getords <- ordinate(fitcbfm)
#' 
#' orddat <- data.frame(dat[,c("x","y")], getords$scores) 
#' 
#' ggplot(orddat, aes(x = x, y = y, color = Axis1)) +
#' geom_point() +
#' scale_color_viridis_c() +
#' labs(x = "x coordinate", y = "y coordinate", color = "Site score 1") +
#' theme_bw()
#' 
#' ggplot(orddat, aes(x = x, y = y, color = Axis2)) +
#' geom_point() +
#' scale_color_viridis_c() +
#' labs(x = "x coordinate", y = "y coordinate", color = "Site score 2") +
#' theme_bw()
#' }
#' 
#' @aliases ordinate ordinate.CBFM 
#' @export
#' @export ordinate.CBFM
#' @md

ordinate.CBFM <- function(object, num_comp = 2, new_B_space = NULL, new_B_time = NULL, new_B_spacetime = NULL, alpha = 0.5, ...) {
        if(!inherits(object, "CBFM")) 
                stop("`object' is not of class \"CBFM\"")

        if(object$family$family[1] %in% c("ztpoisson"))
                warning("For zero-truncated count CBFMs, observational unit-specific scores are produced for all observational units. However these should be treated cautiously if the response matrix itself contains zero counts.")
        
        if(is.null(new_B_space) & is.null(new_B_time) & is.null(new_B_spacetime)) {
                linpred_basisfunctions <- tcrossprod(object$B, object$basis_effects_mat)
                }
        else {
                linpred_basisfunctions <- .predictorforordinate.CBFM(object = object, new_B_space = new_B_space, new_B_time = new_B_time, new_B_spacetime = new_B_spacetime)           
                }

          
        do_svd <- svd(linpred_basisfunctions)
     
        new_lvs <- do_svd$u[, 1:num_comp,drop=FALSE] %*% diag(x = do_svd$d[1:num_comp]^alpha, nrow = num_comp)
        #rownames(new_lvs) <- rownames(object$y)
        new_Loadings <- do_svd$v[, 1:num_comp,drop=FALSE] %*% diag(x = do_svd$d[1:num_comp]^alpha, nrow = num_comp)
        rownames(new_Loadings) <- colnames(object$y)
        colnames(new_lvs) <- colnames(new_Loadings) <- paste0("Axis", 1:num_comp)
     
        out <- list(scores = new_lvs, loadings = new_Loadings)
        return(out)
        }


#' @method ordinate CBFM
#' @export ordinate
ordinate <- function(object, ...) {
        UseMethod("ordinate")
        }    

     
.predictorforordinate.CBFM <- function(object, new_B_space = NULL, new_B_time = NULL, new_B_spacetime = NULL) {
        newB <- NULL
        if(!is.null(new_B_space)) {
                if(object$num_B_space != ncol(new_B_space))
                        stop("The number of columns of new_B_space does not object$num_B_space.")
                newB <- cbind(newB, new_B_space)
                }
        if(!is.null(new_B_time)) {
                if(object$num_B_time != ncol(new_B_time))
                        stop("The number of columns of new_B_time does not object$num_B_time.")
                newB <- cbind(newB, new_B_time)
                }
        if(!is.null(new_B_spacetime)) {
                if(object$num_B_spacetime != ncol(new_B_spacetime))
                        stop("The number of columns of new_B_spacetime does not object$num_B_spacetime.")
                newB <- cbind(newB, new_B_spacetime)
                }
        if(is.null(newB))
                newB <- object$B
        if(!is.null(newB)) {
                if(ncol(newB) != ncol(object$basis_effects_mat))
                        stop("newB does not contain the same number of columns as the number of columns in object$basis_effects_mat.")
                }
          
        num_spp <- nrow(object$betas)
        num_basisfns <- ncol(newB)
     
        ptpred <- tcrossprod(newB, object$basis_effects_mat)
        ptpred <- as.matrix(ptpred)
        return(ptpred)
        }
     
