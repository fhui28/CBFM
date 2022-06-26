#' @title Construct (cross-)correlations due to measured covariates for a CBFM fit
#' 
#' @description 
#' `r lifecycle::badge("experimental")`
#' 
#' Takes a fitted \code{CBFM} object calculates the between-species correlation matrix due to the measured covariates, along with corresponding uncertainty intervals. Both are constructed via a simulation-based approach. Similar to [predict.CBFM()], this correlation matrix can be calculated based on a different sets of covariates to those used to actually fit the model. Additionally, the user can supplied two sets of covariates (data frames), in which case the function calculates cross-correlations (between and within species) between these two sets of covariates. 
#' 
#' 
#' @param object An object of class \code{CBFM}.
#' @param newdata A data frame containing the values of the covariates at which correlations are to be calculated. If this is not provided, then correlations corresponding to the original data are returned. If \code{newdata} is provided then it should contain all the variables needed for constructing correlations, that is, it can construct a model matrix from this as \code{object$formula}.
#' @param newdata2 A second data frame containing the values of the covariates at which cross-correlations are to be calculated. If this is supplied, then \code{newdata} must also be supplied, as the function assumes then the user desires calculation of cross-correlations. 
#' @param coverage The coverage probability of the uncertainty intervals for the correlations. Defaults to 0.95, which corresponds to 95% uncertainty intervals.
#' @param ncores To speed up calculation of the uncertainty estimates, parallelization can be performed, in which case this argument can be used to supply the number of cores to use in the parallelization. Defaults to \code{detectCores()-1}.
#' @param num_sims The number of Monte-Carlo examples to simulate.
#' 
#' 
#' @details 
#' This function is adapted from and behaves somewhat similarly to [boral::get.enviro.cor()], in calculating a between-species correlation matrix due to the measured covariates i.e., shared environmental response, along with corresponding uncertainty intervals. Recall the general form of the mean regression model for the CBFM is given by 
#' 
#' \deqn{g(\mu_{ij}) = \eta_{ij} = x_i^\top\beta_j + b_i^\top a_j,}
#' 
#' where \eqn{x_i} denotes a vector of predictors for unit \eqn{i} i.e., the \eqn{i}-th row from the created model matrix, \eqn{\beta_j} denotes the corresponding regression coefficients for species \eqn{j}.
#' 
#' The covariance and hence correlation between two species \eqn{j} and \eqn{j'} that can be attributed to the measured covariates is then based on examining the components \eqn{x_i^\top\beta_j} and \eqn{x_i^\top\beta_{j'}} across the observational units; see equation 4 in Pollock et al., (2014) for the precise formula, as well as Warton et al., (2015) and Hui (2016) among others. Both the point estimate and the uncertainty intervals for the correlation are constructed by simulation. Specifically, species-specific regression coefficients \eqn{\beta_j} are sampled from their approximate large sample normal distribution (i.e., basically a Gaussian approximation to the posterior distribution of the parameters; see [CBFM()] and the section on estimation and inference), which are then used to calculate the correlations. This sampling and calculation is then performed a large number of times (as governed by \code{num_sims}), after which a point estimate of the correlations is obtained by taking the sample average, while uncertainty intervals are based on taking sample quantiles.
#' 
#' Note that because this function calculates correlations as based on component of the linear predictor \eqn{x_i^\top\beta_j}, then it can not be applied to  \code{CBFM_hurdle} object (which by construction contains two linear predictors, so the user has to decide which component of the hurdle model they are interested in). Analogously, for zero-inflated CBFMs this function currently only calculates the between-species correlation matrix due to the measured covariates in \code{object$formula} i.e., the mean of the count component. It is *not* able to calculate correlations due to measured covariates in \code{object$ziformula} i.e., in modeling the probability of zero-inflation.
#' 
#' With the above definition of the between-species correlation, note that the predictors on which the correlations are constructed need not be the same those used in fitting the original CBFM i.e., the \eqn{x_i}'s can be different to those of \code{object$data}. This is handled via the \code{newdata} argument. Additionally, it is possible to calculate within and between species *cross-correlations* across two different sets of measured predictors. That is, correlations are calculated between \eqn{x_i^\top\beta_j} and \eqn{x_{i2}^\top\beta_{j'}}, where \eqn{x_i} and \eqn{x_{i2}} can be different sets of measured predictors. This is handled by supplying both \code{newdata} and \code{newdata2} arguments simultaneously. Cross-correlations may be useful, say, if the two sets of measurement predictors reflect two different sets of sampling units, and we are interested in how similar (or lack of) the species communities are in terms of their environmental response across these two sets (Ovakainen et al., 2017). Another example if is the same set of observational units are visited at two different time points, and we are interested in how similarity (or lack of) the environmental responses within and between species are between these two time points.
#' 
#' NOTE: A cross-correlation matrix is not going to be a standard correlation matrix in the sense of having the ones alone the diagonal. This is because even for the same species \eqn{j = j'}, the correlation is not guaranteed to be equal to one as the covariates being considered can be different.   
#' 
#'
#' @return A list with the following components is returned:
#' \item{correlation: }{A matrix of (cross-)correlation values.}

#' \item{lower: }{A matrix of the lower bound of the uncertainty intervals for the correlations}

#' \item{upper: }{A matrix of the upper bound of the uncertainty intervals for the correlations}
#' 
#' 
#' @author Francis K.C. Hui <fhui28@gmail.com>, Chris Haak
#' 
#' 
#' @references
#' Hui, F. K. C. (2016). boral-Bayesian ordination and regression analysis of multivariate abundance data in R. Methods in Ecology and Evolution, 7, 744-750.
#' 
#' Ovaskainen, O., Tikhonov, G., Norberg, A., Guillaume Blanchet, F., Duan, L., Dunson, D., and Abrego, N. (2017). How to make more out of community data? A conceptual framework and its implementation as models and software. Ecology letters, 20, 561-576.
#' 
#' Pollock, L. J., Tingley, R., Morris, W. K., Golding, N., O'Hara, R. B., Parris, K. M., Vesk, P. A., and McCarthy, M. A. (2014). Understanding co‐occurrence by modelling species simultaneously with a Joint Species Distribution Model (JSDM). Methods in Ecology and Evolution, 5, 397-406.
#' 
#' Warton, D. I., Blanchet, F. G., O'Hara, R. B., Ovaskainen, O., Taskinen, S., Walker, S. C., and Hui, F. K. C. (2015). So many variables: joint modeling in community ecology. Trends in Ecology and Evolution, 30, 766-779.
#' 
#' @seealso [CBFM()] for fitting CBFMs, and [corB()] for calculating residual between-species (cross-)correlations due to the basis functions.
#'
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
#' library(corrplot)
#' 
#' ##------------------------------
#' ## **Example 1: Fitting a CBFM to spatial multivariate presence-absence data** 
#' ## simulated from a spatial latent variable model
#' ## Please note the data generation process (thus) differs from CBFM.
#' ##------------------------------
#' set.seed(2021)
#' num_sites <- 1000 # 500 (units) sites for training set + 500 sites for external calculation
#' num_spp <- 50 # Number of species
#' num_X <- 4 # Number of regression slopes
#' 
#' spp_slopes <- matrix(runif(num_spp * num_X, -1, 1), nrow = num_spp)
#' spp_intercepts <- runif(num_spp, -2, 0)
#' 
#' # Simulate spatial coordinates and environmental covariate components
#' # We will use this information in later examples as well
#' xy <- data.frame(x = runif(num_sites, 0, 5), y = runif(num_sites, 0, 5))
#' X <- rmvnorm(num_sites, mean = rep(0,4)) 
#' colnames(X) <- c("temp", "depth", "chla", "O2")
#' dat <- data.frame(xy, X)
#' mm <- model.matrix(~ temp + depth + chla + O2 - 1, data = dat) %>% 
#' scale %>% 
#' as.matrix
#' 
#' # Simulate latent variable component
#' # We will use this information in later examples as well
#' true_lvs <- RFsimulate(model = RMexp(var=1, scale=2), 
#' x = xy$x, y = xy$y, n = 2)@data %>% 
#' as.matrix
#' spp_loadings <- matrix(runif(num_spp * 2, -1, 1), nrow = num_spp) 
#' set.seed(NULL)
#' 
#' # Simulate spatial multivariate abundance data (presence-absence)
#' # We will use this information in later examples as well
#' eta <- tcrossprod(cbind(1,mm), cbind(spp_intercepts,spp_slopes)) + 
#' tcrossprod(true_lvs, spp_loadings)
#' simy <- matrix(rbinom(num_sites * num_spp, size = 1, 
#' prob = plogis(eta)), nrow = num_sites)
#' 
#' # Form training and test sets
#' dat_train <- dat[1:500,]
#' dat_external <- dat[501:1000,]
#' simy_train <- simy[1:500,]
#' rm(X, mm, spp_loadings, true_lvs, xy, simy, dat)
#' 
#' 
#' 
#' # Set up spatial basis functions for CBFM -- Most users will start here! 
#' # We will also use this basis functions in some later examples
#' num_basisfunctions <- 25 # Number of spatial basis functions to use
#' # Training set basis functions
#' train_basisfunctions <- mrts(dat_train[,c("x","y")], num_basisfunctions) %>% 
#' as.matrix %>%
#' {.[,-(1)]} # Remove the first intercept column
#' 
#' # Fit CBFM 
#' tic <- proc.time()
#' useformula <- ~ temp + depth + chla + O2
#' fitcbfm <- CBFM(y = simy_train, formula = useformula, data = dat_train, 
#' B_space = train_basisfunctions, family = binomial(), control = list(trace = 1))
#' toc <- proc.time()
#' toc - tic
#' 
#' 
#' # Calculate between-species correlations based on measured covariates in training data
#' getcor <- corX(fitcbfm)
#' corrplot(getcor$corr, method = "square", type = "lower", order = "hclust")
#' 
#' 
#' # Calculate between-species correlations based on measured covariates in external (new) data
#' getcor <- corX(fitcbfm, newdata = dat_external)
#' corrplot(getcor$corr, method = "square", type = "lower", order = "hclust")
#' 
#' 
#' # Calculate species cross-correlations between measured covariates in training and 
#' # external (new) data
#' # This may be useful, for example, if the training and external data were the same sites visited 
#' # at two points in time, and the user is interested in the similarity (or lack of) 
#' # of the environmental responses within and between species are between these two time points. 
#' # Note the resulting cross-correlation matrix is not strictly a correlation matrix 
#' # in the sense of having ones on the diagonals; all elements, including diagonals, 
#' # will lie between -1 and 1.  
#' getcrosscor <- corX(fitcbfm, newdata = dat_train, newdata2 = dat_external)
#' corrplot(getcrosscor$corr, method = "square", type = "lower", is.corr = FALSE)
#' }
#' 
#' @export
#' 
#' @importFrom foreach foreach %dopar%
#' @import Matrix
#' @importFrom abind abind
#' @importFrom doParallel registerDoParallel
#' @importFrom mgcv gam model.matrix.gam predict.gam
#' @importFrom parallel detectCores
#' @importFrom stats qnorm
#' @md

corX <- function(object, newdata = NULL, newdata2 = NULL, coverage = 0.95, ncores = NULL, num_sims = 500) {
        if(!inherits(object, "CBFM")) 
                stop("`object' is not of class \"CBFM\"")

        if(is.null(ncores))
                registerDoParallel(cores = detectCores()-1)
        if(!is.null(ncores))
                registerDoParallel(cores = ncores)
        
        if(!is.null(newdata2) & is.null(newdata)) {
                stop("If newdata2 is supplied, then newdata must also be supplied (for calculation of cross-correlations). Both matrices must have the same dimensions.")
                }
        
        ##--------------------------------
        ## Construct new_X and new_X2, if appropriate.
        ##--------------------------------
        tmp_formula <- as.formula(paste("response", paste(as.character(object$formula),collapse="") ) )
        nullfit <- gam(tmp_formula, data = data.frame(response = runif(nrow(object$y)), object$data), fit = TRUE, control = list(maxit = 1))
        if(is.null(newdata)) 
            new_X <- new_X2 <- predict.gam(nullfit, type = "lpmatrix")
        if(!is.null(newdata))
                new_X <- predict.gam(nullfit, newdata = data.frame(newdata), type = "lpmatrix")
        if(is.null(newdata2))
                new_X2 <- new_X
        if(!is.null(newdata2)) {
                new_X2 <- predict.gam(nullfit, newdata = data.frame(newdata2), type = "lpmatrix")
                }
        rm(tmp_formula, nullfit)
        
        if((ncol(new_X) != ncol(object$betas)) | (ncol(new_X2) != ncol(object$betas)))
                stop("Number of columns in new_X and new_X2 should match the number of columns in object$betas.")
        if(nrow(new_X) != nrow(new_X2))
                stop("Number of rows in new_X and new_X2 should match.")

        
        ##----------------------
        ## Construct (cross)-correlations
        ##----------------------
        num_spp <- nrow(object$betas)
        num_X <- ncol(object$betas)
        num_basisfunctions <- ncol(object$basis_effects_mat)

        if(object$stderrors == FALSE)
                stop("Standard errors can not be calculated since the covariance matrix estimate was not detected to be available in object.")

        eta1 <- as.matrix(tcrossprod(new_X, object$betas))
        eta2 <- as.matrix(tcrossprod(new_X2, object$betas))
        colnames(eta2) <- colnames(eta2) <- rownames(object$betas)          

        # if(!se_fit) {
        #         out <- cor(eta1, eta2)
        #         return(out)
        #         }

        message("Simulation used to calculate for uncertainty intervals for correlations. This could take a while...enjoy a cup of matcha latte while you're waiting uwu")
        ci_alpha <- qnorm((1-coverage)/2, lower.tail = FALSE)

        mu_vec <- as.vector(t(cbind(object$zeroinfl_prob_intercept, object$betas, object$basis_effects_mat)))
        bigcholcovar <- as.matrix(rbind(cbind(object$covar_components$topleft, object$covar_components$topright),
                                          cbind(t(object$covar_components$topright), object$covar_components$bottomright)))
        bigcholcovar <- t(chol(bigcholcovar))
                
                
        innersim_etafn <- function(j) {
                parameters_sim <- matrix(mu_vec + as.vector(bigcholcovar %*% rnorm(length(mu_vec))), nrow = num_spp, byrow = TRUE)
                betas_sim <- parameters_sim[,(ncol(parameters_sim)-num_basisfunctions-num_X+1):(ncol(parameters_sim)-num_basisfunctions), drop=FALSE]
                rm(parameters_sim) # The above could be simplified since I only need to simulate betas, but whatever...

                eta1 <- as.matrix(tcrossprod(new_X, betas_sim))
                eta2 <- as.matrix(tcrossprod(new_X2, betas_sim))
                return(cor(eta1, eta2))
                }
    
        allcors <- foreach(j = 1:num_sims) %dopar% innersim_etafn(j = j)

        rm(bigcholcovar, mu_vec)
        allcors <- abind(allcors, along = 3)
        ptcor <- apply(allcors, c(1,2), mean)
        alllower <- apply(allcors, c(1,2), quantile, prob = (1-coverage)/2)
        allupper <- apply(allcors, c(1,2), quantile, prob = coverage + (1-coverage)/2)
        rownames(ptcor) <- rownames(alllower) <- rownames(alllower) <- colnames(object$y)
        colnames(ptcor) <- colnames(alllower) <- colnames(alllower) <- colnames(object$y)
        rm(allcors)
        gc()
        
        return(list(correlation = ptcor, lower = alllower, upper = allupper))
        }

