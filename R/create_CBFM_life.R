#' @title Simulate data from a CBFM
#' 
#' @description 
#' `r lifecycle::badge("stable")`
#' 
#' Simulates spatio-temporal multivariate abundance data based on a CBFM and given the various parameter values as appropriate.
#'
#' @param family a description of the response distribution to be used in the model, as specified by a family function. Please see details below for more information on the distributions currently permitted.
#' @param formula An object of class "formula", which represents a symbolic description of the model matrix to be created (based on using this argument along with the \code{data} argument). Note there should be nothing on the left hand side of the "~". Formulas based on generalized additive models or GAMs are permitted (at least, for the smoothing terms we have tried so far!); please see [mgcv::formula.gam()] and [mgcv::s()] for more details. 
#' @param ziformula An object of class "formula", which represents a symbolic description of the model matrix to be created for the zero-inflation component (based on using this argument along with the \code{data} argument), if appropriate. Note there should be nothing on the left hand side of the "~". Formulas based on generalized additive models or GAMs are permitted (at least, for the smoothing terms we have tried so far!); please see [mgcv::formula.gam()] and [mgcv::s()] for more details.
#' @param data A data frame containing covariate information, from which the model matrix is to be created (based on this argument along with the \code{formula} argument). 
#' @param B_space An optional matrix of spatial basis functions to be included in the CBFM. One of \code{B_space}, \code{B_time}, or \code{B_spacetime} must be supplied. The basis function matrix may be sparse or dense in form; please see the details and examples later on for illustrations of how they can constructed.
#' @param B_time An optional of matrix of temporal basis functions to be included in the CBFM. One of \code{B_space}, \code{B_time}, or \code{B_spacetime} must be supplied. The basis function matrix may be sparse or dense in form; please see the details and examples later on for illustrations of how they can constructed.
#' @param B_spacetime An optional of matrix of spatio-temporal basis functions to be included in the CBFM e.g., formed from a tensor-product of spatial and temporal basis functions. One of \code{B_space}, \code{B_time}, or \code{B_spacetime} must be supplied. The basis function matrix may be sparse or dense in form; please see the details and examples later on for illustrations of how they can constructed.
#' @param offset A matrix of offset terms to be applied in association with the \code{formula} argument. 
#' @param zioffset A matrix of offset terms to be applied in association with the \code{ziformula} argument. 
#' @param betas A matrix of species-specific regression coefficients corresponding to the model matrix created. The number of rows in \code{betas} is equal to the number of species in the resulting simulated dataset.
#' @param zibetas A matrix of species-specific regression coefficients corresponding to the model matrix created for the zero-inflation component. The number of rows in \code{zibetas} is equal to the number of species in the resulting simulated dataset.
#' @param basis_effects_mat A matrix of species-specific regression coefficients corresponding to the combined matrix of basis functions. If supplied, then number of rows in \code{basis_effects_mat} is equal to the number of species in the resulting simulated dataset. If it is not supplied, then species-specific regression coefficients are simulated based on the \code{Sigma} and \code{G} arguments.   
#' @param Sigma A list containing the covariance matrix of the species-specific regression coefficients, corresponding to the basis functions supplied. This list should contain the elements \code{space}, \code{time} and/or \code{spacetime} as appropriate e.g., if only \code{B_space} is supplied then \code{Sigma$Space} must be supplied.
#' @param G A list containing the baseline between-species correlation matrix, corresponding to the basis functions supplied. This list should contain the elements \code{space}, \code{time} and/or \code{spacetime} as appropriate e.g., if only \code{B_space} is supplied then \code{G$Space} must be supplied.
#' @param trial_size Trial sizes to use for binomial distribution. This can either equal a scalar or a matrix with the same dimension as the simulated response matrix is to be.
#' @param dispparam A vector of species-specific dispersion parameters, to be used for distributions that require one.  
#' @param powerparam A vector of species-specific power parameters, to be used for distributions that require one. 
#' @param zeroinfl_prob A vector of species-specific probabilities of zero-inflation, to be used for distributions that require one. Note \code{ziformula} is supplied, then this argument is ignored.
#' @param max_resp A upper bound to limit the maximum value of responses obtained. This is useful if the user wants, say, all counts to not exceed a particular value. In such case, the function will attempt to simulate counts that do not \code{max_resp}. Note it only \emph{attempts} this: it will give up after 10 unsuccessful attempts and then return whatever is simulated on the 10-th attempt.
#' @param only_y If \code{TRUE}, then only the simulated spatio-temporal multivariate abundance response matrix is returned. Otherwise if \code{FALSE}, then additional information about is returned.
#' 
#' @details 
#' Simulates spatio-temporal multivariate abundance data from a community-level basis function model (CBFM). For the purposes of the package, the CBFM is characterized by the following mean regression model: for observational unit \eqn{i=1,\ldots,N} and species \eqn{j=1,\ldots,m}, we have
#' 
#' \deqn{g(\mu_{ij}) = \eta_{ij} = x_i^\top\beta_j + b_i^\top a_j,}
#'
#' where \eqn{g(.)} is a known link function, \eqn{x_i} denotes a vector of predictors for unit \eqn{i} i.e., the \eqn{i}-th row from the created model matrix, \eqn{\beta_j} denotes the corresponding regression coefficients for species \eqn{j}, \eqn{b_i} denotes a vector of spatial, temporal, and/or spatio-temporal basis functions for unit \eqn{i} , and \eqn{a_j} denotes the corresponding regression coefficients for species \eqn{j}. In the function, \eqn{x_i} is created based on the \code{formula} and \code{data} arguments, \eqn{\beta_j} is supplied as part of the \code{betas} argument, and \eqn{b_i} is formed from the \code{B_space}, \code{B_time} and \code{B_spacetime} arguments. Finally, \eqn{a_j} is either supplied directly as part of \code{basis_effects_mat} argument, or generated based on the \code{Sigma} and \code{G} arguments. 
#' 
#' As an example, suppose we have a CBFM which involves spatial and temporal (but no spatio-temporal) basis functions. Then \eqn{b_i = (b_{i,space}, b_{i,time})} is formed from the \eqn{i}-th rows of \code{B_space} and \code{B_time}, while \eqn{a_j = (a_{j,space}, a_{j,time})} comes from the \eqn{j}-th row \code{basis_effects_mat}. If \code{basis_effects_mat} is not supplied, then it is instead obtain by simulating 
#' 
#' \deqn{(a_{1,space}, \ldots, a_{m,space}) \sim N(0, kronecker(G_{space}, \Sigma_{space})),} 
#' 
#' where \eqn{G_{space}} and \eqn{\Sigma_{space}} are supplied from \code{G$space} and \code{Sigma$space} respectively, and \eqn{kronecker(\cdot)} is the Kroneckker product operator. Similarly, we have \eqn{(a_{1,time}, \ldots, a_{m,time}) \sim N(0, kronecker(G_{time}, \Sigma_{time}))}. 
#' 
#' Based on the mean model given above, responses \eqn{y_{ij}} are then simulated from the assumed distribution, using the additional dispersion and power parameters as appropriate.
#' 
#' \subsection{Distributions}{
#' 
#' Currently the following response distributions are permitted: 
#' \describe{
#' \item{\code{betalogitfam()}: }{Beta distribution using a logit link. The corresponding mean-variance relationship is given by \eqn{V = \mu(1-\mu)/(1+\phi)} where \eqn{\mu} denotes the mean and \eqn{\phi} is the dispersion parameter.}

#' \item{\code{binomial(link = "logit")}: }{Binomial distribution, noting only the logit link is permitted. The corresponding mean-variance relationship is given by \eqn{V = N_{trial}\mu(1-\mu)} where \eqn{\mu} denotes the mean and \eqn{N_{trial}} is the trial size.}

#' \item{\code{Gamma(link = "log")}: }{Gamma distribution, noting only the log link is permitted. The corresponding mean-variance relationship is given by \eqn{V = \phi\mu^2} where \eqn{\mu} denotes the mean and \eqn{\phi} is the dispersion parameter.}

#' \item{\code{gaussian(link = "identity")}: }{Gaussian or normal distribution, noting only the identity link is permitted. The corresponding mean-variance relationship is given by \eqn{V = \phi}, where \eqn{\phi} is the dispersion parameter.}

#' \item{\code{poisson(link = "log")}: }{Poisson distribution, noting only the log link is permitted. The corresponding mean-variance relationship is given by \eqn{V = \mu} where \eqn{\mu} denotes the mean.}

#' \item{\code{nb2()}: }{Negative binomial distribution, noting only the log link is permitted. The corresponding mean-variance relationship is given by \eqn{V = \mu + \phi\mu^2} where \eqn{\mu} denotes the mean and \eqn{\phi} is the dispersion parameter.}

#' \item{\code{tweedielogfam()}: }{Tweedie distribution, noting only the log link is permitted. The corresponding mean-variance relationship is given by \eqn{V = \phi\mu^{\rho}} where \eqn{\mu} denotes the mean, \eqn{\phi} is the dispersion parameter, and \eqn{\rho} is the power parameter.}

#' \item{\code{zipoisson()}: }{Zero-inflated Poisson distribution, noting only the log link for the Poisson part is permitted. This partial mass function of the distribution is given by \eqn{f(y) = \pi I(y=0) + (1-\pi) f_{pois}(y)}, where \eqn{\pi} is the probability of being in the zero-inflation component, while \eqn{f_{pois}(y)} is the usual Poisson distribution. The mean of the Poisson distribution is modeled against covariates and basis functions, while the probability of zero-inflation can either be a single, species-specific probability as given in part of \code{zeroinfl_prob}, or can also be modeled against covariates via \code{ziformula}. In the case of the latter, a logit link function is used.}

#' \item{\code{zinb2()}: }{Zero-inflated negative binomial distribution, noting only the log link for the negative binomial part is permitted. The partial mass function of the distribution is given by \eqn{f(y) = \pi I(y=0) + (1-\pi) f_{NB}(y)}, where \eqn{\pi} is the probability of being in the zero-inflation component, while \eqn{f_{NB}(y)} is the usual negative binomial distribution. The mean of the negative binomial distribution is modeled against covariates and basis functions, while the probability of zero-inflation can either be a single, species-specific probability as given in part of \code{zeroinfl_prob}, or can also be modeled against covariates via \code{ziformula}. In the case of the latter, a logit link function is used.}

#' \item{\code{ztpoisson()}: }{Zero-truncated Poisson distribution, noting only the log link is permitted. The partial mass function of the distribution is given by \eqn{f(y) = f_{pois}(y)/(1-f_{pois}(0)}) where \eqn{f_{pois}(y)} is the usual Poisson distribution. The mean of the Poisson distribution is modeled against covariates and basis functions.}

#' \item{\code{ztnb2()}: }{Zero-truncated negative binomial distribution, noting only the log link is permitted. The partial mass function of the distribution is given by \eqn{f(y) = f_{NB}(y)/(1-f_{NB}(0)}) where \eqn{f_{NB}(y)} is the usual negative binomial distribution. The mean of the negative binomial distribution is modeled against covariates and basis functions.}
#' }
#' 
#' Note with zero truncated distributions being available, generating spatio-temporal multivariate abundance data from a hurdle CBFM is possible by combining it separate mechanisms for generating presence-absence responses and a truncated count responses. Please see the examples below for an illustration.
#' }
#' 
#' @return 
#' If \code{only_y = TRUE}, then the simulated spatio-temporal multivariate abundance response matrix. Otherwise, a list with the following components (if applicable):
#' \describe{
#' \item{y }{The simulated spatio-temporal multivariate abundance response matrix.}
#' \item{basis_effects_mat }{The matrix of species-specific regression coefficients corresponding to the combined matrix of basis functions. This either comes directly from the supplied argument or is a simulated as discussed in Details above.}
#' \item{linear_predictors }{The matrix of linear predictors \eqn{\eta_{ij}}'s.}
#' \item{linear_predictors_B }{The matrix of linear predictors corresponding to the basis functions only i.e., \eqn{b_i^\top a_j}'s.}
#' }
#' 
#' @details # Warning
#' Note **no** checks are made on the arguments \code{Sigma} and \code{G} arguments, if supplied, to see if they are positive definite matrices or not. Please be careful about this!
#'
#' @author Francis K.C. Hui <fhui28@gmail.com>, Chris Haak
#' 
#' @seealso [CBFM()] for fitting CBFMs and [simulate.CBFM()] for simulating spatio-temporal multivariate abundance data from a CBFM fit.
#' 
#' @examples
#' \dontrun{
#' library(autoFRK)
#' library(FRK)
#' library(MASS)
#' library(mvtnorm)
#' library(sp)
#' library(geoR)
#' library(tidyverse)
#' 
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
#' useformula <- ~ temp + depth + chla + O2
#' 
#' # Set up spatial basis functions for CBFM 
#' num_basisfunctions <- 25 # Number of spatial basis functions to use
#' basisfunctions <- mrts(dat[,c("x","y")], num_basisfunctions) %>% 
#' as.matrix %>%
#' {.[,-(1)]} # Remove the first intercept column
#' 
#' true_Sigma_space <- rWishart(1, num_basisfunctions+1, 
#' diag(x = 0.1, nrow = num_basisfunctions-1))[,,1]/10
#' true_G_space <- rWishart(1, num_spp+1, diag(x = 0.1, nrow = num_spp))[,,1] %>% 
#' cov2cor
#' 
#' ##-----------------------------------
#' ## **Example 1: Generate spatial multivariate presence-absence data**
#' ##-----------------------------------
#' # Basis function coefficients are simulated based on the supplied values of Sigma and G 
#' simy <- create_CBFM_life(family = binomial(), formula = useformula, data = dat,
#' B_space = basisfunctions, betas = cbind(spp_intercepts, spp_slopes),
#' Sigma = list(space = true_Sigma_space), G = list(space = true_G_space))
#' 
#' 
#' # Generates spatial multivariate presence-absence data 
#' # Manually supply basis function coefficients 
#' spp_basis_coefs <- matrix(rnorm(num_spp * (num_basisfunctions-1), 0, 0.1), nrow = num_spp)
#' simy <- create_CBFM_life(family = binomial(), formula = useformula, data = dat,
#' betas = cbind(spp_intercepts, spp_slopes), basis_effects_mat = spp_basis_coefs, 
#' B_space = basisfunctions)
#' 
#' 
#' ##-----------------------------------
#' ## **Example 2: Generate spatial multivariate count data from a negative binomial distribution**
#' ##-----------------------------------
#' # Basis function coefficients are simulated based on the supplied values of Sigma and G 
#' spp_dispersion <- runif(num_spp, 0, 5)
#' simy <- create_CBFM_life(family = nb2(), formula = useformula, data = dat,
#' B_space = basisfunctions, betas = cbind(spp_intercepts, spp_slopes),
#' dispparam = spp_dispersion, max_resp = 20000, 
#' Sigma = list(space = true_Sigma_space), G = list(space = true_G_space))
#' 
#' 
#' ##-----------------------------------
#' ## **Example 3a: Generate spatial multivariate count data from a zero-inflated**
#' ## **Poisson distribution with constant species-specific probabilities**
#' ##-----------------------------------
#' # Manually supply basis function coefficients 
#' spp_zeroinfl_prob <- runif(num_spp, 0, 0.5)
#' spp_basis_coefs <- matrix(rnorm(num_spp * (num_basisfunctions-1), 0, 0.1), nrow = num_spp)
#' simy <- create_CBFM_life(family = zipoisson(), formula = useformula, data = dat,
#' betas = cbind(spp_intercepts, spp_slopes), basis_effects_mat = spp_basis_coefs, 
#' B_space = basisfunctions, zeroinfl_prob = spp_zeroinfl_prob)
#' 
#' 
#' ##-----------------------------------
#' ## **Example 3b: Generate spatial multivariate count data from a zero-inflated** 
#' ## **Poisson distribution with species-specific probabilities that also vary with covariates**
#' ##-----------------------------------
#' # Manually supply basis function coefficients 
#' spp_zislopes <- matrix(runif(num_spp * num_X, -1, 1), nrow = num_spp)
#' spp_ziintercepts <- runif(num_spp, -2, 0)
#' spp_basis_coefs <- matrix(rnorm(num_spp * (num_basisfunctions-1), 0, 0.1), nrow = num_spp)
#' simy <- create_CBFM_life(family = zipoisson(), formula = useformula, ziformula = useformula, 
#' data = dat, betas = cbind(spp_intercepts, spp_slopes), 
#' zibetas = cbind(spp_ziintercepts, spp_zislopes), basis_effects_mat = spp_basis_coefs, 
#' B_space = basisfunctions)
#' 
#' 
#' ##-----------------------------------
#' ## **Example 4a: Generate spatial multivariate count data from a zero-inflated** 
#' ## **negative binomial distribution, with constant species-specific probabilities**
#' ##-----------------------------------
#' # Manually supply basis function coefficients 
#' spp_zeroinfl_prob <- runif(num_spp, 0, 0.5)
#' spp_basis_coefs <- matrix(rnorm(num_spp * (num_basisfunctions-1), 0, 0.1), nrow = num_spp)
#' simy <- create_CBFM_life(family = zinb2(), formula = useformula, data = dat,
#' betas = cbind(spp_intercepts, spp_slopes), basis_effects_mat = spp_basis_coefs, 
#' B_space = basisfunctions, dispparam = spp_dispersion, zeroinfl_prob = spp_zeroinfl_prob)
#' 
#' 
#' ##-----------------------------------
#' ## **Example 4b: Generate spatial multivariate count data from a zero-inflated** 
#' ## **negative binomial distribution, with species-specific probabilities that also** 
#' ## **vary with covariates**
#' ##-----------------------------------
#' # Manually supply basis function coefficients 
#' spp_zislopes <- matrix(runif(num_spp * num_X, -1, 1), nrow = num_spp)
#' spp_ziintercepts <- runif(num_spp, -2, 0)
#' spp_basis_coefs <- matrix(rnorm(num_spp * (num_basisfunctions-1), 0, 0.1), nrow = num_spp)
#' simy <- create_CBFM_life(family = zinb2(), formula = useformula, ziformula = useformula, 
#' data = dat, betas = cbind(spp_intercepts, spp_slopes), 
#' zibetas = cbind(spp_ziintercepts, spp_zislopes), basis_effects_mat = spp_basis_coefs, 
#' B_space = basisfunctions, dispparam = spp_dispersion)
#' 
#' 
#' ##-----------------------------------
#' ## **Example 5: Generate spatial multivariate count data from hurdle Poisson distribution**
#' ##-----------------------------------
#' # This can be achieved by combining the mechanisms for presence-absence and zero-trun. Poisson
#' spp_slopes_pa <- matrix(runif(num_spp * num_X, -1, 1), nrow = num_spp)
#' spp_slopes_ztp <- matrix(runif(num_spp * num_X, -1, 1), nrow = num_spp)
#' spp_intercepts_pa <- runif(num_spp, -2, 0)
#' spp_intercepts_ztp <- runif(num_spp, -2, 0)
#' 
#' true_Sigma_space_pa <- rWishart(1, num_basisfunctions+1, 
#' diag(x = 0.1, nrow = num_basisfunctions-1))[,,1]/10
#' true_Sigma_space_ztp <- rWishart(1, num_basisfunctions+1, 
#' diag(x = 0.1, nrow = num_basisfunctions-1))[,,1]/10
#' true_G_space_pa <- rWishart(1, num_spp+1, diag(x = 0.1, nrow = num_spp))[,,1] %>% 
#' cov2cor
#' true_G_space_ztp <- rWishart(1, num_spp+1, diag(x = 0.1, nrow = num_spp))[,,1] %>% 
#' cov2cor
#' 
#' # Generate spatial multivariate presence-absence data first
#' # Basis function coefficients are simulated based on the supplied values of Sigma and G 
#' simy_pa <- create_CBFM_life(family = binomial(), formula = useformula, data = dat,
#' B_space = basisfunctions, betas = cbind(spp_intercepts_pa, spp_slopes_pa),
#' Sigma = list(space = true_Sigma_space_pa), G = list(space = true_G_space_pa))
#' 
#' # Now generate count data from a truncated Poisson distribution
#' simy_ztp <- create_CBFM_life(family = ztpoisson(), formula = useformula, data = dat,
#' B_space = basisfunctions, betas = cbind(spp_intercepts_ztp, spp_slopes_ztp),
#' max_resp = 20000, Sigma = list(space = true_Sigma_space_ztp), 
#' G = list(space = true_G_space_ztp))
#' 
#' # Spatial multivariate count data from a hurdle model is then the product of the two
#' simy_hurdlep <- simy_pa$y *  simy_ztp$y
#' 
#' 
#' ##-----------------------------------
#' ## **Example 6: Generate spatial multivariate count data from hurdle NB distribution**
#' ##-----------------------------------
#' # This can be achieved by combining the mechanisms for presence-absence and zero-trun. NB
#' spp_slopes_ztnb <- matrix(runif(num_spp * num_X, -1, 1), nrow = num_spp)
#' spp_intercepts_ztnb <- runif(num_spp, -2, 0)
#' spp_dispersion <- runif(num_spp, 0, 5)
#' 
#' true_Sigma_space_ztnb <- rWishart(1, num_basisfunctions+1, 
#' diag(x = 0.1, nrow = num_basisfunctions-1))[,,1]/10
#' true_G_space_ztnb <- rWishart(1, num_spp+1, diag(x = 0.1, nrow = num_spp))[,,1] %>% 
#' cov2cor
#' 
#' # Generate spatial multivariate presence-absence data first
#' # Basis function coefficients are simulated based on the supplied values of Sigma and G 
#' simy_pa <- create_CBFM_life(family = binomial(), formula = useformula, data = dat,
#' B_space = basisfunctions, betas = cbind(spp_intercepts_pa, spp_slopes_pa),
#' Sigma = list(space = true_Sigma_space_pa), G = list(space = true_G_space_pa))
#' 
#' # Now generate count data from a truncated NB distribution
#' simy_ztnb <- create_CBFM_life(family = ztnb2(), formula = useformula, data = dat,
#' B_space = basisfunctions, betas = cbind(spp_intercepts_ztp, spp_slopes_ztp),
#' max_resp = 20000, Sigma = list(space = true_Sigma_space_ztp), dispparam = spp_dispersion, 
#' G = list(space = true_G_space_ztp))
#' 
#' # Spatial multivariate count data from a hurdle model is then the product of the two
#' simy_hurdlenb <- simy_pa$y *  simy_ztnb$y
#' }
#' 
#' @export
#' @import Matrix
#' @importFrom actuar rztpois rztnbinom
#' @importFrom mgcv gam model.matrix.gam
#' @importFrom stats model.frame model.offset rbeta rbinom rgamma rnorm rnbinom rpois plogis
#' @importFrom tweedie rtweedie
#' @md

create_CBFM_life <- function(family = binomial(), formula, ziformula = NULL, data,
                             B_space = NULL, B_time = NULL, B_spacetime = NULL,
                             offset = NULL, zioffset = NULL,
                             betas, zibetas = NULL, basis_effects_mat = NULL, 
                             Sigma = list(space = NULL, time = NULL, spacetime = NULL), 
                             G = list(space = NULL, time = NULL, spacetime = NULL), 
                             trial_size = 1, dispparam = NULL, powerparam = NULL, zeroinfl_prob = NULL, 
                             max_resp = Inf, only_y = FALSE) {
     
     formula <- .check_X_formula(formula = formula, data = as.data.frame(data))          
     tmp_formula <- as.formula(paste("response", paste(as.character(formula),collapse = " ") ) )
     nullfit <- gam(tmp_formula, data = data.frame(data, response = rnorm(nrow(data))), fit = TRUE, control = list(maxit = 1))
     X <- model.matrix(nullfit)
     formula_offset <- numeric(nrow(X))
     if(!is.null(model.offset(model.frame(nullfit))))
          formula_offset <- model.offset(model.frame(nullfit))
     
     rm(tmp_formula, nullfit)
     
     if(!is.null(zeroinfl_prob) & !is.null(ziformula))
          warning("Since ziformula has been supplied, then the values supplied to zeroinf_prob will be ignored.")
     
     if(!is.null(ziformula)) {
          if(is.null(zibetas))
               stop("If ziformula is supplied, then zibetas must also be supplied.")
          
          ziformula <- .check_X_formula(formula = ziformula, data = as.data.frame(data))          
          tmp_formula <- as.formula(paste("response", paste(as.character(ziformula),collapse = " ") ) )
          nullfit <- gam(tmp_formula, data = data.frame(data, response = rnorm(nrow(data))), fit = TRUE, control = list(maxit = 1))
          ziX <- model.matrix(nullfit)
          ziformula_offset <- numeric(nrow(ziX))
          if(!is.null(model.offset(model.frame(nullfit))))
               ziformula_offset <- model.offset(model.frame(nullfit))
          rm(tmp_formula, nullfit)
          }

     
     num_units <- nrow(X)
     num_spp <- nrow(betas)
     
     .check_family(family = family, y = matrix(1, nrow = num_units, ncol = num_spp), trial_size = trial_size) 
     
     if(is.null(basis_effects_mat)) {
          .check_B_forms(B_space = B_space, B_time = B_time, B_spacetime = B_spacetime, G = G, Sigma = Sigma, extra_check = TRUE)
          }
     
     B <- cbind(B_space, B_time, B_spacetime)
     B <- Matrix(B, sparse = TRUE)
     if(is.null(rownames(B)))
          rownames(B) <- paste0("units", 1:nrow(B))

     if(!is.null(basis_effects_mat)) {
          message("Because basis_effects_mat is supplied, inputs for Sigma and G are ignored.")
          if(ncol(B) != ncol(basis_effects_mat))
               stop("The number of columns in cbind(B_space, B_time, B_spacetime) must match the number of columns in basis_effects_mat.")
          }
     
     
     ## Generate species coefficients for basis functions, if required
     basismat_notsupplied <- is.null(basis_effects_mat)
     if(basismat_notsupplied) {
          basis_effects_mat <- NULL
     
          if(!is.null(B_space)) {
               true_cholGSigma <- kronecker(t(chol(G$space)), t(chol(Sigma$space)))
               basis_effects_mat <- cbind(basis_effects_mat, matrix(true_cholGSigma %*% rnorm(ncol(B_space)*num_spp), nrow = num_spp, byrow = TRUE))
               rm(true_cholGSigma)
               }
          if(!is.null(B_time)) {
               true_cholGSigma <- kronecker(t(chol(G$time)), t(chol(Sigma$time)))
               basis_effects_mat <- cbind(basis_effects_mat, matrix(true_cholGSigma %*% rnorm(ncol(B_time)*num_spp), nrow = num_spp, byrow = TRUE))
               rm(true_cholGSigma)
               }
          if(!is.null(B_spacetime)) {
               true_cholGSigma <- kronecker(t(chol(G$spacetime)), t(chol(Sigma$spacetime)))
               basis_effects_mat <- cbind(basis_effects_mat, matrix(true_cholGSigma %*% rnorm(ncol(B_spacetime)*num_spp), nrow = num_spp, byrow = TRUE))
               rm(true_cholGSigma)
               }
          }
          
          
     ## Generate response
     true_eta_B <- tcrossprod(B, basis_effects_mat)
     true_eta <- tcrossprod(X, betas) + true_eta_B
     if(!is.null(offset))
          true_eta <- true_eta + offset
     true_eta <- true_eta + formula_offset
     if(!is.null(ziformula)) {
          true_zieta <- tcrossprod(ziX, zibetas)
          if(!is.null(zioffset))
               true_zieta <- true_zieta + zioffset
     
          true_zieta <- true_zieta + ziformula_offset
          }
     
     sim_y <- matrix(NA, nrow = num_units, ncol = num_spp, dimnames = list(units = paste0("unit", 1:num_units), response = paste0("spp", 1:num_spp)))
     for(j in 1:num_spp) {
          if(family$family == "beta")
               sim_y[,j] <- rbeta(num_units, shape1 = plogis(true_eta[,j])*dispparam[j], shape2 = (1-plogis(true_eta[,j]))*dispparam[j])
          if(family$family == "binomial")
               sim_y[,j] <- rbinom(num_units, size = ifelse(length(trial_size) == 1, trial_size, trial_size[,j]), prob = family$linkinv(true_eta[,j]))
          if(family$family == "Gamma")
               sim_y[,j] <- rgamma(num_units, scale = exp(true_eta[,j])*dispparam[j], shape = 1/dispparam[j])
          if(family$family == "gaussian")
               sim_y[,j] <- rnorm(num_units, mean = true_eta[,j], sd = sqrt(dispparam[j]))
          if(family$family == "negative.binomial")
               sim_y[,j] <- rnbinom(num_units, mu = exp(true_eta[,j]), size = 1/dispparam[j])
          if(family$family == "poisson")
               sim_y[,j] <- rpois(num_units, lambda = exp(true_eta[,j]))
          if(family$family == "tweedie")
               sim_y[,j] <- rtweedie(num_units, mu = exp(true_eta[,j]), phi = dispparam[j], power = powerparam[j])
          
          if(family$family == "zipoisson") {
               if(is.null(ziformula)) { 
                    simz <- rbinom(num_units, size = 1, prob = zeroinfl_prob[j]) # Probability of zero-inflation
                    }
               if(!is.null(ziformula)) {
                    simz <- rbinom(num_units, size = 1, prob = plogis(true_zieta[,j])) # Probability of zero-inflation
                    }
               sim_y[,j] <- rpois(num_units, lambda = exp(true_eta[,j]) * (1-simz))
                }
          if(family$family == "zinegative.binomial") {
               if(is.null(ziformula)) { 
                    simz <- rbinom(num_units, size = 1, prob = zeroinfl_prob[j]) # Probability of zero-inflation
                    }
               if(!is.null(ziformula)) {
                    simz <- rbinom(num_units, size = 1, prob = plogis(true_zieta[,j])) # Probability of zero-inflation
                    }
               sim_y[,j] <- rnbinom(num_units, mu = exp(true_eta[,j]) * (1-simz), size = 1/dispparam[j])
                }
           
          if(family$family == "ztpoisson") {
                sim_y[,j] <- actuar::rztpois(num_units, lambda = exp(true_eta[,j])+1e-12) 
                }
           if(family$family == "ztnegative.binomial") {
                make_probs <- 1/(1 + dispparam[j]*exp(true_eta[,j]) + 1e-12)
                sim_y[,j] <- actuar::rztnbinom(num_units, prob = make_probs, size = 1/dispparam[j]) 
                rm(make_probs)
                }
          }

     if(family$family %in% c("poisson", 
                             "negative.binomial", 
                             "tweedie", 
                             "Gamma",
                             "zipoisson", 
                             "zinegative.binomial", 
                             "ztpoisson", 
                             "ztnegative.binomial")) {
          inner_counter <- 0
          while(any(sim_y > max_resp) & inner_counter < 10) {
               if(basismat_notsupplied) {
                    basis_effects_mat <- NULL
                    
                    if(!is.null(B_space)) {
                         true_cholGSigma <- kronecker(t(chol(G$space)), t(chol(Sigma$space)))
                         basis_effects_mat <- cbind(basis_effects_mat, matrix(true_cholGSigma %*% rnorm(ncol(B_space)*num_spp), nrow = num_spp, byrow = TRUE))
                         rm(true_cholGSigma)
                         }
                    if(!is.null(B_time)) {
                         true_cholGSigma <- kronecker(t(chol(G$time)), t(chol(Sigma$time)))
                         basis_effects_mat <- cbind(basis_effects_mat, matrix(true_cholGSigma %*% rnorm(ncol(B_time)*num_spp), nrow = num_spp, byrow = TRUE))
                         rm(true_cholGSigma)
                         }
                    if(!is.null(B_spacetime)) {
                         true_cholGSigma <- kronecker(t(chol(G$spacetime)), t(chol(Sigma$spacetime)))
                         basis_effects_mat <- cbind(basis_effects_mat, matrix(true_cholGSigma %*% rnorm(ncol(B_spacetime)*num_spp), nrow = num_spp, byrow = TRUE))
                         rm(true_cholGSigma)
                         }
                    
                    true_eta_B <- tcrossprod(B, basis_effects_mat)
                    true_eta <- tcrossprod(X, betas) + true_eta_B
                    if(!is.null(offset))
                         true_eta <- true_eta + offset
                    true_eta <- true_eta + formula_offset
                    if(!is.null(ziformula)) {
                         true_zieta <- tcrossprod(ziX, zibetas)
                         if(!is.null(zioffset))
                              true_zieta <- true_zieta + zioffset
                         true_zieta <- true_zieta + ziformula_offset
                         }
                    }
               
               for(j in 1:num_spp) {
                    if(family$family == "beta")
                            sim_y[,j] <- rbeta(num_units, shape1 = plogis(true_eta[,j])*dispparam[j], shape2 = (1-plogis(true_eta[,j]))*dispparam[j])
                    if(family$family == "binomial")
                            sim_y[,j] <- rbinom(num_units, size = ifelse(length(trial_size) == 1, trial_size, trial_size[,j]), prob = family$linkinv(true_eta[,j]))
                    if(family$family == "Gamma")
                            sim_y[,j] <- rgamma(num_units, scale = exp(true_eta[,j])*dispparam[j], shape = 1/dispparam[j])
                    if(family$family == "gaussian")
                            sim_y[,j] <- rnorm(num_units, mean = true_eta[,j], sd = sqrt(dispparam[j]))
                    if(family$family == "negative.binomial")
                            sim_y[,j] <- rnbinom(num_units, mu = exp(true_eta[,j]), size = 1/dispparam[j])
                    if(family$family == "poisson")
                            sim_y[,j] <- rpois(num_units, lambda = exp(true_eta[,j]))
                    if(family$family == "tweedie")
                            sim_y[,j] <- rtweedie(num_units, mu = exp(true_eta[,j]), phi = dispparam[j], power = powerparam[j])
                    
                    if(family$family == "zipoisson") {
                         if(is.null(ziformula)) { 
                              simz <- rbinom(num_units, size = 1, prob = zeroinfl_prob[j]) # Probability of zero-inflation
                              }
                         if(!is.null(ziformula)) {
                              simz <- rbinom(num_units, size = 1, prob = plogis(true_zieta[,j])) # Probability of zero-inflation
                              }
                         sim_y[,j] <- rpois(num_units, lambda = exp(true_eta[,j]) * (1-simz))
                         }
                    if(family$family == "zinegative.binomial") {
                         if(is.null(ziformula)) { 
                              simz <- rbinom(num_units, size = 1, prob = zeroinfl_prob[j]) # Probability of zero-inflation
                              }
                         if(!is.null(ziformula)) {
                              simz <- rbinom(num_units, size = 1, prob = plogis(true_zieta[,j])) # Probability of zero-inflation
                              }
                         sim_y[,j] <- rnbinom(num_units, mu = exp(true_eta[,j]) * (1-simz), size = 1/dispparam[j])
                         }

                    if(family$family == "ztpoisson")
                            sim_y[,j] <- actuar::rztpois(num_units, lambda = exp(true_eta[,j])+1e-12) 
                    if(family$family == "ztnegative.binomial") {
                         make_probs <- 1/(1 + dispparam[j]*exp(true_eta[,j]) + 1e-12)
                         sim_y[,j] <- actuar::rztnbinom(num_units, prob = make_probs, size = 1/dispparam[j]) 
                         rm(make_probs)                         
                         }
                    }
               inner_counter <- inner_counter + 1
               }
          }
          
     out <- list(y = sim_y, 
                 basis_effects_mat = basis_effects_mat, 
                 linear_predictors = true_eta, 
                 linear_predictors_B = true_eta_B)
     
     if(only_y)
          out <- sim_y
     return(out)
     }

