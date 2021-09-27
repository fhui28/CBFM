##----------------------
## Community-level Basis Function Model (CBFM): Estimation done using PQL, plus maximum Laplace approximated (restricted) log-likelihood estimation for the covariance matrices

## See Github issues but also notes for TODO
##----------------------

#' @title Community-level basis function models (CBFMs)
#' 
#' @description 
#' `r lifecycle::badge("experimental")`
#' 
#' Fits CBFMs to spatio-temporal multivariate abundance data, where the basis functions are used to account for spatio-temporal correlation within and between species. Three types of basis functions can supplied and included in conjunction with each other: 1) spatial basis fnuctions; 2) temporal basis functions; 3) spatio-temporal basis functions. For the part of the mean model corresponding to the measured covariates, CBFM currently permits both parametric terms and/or smoothing terms, where the latter makes use are included in a similar manner to [mgcv::gam()]. Estimation and inference for CBFM is based on a maximum penalized quasi-likelihood (PQL) estimation approach.
#'
#' @param y A response matrix, where each row corresponds to an observational unit i.e, a particular space-time coordinate, and each column corresponds to the species.
#' @param formula_X An object of class "formula", which represents a symbolic description of the model matrix to be created (based on using this argument along with the \code{data} argument). Note there should be nothing on the left hand side of the "~". Formulas based on generalized additive models or GAMs are permitted (at least, for the basic smoothing terms we have tried so far!); please see [mgcv::formula.gam()], [mgcv::gam.models()], [mgcv::smooth.terms()], and [mgcv::s()] for more details. 
#' @param data A data frame containing covariate information, from which the model matrix is to be created (based on this argument along with the \code{formula_X} argument). 
#' @param B_space An optional matrix of spatial basis functions to be included in the CBFM. One of \code{B_space}, \code{B_time}, or \code{B_spacetime} must be supplied. The basis function matrix may be sparse or dense in form; please see the details and examples later on for illustrations of how they can constructed.
#' @param B_time An optional of matrix of temporal basis functions to be included in the CBFM. One of \code{B_space}, \code{B_time}, or \code{B_spacetime} must be supplied. The basis function matrix may be sparse or dense in form; please see the details and examples later on for illustrations of how they can constructed.
#' @param B_spacetime An optional of matrix of spatio-temporal basis functions to be included in the CBFM e.g., formed from a tensor-product of spatial and temporal basis functions. One of \code{B_space}, \code{B_time}, or \code{B_spacetime} must be supplied. The basis function matrix may be sparse or dense in form; please see the details and examples later on for illustrations of how they can constructed.
#' @param offset A matrix of offset terms.  
#' @param ncores To speed up fitting, parallelization can be performed, in which case this argument can be used to supply the number of cores to use in the parallelization. Defaults to \code{detectCores()-1}.
#' @param family a description of the response distribution to be used in the model, as specified by a family function. Please see details below for more information on the distributions currently permitted.
#' @param trial_size Trial sizes to use for binomial distribution. This can either equal a scalar or a matrix with the same dimension as \code{y}.
#' @param dofit Should the CBFM be fitted? If set to \code{FALSE}, then the function terminates (and return nothing) immediately after copying the C++ file to the compilation directory; please see the \code{TMB_directories} argument below.
#' @param stderrors Should standard errors of the estimates be calculated? This defaults to \code{TRUE}, but can be set of \code{FALSE} if only point estimations of the regression coefficients for the covariates and basis functions are desired. Please see details later on for more information on how standard errors are constructed. 
#' @param select For cases where \code{formula_X} involves smoothing terms, setting this to \code{TRUE} adds an extra penalty to each smoothing term so that it can be penalized to zero i.e., null space penalization. Please see [mgcv::gam.selection()] and [mgcv::step.gam()] for more details, noting that its implementation for the purposes of CBFM is a *wee bit experimental*. Note this argument has no effect on any parametric terms in the model i.e., it can not shrink parametric terms to zero.  
#' @param gamma For cases where \code{formula_X} involves smoothing terms, setting this to a value greater than one leads to smoother terms i.e., increased penalization. This argument plays exactly the same role as the \code{gamma} argument in [mgcv::gam()], and we refer to the help file for more information. As with the \code{select} argument, its implementation for the purposes of CBFM is a *wee bit experimental*. Note this argument has no effect on any parametric terms or the basis functions part of the CBFM. 
#' @param start_params Starting values for the CBFM. If desired, then a list should be supplied, which must contain at least one the following terms: 
#' \itemize{
#' \item{betas: }{A matrix of starting values for the species-specific regression coefficients related to the covariates, where the number of rows is equal to the number of species.} 

#' \item{basis_effect_mat: }{A matrix of starting values for the species-specific regression coefficients related to the combined matrix of basis functions. Again, the number of rows is equal to the number of species, while the number of columns should equal to \code{ncol(B_space, B_time, B_spacetime)} (or whatever the supplied basis functions are).}

#' \item{dispparam: }{A vector of starting values for the species-specific dispersion parameters, to be used for distributions that require one.}

#' \item{powerparam: }{A vector of starting values for the species-specific power parameters, to be used for distributions that require one.}

#' \item{zeroinfl_prob: }{A vector of species-specific probabilities of zero-inflation, to be used for distributions that require one. }
#' }
#' @param TMB_directories A list with two elements, identifying the directory where TMB C++ file exists (\code{cpp}), and the directory where the corresponding compiled files to be placed (\code{compile}). Unless you really want to do some real mucking around, these should be left at their default i.e., the directory where the packages were installed locally. Please note a version of the C++ file will be copied to the \code{compile} directory.
#' @param control A list of parameters for controlling the fitting process for the "outer" PQL estimation part of the CBFM. This should be a list with the following arguments:
#' \itemize{
#' \item{maxit: }{The maximum number of iterations for the outer algorithm.} 

#' \item{optim_lower/optim_upper: }{Upper and lower box constraints when updating regression coefficients related to the basis functions. Note no constraints are put in place when updating regression coefficients related to the covariates; this are controlled internally by [mgcv::gam.control()] itself.}

#' \item{convergence_type: }{The type of means by which to assess convergence. The current options are "parameters" (default), which assesses convergence based on the mean squared error of the difference between estimated parameters from successive iterations, "linear_predictor" which assesses convergence based on the mean squared error of the difference between estimated linear predictors from successive iterations  and "logLik", which assess convergence based on how close the ratio in the PQL value between successiveiterations is to one.}

#' \item{tol: }{The tolerance value to use when assessing convergence.}

#' \item{initial_betas_dampen: }{A dampening factor which can be used to reduce the magnitudes of the starting  values obtained for the species-specific regression coefficients corresponding to the model matrix i.e., \code{betas}. To elaborate, when starting values are not supplied as part of \code{start_params}, the function will attempt to obtain starting values based on fitting a stacked species distribution model. While this generally works OK, sometimes it can lead to bad starting values for the \code{betas} due to the stacked species distribution model being severely overfitted. An *ad-hoc* fix to this is to dampen/shrink these initial values to be closer to zero, thus allowing the PQL estimation algorithm to actually "work". For instance, setting \code{initial_betas_dampen = 0.5} halves the magnitudes of the staring values for the \code{betas}, including the intercepts.}

#' \item{subsequent_betas_dampen: }{A dampening factor which can be used to reduce the magnitudes of the values obtained for the species-specific regression coefficients corresponding to the model matrix i.e., \code{betas}, during the running of the PQL estimation algorithm. To elaborate, during the PQL algorithm updates are made to the regression coefficients related to the combined matrix of basis functions, conditional on the regression coefficients corresponding to the model matrix. However, sometimes this updating can fails due to the latter producing non-sensible values to condition on e.g., due to severe overfitting in that component. If this occurs, then an *ad-hoc* second attempt is made, but conditioning instead on a dampened/shrunk set of the regression coefficients corresponding to the model matrix, which can often help. This amount of dampening is controlled by this argument. For instance, setting \code{subsequent_betas_dampen = 0.25} sets the magnitudes of the regression coefficients related to the model matrix to a quarter of their original size, including the intercepts. 
#' Note that this argument *only* comes into play when the first attempt, which can be thought of as updating with \code{subsequent_betas_dampen = 1}, to update the regression coefficients associated with the combined matrix of basis functions fails. } 

#' \item{seed: }{The seed to use for the PQL algorithm. This is only applicable when the starting values are randomly generated, which be default should not be the case.}

#' \item{trace: }{If set to \code{TRUE} or \code{1}, then information at each iteration step of the outer algorithm will be printed. }

#' \item{ridge: }{A additional ridge parameter that can be included to act as a ridge penalty when estimating the regression coefficients related to the covariates.}
#' }
#' @param Sigma_control A list of parameters for controlling the fitting process for the "inner" estimation part of the CBFM pertaining to the community-level covariance matrices of the basis function regression coefficients. This should be a list with the following arguments:
#' \itemize{
#' \item{rank: }{The rank of the community-level covariance matrices of the basis function regression coefficients. This either equals to a scalar, or a vector with length equal to how many of \code{B_space/B_time/B_spacetime} are supplied. If it is a scalar, then it is assumed that the same rank is used for all the community-level covariance matrices. The ranks should be at least equal to 2, and not larger than the number of species. Please see details below for more information.} 

#' \item{maxit: }{The maximum number of iterations for inner update of the community-level covariance matrices.} 

#' \item{tol: }{The tolerance value to use when assessing convergence. Convergence for the inner algorithm is assessed based on the norm of the difference between estimated parameters from successive iterations.} 

#' \item{method: }{The method by which to update the community-level covariance matrices. The current options are "LA" (default) which uses optimizing the Laplace approximated restricted maximum likelihood (REML), and "simple" which uses a fast large sample covariance update. *The latter is \emph{much} faster than the former, but is much less accurate and we only recommend using it for pilot testing.*} 

#' \item{trace: }{If set to \code{TRUE} or \code{1}, then information at each iteration step of the inner algorithm will be printed.}
#' }
#' @param G_control A list of parameters for controlling the fitting process for the "inner" estimation part of the CBFM pertaining to the so-called baseline between-species correlation matrices of the basis function regression coefficients. This should be a list with the following arguments:
#' \itemize{
#' \item{rank}{The rank of the between-species correlation matrices of the basis function regression coefficients. This either equals to a scalar, or a vector with length equal to how many of \code{B_space/B_time/B_spacetime} are supplied. If it is a scalar, then it is assumed that the same rank is used for all the correlation matrices. The ranks should be at least equal to 2, and not larger than the number of species. Please see details below for more information.} 

#' \item{nugget_profile: }{The sequence of values to try for calculating the nugget effect in each between-species correlation matrix. Please see details below for more information.} 

#' \item{maxit: }{The maximum number of iterations for inner update of the community-level covariance matrices.} 

#' \item{tol: }{The tolerance value to use when assessing convergence. Convergence for the inner algorithm is assessed based on the norm of the difference between estimated parameters from successive iterations.} 

#' \item{method: }{The method by which to update the correlation matrices. The current options are "LA" (default) which uses optimizing the Laplace approximated restricted maximum likelihood (REML), and "simple" which uses a fast large sample covariance update. *The latter is \emph{much} faster than the former, but is much less accurate and we only recommend using it for pilot testing.*} 

#' \item{trace: }{If set to \code{TRUE} or \code{1}, then information at each iteration step of the inner algorithm will be printed.}
#' }
#' @param k_check_control A list of parameters for controlling [mgcv::k.check()] when it is applied to CBFMs involving smoothing terms for the measured covariates i.e., when smoothing terms are involved in \code{formula_X}. Please see [mgcv::k.check()] for more details on how this test works. This should be a list with the following two arguments:
#' \itemize{
#' \item{subsample: }{If the number of observational units i.e., \code{nrow(y)} exceeds this number, then testing is done using a random sub-sample of units of this size.} 

#' \item{n.rep: }{How many re-shuffles of the residuals should be done in order to a P-value for testing. } 
#' }
#'
#'
#' @details 
#'
#' Community-level basis function models (CBFMs) are a class of spatio-temporal joint species distribution models for multivariate abundance data, which builds on the ideas of fixed rank kriging (FRK, Cressie and Johannesson, 2008; Zammit-Mangion and Cressie, 2017) and multivariate spatio-temporoal mixed models (Bradley et al., 2018) and adapts them specifically for multivariate abundance data in community ecology. CBFMs provide an alternative and not necessarily superior approach to the increasingly popular latent variable models (LVMs) approach for joint species distribution modeling, as available in a number of packages such as [Hmsc::Hmsc-package()] (Tikhonov et al., 2020), [gllvm::gllvm()] (Niku et al., 2019), and [boral::boral()] (Hui, 2016); see also Warton et al., (2015a,b), Thorson et al. (2016) and Ovaskainen and Abrego (2020) among others for introductions to the use of LVMs in community ecology.  The key difference between LVMs and CBFMs is that rather than using a small number of latent variables which are assumed to be random across observational units to induce spatio-temporal correlations within and between species, CBFMs use a larger number of spatially, temporal, and/or spatio-temporal basis functions that are specified \emph{a-priori} and remain fixed in the model. The randomness instead comes from species-specific regression coefficients related to these basis functions, which in turn induce spatio-temporal correlations within and between species. 
#' 
#' In using a basis function approach, CBFMs can thus be considered both as a type of generalized additive model (GAM, Guisan et al., 2002; Wood, 2017) and a generalized linear mixed model (GLMM, Bolker et al., 2009; Brooks et al., 2017). This in turn means CBFMs can leverage from the plethora of techniques that have been already developed for such methods, with one notable benefit being that computationally, CBFMs tend to more efficient and scale better than many existing implementations of LVMs. 
#' 
#'
#' ## Some more mathematics
#' 
#' Turning to more mathematical details, for the purposes of the package the CBFM is characterized by the following mean regression model: for observational unit \eqn{i=1,\ldots,N} and species \eqn{j=1,\ldots,m}, we have
#' 
#' \deqn{g(\mu_{ij}) = \eta_{ij} = x_i^\top\beta_j + b_i^\top a_j,}
#'
#' where \eqn{g(.)} is a known link function, \eqn{x_i} denotes a vector of predictors for unit i i.e., the i-th row from the created model matrix, \eqn{\beta_j} denotes the corresponding regression coefficients for species j, \eqn{b_i} denotes a vector of spatial and/or temporal basis functions for unit i , and \eqn{a_j} denotes the corresponding regression coefficients for species j. 
#' 
#' In the function, the vector of predictors \eqn{x_i} is created based on the \code{formula_X} and \code{data} arguments. Smoothing terms are permitted in \code{formula_X}, and these can be included in the same way as in [mgcv::gam.models()]; see also [mgcv::smooth.terms()]. Note smoothing terms in this context also permits the inclusion of (species-specific) random intercepts and slopes, through the use of the \code{s(..., bs = "re")}; please see [mgcv::random.effects()] and [mgcv::gam.vcomp()] for more details. These may be included, say, as a simple approach to account for nested sampling designs, multiple data sources/surveys etc..., although please note these random effects are specific to each species i.e., they are *not* random row effects as found in packages such as [boral::boral()] and [gllvm::gllvm()]. 
#' 
#' When smoothing terms are included in the CBFM, a check of the smooth basis dimension and whether it is adequate is also automatically performed, courtesy of the [mgcv::k.check()] function; see that function's help file as well as [mgcv::choose.k()] for more general details. Furthermore, selection of smoothing terms is also possible, using either shrinkage smoothers or null space penalization; please see [mgcv::gam.selection()] and [mgcv::step.gam()] for more details. However, we must warn the practitioner that **some of the smoothers that \code{mgcv} e.g., [mgcv::linear.functional.terms()] has available have not been fully tested for CBFM**, so some make may not work. If you encounter any problems, please post a Github issue on the CBFM repository!  
#' 
#' Next, the vector basis functions \eqn{b_i} is formed from the \code{B_space}, \code{B_time} and \code{B_spacetime} arguments. At least one of these arguments must be supplied. As an example, suppose we wish to fit a CBFM with spatial and temporal basis functions which are included in an additive manner. Then only \code{B_space} and \code{B_time} should be supplied, in which case the mean regression model for the CBFM can be rewritten as:
#' 
#' \deqn{g(\mu_{ij}) = \eta_{ij} = x_i^\top\beta_j + b_{i,space}^\top a_{j,space} + b_{i,time}^\top a_{j,time},}
#'
#' where \eqn{b_i = (b_{i,space}, b_{i,time})} and \eqn{a_j = (a_{j,space}, a_{j,time})}. If purely spatial or temporal multivariate abundance data is recorded, then one should only supply \code{B_space} and \code{B_time}, respectively, ahd the mean regression model is simplified accordingly. 
#' 
#' As another example, suppose we wish to include spatio-temporal basis functions (which are formed from a tensor-product). Then only \code{B_spacetime} should be supplied, in which case the mean regression model for the CBFM can be rewritten as:
#' 
#' \deqn{g(\mu_{ij}) = \eta_{ij} = x_i^\top\beta_j + b_{i,spacetime}^\top a_{j,spacetime},}
#'
#' where \eqn{b_i = b_{i,spacetime}} and \eqn{a_j = a_{j,spacetime}}. More details and recommendations on how to construct this basis functions (including the tensor-product mentioned above) are provided later on. 
#' 
#' Note that for zero-inflated distributions, it is the *mean of the non-zero inflated component that is modeled and not the mean of the entire distribution.*
#' 
#' ** Remark on flavors and choices of CBFMs:** It is up to the practitioner as to what 'flavor' of CBFM that wish to fit, depending on interpretation and question of interests. For instance, with spatio-temporal multivariate abundance data, one may want the separate sets of spatial and temporal basis functions to be included in an additive manner (as seen above, which in analogous to an LVM where separate spatial LVs and temporal LVMs are added together), or have a single set of spatio-temporal basis functions formed from a tensor-product say (also as seen above, which in analogous to an LVM with a set of spatio-temporal LVs), or have a combination of the two where (say) the basis functions included in \code{B_space} and \code{B_time} are accounting for correlations on a course scale while the basis functions included in \code{B_spacetime} are accounting for resolutions on a fine scale. We refer the interested reader to Thorson et al., (2016) and Thorson (2019) for examples of similar kinds of constructs and flavors within the LVM framework. 
#' 
#' We also point out that this package only implements one possible version of a wider class of CBFMs; other potentially superior versions e.g., spatial basis functions with temporally varying corresponding regressions coefficients, are possible under the CBFM framework, but are far outside the scope of this package (sorry!).      
#'    
#' In principle, it is also possible to employ a more data-driven approach such as cross-validation or information criteria to choose the "flavor" of CBFM for a particular data set, although this is not currently not explicitly implemented in the package (sorry again!). The same discourse also applies to choosing the number of basis functions to include in the arguments \code{B_space/B_time/B_spacetime}, similar to choosing the number of latent variables in a LVM, although this choice is also heavily dependent on the type of basis functions used. We refer the reader to [mgcv::choose.k()] as some of the advice provided there may be applicable to CBFMs e.g., using residual analysis to informally check whether an increase the number of spatial and/or temporal basis functions is required. Furthermore, we echo a sentiment written there (while acknowledging things are more tricky with spatial and/or temporal basis functions, as well as for discrete responses!): 
#' 
#' *"So, exact choice of \eqn{k} (the number of basis functions in our situation) is not generally critical: it should be chosen to be large enough that you are reasonably sure of having enough degrees of freedom to represent the underlying 'truth' reasonably well, but small enough to maintain reasonable computational efficiency. Clearly 'large' and 'small' are dependent on the particular problem being addressed."* 
#' 
#' 
#' In the CBFM, basis functions \eqn{b_i} are specified \emph{a-priori} and remain fixed throughout. Instead, it is the associated species-specific regression coefficients \eqn{a_j} which are assumed to be random. Specifically, in this package we assume follow a multivariate normal distribution as follows:
#'   
#' \deqn{(a_1,\ldots,a_m) \sim N(0, kronecker(G, \Sigma)),} 
#' 
#' where \eqn{G} is a so-called baseline between-species correlation matrix, \eqn{\Sigma} is the community-level covariance matrix for the basis function regression coefficients, and \eqn{kronecker(\cdot)} is the Kroneckker product operator. When multiple sets of basis functions are included, then this carries over. For instance, in the example above involving a CBFM with spatial and temporal basis functions, with only \code{B_space} and \code{B_time} supplied, then we have
#'   
#' \deqn{(a_{1,space},\ldots,a_{m,space}) \sim N(0, kronecker(G_{space}, \Sigma_{space})),} 
#' and
#' \deqn{(a_{1,time},\ldots,a_{m,time}) \sim N(0, kronecker(G_{time}, \Sigma_{time})).} 
#' 
#' Furthermore, to reduce the number of parameters needed to be estimated in both the \eqn{G}'s and \eqn{\Sigma}'s, a rank-reduced structure is adopted in both (inspired and similar to that of LVMs). Specifically, we assume \eqn{G = \Lambda_{G}\Lambda_{G}^top + \kappa_G I_m} where \eqn{\Lambda_{G}} is an \eqn{m \times d_G} loading matrix and \eqn{\kappa_G > 0} is a nugget effect, with \eqn{I_m} being an identity matrix with dimenson \eqn{m}. The quantity \eqn{d_G << m} is tha rank, and similar to LVMs the larger the rank the more flexible this structure is at capturing the baseline-between species correlations (at the cost of more parameters). Similarly, we have \eqn{\Sigma = \Lambda_{\Sigma}\Lambda_{\Sigma}^top + \kappa_{\Sigma} I_{q}}, where \eqn{\Lambda_{\Sigma}} is an \eqn{q \times d_{\Sigma}} loading matrix, and \eqn{q} is the number of basis functions included in the model. When multiple sets of basis functions are included e.g., both \code{B_space} and \code{B_time}, then rank-reduced structures are used accordingly. 
#'
#' The ranks \eqn{d_G} and \eqn{d_{\Sigma}} should be smaller than the number of species and basis functions respecitvely included in the model, and generally speaking provided the rank/s is large enough then results should not depend much on their choice. The nugget effect is included to ensure that resulting rank-reduced forms of \eqn{G} and \eqn{\Sigma} remain positive definite. Moreover they have the interpretation of adjusting for the relative strength of correlation between species, say. 
#' 
#'
#' ## Distributions
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

#' \item{\code{zipoisson()}: }{Zero-inflated Poisson distribution, noting only the log link for the Poisson part is permitted. This partial mass function of the distribution is given by \eqn{f(y) = \pi I(y=0) + (1-pi) f_{pois}(y)}, where \eqn{\pi} is the probability of being in the zero-inflation component, while \eqn{f_{pois}(y)} is the usual Poisson distribution. The mean of the Poisson distribution is modeled against covariates and basis functions, while the probability of zero-inflation is a single, species-specific quantity that is estimated.}

#' \item{\code{zinb2()}: }{Zero-inflated negative binomial distribution, noting only the log link for the negative binomial part is permitted. The partial mass function of the distribution is given by \eqn{f(y) = \pi I(y=0) + (1-pi) f_{NB}(y)}, where \eqn{\pi} is the probability of being in the zero-inflation component, while \eqn{f_{NB}(y)} is the usual negative binomial distribution. The mean of the negative binomial distribution is modeled against covariates and basis functions, while the probability of zero-inflation is a single, species-specific quantity that is estimated.}
#' }
#' 
#' 
#' ## Constructing basis functions
#' 
#' The CBFM approach relies on the inclusion of the \emph{pre-specified} spatial, temporal, and/or spatio-temporal basis functions to account for spatio-temporal correlations within and between species (see Hefley et al., 2017, for a general overview of using basis functions to model autocorrelation in ecological data). Currently, the package does not provide default arguments to use for this, and this is deliberately the case as we wish to compel the practitioner to work and think a bit harder on designing the right basis functions for use when CBFMs to their particular analysis.
#' 
#' At the same time, it would be remiss not to at least provide some brief recommendations based on previous experience, and we do so below. Please also see the examples later on for some more concrete applications.
#' \describe{
#' \item{\code{B_space}: }{We have found that the multi-resolution thin-plate spline basis functions (Tzeng and Huang, 2018), as implemented in [autoFRK::mrts()], work fairly well here as spatial basis functions. They are simple to use and require the user to only supply the number of basis functions, which itself is tied to the resolution at which the practitioner wants to model their spatial correlations. For spatial multivariate abundance data we have usually found that 50 or less spatial basis functions of such type are required.
#' 
#' Another option for spatial basis functions is to use [FRK::auto_basis()], which produces basis functions that are sparse in design but consequently require \emph{many} more in number compared to the multi-resolution thin-plate splines mentioned above. This approach is more customizable however, with the choice of resolutions, basis function centers, and aperture among other choices; please see Zammit-Mangion and Cressie (2017) and Wilke et al. (2019) for more details.}
#' \item{\code{B_time}: }{Both of the approaches mentioned above for \code{B_space} can also be applied here, although with temporal basis functions we have generally found the approach implemented in [FRK::auto_basis()] to work satisfactorily in many cases, given their customizability and sparsity (local support). The multi-resolution thin-plate spline basis functions approach, when applied solely in the 1-D temporal dimension, can produce long-term temporal trends that may be undesirable.}
#' \item{\code{B_spacetime}: }{A general and most parsimonious starting point for constructing spatio-temporal basis functions is to make use of a tensor-product form (analogous to [mgcv::te()]). That is, after constructing a set of spatial and a set of temporal basis functions, we can use the [tensorproduct()] function to construct the tensor-product and include them \code{B_spacetime}. 
#' 
#' It is recommended that both 'ingredient' basis functions used in the tensor-product are sparse in design to facilitate computation e.g., as implemented in [FRK::auto_basis()]; see the [FRK::FRK-package()] package as well as Wilke et al. (2019) for some examples. Also, we recommend you do not use these same 'ingredient' basis functions in the \code{B_space} and \code{B_time} arguments, as this may lead to overfitting. Put another way, and as hinted at previously, if \code{B_spacetime} is supplied at the same time as either \code{B_space} and/or \code{B_time} is supplied, then they should generally be constructed to act different resolutions of the spatio-temporal correlation.} 
#' }
#' 
#' ## A note on estimation and inference
#' 
#' As mentioned above, because CBFMs uses a basis function approach to model spatio-temporal correlations between and within species, then they can be thought of as a type of GAM. Similar to a common implementation of GAMs then, this package uses a maximized penalized quasi-likelihood (PQL) approach for estimation and inference (Breslow and Clayton, 1993; Wood, 2017), while the baseline between-response correlation and community-level covariance matrices are estimated by maximum Laplace approximated residual maximum likelihood (REML) estimation (Wood, 2011). Currently, CBFM makes use of both the machinery available in the [mgcv] package (Wood, 2017) as well as that of Template Model Builder (TMB, Kristensen et al., 2016) to facilitate this. 
#' 
#' If \code{start_params} is not supplied, then CBFM attempts to obtain starting values based on fitting an appropriate stacked species distribution model. This generally works OK, but can sometimes fail badly e.g., if the stacked species distribution model severely overfits for one or more species. A tell-tale sign of when it occurs is if from the returned CBFM fit, the estimates of regression coefficients corresponding to the spatial and/or temporal basis functions i.e., \code{basis_effects_mat}, are extremely close to zero for these problematic species. There are no easy, principled solutions for such situations (as it may reflect an underlying intriguing feature of the proposed model for the predictors, data, or it may genuinely be that the stacked species distribution model is already fitting incredibly well!). One *ad-hoc* fix is available through \code{control$initial_betas_dampen}, but it is not guaranteed to work.  
#' 
#' Standard errors and resulting techniques like confidence intervals are based on the approximate large sample distribution of the regression coefficients, and use the so-called Bayesian posterior covariance matrix for the coefficients, similar to (but not as sophisticated as!) what is provided  [mgcv::summary.gam()]. Please note that **all standard errors and thus inference are currently computed without considering uncertainty in estimation of covariance \eqn{\Sigma} and correlation matrices \eqn{G}. They can lead to standard errors that are potentially too small, so please keep this in mind.** 
#' 
#' Also, the current estimation approach **does not provide uncertainty quantification of \eqn{\Sigma} and \eqn{G}**, does not provide uncertainty estimates in the smoothing parameter. This is in line with the current main aims of this CBFM package, which are tailored more towards estimation and inference of regression coefficients and spatio-temporal prediction (in a relatively computationally efficient and scalable manner). Future versions of package may seek to rectify this, but for now apologies!  
#' 
#' 
#' @return An object of class "CBFM" which includes the following components, not necessarily in the order below (and as appropriate):
#' \item{call: }{The matched function call.}

#' \item{family: }{The supplied response distribution i.e., family function, to be used in the model.}

#' \item{y, data, trial_size: }{The supplied response matrix, covariate information data frame, and trial size(s).}

#' \item{formula_X: }{The supplied symbolic description of the model matrix to be created.}

#' \item{B: }{The full matrix basis functions i.e., basically the result of \code{cbind(B_space, B_time, B_spacetime)}.}

#' \item{which_B_used: }{A vector of length three, indicating which of \code{B_space, B_time, B_spacetime} was supplied. For example \code{which_B_bused = c(1,0,0)} implies only \code{B_space} was supplied.}

#' \item{num_B_space: }{The number of spatial basis functions supplied i.e., \code{ncol(B_space)}.} 

#' \item{num_B_time: }{The number of temporal basis functions supplied i.e., \code{ncol(B_time)}.} 

#' \item{num_B_spacetime: }{The number of spatio-temporal basis functions supplied i.e., \code{ncol(B_spacetime)}.} 

#' \item{num_B: }{The total number of basis functions supplied i.e., \code{ncol(cbind(B_space, B_time, B_spacetime))}.}

#' \item{converged: }{Indicates whether or not the PQL estimation algorithm converged. Note results may still be outputed even if this is \code{FALSE}.}

#' \item{logLik: }{The value of the likelihood (excluding the quadratic penalty term in the PQL) upon convergence.}

#' \item{pql_logLik: }{The value of the PQL i.e., the likelihood plus the quadratic penalty term, upon convergence.}

#' \item{deviance: }{The deviance for the fitted model. Note the deviance calculation here does *not* include the quadratic term of the PQL.}

#' \item{null_deviance: }{The null deviance i.e., deviance of a stacked model (GLM) where each species model contains only an intercept. Note the deviance calculation here does *not* include the quadratic term of the PQL}

#' \item{deviance_explained: }{The *percentage* of null deviance explained by the model. In community ecology this is typically not very high (haha!); please see [varpart()] for more capacity to perform variance partitioning in a CBFM.}

#' \item{edf/edf1: }{Matrix of estimated degrees of freedom for each model parameter in \code{formula_X}. The number of columns of the matrix should be equal to the number of species i.e., \code{ncol(y)}. Penalization means that many of these are less than one. \code{edf1} is an alternative estimate of EDF. Note these values are pulled straight from the GAM part of the estimation algorithm, and consequently may only be *very* approximate. }

#' \item{pen_edf: }{A list with each element containing a vector of the estimated degrees of freedom associated with each smoothing term in \code{formula_X}. The length of the list should be equal to the number of species i.e., \code{ncol(y)}. Note these values are pulled straight from the GAM part of the estimation algorithm, and consequently may only be *very* approximate.}

#' \item{k_check: }{A list resulting from the application of [mgcv::k.check()], used as a diagnostic test of whether the smooth basis dimension is adequate for smoothing terms included in \code{formula_X}, on a per-species basis. Please see [mgcv::k.check()] for more details on the test and the output. Note that if no smoothing terms are included in \code{formula_X}, then this will be a list of \code{NULL} elements.}

#' \item{vcomp: }{A list with length equal to \code{ncol(y)}, where each element contains a vector of the estimated variance components (as standard deviations) associated with the smoothing terms included in \code{formula_X}. This output is only really useful when one or more of the smoothing terms were included in the CBFM as species-specific intercepts/slopes (see [mgcv::random.effects()] for more details), in which case the corresponding values in \code{vcomp} are the estimated variance components (estimated standard deviations to be precise) associated with these random effects; see [mgcv::random.effects()] and [mgcv::gam.vcomp()] for more details on the one-to-one relationship between smoothing parameters in GAMs and variance components in mixed models. Note that if no smoothing terms are included in \code{formula_X}, then this will be a list of \code{NULL} elements.}

#' \item{betas: }{The estimated matrix of species-specific regression coefficients corresponding to the model matrix created. The number of rows in \code{betas} is equal to the number of species i.e., \code{ncol(y)}.}

#' \item{basis_effects_mat: }{The estimated matrix of species-specific regression coefficients corresponding to the combined matrix of basis functions. The number of rows in \code{basis_effects_mat} is equal to the number of species i.e., \code{ncol(y)}.}

#' \item{dispparam: }{The estimated vector of species-specific dispersion parameters, for distributions which require one. }

#' \item{powerparam: }{The estimated vector of species-specific power parameters, for distributions which require one. }

#' \item{zeroinfl_prob_intercept: }{The estimated vector of species-specific probabilities of zero-inflation, for distributions which require one. *Note this is presented on the logit scale*, that is the model returns \eqn{log(\pi_j/(1-\pi_j))} where \eqn{\pi_j} is the probability of zero-inflation. This is the same as the intercept term of a logistic regression model for the probabilities of zero-inflation, hence the name. }

#' \item{linear_predictor: }{The estimated matrix of linear predictors. Note that for zero-inflated distributions, the mean of the non-zero-inflated component is modeled in CBFM, and the function returns the linear predictors corresponding to this non-zero-inflated component in the CBFM. }

#' \item{fitted: }{The estimated matrix of fitted mean values. Note that for zero-inflated distributions, while the mean of the non-zero-inflated component is modeled in CBFM, the fitted values are the *actual expected mean values* i.e., it returns estimated values of \eqn{(1-\pi_j)*\mu_{ij}} where \eqn{\pi_j} is the species-specific probability of zero inflation and \eqn{\mu_{ij}} is the mean of the non-zero-inflated component.}

#' \item{Sigma_space/Loading_Sigma_space/nugget_Sigma_space: }{The estimated community-level covariance matrix/loadings/nugget effect associated with the spatial basis functions, if \code{B_space} is supplied.}

#' \item{G_space/Loading_G_space/nugget_G_space: }{The estimated baseline between species correlation matrix/loadings/nugget effect associated with the spatial basis functions, if \code{B_space} is supplied.}

#' \item{Sigma_time/Loading_Sigma_time/nugget_Sigma_time: }{The estimated community-level covariance matrix/loadings/nugget effect associated with the temporal basis functions, if \code{B_time} is supplied.}

#' \item{G_time/Loading_G_time/nugget_G_time: }{The estimated baseline between species correlation matrix/loadings/nugget effect associated with the temporal basis functions, if \code{B_time} is supplied.}

#' \item{Sigma_spacetime/Loading_Sigma_spacetime/nugget_Sigma_spacetime: }{The estimated community-level covariance matrix/loadings/nugget effect associated with the spatio-temporal basis functions, if \code{B_spacetime} is supplied.}

#' \item{G_spacetime/Loading_G_spacetime/nugget_G_spacetime: }{The estimated baseline between species correlation matrix/loadings/nugget effect associated with the spatio-temporal basis functions, if \code{B_spacetime} is supplied.}

#' \item{stderrors: }{The supplied argument for \code{stderrors} i.e., whether standard errors were calculated.}

#' \item{covar_components: }{If \code{stderrors = TRUE}, then a list containing with the following components: 
#' 1) \code{topleft}, which is a matrix corresponding to the top-left block of the full Bayesian posterior covariance matrix. The top-left block specifically relates to the regression coefficients associated with the measured predictors i.e., the covariance matrix associated with \code{object$betas}, and the species-specific zero-inflated probabilities on the logit scale if the response distribution involved one;
#' 2) \code{topright}, which is a matrix of the top-right block of the full Bayesian posterior covariance matrix. The top-right block specifically relates to the cross-covariance of the regression coefficients associated with the measured predictors (plus the species-specific zero-inflated probabilities on the logit scale) and the basis functions i.e., the cross-covariance matrix between \code{object$betas} and \code{object$basis_effects_mat}; 
#' 3) \code{bottomright}, which is a matrix containing components of the bottom-right block of the full Bayesian posterior covariance matrix. The bottom-left block specifically relates to the regression coefficients associated with the basis functions i.e., the covariance matrix associated with \code{object$basis_effects_mat}.
#' 
#' Please use the [summary.CBFM()] function to obtain standard errors and confidence interval limits in a (slightly) more user-friendly form.}
#' \item{time_taken: }{The time taken to run the PQL estimation algorithm, in seconds. This is calculated simply using differences in calls of [base::proc.time()].}
#' 
#'
#' @details # Warning
#' CBFMs are designed for \emph{spatio-temporal} multivariate abundance data, such that you can sensibly construct basis functions from the space-time coordinate of each observational unit. Please do not use them for data that are **not** spatially or temporally indexed. We recommend you fit standard LVMs in those scenarios, such that made available in [gllvm::gllvm()] and [Hmsc::sampleMcmc()].
#' 
#' Not for some distributions it is not the mean of the entire distribution which is modeled. For example, in zero-inflated distributions it is the mean of the non-zero-inflated component that is modeled with the regression model described above.
#' 
#' Not all (in fact, not many) of the smoothing available that are available in [mgcv::gam.models()] have been fully tested out, so please be aware that some make not work well if at all! 
#' 
#' Please note that all standard errors and thus inference are currently computed without considering uncertainty in estimation of covariance \eqn{\Sigma} and correlation matrices \eqn{G}, as well as the any dispersion/power paameters, analogous to [mgcv::summary.gam()]. This can lead to standard errors that are potentially too small, so please keep this in mind. Also, the current estimation approach does not provide uncertainty quantification of \eqn{\Sigma} and \eqn{G}. Indeed, the "strength" of the CBFM approach (especially with the current approach to estimation) is its competitive predictive performance relative to computation efficiency and scalability; **estimates of \eqn{\Sigma} and \eqn{G} may not be too reliable.**
#'
#' Missing values are currently not handled in any manner or form in this package (sorry!). If you do have any missing values, then the standard course of action (and which is what functions such as [stats::lm()] and [mgcv::gam()] do as a default; see also \code{options("na.action")}) is to apply [stats::na.omit()] and remove all observational units from your data with one or more missing values. This can of course also be done manually by the practitioner.  
#'
#' @author Francis K.C. Hui <fhui28@gmail.com>, Chris Haak
#'
#'
#' @references
#' Bolker, B. M., Brooks, M. E., Clark, C. J., Geange, S. W., Poulsen, J. R., Stevens, M. H. H., and White, J. S. S. (2009). Generalized linear mixed models: a practical guide for ecology and evolution. Trends in Ecology & Evolution, 24, 127-135.
#'  
#' Bradley, J. R., Holan, S. H., and Wikle, C. K. (2018). Computationally efficient multivariate spatio-temporal models for high-dimensional count-valued data (with discussion). Bayesian Analysis, 13, 253-310.
#' 
#' Breslow, N. E., and Clayton, D. G. (1993). Approximate inference in generalized linear mixed models. Journal of the American statistical Association, 88, 9-25.
#' 
#' Brooks, M. E., Kristensen, K., Van Benthem, K. J., Magnusson, A., Berg, C. W., Nielsen, A., and Bolker, B. M. (2017). glmmTMB balances speed and flexibility among packages for zero-inflated generalized linear mixed modeling. The R journal, 9, 378-400.
#' 
#' Cressie, N., and Johannesson, G. (2008). Fixed rank kriging for very large spatial data sets. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 70, 209-226.
#' 
#' Guisan, A., Edwards Jr, T. C., and Hastie, T. (2002). Generalized linear and generalized additive models in studies of species distributions: setting the scene. Ecological modelling, 157, 89-100.
#' 
#' Hefley, T. J., Broms, K. M., Brost, B. M., Buderman, F. E., Kay, S. L., Scharf, H. R., and Hooten, M. B. (2017). The basis function approach for modeling autocorrelation in ecological data. Ecology, 98, 632-646.
#' 
#' Hui, F. K. C. (2016). boral-Bayesian ordination and regression analysis of multivariate abundance data in R. Methods in Ecology and Evolution, 7, 744-750.
#' 
#' Kristensen, K., Nielsen, A., Berg, C. W., Skaug, H., and Bell, B. M. (2016). TMB: Automatic Differentiation and Laplace Approximation. Journal of Statistical Software, 70, 1-21.
#' 
#' Niku, J., Hui, F. K.C., Taskinen, S., and Warton, D. I. (2019). gllvm: Fast analysis of multivariate abundance data with generalized linear latent variable models in R. Methods in Ecology and Evolution, 10, 2173-2182.
#' 
#' Ovaskainen, O., and Abrego, N. (2020). Joint species distribution modelling: with applications in R. Cambridge University Press.
#'
#' Thorson, J. T., Ianelli, J. N., Larsen, E. A., Ries, L., Scheuerell, M. D., Szuwalski, C., and Zipkin, E. F. (2016). Joint dynamic species distribution models: a tool for community ordination and spatio-temporal monitoring. Global Ecology and Biogeography, 25, 1144-1158.
#'
#' Thorson, J. T. (2019). Guidance for decisions using the Vector Autoregressive Spatio-Temporal (VAST) package in stock, ecosystem, habitat and climate assessments. Fisheries Research, 210, 143-161.
#'  
#' Tikhonov, G., Opedal, O. H., Abrego, N., Lehikoinen, A., de Jonge, M. M., Oksanen, J., and Ovaskainen, O. (2020). Joint species distribution modelling with the R-package Hmsc. Methods in ecology and evolution, 11, 442-447.
#' 
#' Tzeng, S., and Huang, H. C. (2018). Resolution adaptive fixed rank kriging. Technometrics, 60, 198-208.
#' 
#' Wikle, C. K., Zammit-Mangion, A., and Cressie, N. (2019). Spatio-temporal Statistics with R. CRC Press.
#' 
#' Warton, D. I., Blanchet, F. G., O'Hara, R. B., Ovaskainen, O., Taskinen, S., Walker, S. C., and Hui, F. K. C. (2015). So many variables: joint modeling in community ecology. Trends in Ecology and Evolution, 30, 766-779.
#'
#' Warton, D. I., Blanchet, F. G., O'Hara, R., Ovaskainen, O., Taskinen, S., Walker, S. C., and Hui, F. K. C. (2016). Extending Joint Models in Community Ecology: A Response to Beissinger et al. Trends in ecology & evolution, 31, 737-738.
#'
#' Wood, S. N. (2011). Fast stable restricted maximum likelihood and marginal likelihood estimation of semiparametric generalized linear models. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 73, 3-36.
#' 
#' Wood, S. N. (2017). Generalized additive models: An introduction with R. CRC press.
#' 
#' Zammit-Mangion, A., and Cressie, N. (2017). FRK: An R package for spatial and spatio-temporal prediction with large datasets. arXiv preprint arXiv:1705.08105.
#'
#'
#' @seealso [fitted.CBFM()] for extracting the fitted values from a CBFM fit, [influence_CBFM()] for calculating some basic influence measures from a CBFM fit, [ordinate.CBFM()] for an *ad-hoc* approach to constructing spatio-temporal ordinations from a CBFM fit, [plot.CBFM()] for basic residual diagnostics from a CBFM fit, [predict.CBFM()] for constructing predictions from a CBFM fit, [residuals.CBFM()] for calculating residuals from a CBFM fit, [simulate.CBFM()] for simulating spatio-temporal multivariate abundance data from a CBFM fit, [summary.CBFM()] for summaries including standard errors and confidence intervals, and [varpart()] for variance partitioning of a CBFM fit.
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
#' 
#' ##------------------------------
#' ## **Example 1a: Fitting a CBFM to spatial multivariate presence-absence data**
#' ## simulated from a spatial latent variable model
#' ## Please note the data generation process (thus) differs from CBFM.
#' ##------------------------------
#' set.seed(2021)
#' num_sites <- 1000 # 500 (units) sites for training set + 500 sites for testing.
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
#' dat_test <- dat[501:1000,]
#' simy_train <- simy[1:500,]
#' simy_test <- simy[501:1000,]
#' rm(X, mm, spp_loadings, true_lvs, xy, simy, dat)
#' 
#' 
#' # Fit stacked GLM as a baseline
#' fitstacked <- manyglm(simy_train ~ temp + depth + chla + O2, family = binomial(), data = dat_train)
#' 
#' 
#' # Set up spatial basis functions for CBFM -- Most practitioners will start here! 
#' # We will also use this basis functions in some later examples
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
#' # Fit CBFM 
#' tic <- proc.time()
#' useformula <- ~ temp + depth + chla + O2
#' fitcbfm <- CBFM(y = simy_train, formula_X = useformula, data = dat_train, 
#' B_space = train_basisfunctions, family = binomial(), control = list(trace = 1))
#' toc <- proc.time()
#' toc - tic
#' 
#' summary(fitcbfm) %>% 
#' str
#' 
#' 
#' # Calculate predictions onto test dataset
#' predictions_stacked <- predict(fitstacked, newdata = dat_test, type = "response")
#' predictions_cbfm <- predict(fitcbfm, newdata = dat_test, type = "response", 
#' new_B_space = test_basisfunctions)
#' 
#' # Evaluation predictions
#' # Tjur R-squared across species
#' tjurR2 <- data.frame(
#' stacked = sapply(1:num_spp, function(j) { 
#' m1 <- predictions_stacked[which(simy_test[,j] > 0),j] %>%
#' mean(na.rm = TRUE)
#' m0 <- predictions_stacked[which(simy_test[,j] == 0),j] %>%
#' mean(na.rm = TRUE)
#' m1 - m0     
#' }),
#' cbfm = sapply(1:num_spp, function(j) { 
#' m1 <- predictions_cbfm[which(simy_test[,j] > 0),j] %>%
#' mean(na.rm = TRUE)
#' m0 <- predictions_cbfm[which(simy_test[,j] == 0),j] %>%
#' mean(na.rm = TRUE)
#' m1 - m0     
#' })
#' )
#' 
#' boxplot(tjurR2, main = "Tjur-R2", names = c("Stacked GLM", "CBFM"))
#' 
#' ggplot(tjurR2, aes(x = stacked, y = cbfm)) +
#' geom_point() +
#' geom_abline(intercept = 0, slope = 1, linetype = 2) +
#' labs(x = "Stacked SDM", y = "CBFM", main = "Tjur-R2") +
#' theme_bw()
#' 
#' # AUC across species
#' aucs <- data.frame(
#' stacked = sapply(1:num_spp, function(j) { 
#' pred <- prediction(predictions_stacked[,j], labels = simy_test[,j]) %>%
#' performance(measure = "auc")  
#' pred@y.values[[1]]
#' }),
#' cbfm = sapply(1:num_spp, function(j) { 
#' pred <- prediction(predictions_cbfm[,j], labels = simy_test[,j]) %>%
#' performance(measure = "auc") 
#' pred@y.values[[1]]
#' })
#' )
#' 
#' boxplot(aucs, main = "AUC", names = c("Stacked GLM", "CBFM"))
#' 
#' ggplot(aucs, aes(x = stacked, y = cbfm)) +
#' geom_point() +
#' geom_abline(intercept = 0, slope = 1, linetype = 2) +
#' labs(x = "Stacked SDM", y = "CBFM", main = "AUC") +
#' theme_bw()
#' 
#' 
#' 
#' ##------------------------------
#' ## **Example 1b: Repeat Example 1a but illustrate the use of smoothing terms in CBFM**
#' ## Since the true model only involves parametric terms, then we do not expect its performance
#' ## to be as good as assuming the right form for the mean model.
#' ## It is purely for illustration purposes. 
#' ## Please note this will take a while to run...get a cup of tea and stretch your legs! 
#' ##------------------------------
#' # Set up spatial basis functions for CBFM -- Most practitioners will start here! 
#' # This is the same set up as Example 1a
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
#' # Fit CBFM 
#' tic <- proc.time()
#' useformula <- ~ s(temp) + s(depth) + s(chla) + s(O2)
#' fitcbfm_gam <- CBFM(y = simy_train, formula_X = useformula, 
#' data = dat_train, B_space = train_basisfunctions, family = binomial(), control = list(trace = 1))
#' toc <- proc.time()
#' toc-tic
#' 
#' summary(fitcbfm_gam) %>% 
#' str
#' 
#' 
#' # Calculate predictions onto test dataset
#' predictions_cbfm_gam <- predict(fitcbfm_gam, newdata = dat_test, type = "response", 
#' new_B_space = test_basisfunctions)
#' 
#' # Evaluation predictions
#' # Tjur R-squared across species
#' tjurR2$cbfm_gam = sapply(1:num_spp, function(j) { 
#' m1 <- predictions_cbfm_gam[which(simy_test[,j] > 0),j] %>%
#' mean(na.rm = TRUE)
#' m0 <- predictions_cbfm_gam[which(simy_test[,j] == 0),j] %>%
#' mean(na.rm = TRUE)
#' m1 - m0     
#' })
#' 
#' boxplot(tjurR2, main = "Tjur-R2", names = c("Stacked GLM", "CBFM", "CBFM (GAM)"))
#' 
#' ggplot(tjurR2, aes(x = stacked, y = cbfm_gam)) +
#' geom_point() +
#' geom_abline(intercept = 0, slope = 1, linetype = 2) +
#' labs(x = "Stacked SDM", y = "CBFM (GAM)", main = "Tjur-R2") +
#' theme_bw()
#' 
#' # AUC across species
#' aucs$cbfm_gam <- sapply(1:num_spp, function(j) { 
#' pred <- prediction(predictions_cbfm_gam[,j], labels = simy_test[,j]) %>%
#' performance(measure = "auc") 
#' pred@y.values[[1]]
#' })
#' 
#' boxplot(aucs, main = "AUC", names = c("Stacked GLM", "CBFM", "CBFM (GAM)"))
#' 
#' ggplot(aucs, aes(x = stacked, y = cbfm_gam)) +
#' geom_point() +
#' geom_abline(intercept = 0, slope = 1, linetype = 2) +
#' labs(x = "Stacked SDM", y = "CBFM (GAM)", main = "AUC") +
#' theme_bw()
#'
#' 
#' 
#' ##------------------------------
#' ## **Example 1c: Repeat Example 1a but illustrate applications to Poisson count data**
#' ##------------------------------
#' # Simulate spatial multivariate abundance data
#' simy <- matrix(rpois(num_sites * num_spp, lambda = exp(eta)), nrow = num_sites)
#' 
#' # Form training and test sets
#' simy_train <- simy[1:500,]
#' simy_test <- simy[501:1000,]
#' 
#' 
#' # Fit stacked GLM as a baseline
#' fitstacked <- manyglm(simy_train ~ temp + depth + chla + O2, family = "poisson", 
#' data = dat_train)
#' 
#' 
#' # Set up spatial basis functions for CBFM -- Most practitioners will start here! 
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
#' # Fit Poisson CBFM 
#' tic <- proc.time()
#' useformula <- ~ temp + depth + chla + O2
#' fitcbfm <- CBFM(y = simy_train, formula_X = useformula, data = dat_train, 
#' B_space = train_basisfunctions, family = poisson(), control = list(trace = 1))
#' toc <- proc.time()
#' toc - tic
#' 
#' summary(fitcbfm) %>% 
#' str
#' 
#' 
#' # Calculate predictions onto test dataset
#' predictions_stacked <- predict(fitstacked, newdata = dat_test, type = "response")
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
#' preddeviance <- data.frame(
#' stacked = sapply(1:num_spp, function(j) { 
#' -2*sum(dpois(simy_test[,j], lambda = predictions_stacked[,j], log = TRUE))
#' }),
#' cbfm = sapply(1:num_spp, function(j) { 
#' -2*sum(dpois(simy_test[,j], lambda = predictions_cbfm[,j], log = TRUE))
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
#' ##------------------------------
#' ## **Example 1d: Repeat Example 1a but illustrate applications to negative binomial count data** 
#' ##------------------------------
#' # Simulate spatial multivariate abundance data
#' spp_dispersion <- runif(num_spp)
#' simy <- matrix(rnbinom(num_sites * num_spp, mu = exp(eta), 
#' size = matrix(1/spp_dispersion, nrow = num_sites, ncol = num_spp, byrow = TRUE)),
#' nrow = num_sites)
#' 
#' # Form training and test sets
#' simy_train <- simy[1:500,]
#' simy_test <- simy[501:1000,]
#' 

#' # Fit stacked GLM as a baseline
#' fitstacked <- manyglm(simy_train ~ temp + depth + chla + O2, family = "negative.binomial",
#' data = dat_train)
#' 
#' 
#' # Set up spatial basis functions for CBFM -- Most practitioners will start here! 
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
#' # Fit negative binomial CBFM
#' # Please make sure you set up spatial basis functions for CBFM, we do so in Example 1a
#' # Please note the negative binomial distribution is not necessarily required for 
#' # overdispersed count data, since the latent variables can to some degree account for
#' # overdispersion. 
#' tic <- proc.time()
#' useformula <- ~ temp + depth + chla + O2
#' fitcbfm <- CBFM(y = simy_train, formula_X = useformula, data = dat_train, 
#' B_space = train_basisfunctions, family = nb2(), control = list(trace = 1))
#' toc <- proc.time()
#' toc - tic
#' 
#' summary(fitcbfm) %>% 
#' str
#' 
#' 
#' # Calculate predictions onto test dataset
#' predictions_stacked <- predict(fitstacked, newdata = dat_test, type = "response")
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
#' preddeviance <- data.frame(
#' stacked = sapply(1:num_spp, function(j) { 
#' -2*sum(dnbinom(simy_test[,j], mu = predictions_stacked[,j], size = 1/fitstacked$phi[j], 
#' log = TRUE))
#' }),
#' cbfm = sapply(1:num_spp, function(j) { 
#' -2*sum(dnbinom(simy_test[,j], mu = predictions_cbfm[,j], size = 1/fitcbfm$dispparam[j], 
#' log = TRUE))
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
#' ##------------------------------
#' ## **Example 1e: Repeat Example 1a but illustrate applications to ZIP count data** 
#' ##------------------------------
#' library(pscl)
#' 
#' # Probability of zero-inflation 
#' spp_zeroinfl_prob <- runif(num_spp, 0.1, 0.5) 
#' 
#' # Simulate spatial multivariate abundance data
#' component_ind <- matrix(rbinom(num_sites * num_spp, size = 1, 
#' prob = matrix(spp_zeroinfl_prob, num_sites, num_spp, byrow = TRUE)), num_sites, num_spp)
#' simy <- matrix(rpois(num_sites * num_spp, lambda = exp(eta) * (1-component_ind)), 
#' num_sites, num_spp)
#' rm(component_ind)
#' 
#' # Form training and test sets
#' simy_train <- simy[1:500,]
#' simy_test <- simy[501:1000,]
#' 
#' 
#' # Fit stacked ZIP regression models as a baseline
#' fitstacked <- NULL 
#' for(j in 1:num_spp) {
#' fitstacked[[j]] <- zeroinfl(resp ~ temp + depth + chla + O2 | 1, 
#' data = data.frame(resp = simy_train[,j], dat_train))
#' }
#' 
#' 
#' # Set up spatial basis functions for CBFM -- Most practitioners will start here! 
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
#' # Fit zero-inflated Poisson CBFM
#' tic <- proc.time()
#' useformula <- ~ temp + depth + chla + O2
#' fitcbfm <- CBFM(y = simy_train, formula_X = useformula, data = dat_train, 
#' B_space = train_basisfunctions, family = zipoisson(), control = list(trace = 1))
#' toc <- proc.time()
#' toc - tic
#' 
#' summary(fitcbfm) %>% 
#' str
#' 
#' 
#' # Calculate predictions onto test dataset
#' predictions_stacked <- sapply(1:num_spp, function(j) predict(fitstacked[[j]], 
#' newdata = dat_test, type = "response"))
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
#' # Need to define density of zero-inflated Poisson distribution first (or get it from a package)
#' dlzipois <- function(y, lambda, p0) {
#' logp <- log(1-p0) + dpois(y, lambda = lambda, log=TRUE)
#' logp[y == 0] <- log(exp(logp[y == 0]) + p0) 
#' return(logp)
#' }
#' preddeviance <- data.frame(
#' stacked = sapply(1:num_spp, function(j) { 
#' -2*sum(dlzipois(simy_test[,j], lambda = predictions_stacked[,j], 
#' p0 = plogis(fitstacked[[j]]$coefficients$zero)))
#' }),
#' cbfm = sapply(1:num_spp, function(j) { 
#' -2*sum(dlzipois(simy_test[,j], lambda = predictions_cbfm[,j], 
#' p0 = plogis(fitcbfm$zeroinfl_prob_intercept[j])))
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
#' 
#' ##------------------------------
#' ## **Example 1f: Repeat Example 1a but illustrate applications to ZINB count data**
#' ##------------------------------
#' # Probability of zero-inflation 
#' spp_zeroinfl_prob <- runif(num_spp, 0.1, 0.5) 
#' spp_dispersion <- runif(num_spp)
#' 
#' # Simulate spatial multivariate abundance data
#' # Note the deliberate "+2" on the linear predictor: This creates data that is a bit more 
#' # clearly overdispersed and zero-inflated...
#' component_ind <- matrix(rbinom(num_sites * num_spp, size = 1, 
#' prob = matrix(spp_zeroinfl_prob, num_sites, num_spp, byrow = TRUE)), num_sites, num_spp)
#' simy <- matrix(rnbinom(num_sites * num_spp, mu = exp(eta+2) * (1-component_ind),
#' size = matrix(1/spp_dispersion, nrow = num_sites, ncol = num_spp, byrow = TRUE)),
#' num_sites, num_spp)
#' rm(component_ind)
#' 
#' # Form training and test sets
#' simy_train <- simy[1:500,]
#' simy_test <- simy[501:1000,]
#' 
#' 
#' # Fit stacked ZIP regression models as a baseline
#' fitstacked <- NULL 
#' for(j in 1:num_spp) {
#' fitstacked[[j]] <- zeroinfl(resp ~ temp + depth + chla + O2 | 1, 
#' dist = "negbin", data = data.frame(resp = simy_train[,j], dat_train))
#' }
#' 
#' 
#' # Set up spatial basis functions for CBFM -- Most practitioners will start here! 
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
#' # Fit zero-inflated negative binomial CBFM
#' tic <- proc.time()
#' useformula <- ~ temp + depth + chla + O2
#' fitcbfm <- CBFM(y = simy_train, formula_X = useformula, data = dat_train, 
#' B_space = train_basisfunctions, family = zinb2(), control = list(trace = 1))
#' toc <- proc.time()
#' toc - tic
#' 
#' summary(fitcbfm) %>% 
#' str
#' 
#' 
#' # Calculate predictions onto test dataset
#' predictions_stacked <- sapply(1:num_spp, function(j) predict(fitstacked[[j]], 
#' newdata = dat_test, type = "response"))
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
#' # Need to define density of zero-inflated Poisson distribution first (or get it from a package)
#' dlzinb <- function(y, lambda, p0, phi) {
#' logp <- log(1-p0) + dnbinom(y, mu = lambda, size = 1/phi, log=TRUE)
#' logp[y == 0] <- log(exp(logp[y == 0]) + p0) 
#' return(logp)
#' }
#' preddeviance <- data.frame(
#' stacked = sapply(1:num_spp, function(j) { 
#' -2*sum(dlzinb(simy_test[,j], lambda = predictions_stacked[,j], 
#' p0 = plogis(fitstacked[[j]]$coefficients$zero),
#' phi = 1/fitstacked[[j]]$theta))
#' }),
#' cbfm = sapply(1:num_spp, function(j) { 
#' -2*sum(dlzinb(simy_test[,j], lambda = predictions_cbfm[,j], 
#' p0 = plogis(fitcbfm$zeroinfl_prob_intercept[j]),
#' phi = fitcbfm$dispparam[j]))
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
#' 
#' ##------------------------------
#' ## **Example 1g Repeat Example 1a but illustrate applications to biomass**
#' ##------------------------------
#' library(tweedie)
#' library(statmod)
#' 
#' spp_dispersion <- runif(num_spp, 0, 5)
#' simy <- matrix(rtweedie(num_sites * num_spp, mu = exp(eta), 
#' phi = matrix(spp_dispersion, nrow = num_sites, ncol = num_spp, byrow = TRUE),
#' power = 1.6),
#' nrow = num_sites)
#' 
#' # Form training and test sets
#' simy_train <- simy[1:500,]
#' simy_test <- simy[501:1000,]
#' 
#' # Fit stacked GLM as a baseline
#' fitstacked <- lapply(1:num_spp, function(j) {
#' # Note power parameter is assumed to be known for stacked model
#' glm(simy_train[,j] ~ temp + depth + chla + O2, 
#' family = tweedie(var.power = 1.6, link.power = 0), data = dat_train)
#' })
#' 
#' 
#' # Set up spatial basis functions for CBFM -- Most practitioners will start here! 
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
#' # Fit Tweedie CBFM
#' tic <- proc.time()
#' useformula <- ~ temp + depth + chla + O2
#' fitcbfm <- CBFM(y = simy_train, formula_X = useformula, data = dat_train, 
#' B_space = train_basisfunctions, family = tweedielogfam(), control = list(trace = 1))
#' toc <- proc.time()
#' toc - tic
#' 
#' summary(fitcbfm) %>% 
#' str
#' 
#' 
#' # Calculate predictions onto test dataset
#' predictions_stacked <- sapply(1:num_spp, function(j) {
#' predict(fitstacked[[j]], newdata = dat_test, type = "response")
#' })
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
#' preddeviance <- data.frame(
#' stacked = sapply(1:num_spp, function(j) { 
#' -2*sum(log(dtweedie(simy_test[,j], mu = predictions_stacked[,j], 
#' phi = mean(fitstacked[[j]]$residuals^2), power = 1.6)))
#' }),
#' cbfm = sapply(1:num_spp, function(j) { 
#' -2*sum(log(dtweedie(simy_test[,j], mu = predictions_cbfm[,j], 
#' phi = fitcbfm$dispparam[j], power = fitcbfm$powerparam[j])))
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
#' rm(list = ls())
#' 
#' 
#' 
#' ##------------------------------
#' ## **Example 2a: Fitting a CBFM to multivariate spatial continuous data**
#' ## with a data generation pulled from the Hmsc vignette (also spatial latent variable model)
#' ## Please see Case 6 in
#' ## <https://cran.r-project.org/web/packages/Hmsc/vignettes/vignette_5_performance.pdf>
#' ##------------------------------
#' library(Hmsc)
#' library(devtools)
#' 
#' # The makedata function included in the vignette folder produces datasets 
#' # based on the Hmsc model
#' source_url("https://raw.githubusercontent.com/hmsc-r/HMSC/master/vignettes/makedata.R")
#' 
#' tmp = makedata(ns=50, ny=200, spatial=TRUE)
#' all.data=tmp[[1]]
#' all.parameters=tmp[[2]]
#' 
#' # Fit stacked GLM as a baseline
#' fitstacked <- manylm(all.data$Y ~ X.categorical + X.covariate, data = all.data$X.data)
#' 
#' 
#' # Set up spatial basis functions for CBFM -- Most practitioners will start here! 
#' num_basisfunctions <- 20 # Number of spatial basis functions to use
#' basisfunctions <- mrts(all.data$xy, num_basisfunctions) %>% 
#' as.matrix %>%
#' {.[,-(1)]} # Remove the first intercept column
#' 
#' # Fit CBFM 
#' tic <- proc.time()
#' fitcbfm <- CBFM(y = all.data$Y, formula_X = all.data$X.formula, data = all.data$X.data, 
#' B_space = basisfunctions, family = gaussian(), control = list(trace = 1))
#' toc <- proc.time()
#' toc - tic
#' 
#' summary(fitcbfm) %>% 
#' str
#' 
#' 
#' # Evaluate in-sample performance (similar to what was done in the Hmsc vignette)
#' # with the = evaluateModelFit() function 
#' predictions_stacked <- fitstacked$fitted
#' predictions_cbfm <- fitcbfm$fitted 
#' 
#' num_spp <- ncol(all.data$Y)
#' 
#' # Root mean-squared error (RMSE)
#' RMSE <- data.frame(stacked = sapply(1:num_spp, function(j) {
#' sqrt(mean((all.data$Y[,j]-predictions_stacked[,j])^2))
#' }),
#' cbfm = sapply(1:num_spp, function(j) {
#' sqrt(mean((all.data$Y[,j]-predictions_cbfm[,j])^2))
#' })
#' )
#' boxplot(RMSE, main = "RMSE", names = c("Stacked GLM", "CBFM"))
#' 
#' pearsonR2 <- data.frame(stacked = sapply(1:num_spp, function(j) {
#' out <- cor(all.data$Y[,j], predictions_stacked[,j])
#' out^2 * sign(out)
#' }),
#' cbfm = sapply(1:num_spp, function(j) {
#' out <- cor(all.data$Y[,j],predictions_cbfm[,j])
#' out^2 * sign(out)
#' })
#' )
#' boxplot(pearsonR2, main = "Pearson R-squared", names = c("Stacked GLM", "CBFM"))
#' 
#' pearsonR2 <- data.frame(stacked = sapply(1:num_spp, function(j) {
#' out <- cor(all.data$Y[,j], predictions_stacked[,j])
#' out^2 * sign(out)
#' }),
#' cbfm = sapply(1:num_spp, function(j) {
#' out <- cor(all.data$Y[,j],predictions_cbfm[,j])
#' out^2 * sign(out)
#' })
#' )
#' boxplot(pearsonR2, main = "Pearson R-squared", names = c("Stacked GLM", "CBFM"))
#' 
#' 
#' 
#' ##------------------------------
#' ## **Example 2b: Repeat example 2a with presence-absence data**
#' ## Similar to a combination of Cases 2 and 6 in
#' ## <https://cran.r-project.org/web/packages/Hmsc/vignettes/vignette_5_performance.pdf>
#' ##------------------------------
#' # Generate data
#' L1 = all.parameters$L
#' Y2 = 1*(L1 + matrix(rnorm(n = nrow(all.data$Y)*ncol(all.data$Y)), ncol = ncol(all.data$Y)) > 0)
#' all.data$Y = Y2
#' 
#' 
#' # Fit stacked GLM as a baseline
#' fitstacked <- manyglm(all.data$Y ~ X.categorical + X.covariate, data = all.data$X.data,
#' family = binomial())
#' 
#' 
#' # Set up spatial basis functions for CBFM -- Most practitioners will start here! 
#' num_basisfunctions <- 20 # Number of spatial basis functions to use
#' basisfunctions <- mrts(all.data$xy, num_basisfunctions) %>% 
#' as.matrix %>%
#' {.[,-(1)]} # Remove the first intercept column
#' 
#' # Fit CBFM 
#' # Note also that Hmsc generates and fits models assuming a probit link, 
#' # but CBFM uses a logit link
#' tic <- proc.time()
#' fitcbfm <- CBFM(y = all.data$Y, formula_X = all.data$X.formula, data = all.data$X.data, 
#' B_space = basisfunctions, family = binomial(), control = list(trace = 1))
#' toc <- proc.time()
#' toc - tic
#' 
#' 
#' # Evaluate in-sample performance (similar to what was done in the Hmsc vignette)
#' # with the = evaluateModelFit() function 
#' predictions_stacked <- fitstacked$fitted
#' predictions_cbfm <- fitcbfm$fitted 
#' 
#' num_spp <- ncol(all.data$Y)
#' 
#' # Root mean-squared error (RMSE)
#' RMSE <- data.frame(stacked = sapply(1:num_spp, function(j) {
#' sqrt(mean((all.data$Y[,j]-predictions_stacked[,j])^2))
#' }),
#' cbfm = sapply(1:num_spp, function(j) {
#' sqrt(mean((all.data$Y[,j]-predictions_cbfm[,j])^2))
#' })
#' )
#' boxplot(RMSE, main = "RMSE", names = c("Stacked GLM", "CBFM"))
#' 
#' # Tjur R-squared across species
#' tjurR2 <- data.frame(stacked = sapply(1:num_spp, function(j) { 
#' m1 <- predictions_stacked[which(all.data$Y[,j] > 0),j] %>%
#' mean(na.rm = TRUE)
#' m0 <- predictions_stacked[which(all.data$Y[,j] == 0),j] %>%
#' mean(na.rm = TRUE)
#' m1 - m0     
#' }),
#' cbfm = sapply(1:num_spp, function(j) { 
#' m1 <- predictions_cbfm[which(all.data$Y[,j] > 0),j] %>%
#' mean(na.rm = TRUE)
#' m0 <- predictions_cbfm[which(all.data$Y[,j] == 0),j] %>%
#' mean(na.rm = TRUE)
#' m1 - m0     
#' })
#' )
#' boxplot(tjurR2, main = "Tjur-R2", names = c("Stacked GLM", "CBFM"))
#' 
#' # AUC across species
#' aucs <- data.frame(
#' stacked = sapply(1:num_spp, function(j) { 
#' pred <- prediction(predictions_stacked[,j], labels = all.data$Y[,j]) %>%
#' performance(measure = "auc")  
#' pred@y.values[[1]]
#' }),
#' cbfm = sapply(1:num_spp, function(j) { 
#' pred <- prediction(predictions_cbfm[,j], labels = all.data$Y[,j]) %>%
#' performance(measure = "auc") 
#' pred@y.values[[1]]
#' })
#' )
#' boxplot(aucs, main = "AUC", names = c("Stacked GLM", "CBFM"))
#' 
#' set.seed(NULL)
#' rm(list = ls())
#' 
#' 
#' 
#' ##------------------------------
#' ## **Example 3a: Fitting an additive CBFM to spatio-temporal multivariate presence-absence data**
#' ## simulated from a spatio-temporal latent variable model
#' ## Please note the data generation process (thus) differs from CBFM.
#' ## Also, the additive CBFM might take a while to fit...grab a cup of tea or coffee!
#' ##------------------------------
#' set.seed(2021)
#' num_sites <- 1000 # 500 (units) sites for training set + 500 sites for testing.
#' num_spp <- 50 # Number of species
#' num_X <- 4 # Number of regression slopes
#' 
#' spp_slopes <- matrix(runif(num_spp * num_X, -1, 1), nrow = num_spp)
#' spp_intercepts <- runif(num_spp, -2, 0)
#' 
#' # Simulate spatio-temporal coordinates and environmental covariate components
#' # Note we assume that each site is only visited once, but the code below can be adapted to 
#' # when the same sites are repeatedly visited
#' # We will also use this information in examples below
#' xy <- data.frame(x = runif(num_sites, 0, 5), y = runif(num_sites, 0, 5))
#' X <- rmvnorm(num_sites, mean = rep(0,4)) 
#' colnames(X) <- c("temp", "depth", "chla", "O2")
#' dat <- data.frame(xy, time = sort(runif(1000, 0, 10)) , X)
#' mm <- model.matrix(~ temp + depth + chla + O2 - 1, data = dat) %>% 
#' scale %>% 
#' as.matrix
#' 
#' # Simulate latent variable component
#' # We will also use this information in examples below
#' true_space_lvs <- RFsimulate(model = RMexp(var = 1, scale = 2), x = xy$x, y = xy$y, 
#' n = 2)@data %>% 
#' as.matrix
#' true_time_lvs <- RFsimulate(model = RMgauss(var = 1, scale = 1), x = dat$time, 
#' n = 2)@data %>% 
#' as.matrix
#' spp_space_loadings <- matrix(runif(num_spp * 2, -1, 1), nrow = num_spp) 
#' spp_time_loadings <- matrix(runif(num_spp * 2, -0.5, 0.5), nrow = num_spp) 
#' 
#' # Simulate spatial multivariate abundance data (presence-absence)
#' # We will also use this information in examples below
#' eta <- tcrossprod(cbind(1,mm), cbind(spp_intercepts,spp_slopes)) + 
#' tcrossprod(true_space_lvs, spp_space_loadings) +
#' tcrossprod(true_time_lvs, spp_time_loadings)
#' simy <- matrix(rbinom(num_sites * num_spp, size = 1, 
#' prob = plogis(eta)), nrow = num_sites)
#' 
#' # Form training and test sets
#' simy_train <- simy[1:500,]
#' simy_test <- simy[501:1000,]
#' dat_train <- dat[1:500,]
#' dat_test <- dat[501:1000,]
#' rm(X, eta, mm, spp_space_loadings, spp_time_loadings, true_space_lvs, true_time_lvs,
#' xy, simy)
#' 
#' 
#' # Fit stacked GLM as a baseline
#' fitstacked <- manyglm(simy_train ~ temp + depth + chla + O2, family = binomial(), data = dat_train)
#' 
#' 
#' # Set up spatial basis functions for CBFM -- Most practitioners will start here! 
#' num_space_basisfunctions <- 20 # Number of spatial basis functions to use
#' # Training set basis functions
#' train_space_basisfunctions <- mrts(dat_train[,c("x","y")], num_space_basisfunctions) %>% 
#' as.matrix %>%
#' {.[,-(1)]} # Remove the first intercept column
#' # Testing set basis functions
#' test_space_basisfunctions <- mrts(dat_train[,c("x","y")], num_space_basisfunctions) %>% 
#' predict(newx = dat_test[,c("x","y")]) %>% 
#' as.matrix %>%
#' {.[,-c(1)]} 
#' 
#' # Training and test temporal basis functions
#' num_time_basisfunctions <- 10 # Number of temporal basis functions to use
#' time_knots <- seq(0, 10, length = num_time_basisfunctions)
#' time_basisfunctions <- local_basis(manifold = real_line(), loc = as.matrix(time_knots),
#' scale = rep(2, length(time_knots)), type = "bisquare")
#' time_basisfunctions <- eval_basis(time_basisfunctions, s = as.matrix(dat$time)) 
#' train_time_basisfunctions <- time_basisfunctions[1:500,] 
#' test_time_basisfunctions <- time_basisfunctions[501:1000,] 
#' rm(time_basisfunctions, time_knots)
#' 
#' # Fit CBFM with Additive spatial and temporal basis functions 
#' tic <- proc.time()
#' useformula <- ~ temp + depth + chla + O2
#' fitcbfm_additive <- CBFM(y = simy_train, formula_X = useformula, data = dat_train, 
#' B_space = train_space_basisfunctions, B_time = train_time_basisfunctions, family = binomial(), 
#' G_control = list(rank = c(5,5)), Sigma_control = list(rank = c(5,5)), control = list(trace = 1))
#' toc <- proc.time()
#' toc - tic
#' 
#'  
#' # Calculate predictions onto test dataset
#' predictions_stacked <- predict(fitstacked, newdata = dat_test, type = "response")
#' predictions_cbfm_additive <- predict(fitcbfm_additive, newdata = dat_test, type = "response", 
#' new_B_space = test_space_basisfunctions, new_B_time = test_time_basisfunctions)
#' 
#' # Evaluation predictions
#' # Tjur R-squared across species
#' tjurR2 <- data.frame(
#' stacked = sapply(1:num_spp, function(j) { 
#' m1 <- predictions_stacked[which(simy_test[,j] > 0),j] %>%
#' mean(na.rm = TRUE)
#' m0 <- predictions_stacked[which(simy_test[,j] == 0),j] %>%
#' mean(na.rm = TRUE)
#' m1 - m0     
#' }),
#' cbfm = sapply(1:num_spp, function(j) { 
#' m1 <- predictions_cbfm_additive[which(simy_test[,j] > 0),j] %>%
#' mean(na.rm = TRUE)
#' m0 <- predictions_cbfm_additive[which(simy_test[,j] == 0),j] %>%
#' mean(na.rm = TRUE)
#' m1 - m0     
#' })
#' )
#' 
#' boxplot(tjurR2, main = "Tjur-R2", names = c("Stacked GLM", "CBFM (additive)"))
#' 
#' ggplot(tjurR2, aes(x = stacked, y = cbfm)) +
#' geom_point() +
#' geom_abline(intercept = 0, slope = 1, linetype = 2) +
#' labs(x = "Stacked SDM", y = "CBFM (additive)", main = "Tjur-R2") +
#' theme_bw()
#' 
#' # AUC across species
#' aucs <- data.frame(
#' stacked = sapply(1:num_spp, function(j) { 
#' pred <- prediction(predictions_stacked[,j], labels = simy_test[,j]) %>%
#' performance(measure = "auc")  
#' pred@y.values[[1]]
#' }),
#' cbfm = sapply(1:num_spp, function(j) { 
#' pred <- prediction(predictions_cbfm_additive[,j], labels = simy_test[,j]) %>%
#' performance(measure = "auc") 
#' pred@y.values[[1]]
#' })
#' )
#' 
#' boxplot(aucs, main = "AUC", names = c("Stacked GLM", "CBFM (additive)"))
#' 
#' ggplot(aucs, aes(x = stacked, y = cbfm)) +
#' geom_point() +
#' geom_abline(intercept = 0, slope = 1, linetype = 2) +
#' labs(x = "Stacked SDM", y = "CBFM (additive)", main = "AUC") +
#' theme_bw()
#' 
#' 
#' ##------------------------------
#' ## **Example 3b: Repeat Example 3a but with tensor product basis functions**
#' ## Please note the data generation process (thus) differs from CBFM.
#' ## To save some time, and for illustrative purposes, we wil use the fast method estimating
#' ## the covariance matrices
#' # Nevertheless, please note this might take quite a while...grab some senbei or a melonpan!
#' ##------------------------------
#' tic <- proc.time()
#' useformula <- ~ temp + depth + chla + O2
#' train_st_basisfunctions <- tensorproduct(train_space_basisfunctions, train_time_basisfunctions)
#' dim(train_st_basisfunctions)
#' 
#' fitcbfm_tensor <- CBFM(y = simy_train, formula_X = useformula, data = dat_train, 
#' B_spacetime = train_st_basisfunctions, family = binomial(), 
#' G_control = list(rank = 10, method = "simple"), 
#' Sigma_control = list(rank = 10, method = "simple"), control = list(trace = 1))
#' toc <- proc.time()
#' toc - tic
#' 
#' 
#' test_st_basisfunctions <- tensorproduct(test_space_basisfunctions, test_time_basisfunctions)
#' predictions_cbfm_tensor <- predict(fitcbfm_tensor, newdata = dat_test, type = "response", 
#' new_B_spacetime = test_st_basisfunctions)
#' 
#' # Tjur-R2 across species
#' tjurR2$cbfm_tensor = sapply(1:num_spp, function(j) { 
#' m1 <- predictions_cbfm_tensor[which(simy_test[,j] > 0),j] %>%
#' mean(na.rm = TRUE)
#' m0 <- predictions_cbfm_tensor[which(simy_test[,j] == 0),j] %>%
#' mean(na.rm = TRUE)
#' m1 - m0     
#' })
#' 
#' boxplot(tjurR2, main = "Tjur-R2", 
#' names = c("Stacked GLM", "CBFM (additive)", "CBFM (tensor)"))
#' 
#' ggplot(tjurR2, aes(x = stacked, y = cbfm_tensor)) +
#' geom_point() +
#' geom_abline(intercept = 0, slope = 1, linetype = 2) +
#' labs(x = "Stacked SDM", y = "CBFM (tensor)", main = "Tjur-R2") +
#' theme_bw()
#' 
#' # AUC across species
#' aucs$cbfm_tensor = sapply(1:num_spp, function(j) { 
#' pred <- prediction(predictions_cbfm_tensor[,j], labels = simy_test[,j]) %>%
#' performance(measure = "auc") 
#' pred@y.values[[1]]
#' })
#' 
#' boxplot(aucs, main = "AUC", 
#' names = c("Stacked GLM", "CBFM (additive)", "CBFM (tensor)"))
#' 
#' ggplot(aucs, aes(x = stacked, y = cbfm_tensor)) +
#' geom_point() +
#' geom_abline(intercept = 0, slope = 1, linetype = 2) +
#' labs(x = "Stacked SDM", y = "CBFM (tensor)", main = "AUC") +
#' theme_bw()
#' 
#' # As one can see, the tensor product CBFM here performs a lot worse than the additive CBFM.
#' # This is partly not surprising given that we used a more approximate fitting method for 
#' # the tensor product CBFM but also because of the additive latent variables underlying the data
#' # generation process.
#' 
#' # Indeed, both the authors of this package have found that in practice, the additive CBFM 
#' # oftens tends to perform better than the tensor product CBFM for 
#' # spatio-temporal multivariate abundance data.
#' }
#' 
#' @export
#'
#' @import foreach  
#' @import Matrix
#' @importFrom compiler cmpfun
#' @importFrom gamlss gamlss gamlss.control getSmo
#' @importFrom gamlss.add ga
#' @importFrom doParallel registerDoParallel
#' @importFrom MASS theta.mm
#' @importFrom methods as
#' @importFrom mgcv betar gam gam.vcomp k.check ldTweedie logLik.gam model.matrix.gam pen.edf nb rTweedie Tweedie tw ziplss
#' @importFrom numDeriv grad
#' @importFrom parallel detectCores
#' @importFrom stats fitted dnorm pnorm qnorm rnorm dbinom pbinom rbinom dnbinom pnbinom rnbinom dbeta pbeta rbeta dexp pexp rexp dgamma pgamma rgamma dlogis plogis qlogis dpois ppois rpois runif dchisq pchisq qchisq qqnorm as.formula binomial formula Gamma logLik model.matrix optim nlminb residuals 
#' @importFrom TMB MakeADFun
#' @importFrom utils capture.output
#' @md


## COMMENCE FUNCTION!
CBFM <- function(y, formula_X, data, B_space = NULL, B_time = NULL, B_spacetime = NULL, 
     offset = NULL, ncores = NULL, family = stats::gaussian(), trial_size = 1, dofit = TRUE, stderrors = TRUE, select = FALSE, gamma = 1,
     start_params = list(betas = NULL, basis_effects_mat = NULL, dispparam = NULL, powerparam = NULL, zeroinfl_prob = NULL),
     TMB_directories = list(cpp = system.file("executables", package = "CBFM"), compile = system.file("executables", package = "CBFM")),
     control = list(maxit = 100, optim_lower = -10, optim_upper = 10, convergence_type = "parameters", tol = 1e-4, initial_beta_dampen = 1, 
                    subsequent_betas_dampen = 0.25, seed = NULL, trace = 0, ridge = 0), 
     Sigma_control = list(rank = 5, maxit = 100, tol = 1e-4, method = "LA", trace = 0), 
     G_control = list(rank = 5, nugget_profile = seq(0.05, 0.95, by = 0.05), maxit = 100, tol = 1e-4, method = "LA", trace = 0),
     k_check_control = list(subsample = 5000, n.rep = 400)
     ) {
          
     ##----------------
     ## Opening checks and all that jazz
     ##----------------
     if(is.null(ncores))
          registerDoParallel(cores = detectCores()-1)
     if(!is.null(ncores))
          registerDoParallel(cores = ncores)

     .check_family(family = family, y = y, trial_size = trial_size) 
     trial_size_length <- as.numeric(length(trial_size) == 1)
          
     if(!is.matrix(y))
          stop("y should be a matrix.")
     if(is.null(colnames(y)))
          colnames(y) <- paste0("response", 1:ncol(y))
     if(is.null(rownames(y)))
          rownames(y) <- paste0("units", 1:nrow(y))
          

     ## Form full basis function matrix B
     .check_B_forms(B_space = B_space, B_time = B_time, B_spacetime = B_spacetime)
     which_B_used <- c(0,0,0)
     num_spacebasisfns <- num_timebasisfns <- num_spacetimebasisfns <- 0
     if(!is.null(B_space)) {
          which_B_used[1] <- 1
          num_spacebasisfns <- ncol(B_space)
          if(is.null(colnames(B_space)))
               colnames(B_space) <- paste0("B_space_",1:num_spacebasisfns)
          }
     if(!is.null(B_time)) {
          which_B_used[2] <- 1
          num_timebasisfns <- ncol(B_time)
          if(is.null(colnames(B_time)))
               colnames(B_time) <- paste0("B_time_",1:num_timebasisfns)
          }
     if(!is.null(B_spacetime)) {
          which_B_used[3] <- 1
          num_spacetimebasisfns <- ncol(B_spacetime)
          if(is.null(colnames(B_spacetime)))
               colnames(B_spacetime) <- paste0("B_spacetime_",1:num_spacetimebasisfns)
          }

     B <- cbind(B_space, B_time, B_spacetime)
     B <- Matrix(B, sparse = TRUE)
     #if(is.null(rownames(B)))
     #     rownames(B) <- paste0("units", 1:nrow(B_space))
          
          
     control <- .fill_control(control = control)
     Sigma_control <- .fill_Sigma_control(control = Sigma_control, which_B_used = which_B_used)
     G_control <- .fill_G_control(control = G_control, which_B_used = which_B_used)

     ## Form covariate model matrix B
     formula_X <- .check_X_formula(formula_X = formula_X, data = as.data.frame(data))          
     tmp_formula <- as.formula(paste("response", paste(as.character(formula_X),collapse="") ) )
     nullfit <- gam(tmp_formula, data = data.frame(data, response = y[,1]), fit = TRUE, control = list(maxit = 1))
     X <- model.matrix(nullfit)
     rm(tmp_formula, nullfit)
     rownames(X) <- rownames(y)

     .check_BX(B = B, X = X) 
     .check_offset(offset = offset, y = y) 
     
     if(Sigma_control$inv_method == "simple" || G_control$inv_method == "simple")
          warning("The simple method is fast but not very reliable for producing good estimates of covariance matrices!")
     
     if(!is.null(control$seed)) {
          set.seed(control$seed)
          on.exit(set.seed(NULL))
          }
               
     num_units <- nrow(y)
     num_spp <- ncol(y)
     num_X <- ncol(X)
     num_basisfns <- ncol(B)
     .check_ranks2(num_spp = num_spp, which_B_used = which_B_used, G_control = G_control, 
                   vec_num_basisfns = c(num_spacebasisfns,num_timebasisfns,num_spacetimebasisfns), Sigma_control = Sigma_control)
     
     ##----------------
     ## Compile TMB C++ files
     ## Adapted from how JT does it with VAST, which is simple and I rather like!
     ##----------------
     if(sum(which_B_used) == 1)
          getDLL <- "cbfm_oneB"
     if(sum(which_B_used) == 2)
          getDLL <- "cbfm_twoB"
     if(sum(which_B_used) == 3)
          getDLL <- "cbfm_threeB"
     if(control$trace > 0) {
          message("Compiling TMB C++ file...")
          }
     
     file.copy(from = paste0(TMB_directories$cpp, "/", getDLL, ".cpp"), to = paste0(TMB_directories$compile, "/", getDLL, ".cpp"), overwrite = FALSE)
     if(!dofit) {
         message("CBFM not fitted. Function is terminated after the C++ file in copied into ", TMB_directories$compile)
         message("Otsukaresama deshita uwu")
         return()
         }
          
     origwd <- getwd()          
     # Please see https://github.com/kaskr/adcomp/issues/321 for flags argument
     setwd(TMB_directories$compile)
     TMB::compile(paste0(getDLL, ".cpp"), flags = "-Wno-ignored-attributes -O2 -mfpmath=sse -msse2 -mstackrealign" )
     dyn.load(paste0(TMB_directories$compile, "/", TMB::dynlib(getDLL)))
     setwd(origwd)
     
     
     ##----------------
     ## Obtain starting values -- Note that no selection is attempted here  
     ##----------------
     .check_start_params(start_params = start_params, num_spp = num_spp, num_basisfns = num_basisfns, num_X = num_X)
     initfit_fn <- function(j, formula_X, fornulldeviance = FALSE) {
          tmp_formula <- as.formula(paste("response", paste(as.character(formula_X),collapse="") ) )
          
          if(family$family %in% c("gaussian","poisson","Gamma")) {
               fit0 <- try(gam(tmp_formula, data = data.frame(response = y[,j], data), offset = offset[,j], method = "ML", family = family, gamma = gamma), silent = TRUE)
               fit0$logLik <- try(logLik(fit0), silent = TRUE)
               }
          if(family$family %in% c("binomial")) {
               tmp_formula <- as.formula(paste("cbind(response, size - response)", paste(as.character(formula_X),collapse="") ) )
               use_size <- .ifelse_size(trial_size = trial_size, trial_size_length = trial_size_length, j = j, num_units = num_units)
               
               fit0 <- try(gam(tmp_formula, data = data.frame(response = y[,j], data, size = use_size), offset = offset[,j], method = "ML", 
                    family = family, gamma = gamma), silent = TRUE)
               fit0$logLik <- try(logLik(fit0), silent = TRUE)
               }
          if(family$family %in% c("negative.binomial")) {
               fit0 <- try(gam(tmp_formula, data = data.frame(response = y[,j], data), offset = offset[,j], method = "ML", family = nb(), gamma = gamma), silent = TRUE)
               fit0$logLik <- try(logLik(fit0), silent = TRUE)
               }
          if(family$family %in% c("Beta")) {
               fit0 <- try(gam(tmp_formula, data = data.frame(response = y[,j], data), offset = offset[,j], method = "ML", family = betar(link = "logit"), gamma = gamma), silent = TRUE)
               fit0$logLik <- try(logLik(fit0), silent = TRUE)
               }
          if(family$family == "tweedie") {
               fit0 <- try(gam(tmp_formula, data = data.frame(response = y[,j], data), offset = offset[,j], method = "ML", family = Tweedie(p = 1.6, link = "log"), gamma = gamma), silent = TRUE)
               fit0$logLik <- try(logLik(fit0), silent = TRUE)
               }
          if(family$family == "zipoisson") {
               # Initial weights/posterior probabilities of being in zero-inflation component
               init_pi <- mean(y[,j] == 0)
               init_lambda <- mean(y[,j])
               w <- ifelse(y[,j] == 0, init_pi / (init_pi + (1-init_pi) * dpois(0, init_lambda)), 0)
               fit0 <- try(gam(tmp_formula, data = data.frame(response = y[,j], data), weights = 1-w, offset = offset[,j], method = "ML", family = "poisson", gamma = gamma), silent = TRUE)
               
               if(fornulldeviance) {
                    inner_err <- Inf
                    cw_inner_logL <- logLik(fit0)
                    while(inner_err > 0.25) {
                         w <- ifelse(y[,j] == 0, init_pi / (init_pi + (1-init_pi) * dpois(0, lambda = fitted(fit0))), 0) # Posterior probabilities of being in zero-inflation component
                         init_pi <- mean(w + 1e-3)
                         fit0 <- try(gam(tmp_formula, data = data.frame(response = y[,j], data), weights = 1-w, offset = offset[,j], method = "ML", family = "poisson", gamma = gamma), silent = TRUE)
                         inner_err <- abs(logLik(fit0)/cw_inner_logL-1)
                         cw_inner_logL <- logLik(fit0)
                         }
                    }
               
               fit0$logLik <- try(sum(log(init_pi*as.numeric(y[,j]==0) + (1-init_pi)*dpois(y[,j], lambda = fitted(fit0)) )), silent = TRUE)
               }
          if(family$family[1] == "zinegative.binomial") {
               # Initial weights/posterior probabilities of being in zero-inflation component
               init_pi <- mean(y[,j] == 0)
               init_lambda <- mean(y[,j])
               w <- ifelse(y[,j] == 0, init_pi / (init_pi + (1-init_pi) * dnbinom(0, mu = init_lambda, size = 1/0.2)), 0)
               fit0 <- try(gam(tmp_formula, data = data.frame(response = y[,j], data), weights = 1-w, offset = offset[,j], method = "ML", family = nb(), gamma = gamma), silent = TRUE)
               
               if(fornulldeviance) {
                    inner_err <- Inf
                    cw_inner_logL <- logLik(fit0)
                    while(inner_err > 0.25) {
                         w <- ifelse(y[,j] == 0, init_pi / (init_pi + (1-init_pi) * dnbinom(0, mu = fitted(fit0), size = fit0$family$getTheta(TRUE))), 0) # Posterior probabilities of being in zero-inflation component
                         init_pi <- mean(w + 1e-3)
                         fit0 <- try(gam(tmp_formula, data = data.frame(response = y[,j], data), weights = 1-w, offset = offset[,j], method = "ML", family = nb(), gamma = gamma), silent = TRUE)
                         inner_err <- abs(logLik(fit0)/cw_inner_logL-1)
                         cw_inner_logL <- logLik(fit0)
                         }
                    }
               
               fit0$logLik <- try(sum(log(init_pi*as.numeric(y[,j]==0) + (1-init_pi)*dnbinom(y[,j], mu = fitted(fit0), size = fit0$family$getTheta(TRUE)) )), silent = TRUE)
               }
          if(family$family %in% c("ztpoisson")) {
               tmp_formula <- as.formula(paste("response", paste(as.character(formula_X),collapse=""), "+ offset(off)" ) )
               tmp_dat <- data.frame(response = c(y[,j], numeric(20)), data[c(1:nrow(data), 1:20),],
                                     off = c(offset[,j], numeric(20))) # Append some zeros speeds ziplss up a heck of a lot!
               fit0 <- try(gam(list(tmp_formula, ~1), data = tmp_dat, method = "ML", family = ziplss(), gamma = gamma), silent = TRUE)
               if(!inherits(fit0, "try-error")) {
                    fit0$coefficients <- fit0$coefficients[1:num_X] 
                    fit0$logLik <- sum(.dztpois(y[,j], lambda = exp(fitted(fit0)[1:num_units,1]), log = TRUE))
                    }
               }
          # if(family$family %in% c("ztnegative.binomial")) {
          #      gac <- ga.control(method = "ML", gamma = gamma)
          #      tmp_formula <- as.formula(paste("response ~ ga(", paste(as.character(formula_X),collapse=""), ", control = gac)"))
          #      # Due to the way ga works, you need to attach things to the global environment. But then delete them immediately afterwards to be save 
          #      tmp_dat <<- data.frame(response = y[,j], data)
          #      tmp_offset <<- offset[,j]
          #      fit0 <- try(gamlss(tmp_formula, data = tmp_dat, offset = tmp_offset, family = ztnb2(), 
          #                         control = gamlss.control(trace = TRUE, n.cyc = control$maxit)), silent = TRUE)
          #      if(inherits(fit0, "try-error")) {
          #           fit0 <- try(gamlss(tmp_formula, data = tmp_dat, offset = tmp_offset, family = NBI(), 
          #                              control = gamlss.control(trace = FALSE, n.cyc = control$maxit)), silent = TRUE)
          #           fit0 <- try(gamlss(tmp_formula, data = tmp_dat, offset = tmp_offset, family = ztnb2(), start.from = fit0, 
          #                              method = CG(), control = gamlss.control(trace = FALSE, n.cyc = control$maxit)), silent = TRUE) # Using CG tends to be more stable but comes at the risk of non-convergence
          #           }
          #      if(fornulldeviance) { # This should always work given it is a null model. We run it separately as interfacing to ga can run into trouble as well 
          #           fit0 <- try(gamlss(response ~ 1, data = tmp_dat, offset = tmp_offset, family = ztnb2(), 
          #                              control = gamlss.control(trace = FALSE, n.cyc = control$maxit)), silent = TRUE)
          #           }
          #      fit1 <- try(getSmo(fit0), silent = TRUE)
          #      fit1$logLik <- try(logLik(fit0), silent = TRUE)
          #      fit0 <- fit1
          #      rm(tmp_dat, tmp_offset, pos = 1)
          #      }
                
                   
          if(inherits(fit0, "try-error"))
               fit0 <- list(coefficients = runif(num_X), dispparam = 1)
          
          return(fit0)
          }
     
     if(is.null(start_params$betas)) { 
          if(control$trace > 0)
               message("Calculating starting values...")
          
          all_start_fits <- foreach(j = 1:num_spp) %dopar% initfit_fn(j = j, formula_X = formula_X)              
          start_params$betas <- do.call(rbind, lapply(all_start_fits, function(x) x$coefficients))
          start_params$betas <- start_params$betas * control$initial_betas_dampen
          rm(all_start_fits)
          gc()
          }
     if(is.null(start_params$basis_effects_mat))
          start_params$basis_effects_mat <- matrix(0, nrow = num_spp, ncol = num_basisfns)
     if(is.null(start_params$dispparam)) {
          if(family$family[1] %in% c("poisson","binomial","zipoisson","ztpoisson"))                        
               start_params$dispparam <- rep(1, num_spp)          
          if(family$family[1] %in% c("gaussian","Gamma","tweedie"))                        
                    start_params$dispparam <- 0.5 * colMeans((y - tcrossprod(X, start_params$betas) - tcrossprod(B, start_params$basis_effects_mat))^2)
          if(family$family[1] %in% c("negative.binomial")) 
               start_params$dispparam <- sapply(1:num_spp, function(j) 0.5/theta.mm(y = y[,j], mu = exp(X%*%start_params$betas[j,] + B%*%start_params$basis_effects_mat[j,]), dfr = num_units - num_X - num_basisfns))
          if(family$family[1] %in% c("zinegative.binomial","ztnegative.binomial"))                        
               start_params$dispparam <- rep(0.2, num_spp)        
          if(family$family[1] %in% c("Beta")) 
               start_params$dispparam <- rep(0.5, num_spp)          
          }
     if(is.null(start_params$powerparam))
          start_params$powerparam <- rep(ifelse(family$family[1] == "tweedie", 1.6, 0), num_spp)
     if(is.null(start_params$zeroinfl_prob)) {
          if(family$family[1] %in% c("zipoisson","zinegative.binomial"))
               start_params$zeroinfl_prob_intercept <- binomial()$linkfun(colMeans(y == 0)+1e-3) 
          if(!(family$family[1] %in% c("zipoisson","zinegative.binomial")))
               start_params$zeroinfl_prob_intercept <- rep(0, num_spp) 
          }
     
     if(which_B_used[1]) {
          start_params$Sigma_space <- diag(1, nrow = num_spacebasisfns)
          start_params$G_space <- diag(1, nrow = num_spp)
          new_LoadingnuggetSigma_space <- list(covinv = start_params$Sigma_space)
          new_LoadingnuggetG_space <- list(covinv = start_params$G_space)
          }
     if(which_B_used[2]) {
          start_params$Sigma_time <- diag(1, nrow = num_timebasisfns)
          start_params$G_time <- diag(1, nrow = num_spp)
          new_LoadingnuggetSigma_time <- list(covinv = start_params$Sigma_time)
          new_LoadingnuggetG_time <- list(covinv = start_params$G_time)
          }
     if(which_B_used[3]) {
          start_params$Sigma_spacetime <- diag(1, nrow = num_spacetimebasisfns)
          start_params$G_spacetime <- diag(1, nrow = num_spp)
          new_LoadingnuggetSigma_spacetime <- list(covinv = start_params$Sigma_spacetime)
          new_LoadingnuggetG_spacetime <- list(covinv = start_params$G_spacetime)
          }
     start_params$logLik <- -Inf
     
     
     ##-------------------------
     ## Run PQL algorithm
     ##-------------------------          
     tic <- proc.time()
     counter <- 0
     diff <- 10
     converged <- FALSE
     if(is.null(offset)) { ## Only make offset here, after initial fitting done
          offset <- Matrix(0, nrow = nrow(y), ncol = ncol(y), sparse = TRUE)
          }
     if(control$trace > 0)
          message("Commencing model fitting...")
     
     while(diff > control$tol & counter < control$maxit) {
          make_tidibits_data <- function() {
               out <- list(B = B)
               if(identical(which_B_used, c(1,0,0))) {
                    out$Sigmainv <- as.matrix(new_LoadingnuggetSigma_space$covinv)
                    out$Ginv <- as.matrix(new_LoadingnuggetG_space$covinv)
                    }
               if(identical(which_B_used, c(0,1,0))) {
                    out$Sigmainv <- as.matrix(new_LoadingnuggetSigma_time$covinv)
                    out$Ginv <- as.matrix(new_LoadingnuggetG_time$covinv)
                    }
               if(identical(which_B_used, c(0,0,1))) {
                    out$Sigmainv <- as.matrix(new_LoadingnuggetSigma_spacetime$covinv)
                    out$Ginv <- as.matrix(new_LoadingnuggetG_spacetime$covinv)
                    }
                    
               if(identical(which_B_used, c(1,1,0))) {
                    out$Sigmainv_B1 <- as.matrix(new_LoadingnuggetSigma_space$covinv)
                    out$Ginv_B1 <- as.matrix(new_LoadingnuggetG_space$covinv)
                    out$Sigmainv_B2 <- as.matrix(new_LoadingnuggetSigma_time$covinv)
                    out$Ginv_B2 <- as.matrix(new_LoadingnuggetG_time$covinv)
                    }
               if(identical(which_B_used, c(1,0,1))) {
                    out$Sigmainv_B1 <- as.matrix(new_LoadingnuggetSigma_space$covinv)
                    out$Ginv_B1 <- as.matrix(new_LoadingnuggetG_space$covinv)
                    out$Sigmainv_B2 <- as.matrix(new_LoadingnuggetSigma_spacetime$covinv)
                    out$Ginv_B2 <- as.matrix(new_LoadingnuggetG_spacetime$covinv)
                    }
               if(identical(which_B_used, c(0,1,1))) {
                    out$Sigmainv_B1 <- as.matrix(new_LoadingnuggetSigma_time$covinv)
                    out$Ginv_B1 <- as.matrix(new_LoadingnuggetG_time$covinv)
                    out$Sigmainv_B2 <- as.matrix(new_LoadingnuggetSigma_spacetime$covinv)
                    out$Ginv_B2 <- as.matrix(new_LoadingnuggetG_spacetime$covinv)
                    }
               if(identical(which_B_used, c(1,1,1))) {
                    out$Sigmainv_B1 <- as.matrix(new_LoadingnuggetSigma_space$covinv)
                    out$Ginv_B1 <- as.matrix(new_LoadingnuggetG_space$covinv)
                    out$Sigmainv_B2 <- as.matrix(new_LoadingnuggetSigma_time$covinv)
                    out$Ginv_B2 <- as.matrix(new_LoadingnuggetG_time$covinv)
                    out$Sigmainv_B3 <- as.matrix(new_LoadingnuggetSigma_spacetime$covinv)
                    out$Ginv_B3 <- as.matrix(new_LoadingnuggetG_spacetime$covinv)
                    }                    

               return(out)
               }

          if(counter == 0) {
               tidbits_data <- make_tidibits_data()
                                   
               new_fit_CBFM_ptest <- list(betas = start_params$betas,
                    basis_effects_mat = start_params$basis_effects_mat,
                    dispparam = start_params$dispparam,
                    powerparam = start_params$powerparam,
                    zeroinfl_prob_intercept = start_params$zeroinfl_prob_intercept,
                    logLik = start_params$logLik
                    )
                    
               cw_params <- as.vector(unlist(new_fit_CBFM_ptest))
               cw_params <- cw_params[-length(cw_params)] # Get rid of the logLik term
               cw_logLik <- new_fit_CBFM_ptest$logLik
               cw_linpred <- Matrix(0, nrow = num_units, ncol = num_spp, sparse = TRUE)
               rm(start_params)
               }
          if(counter > 0) {          
               tidbits_data <- make_tidibits_data()
               }

          
          inner_err <- Inf
          cw_inner_logL <- -Inf
          if(control$trace > 0)
               message("Updating all coefficients and dispersion/power parameters (also running inner EM algorithm if appropriate)...")         
          while(inner_err > 0.01) {
               ##-------------------------
               ## Updating zero-inflation probabilities for distributions that need it. Note it is parameterized in terms of an intercept on the logit scale
               ## Also E-step before that if need be
               ##-------------------------
               getweights <- .estep_fn(family = family, cwfit = new_fit_CBFM_ptest, y = y, X = X, B = B) # Posterior probabilities of zero-inflation
               #if(family$family[1] == "zipoisson") {
               #     for(j in 1:num_spp) {
               #          cwfit <- glm(resp ~ 1, data = data.frame(resp = getweights[,j]), family = binomial())
               #          new_fit_CBFM_ptest$zeroinfl_prob_intercept[j] <- coef(cwfit)[1]
               #          }
               #     }
               #rm(cwfit,j)
               new_fit_CBFM_ptest$zeroinfl_prob_intercept <- binomial()$linkfun(colMeans(getweights + 1e-3))
               
               
               ##-------------------------
               ## Update smoothing coefficients for all basis functions, one species at a time
               ##-------------------------
               update_basiscoefsspp_fn <- function(j) {
                    tidbits_parameters <- list(basis_effects = new_fit_CBFM_ptest$basis_effects_mat[j,])
                    
                    tidbits_data$y <- y[,j]
                    tidbits_data$Xbeta <- as.vector(X %*% new_fit_CBFM_ptest$betas[j,])
                    tidbits_data$dispparam <- new_fit_CBFM_ptest$dispparam[j]
                    tidbits_data$powerparam <- new_fit_CBFM_ptest$powerparam[j]
                    tidbits_data$estep_weights <- as.vector(getweights[,j])
                    tidbits_data$offset <- offset[,j]
                    if(sum(which_B_used) == 1)
                         tidbits_data$other_basis_effects_mat <- new_fit_CBFM_ptest$basis_effects_mat
                    if(sum(which_B_used) == 2) {
                         if(identical(which_B_used, c(1,1,0))) { 
                              tidbits_data$other_basis_effects_mat_B1 <- new_fit_CBFM_ptest$basis_effects_mat[,1:num_spacebasisfns,drop=FALSE]
                              tidbits_data$other_basis_effects_mat_B2 <- new_fit_CBFM_ptest$basis_effects_mat[,num_spacebasisfns+(1:num_timebasisfns),drop=FALSE]
                              }
                         if(identical(which_B_used, c(1,0,1))) { 
                              tidbits_data$other_basis_effects_mat_B1 <- new_fit_CBFM_ptest$basis_effects_mat[,1:num_spacebasisfns,drop=FALSE]
                              tidbits_data$other_basis_effects_mat_B2 <- new_fit_CBFM_ptest$basis_effects_mat[,num_spacebasisfns+num_timebasisfns+(1:num_spacetimebasisfns),drop=FALSE]
                              }
                         if(identical(which_B_used, c(0,1,1))) { 
                              tidbits_data$other_basis_effects_mat_B1 <- new_fit_CBFM_ptest$basis_effects_mat[,num_spacebasisfns+(1:num_timebasisfns),drop=FALSE]
                              tidbits_data$other_basis_effects_mat_B2 <- new_fit_CBFM_ptest$basis_effects_mat[,num_spacebasisfns+num_timebasisfns+(1:num_spacetimebasisfns),drop=FALSE]
                              }
                         }
                    if(sum(which_B_used) == 3) {
                         tidbits_data$other_basis_effects_mat_B1 <- new_fit_CBFM_ptest$basis_effects_mat[,1:num_spacebasisfns,drop=FALSE]
                         tidbits_data$other_basis_effects_mat_B2 <- new_fit_CBFM_ptest$basis_effects_mat[,num_spacebasisfns+(1:num_timebasisfns),drop=FALSE]
                         tidbits_data$other_basis_effects_mat_B3 <- new_fit_CBFM_ptest$basis_effects_mat[,num_spacebasisfns+num_timebasisfns+(1:num_spacetimebasisfns),drop=FALSE]
                         }
                    
                    tidbits_data$spp_ind <- j
                    tidbits_data$family <- .family_to_counter(family = family)
                    tidbits_data$trial_size <- .ifelse_size(trial_size = trial_size, trial_size_length = trial_size_length, j = j, num_units = num_units, family = family)
     
                    tidbits_constraints <- list(lower = rep(control$optim_lower, sum(sapply(tidbits_parameters, length))),
                         upper = rep(control$optim_upper, sum(sapply(tidbits_parameters, length)))
                         )
                         
                    CBFM_objs <- TMB::MakeADFun(data = tidbits_data, parameters = tidbits_parameters, DLL = getDLL, hessian = FALSE, silent = TRUE)
                    
                    new_fit_CBFM <- try(nlminb(start = CBFM_objs$par, objective = CBFM_objs$fn, gradient = CBFM_objs$gr,
                         lower = tidbits_constraints$lower, upper = tidbits_constraints$upper), silent = TRUE)
                    # Dampen Xbeta component and run it again...it is kind of ad-hoc but has been shown to be helpful especially with GAM fits to extremely overdispersed counts 
                    if(inherits(new_fit_CBFM, "try-error")) {
                         tidbits_data$Xbeta <- tidbits_data$Xbeta * control$subsequent_betas_dampen 
                         CBFM_objs <- TMB::MakeADFun(data = tidbits_data, parameters = tidbits_parameters, DLL = getDLL, hessian = FALSE, silent = TRUE)
                         new_fit_CBFM <- try(nlminb(start = CBFM_objs$par, objective = CBFM_objs$fn, gradient = CBFM_objs$gr,
                              lower = tidbits_constraints$lower, upper = tidbits_constraints$upper), silent = TRUE)
                         }
                    if(inherits(new_fit_CBFM, "try-error")) {
                         new_fit_CBFM <- list(par = tidbits_parameters$basis_effects)
                         }
                    #new_fit_CBFM <- optim(par = CBFM_objs$par, fn = CBFM_objs$fn, gr = CBFM_objs$gr, method = "BFGS")
                    
                    return(new_fit_CBFM)
                    }
               update_basiscoefsspp_cmpfn <- compiler::cmpfun(update_basiscoefsspp_fn)
                    
               all_update_coefs <- foreach(j = 1:num_spp) %dopar% update_basiscoefsspp_cmpfn(j = j)
               for(j in 1:num_spp) {
                    new_fit_CBFM_ptest$basis_effects_mat[j,] <- all_update_coefs[[j]]$par #[grep("basis_effects", names(all_update_coefs[[j]]$par))]
                    }
               rm(all_update_coefs, update_basiscoefsspp_fn)

          
               ##-------------------------
               ## Update coefficients related to covariate model matrix X, and other nuisance parameters, one response at a time
               ##-------------------------
               update_Xcoefsspp_fn <- function(j) {
                    tmp_formula <- as.formula(paste("response", paste(as.character(formula_X),collapse="") ) )
                    new_offset <- offset[,j] + as.vector(B %*% new_fit_CBFM_ptest$basis_effects_mat[j,])
                    Hmat <- diag(control$ridge+1e-15, nrow = num_X)
                    
                    if(family$family %in% c("gaussian","poisson","Gamma")) {
                         fit0 <- try(gam(tmp_formula, data = data.frame(response = y[,j], data), offset = new_offset, method = "ML", 
                                         H = Hmat, family = family, select = select, gamma = gamma), silent = TRUE)
                         if(inherits(fit0, "try-error"))
                              fit0 <- gam(tmp_formula, data = data.frame(response = y[,j], data), offset = new_offset, method = "ML", 
                                          family = family, select = select, gamma = gamma)
                         }
                    if(family$family %in% c("binomial")) {
                         tmp_formula <- as.formula(paste("cbind(response, size - response)", paste(as.character(formula_X),collapse="") ) )
                         use_size <- .ifelse_size(trial_size = trial_size, trial_size_length = trial_size_length, j = j, num_units = num_units)
     
                         fit0 <- try(gam(tmp_formula, data = data.frame(response = y[,j], data, size = use_size), offset = new_offset, method = "ML", 
                                         H = Hmat, family = family, select = select, gamma = gamma), silent = TRUE)
                         if(inherits(fit0, "try-error"))
                              fit0 <- gam(tmp_formula, data = data.frame(response = y[,j], data), offset = new_offset, method = "ML", 
                                          family = family, select = select, gamma = gamma)
                         }
                    if(family$family %in% c("negative.binomial")) {
                         fit0 <- try(gam(tmp_formula, data = data.frame(response = y[,j], data), offset = new_offset, method = "ML",
                                         H = Hmat, family = nb(), select = select, gamma = gamma), silent = TRUE)
                         if(inherits(fit0, "try-error"))
                              fit0 <- gam(tmp_formula, data = data.frame(response = y[,j], data), offset = new_offset, method = "ML", 
                                          family = nb(), select = select, gamma = gamma)
                         }
                    if(family$family %in% c("Beta")) {
                         fit0 <- try(gam(tmp_formula, data = data.frame(response = y[,j], data), offset = new_offset, method = "ML", 
                                         H = Hmat, family = betar(link = "logit"), select = select, gamma = gamma), silent = TRUE)
                         if(inherits(fit0, "try-error"))
                              fit0 <- gam(tmp_formula, data = data.frame(response = y[,j], data), offset = new_offset, method = "ML", 
                                          family = betar(link = "logit"), select = select, gamma = gamma)
                         }
                    if(family$family %in% c("tweedie")) {
                         fit0 <- try(gam(tmp_formula, data = data.frame(response = y[,j], data), offset = new_offset, method = "ML",
                                         H = Hmat, family = tw(link = "log"), select = select, gamma = gamma), silent = TRUE)
                         if(inherits(fit0, "try-error"))
                              fit0 <- gam(tmp_formula, data = data.frame(response = y[,j], data), offset = new_offset, method = "ML", 
                                          family = tw(link = "log"), select = select, gamma = gamma)
                         }
                    if(family$family %in% c("zipoisson")) {
                         fit0 <- try(gam(tmp_formula, data = data.frame(response = y[,j], data), offset = new_offset, method = "ML", 
                                         weights = 1-getweights[,j], H = Hmat, family = "poisson", select = select, gamma = gamma), silent = TRUE)
                         if(inherits(fit0, "try-error"))
                              fit0 <- gam(tmp_formula, data = data.frame(response = y[,j], data), offset = new_offset, method = "ML", 
                                          weights = 1-getweights[,j], family = "poisson", select = select, gamma = gamma)
                         }
                    if(family$family %in% c("zinegative.binomial")) {
                         fit0 <- try(gam(tmp_formula, data = data.frame(response = y[,j], data), offset = new_offset, method = "ML", 
                                         weights = 1-getweights[,j], H = Hmat, family = nb(), select = select, gamma = gamma), silent = TRUE)
                         if(inherits(fit0, "try-error"))
                              fit0 <- gam(tmp_formula, data = data.frame(response = y[,j], data), offset = new_offset, method = "ML", 
                                          weights = 1-getweights[,j], family = nb(), select = select, gamma = gamma)
                         }
                    if(family$family %in% c("ztpoisson")) {
                         Hmat <- diag(control$ridge+1e-15, nrow = num_X+1)
                         tmp_dat <- data.frame(response = c(y[,j], numeric(20)), data[c(1:nrow(data), 1:20),], 
                                               off = c(offset[,j] + as.vector(B %*% new_fit_CBFM_ptest$basis_effects_mat[j,]), runif(20))) # Append some zeros speeds ziplss up a heck of a lot!
                         tmp_formula <- as.formula(paste("response", paste(as.character(formula_X),collapse=""), "+ offset(off)" ) )
                         fit0 <- try(gam(list(tmp_formula, ~1), data = tmp_dat, offset = NULL, method = "ML", H = Hmat, family = ziplss(), 
                                            select = select, gamma = gamma), silent = TRUE)
                         if(inherits(fit0, "try-error")) {
                              fit0 <- try(gam(tmp_formula, data = tmp_dat, offset = new_offset, method = "ML", family = ziplss(), select = select, gamma = gamma), silent = TRUE)
                              }
                         if(!inherits(fit0, "try-error")) {
                              fit0$coefficients <- fit0$coefficients[1:num_X] 
                              fit0$linear.predictors <- fit0$linear.predictors[1:num_units,1]
                              fit0$aic <- -2*sum(.dztpois(y[,j], lambda = exp(fitted(fit0)[1:num_units,1]), log = TRUE)) + 2*sum(fit0$edf) # Hack this so that logLik returns the right value
                              }
                         }
                    # if(family$family %in% c("ztnegative.binomial")) {
                    #      gac <- ga.control(method = "ML", gamma = gamma)
                    #      tmp_formula <- as.formula(paste("response ~ ga(", paste(as.character(formula_X),collapse=""), ", control = gac)"))
                    #      tmp_dat <<- data.frame(response = y[,j], data)
                    #      tmp_offset <<- offset[,j]
                    #      fit0 <- try(gamlss(tmp_formula, data = tmp_dat, offset = tmp_offset, family = ztnb2(), 
                    #                         control = gamlss.control(trace = FALSE, n.cyc = control$maxit)), silent = TRUE)
                    #      if(inherits(fit0, "try-error")) {
                    #           fit0 <- try(gamlss(tmp_formula, data = tmp_dat, offset = tmp_offset, family = NBI(), 
                    #                              control = gamlss.control(trace = FALSE, n.cyc = control$maxit)), silent = TRUE)
                    #           fit0 <- try(gamlss(tmp_formula, data = tmp_dat, offset = tmp_offset, family = ztnb2(), start.from = fit0, 
                    #                              method = CG(), control = gamlss.control(trace = FALSE, n.cyc = control$maxit)), silent = TRUE) 
                    #           }
                    #      fit1 <- getSmo(fit0)
                    #      fit0$coefficients <- fit1$coefficients
                    #      fit0$linear.predictors <- fit1$linear.predictors
                    #      fit0$edf <- fit1$edf
                    #      fit0$edf1 <- fit1$edf1
                    #      rm(tmp_dat, tmp_offset, pos = 1)
                    #      rm(fit1)
                    #      }
                    
                    out <- list(coefficients = fit0$coefficients, linear.predictors = fit0$linear.predictors, logLik = as.vector(logLik(fit0)), 
                         fit = fit0, S = .get_bigS(fit_gam = fit0, num_X = num_X))
                    if(family$family %in% c("gaussian","Gamma"))                        
                         out$dispparam <- fit0$sig2
                    if(family$family %in% c("negative.binomial","zinegative.binomial"))                        
                         out$dispparam <- 1/fit0$family$getTheta(TRUE)
                    if(family$family %in% c("ztnegative.binomial"))                        
                         out$dispparam <- as.vector(exp(fit0$sigma.coefficients))
                    if(family$family == "Beta")                        
                         out$dispparam <- exp(fit0$family$getTheta(TRUE))
                    if(family$family == "tweedie") {                        
                         out$dispparam <- fit0$sig2
                         out$powerparam <- fit0$family$getTheta(TRUE)
                         }
               
                    return(out)
                    }               
               update_Xcoefsspp_cmpfn <- compiler::cmpfun(update_Xcoefsspp_fn)
               
               all_update_coefs <- foreach(j = 1:num_spp) %dopar% update_Xcoefsspp_cmpfn(j = j)
               new_fit_CBFM_ptest$betas <- do.call(rbind, lapply(all_update_coefs, function(x) x$coefficients))
               new_fit_CBFM_ptest$linear_predictor <- sapply(all_update_coefs, function(x) x$linear.predictors)         
               for(j in 1:num_spp) {
                    if(family$family[1] %in% c("gaussian","Gamma","negative.binomial","tweedie","Beta","zinegative.binomial","ztnegative.binomial"))                        
                         new_fit_CBFM_ptest$dispparam[j] <- all_update_coefs[[j]]$dispparam
                    if(family$family[1] == "tweedie") {                        
                         new_fit_CBFM_ptest$powerparam[j] <- all_update_coefs[[j]]$powerparam
                         }
                    }
               new_fit_CBFM_ptest$logLik <- sum(sapply(all_update_coefs, function(x) x$logLik))          
               if(family$family[1] == "zipoisson") {
                    cw_logL <- 0
                    for(j in 1:num_spp) {
                         cw_logL <- cw_logL + sum(.dzipoisson_log(y = y[,j], eta = new_fit_CBFM_ptest$linear_predictor[,j], 
                                                                  zeroinfl_prob = plogis(new_fit_CBFM_ptest$zeroinfl_prob_intercept[j]))
                                                  )
                         }
                    new_fit_CBFM_ptest$logLik <- cw_logL
                    rm(cw_logL)
                    }
               if(family$family[1] == "zinegative.binomial") {
                    cw_logL <- 0
                    for(j in 1:num_spp) {
                         cw_logL <- cw_logL + sum(.dzinegativebinomial_log(y = y[,j], eta = new_fit_CBFM_ptest$linear_predictor[,j], 
                                                                  zeroinfl_prob = plogis(new_fit_CBFM_ptest$zeroinfl_prob_intercept[j]),
                                                                  phi = new_fit_CBFM_ptest$dispparam[j])
                                                  )
                         }
                    new_fit_CBFM_ptest$logLik <- cw_logL
                    rm(cw_logL)
                    }
               
               inner_err <- abs(new_fit_CBFM_ptest$logLik/cw_inner_logL-1)
               cw_inner_logL <- new_fit_CBFM_ptest$logLik
               #print(new_fit_CBFM_ptest$logLik)
               rm(all_update_coefs, update_Xcoefsspp_fn)

               if(!(family$family[1] %in% c("zipoisson","zinegative.binomial")))
                    break;
               }          
          
          
          ##-------------------------
          ## Update between spp correlation matrix G. First assume unstructured, then cov2cor, then update loading and nugget. 
          ##-------------------------
          if(control$trace == 1)
               message("Updating between response correlaton matrices, G")
          if(which_B_used[1]) {
               new_G_space <- update_G_fn(Ginv = new_LoadingnuggetG_space$covinv, 
                    basis_effects_mat = new_fit_CBFM_ptest$basis_effects_mat[,1:num_spacebasisfns,drop=FALSE]+G_control$tol, 
                    Sigmainv = new_LoadingnuggetSigma_space$covinv, B = B_space, X = X, y_vec = as.vector(y), 
                    linpred_vec = c(new_fit_CBFM_ptest$linear_predictor),  dispparam = new_fit_CBFM_ptest$dispparam, 
                    powerparam = new_fit_CBFM_ptest$powerparam, zeroinfl_prob_intercept = new_fit_CBFM_ptest$zeroinfl_prob_intercept, 
                    trial_size = trial_size, family = family, G_control = G_control)
               new_LoadingnuggetG_space <- update_LoadingG_fn(G = new_G_space, G_control = G_control, use_rank_element = 1)
               rm(new_G_space)
               }
          if(which_B_used[2]) {
               new_G_time <- update_G_fn(Ginv = new_LoadingnuggetG_time$covinv, 
                    basis_effects_mat = new_fit_CBFM_ptest$basis_effects_mat[,num_spacebasisfns + (1:num_timebasisfns),drop=FALSE]+G_control$tol, 
                    Sigmainv = new_LoadingnuggetSigma_time$covinv, B = B_time, X = X, y_vec = as.vector(y), 
                    linpred_vec = c(new_fit_CBFM_ptest$linear_predictor), dispparam = new_fit_CBFM_ptest$dispparam, 
                    powerparam = new_fit_CBFM_ptest$powerparam, zeroinfl_prob_intercept = new_fit_CBFM_ptest$zeroinfl_prob_intercept, 
                    trial_size = trial_size, family = family, G_control = G_control)
               new_LoadingnuggetG_time <- update_LoadingG_fn(G = new_G_time, G_control = G_control, use_rank_element = sum(which_B_used[1:2]))
               rm(new_G_time)
               }
          if(which_B_used[3]) {
               new_G_spacetime <- update_G_fn(Ginv = new_LoadingnuggetG_spacetime$covinv, 
                    basis_effects_mat = new_fit_CBFM_ptest$basis_effects_mat[,num_spacebasisfns + num_timebasisfns + (1:num_spacetimebasisfns),drop=FALSE]+G_control$tol, 
                    Sigmainv = new_LoadingnuggetSigma_spacetime$covinv, B = B_spacetime, X = X, y_vec = as.vector(y), 
                    linpred_vec = c(new_fit_CBFM_ptest$linear_predictor), dispparam = new_fit_CBFM_ptest$dispparam, 
                    powerparam = new_fit_CBFM_ptest$powerparam, zeroinfl_prob_intercept = new_fit_CBFM_ptest$zeroinfl_prob_intercept, 
                    trial_size = trial_size, family = family, G_control = G_control)
               new_LoadingnuggetG_spacetime <- update_LoadingG_fn(G = new_G_spacetime, G_control = G_control, use_rank_element = sum(which_B_used[1:3]))
               rm(new_G_spacetime)
               }
               
          

          ##-------------------------
          ## Update covariance function for basis functions Sigma. First assuming unstructured, then update loading and nugget 
          ##-------------------------
          if(control$trace == 1)
               message("Updating covariance matrices for basis functions, Sigma")
          if(which_B_used[1]) {
               new_Sigma_space <- update_Sigma_fn(Sigmainv = new_LoadingnuggetSigma_space$covinv, 
                    basis_effects_mat = new_fit_CBFM_ptest$basis_effects_mat[,1:num_spacebasisfns,drop=FALSE]+Sigma_control$tol, 
                    Ginv = new_LoadingnuggetG_space$covinv, B = B_space, X = X, y_vec = as.vector(y), 
                    linpred_vec = c(new_fit_CBFM_ptest$linear_predictor), dispparam = new_fit_CBFM_ptest$dispparam, 
                    powerparam = new_fit_CBFM_ptest$powerparam, zeroinfl_prob_intercept = new_fit_CBFM_ptest$zeroinfl_prob_intercept, 
                    trial_size = trial_size, family = family, Sigma_control = Sigma_control)
               new_LoadingnuggetSigma_space <- update_LoadingSigma_fn(Sigma = new_Sigma_space, Sigma_control = Sigma_control, use_rank_element = 1)
               rm(new_Sigma_space)
               }
          if(which_B_used[2]) {
               new_Sigma_time <- update_Sigma_fn(Sigmainv = new_LoadingnuggetSigma_time$covinv, 
                    basis_effects_mat = new_fit_CBFM_ptest$basis_effects_mat[,num_spacebasisfns + (1:num_timebasisfns),drop=FALSE]+Sigma_control$tol, 
                    Ginv = new_LoadingnuggetG_time$covinv, B = B_time, X = X, y_vec = as.vector(y), 
                    linpred_vec = c(new_fit_CBFM_ptest$linear_predictor),  dispparam = new_fit_CBFM_ptest$dispparam, 
                    powerparam = new_fit_CBFM_ptest$powerparam, zeroinfl_prob_intercept = new_fit_CBFM_ptest$zeroinfl_prob_intercept, 
                    trial_size = trial_size, family = family, Sigma_control = Sigma_control)
               new_LoadingnuggetSigma_time <- update_LoadingSigma_fn(Sigma = new_Sigma_time, Sigma_control = Sigma_control, use_rank_element = sum(which_B_used[1:2]))
               rm(new_Sigma_time)
               }
          if(which_B_used[3]) {
               new_Sigma_spacetime <- update_Sigma_fn(Sigmainv = new_LoadingnuggetSigma_spacetime$covinv, 
                    basis_effects_mat = new_fit_CBFM_ptest$basis_effects_mat[,num_spacebasisfns + num_timebasisfns + (1:num_spacetimebasisfns),drop=FALSE]+Sigma_control$tol, 
                    Ginv = new_LoadingnuggetG_spacetime$covinv, B = B_spacetime, X = X, y_vec = as.vector(y), 
                    linpred_vec = c(new_fit_CBFM_ptest$linear_predictor), dispparam = new_fit_CBFM_ptest$dispparam, 
                    powerparam = new_fit_CBFM_ptest$powerparam, zeroinfl_prob_intercept = new_fit_CBFM_ptest$zeroinfl_prob_intercept, 
                    trial_size = trial_size, family = family, Sigma_control = Sigma_control)
               new_LoadingnuggetSigma_spacetime <- update_LoadingSigma_fn(Sigma = new_Sigma_spacetime, Sigma_control = Sigma_control, 
                    use_rank_element = sum(which_B_used[1:3]))
               rm(new_Sigma_spacetime)
               }
     

          ##-------------------------
          ## Finish iteration
          ##-------------------------
          new_params <- c(c(new_fit_CBFM_ptest$betas), c(new_fit_CBFM_ptest$basis_effects_mat), 
                          log(new_fit_CBFM_ptest$dispparam), new_fit_CBFM_ptest$powerparam, new_fit_CBFM_ptest$zeroinfl_prob_intercept)
          new_logLik <- new_fit_CBFM_ptest$logLik
          if(which_B_used[1])
               new_logLik <- new_logLik + .calc_pqlquadraticterm_basiseffects(
                    basis_effects_mat = new_fit_CBFM_ptest$basis_effects_mat[,1:num_spacebasisfns,drop=FALSE], 
                    Ginv = new_LoadingnuggetG_space$covinv, Sigmainv = new_LoadingnuggetSigma_space$covinv)
          if(which_B_used[2])
               new_logLik <- new_logLik + .calc_pqlquadraticterm_basiseffects(
                    basis_effects_mat = new_fit_CBFM_ptest$basis_effects_mat[,num_spacebasisfns + (1:num_timebasisfns),drop=FALSE], 
                    Ginv = new_LoadingnuggetG_time$covinv, Sigmainv = new_LoadingnuggetSigma_time$covinv)
          if(which_B_used[3])
               new_logLik <- new_logLik + .calc_pqlquadraticterm_basiseffects(
                    basis_effects_mat = new_fit_CBFM_ptest$basis_effects_mat[,num_spacebasisfns + num_timebasisfns + (1:num_spacetimebasisfns),drop=FALSE], 
                    Ginv = new_LoadingnuggetG_spacetime$cov, Sigmainv = new_LoadingnuggetSigma_spacetime$covinv)

          
          if(control$convergence_type == "parameters")
               diff <- mean((new_params - cw_params)^2) 
          if(control$convergence_type == "linear_predictor")
               diff <- mean((new_fit_CBFM_ptest$linear_predictor - cw_linpred)^2) 
          if(control$convergence_type == "logLik")
               diff <- abs(new_logLik/cw_logLik-1) 
          if(control$trace > 0) {
               if(control$convergence_type == "parameters")
                    message("Iteration: ", counter, "\t Difference in parameter estimates: ", round(diff,5))
               if(control$convergence_type == "linear_predictor")
                    message("Iteration: ", counter, "\t MSE of linear predictors: ", round(mean((new_fit_CBFM_ptest$linear_predictor - cw_linpred)^2),5))
               if(control$convergence_type == "logLik")
                    message("Iteration: ", counter, "\t Ratio in PQL log-likelihoods: ", round(diff,5))
               }
               
          cw_params <- new_params
          cw_logLik <- new_logLik
          cw_linpred <- new_fit_CBFM_ptest$linear_predictor
          counter <- counter + 1
          }
     toc <- proc.time()
     gc()

     
     ##-------------------------
     ## Do final fit -- just of coefficients and dispersion/power parameters only. 
     ##-------------------------
     if(diff < control$tol)
          converged <- TRUE
     tidbits_data <- make_tidibits_data()
     inner_err <- Inf
     cw_inner_logL <- -Inf
     while(inner_err > 0.01) {
          getweights <- .estep_fn(family = family, cwfit = new_fit_CBFM_ptest, y = y, X = X, B = B) # Posterior probabilities of zero-inflation
          new_fit_CBFM_ptest$zeroinfl_prob_intercept <- binomial()$linkfun(colMeans(getweights + 1e-5))
     
          all_update_coefs <- foreach(j = 1:num_spp) %dopar% update_basiscoefsspp_cmpfn(j = j)
          for(j in 1:num_spp) {
               new_fit_CBFM_ptest$basis_effects_mat[j,] <- all_update_coefs[[j]]$par#[grep("basis_effects", names(all_update_coefs[[j]]$par))]
               }
          rm(all_update_coefs)
          
          all_update_coefs <- foreach(j = 1:num_spp) %dopar% update_Xcoefsspp_cmpfn(j = j)
          new_fit_CBFM_ptest$betas <- do.call(rbind, lapply(all_update_coefs, function(x) x$coefficients))
          new_fit_CBFM_ptest$linear_predictor <- sapply(all_update_coefs, function(x) x$linear.predictors)          
          new_fit_CBFM_ptest$edf <- sapply(all_update_coefs, function(x) x$fit$edf) # Maybe shoddy for ziplss     
          new_fit_CBFM_ptest$edf1 <- sapply(all_update_coefs, function(x) x$fit$edf1) # Maybe shoddy for ziplss          
          new_fit_CBFM_ptest$pen_edf <- lapply(all_update_coefs, function(x) pen.edf(x$fit)) # Maybe shoddy for ziplss     
          for(j in 1:num_spp) {
               if(family$family %in% c("gaussian","Gamma","negative.binomial","tweedie","Beta", "zinegative.binomial","ztnegative.binomial"))
                    new_fit_CBFM_ptest$dispparam[j] <- all_update_coefs[[j]]$dispparam
               if(family$family == "tweedie") {                        
                    new_fit_CBFM_ptest$powerparam[j] <- all_update_coefs[[j]]$powerparam
                    }
               }
          new_fit_CBFM_ptest$logLik <- sum(sapply(all_update_coefs, function(x) x$logLik))          
          if(family$family[1] == "zipoisson") {
               cw_logL <- 0
               for(j in 1:num_spp) {
                    cw_logL <- cw_logL + sum(.dzipoisson_log(y = y[,j], eta = new_fit_CBFM_ptest$linear_predictor[,j], 
                                                        zeroinfl_prob = plogis(new_fit_CBFM_ptest$zeroinfl_prob_intercept[j]))
                                             )
                    }
               new_fit_CBFM_ptest$logLik <- cw_logL
               rm(cw_logL)
               }
          if(family$family[1] == "zinegative.binomial") {
               cw_logL <- 0
               for(j in 1:num_spp) {
                    cw_logL <- cw_logL + sum(.dzinegativebinomial_log(y = y[,j], eta = new_fit_CBFM_ptest$linear_predictor[,j], 
                                                        zeroinfl_prob = plogis(new_fit_CBFM_ptest$zeroinfl_prob_intercept[j]),
                                                        phi = new_fit_CBFM_ptest$dispparam[j])
                                             )
                    }
               new_fit_CBFM_ptest$logLik <- cw_logL
               rm(cw_logL)
               }
          
          inner_err <- abs(new_fit_CBFM_ptest$logLik/cw_inner_logL-1)
          cw_inner_logL <- new_fit_CBFM_ptest$logLik
          #print(new_fit_CBFM_ptest$logLik)

          if(!(family$family[1] %in% c("zipoisson","zinegative.binomial")))
               break;
          }
     all_S <- sapply(all_update_coefs, function(x) x$S)
     all_k_check <- foreach(j = 1:num_spp) %dopar% k.check(all_update_coefs[[j]]$fit, subsample = k_check_control$subsample, n.rep = k_check_control$n.rep)
     names(all_k_check) <- colnames(y)
     if(family$family[1] %in% c("ztpoisson","ztnegative.binomial"))
          warning("k_check may not be terrible or not available for zero-truncated distributions. Take any results given here with a big grain of salt!")
     invisible(capture.output( all_vcomp <- lapply(1:num_spp, function(j) {
          if(class(all_update_coefs[[j]]$fit)[1] == "gam")
               out <- gam.vcomp(all_update_coefs[[j]]$fit)
          if(class(all_update_coefs[[j]]$fit)[1] == "gamlss")
               out <- gam.vcomp(getSmo(all_update_coefs[[j]]$fit))
          if(is.matrix(out))
               return(out[,1])
          if(is.list(out))
               return(out$all)
          if(is.numeric(out))
               return(out)
          }
          )))
     names(all_vcomp) <- colnames(y)
     rm(all_update_coefs, tidbits_data, inner_err, cw_inner_logL, cw_logLik, cw_params, new_params, diff, counter)
     gc()
     
     
     # Calculate deviance, null deviance etc...note deviance calculation excludes the quadratic term in the PQL
     nulldeviance <- foreach(j = 1:num_spp) %dopar% initfit_fn(j = j, formula_X = ~ 1, fornulldeviance = TRUE)
     nulldeviance <- sum(sapply(nulldeviance, function(x) -2*x$logLik))
     rm(initfit_fn)
     
     new_logLik <- new_fit_CBFM_ptest$logLik
     if(which_B_used[1])
          new_logLik_pql <- new_logLik + .calc_pqlquadraticterm_basiseffects(
               basis_effects_mat = new_fit_CBFM_ptest$basis_effects_mat[,1:num_spacebasisfns,drop=FALSE], 
               Ginv = new_LoadingnuggetG_space$covinv, Sigmainv = new_LoadingnuggetSigma_space$covinv)
     if(which_B_used[2])
          new_logLik_pql <- new_logLik + .calc_pqlquadraticterm_basiseffects(
               basis_effects_mat = new_fit_CBFM_ptest$basis_effects_mat[,num_spacebasisfns + (1:num_timebasisfns),drop=FALSE], 
               Ginv = new_LoadingnuggetG_time$covinv, Sigmainv = new_LoadingnuggetSigma_time$covinv)
     if(which_B_used[3])
          new_logLik_pql <- new_logLik + .calc_pqlquadraticterm_basiseffects(
               basis_effects_mat = new_fit_CBFM_ptest$basis_effects_mat[,num_spacebasisfns + num_timebasisfns + (1:num_spacetimebasisfns),drop=FALSE],
               Ginv = new_LoadingnuggetG_spacetime$cov, Sigmainv = new_LoadingnuggetSigma_spacetime$covinv)
     
     
     ##-----------------
     ## Format output
     ##-----------------
     out_CBFM <- list(call = match.call())
     out_CBFM$family <- family
     out_CBFM$y <- y
     out_CBFM$data <- data
     out_CBFM$trial_size <- trial_size
     out_CBFM$formula_X <- formula_X
     out_CBFM$select <- select
     out_CBFM$gamma <- gamma
     out_CBFM$B <- B
     out_CBFM$which_B_used <- which_B_used
     out_CBFM$num_B_space <- num_spacebasisfns
     out_CBFM$num_B_time <- num_timebasisfns
     out_CBFM$num_B_spacetime <- num_spacetimebasisfns
     out_CBFM$num_B <- num_basisfns
     out_CBFM$converged <- converged
     out_CBFM$logLik <- new_logLik
     out_CBFM$pql_logLik <- new_logLik_pql
     out_CBFM$deviance <- -2*out_CBFM$logLik
     out_CBFM$null_deviance <- nulldeviance
     out_CBFM$deviance_explained <- 100*(out_CBFM$null_deviance - out_CBFM$deviance)/out_CBFM$null_deviance
     out_CBFM$edf <- new_fit_CBFM_ptest$edf
     out_CBFM$edf1 <- new_fit_CBFM_ptest$edf1
     out_CBFM$pen_edf <- new_fit_CBFM_ptest$pen_edf
     out_CBFM$k_check <- all_k_check
     out_CBFM$vcomp <- all_vcomp
     
     out_CBFM$betas <- new_fit_CBFM_ptest$betas
     out_CBFM$basis_effects_mat <- new_fit_CBFM_ptest$basis_effects_mat
     out_CBFM$dispparam <- new_fit_CBFM_ptest$dispparam
     out_CBFM$powerparam <- new_fit_CBFM_ptest$powerparam
     out_CBFM$zeroinfl_prob_intercept <- new_fit_CBFM_ptest$zeroinfl_prob_intercept
     out_CBFM$linear_predictor <- new_fit_CBFM_ptest$linear_predictor
     rm(new_fit_CBFM_ptest, getweights, converged, all_k_check)

     # ## Get restricted CBFM estimates of the coefficients (really more for exploration at this point in time)
     #OLSmatrix_transpose <- X %*% solve(crossprod(X))
     #out_CBFM$partitioned_beta <- out_CBFM$betas + tcrossprod(out_CBFM$basis_effects_mat, B) %*% OLSmatrix_transpose
     #rm(OLSmatrix_transpose)

     if(!(family$family %in% c("zipoisson","zinegative.binomial","ztpoisson","ztnegative.binomial"))) 
          out_CBFM$fitted <- family$linkinv(out_CBFM$linear_predictor)
     if(family$family %in% c("zipoisson","zinegative.binomial"))
          out_CBFM$fitted <- family$linkinv(out_CBFM$linear_predictor) * matrix(1-plogis(out_CBFM$zeroinfl_prob_intercept), nrow = num_units, ncol = num_spp, byrow = TRUE)
     if(family$family == "ztpoisson")
          out_CBFM$fitted <- exp(out_CBFM$linear_predictor) / (1 - dpois(0, lambda = exp(out_CBFM$linear_predictor)))
     if(family$family == "ztnegative.binomial")
          out_CBFM$fitted <- exp(out_CBFM$linear_predictor) / (1 - dnbinom(0, mu = exp(out_CBFM$linear_predictor), size = matrix(1/out_CBFM$dispparam, nrow = num_units, ncol = num_spp, byrow = TRUE)))
     
     names(out_CBFM$dispparam) <- names(out_CBFM$powerparam) <- names(out_CBFM$zeroinfl_prob_intercept) <- 
          colnames(out_CBFM$edf) <- colnames(out_CBFM$edf1) <- colnames(y)
     if(!is.null(all_S[[1]])) {
          names(out_CBFM$pen_edf) <- colnames(y)
          }
          
     
     if(which_B_used[1]) {
          out_CBFM$Sigma_space <- new_LoadingnuggetSigma_space$cov
          out_CBFM$Loading_Sigma_space <- new_LoadingnuggetSigma_space$Loading
          out_CBFM$nugget_Sigma_space <- new_LoadingnuggetSigma_space$nugget
          out_CBFM$G_space <- new_LoadingnuggetG_space$cov
          out_CBFM$Loading_G_space <- new_LoadingnuggetG_space$Loading
          out_CBFM$nugget_G_space <- new_LoadingnuggetG_space$nugget
          rm(new_LoadingnuggetG_space)
          
          rownames(out_CBFM$G_space) <- colnames(out_CBFM$G_space) <- colnames(y)
          colnames(out_CBFM$Sigma_space) <- rownames(out_CBFM$Loading_Sigma_space) <- colnames(B_space)
          if(num_spp > 2)
               colnames(out_CBFM$Loading_G_space) <- paste0("Loading", 1:G_control$rank[1])
          if(num_spacebasisfns > 2) {
               colnames(out_CBFM$Loading_Sigma_space) <- paste0("Loading", 1:Sigma_control$rank[1])
               }
          }          
     if(which_B_used[2]) {
          out_CBFM$Sigma_time <- new_LoadingnuggetSigma_time$cov
          out_CBFM$Loading_Sigma_time <- new_LoadingnuggetSigma_time$Loading
          out_CBFM$nugget_Sigma_time <- new_LoadingnuggetSigma_time$nugget

          out_CBFM$G_time <- new_LoadingnuggetG_time$cov
          out_CBFM$Loading_G_time <- new_LoadingnuggetG_time$Loading
          out_CBFM$nugget_G_time <- new_LoadingnuggetG_time$nugget
          rm(new_LoadingnuggetSigma_time)
          
          rownames(out_CBFM$G_time) <- colnames(out_CBFM$G_time) <- colnames(y)
          colnames(out_CBFM$Sigma_time) <- rownames(out_CBFM$Loading_Sigma_time) <- colnames(B_time)
          if(num_spp > 2)
               colnames(out_CBFM$Loading_G_time) <- paste0("Loading", 1:G_control$rank[sum(which_B_used[1:2])])
          if(num_timebasisfns > 2)
               colnames(out_CBFM$Loading_Sigma_time) <- paste0("Loading", 1:Sigma_control$rank[sum(which_B_used[1:2])])
          }
     if(which_B_used[3]) {
          out_CBFM$Sigma_spacetime <- new_LoadingnuggetSigma_spacetime$cov
          out_CBFM$Loading_Sigma_spacetime <- new_LoadingnuggetSigma_spacetime$Loading
          out_CBFM$nugget_Sigma_spacetime <- new_LoadingnuggetSigma_spacetime$nugget
          out_CBFM$G_spacetime <- new_LoadingnuggetG_spacetime$cov
          out_CBFM$Loading_G_spacetime <- new_LoadingnuggetG_spacetime$Loading
          out_CBFM$nugget_G_spacetime <- new_LoadingnuggetG_spacetime$nugget
          rm(new_LoadingnuggetSigma_spacetime)

          rownames(out_CBFM$G_spacetime) <- colnames(out_CBFM$G_spacetime) <- colnames(y)
          colnames(out_CBFM$Sigma_spacetime) <- rownames(out_CBFM$Loading_Sigma_spacetime) <- colnames(B_spacetime)          
          if(num_spp > 2)
               colnames(out_CBFM$Loading_G_spacetime) <- paste0("Loading", 1:G_control$rank[sum(which_B_used[1:3])])
          if(num_spacetimebasisfns > 2)
               colnames(out_CBFM$Loading_Sigma_spacetime) <- paste0("Loading", 1:Sigma_control$rank[sum(which_B_used[1:3])])
          }

     rownames(out_CBFM$betas) <- rownames(out_CBFM$basis_effects_mat) <- colnames(out_CBFM$linear_predictor) <- colnames(out_CBFM$fitted) <- colnames(y)
     rownames(out_CBFM$linear_predictor) <- rownames(out_CBFM$fitted) <- rownames(X)
     colnames(out_CBFM$betas) <- colnames(X)
     colnames(out_CBFM$basis_effects_mat) <- colnames(B)     
     

     ##-------------------------
     ## Calculate structures needed for producing standard errors for coefficients
     ## The Bayesian posterior covariance matrix used, as opposed to the frequentist sandwich form. This is consistent with the default available in mgcv     
     ## Similar to the default with summary.gam in mgcv, the uncertainty in the nuisance parameters or the covariance matrix is not accounted for! 
     ## Make use of blockwise inversion 
     ##-------------------------
     out_CBFM$stderrors <- stderrors
     if(stderrors) {          
          if(control$trace)
               message("Calculating (components of) the covariance (standard error) matrix...")
          
          weights_mat <- .neghessfamily(family = family, eta = out_CBFM$linear_predictor, y = y, 
                                        phi = matrix(out_CBFM$dispparam, num_units, num_spp, byrow = TRUE), 
                                        powerparam = matrix(out_CBFM$powerparam, num_units, num_spp, byrow = TRUE),
                                        zeroinfl_prob_intercept = matrix(out_CBFM$zeroinfl_prob_intercept, num_units, num_spp, byrow = TRUE), 
                                        trial_size = trial_size, domore = TRUE)
          if(!(family$family[1] %in% c("zipoisson","zinegative.binomial")))
               weights_mat <- matrix(weights_mat$out, nrow = num_units, ncol = num_spp) # Overwrite weights_mat since only one quantity needed
          if(family$family[1] %in% c("zipoisson","zinegative.binomial"))
               weights_mat_betabeta <- matrix(weights_mat$out, nrow = num_units, ncol = num_spp)
          
          # Bottom right of covariance matrix
          D1minusCAinvB_fn <- function(j) {                
               if(!(family$family[1] %in% c("zipoisson","zinegative.binomial"))) {
                    XTWX_inv <- crossprod(X*sqrt(weights_mat[,j])) + Diagonal(x = control$ridge+1e-15, n = num_X) + all_S[[j]]
                    XTWX_inv <- chol2inv(chol( 0.5*(XTWX_inv + t(XTWX_inv)) ))
                    BTWX <- crossprod(B, X*weights_mat[,j])               
                    return(crossprod(B*sqrt(weights_mat[,j])) - BTWX %*% tcrossprod(XTWX_inv, BTWX))
                    }

               if(family$family[1] %in% c("zipoisson","zinegative.binomial")) {
                    Xi <- bdiag(matrix(1, num_units, 1), X)
                    bigW <- cbind(Diagonal(x = weights_mat$out_zeroinflzeroinfl[,j]),  Diagonal(x = weights_mat$out_zeroinflbetas[,j]))
                    bigW <- rbind(bigW, cbind(Diagonal(x = weights_mat$out_zeroinflbetas[,j]),  Diagonal(x = weights_mat_betabeta[,j])))
                    XTWX_inv <- crossprod(Xi, bigW) %*% Xi + Diagonal(x = control$ridge+1e-15, n = num_X+1) + bdiag(matrix(0,1,1), all_S[[j]])
                    XTWX_inv <- chol2inv(chol( 0.5*(XTWX_inv + t(XTWX_inv)) ))
                    BTWX <- crossprod(B, cbind(Diagonal(x = weights_mat$out_zeroinflbetas[,j]), Diagonal(x = weights_mat_betabeta[,j]))) %*% Xi               
                    return(crossprod(B*sqrt(weights_mat_betabeta[,j])) - BTWX %*% tcrossprod(XTWX_inv, BTWX))
                    }
               }
          all_D1minusCAinvB <- foreach(j = 1:num_spp) %dopar% D1minusCAinvB_fn(j = j)
          all_D1minusCAinvB <- bdiag(all_D1minusCAinvB)
          if(identical(which_B_used, c(1,0,0)))
               DminusCAinvB_inv <- forceSymmetric(all_D1minusCAinvB + kronecker(chol2inv(chol(out_CBFM$G_space)), chol2inv(chol(out_CBFM$Sigma_space))))
          if(identical(which_B_used, c(0,1,0)))
               DminusCAinvB_inv <- forceSymmetric(all_D1minusCAinvB + kronecker(chol2inv(chol(out_CBFM$G_time)), chol2inv(chol(out_CBFM$Sigma_time))))
          if(identical(which_B_used, c(0,0,1)))
               DminusCAinvB_inv <- forceSymmetric(all_D1minusCAinvB + kronecker(chol2inv(chol(out_CBFM$G_spacetime)), chol2inv(chol(out_CBFM$Sigma_spacetime))))
          if(identical(which_B_used, c(1,1,0)))
               DminusCAinvB_inv <- forceSymmetric(all_D1minusCAinvB + 
                    .kkproduct(G1 = out_CBFM$G_space, G2 = out_CBFM$G_time, Sigma1 = out_CBFM$Sigma_space, Sigma2 = out_CBFM$Sigma_time)) 
          if(identical(which_B_used, c(1,0,1)))
               DminusCAinvB_inv <- forceSymmetric(all_D1minusCAinvB + 
                    .kkproduct(G1 = out_CBFM$G_space, G2 = out_CBFM$G_spacetime, Sigma1 = out_CBFM$Sigma_space, Sigma2 = out_CBFM$Sigma_spacetime)) 
          if(identical(which_B_used, c(0,1,1)))
               DminusCAinvB_inv <- forceSymmetric(all_D1minusCAinvB + 
                    .kkproduct(G1 = out_CBFM$G_time, G2 = out_CBFM$G_spacetime, Sigma1 = out_CBFM$Sigma_time, Sigma2 = out_CBFM$Sigma_spacetime)) 
          if(identical(which_B_used, c(1,1,1)))
               DminusCAinvB_inv <- forceSymmetric(all_D1minusCAinvB + 
                    .kkproduct(G1 = out_CBFM$G_space, G2 = out_CBFM$G_time, G3 = out_CBFM$G_spacetime, 
                         Sigma1 = out_CBFM$Sigma_space, Sigma2 = out_CBFM$Sigma_time, Sigma3 = out_CBFM$Sigma_spacetime)) 
          
          rm(all_D1minusCAinvB, D1minusCAinvB_fn)
          DminusCAinvB_inv <- chol2inv(chol(DminusCAinvB_inv)) ## Bottleneck! 

          # Top right of covariance matrix -- Could probably remove this if you use the same calculations above
          AinvandB_fn <- function(j) {
               if(!(family$family[1] %in% c("zipoisson","zinegative.binomial"))) {
                    XTWX_inv <- crossprod(X*sqrt(weights_mat[,j])) + Diagonal(x=control$ridge+1e-15, n = num_X) + all_S[[j]]
                    XTWX_inv <- chol2inv(chol( 0.5*(XTWX_inv + t(XTWX_inv)) ))
                    XTWB <- crossprod(X*weights_mat[,j], B)
                    }
               if(family$family[1] %in% c("zipoisson","zinegative.binomial")) {
                    Xi <- bdiag(matrix(1, num_units, 1), X)
                    bigW <- cbind(Diagonal(x = weights_mat$out_zeroinflzeroinfl[,j]),  Diagonal(x = weights_mat$out_zeroinflbetas[,j]))
                    bigW <- rbind(bigW, cbind(Diagonal(x = weights_mat$out_zeroinflbetas[,j]),  Diagonal(x = weights_mat_betabeta[,j])))
                    XTWX_inv <- crossprod(Xi, bigW) %*% Xi + Diagonal(x = control$ridge+1e-15, n = num_X+1) + bdiag(matrix(0,1,1), all_S[[j]])
                    XTWX_inv <- chol2inv(chol( 0.5*(XTWX_inv + t(XTWX_inv)) ))
                    XTWB <- crossprod(Xi, rbind(Diagonal(x = weights_mat$out_zeroinflbetas[,j]), Diagonal(x = weights_mat_betabeta[,j]))) %*% B
                    }
               return(list(Ainv = XTWX_inv, B = XTWB))
               }
          all_AinvandB <- foreach(j = 1:num_spp) %dopar% AinvandB_fn(j = j)
          
          AinvBDminusCAinvB_inv <- bdiag(lapply(all_AinvandB, function(x) x$Ainv)) %*% bdiag(lapply(all_AinvandB, function(x) x$B)) %*% DminusCAinvB_inv 
               
          # Collect parts of covariance matrix
          out_CBFM$covar_components <- list(
               topleft = bdiag(lapply(all_AinvandB, function(x) x$Ainv)) + AinvBDminusCAinvB_inv %*% bdiag(lapply(all_AinvandB, function(x) t(x$B))) %*% bdiag(lapply(all_AinvandB, function(x) x$Ainv)),
               topright = -AinvBDminusCAinvB_inv,
               bottomright = DminusCAinvB_inv
               )
          rm(all_AinvandB, AinvBDminusCAinvB_inv, AinvandB_fn, DminusCAinvB_inv)
                    
                    
          # ## All matrices are dense. But to save memory and because only certain components of the matrices are needed later on, extract relevant principle submatrices only. This is now abandoned
          #out_CBFM$covar_components$topright <- foreach(j = 1:num_spp) %dopar% .extractcovarblocks_topright(j = j, Q = out_CBFM$covar_components$topright)               
          #out_CBFM$covar_components$bottomright <- foreach(j = 1:num_spp) %dopar% .extractcovarblocks_bottomright(j = j, Q = out_CBFM$covar_components$bottomright)
               
          if(!(family$family[1] %in% c("zipoisson","zinegative.binomial")))
               rownames(out_CBFM$covar_components$topleft) <- colnames(out_CBFM$covar_components$topleft) <- rownames(out_CBFM$covar_components$topright) <-
                    apply(as.data.frame.table(t(out_CBFM$betas))[,1:2],1,function(x) paste(x, collapse = ":"))
          if(family$family[1] %in% c("zipoisson","zinegative.binomial")) {
               make_tab <- cbind(out_CBFM$zeroinfl_prob_intercept, out_CBFM$betas)
               colnames(make_tab)[1] <- "ZeroInfl(Intercept)"
               rownames(out_CBFM$covar_components$topleft) <- colnames(out_CBFM$covar_components$topleft) <- rownames(out_CBFM$covar_components$topright) <-
                    apply(as.data.frame.table(t(make_tab))[,1:2],1,function(x) paste(x, collapse = ":"))
               rm(make_tab)
               }
          
          colnames(out_CBFM$covar_components$topright) <- rownames(out_CBFM$covar_components$bottomright) <- colnames(out_CBFM$covar_components$bottomright) <-
               apply(as.data.frame.table(t(out_CBFM$basis_effects_mat))[,1:2],1,function(x) paste(x, collapse = ":"))
          } 
 
 
     if(!(family$family[1] %in% c("Beta","gaussian","Gamma","negative.binomial","tweedie","zinegative.binomial","ztnegative.binomial")))
          out_CBFM$dispparam <- NULL
     if(!(family$family %in% c("tweedie")))                        
          out_CBFM$powerparam <- NULL
     if(!(family$family %in% c("zipoisson","zinegative.binomial")))                        
          out_CBFM$zeroinfl_prob_intercept <- NULL
     
     out_CBFM$time_taken <- toc[3] - tic[3] 
     class(out_CBFM) <- "CBFM"
     return(out_CBFM)
     }
     
