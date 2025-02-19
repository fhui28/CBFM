#' @title Community-level basis function models (CBFMs)
#' 
#' @description 
#' `r lifecycle::badge("experimental")`
#' 
#' Fits CBFMs to spatio-temporal multivariate abundance data, where the basis functions are used to account for spatio-temporal correlation within and between-species. Three types of basis functions can supplied and included in conjunction with each other: 1) spatial basis functions; 2) temporal basis functions; 3) spatio-temporal basis functions. For the part of the mean model corresponding to the measured covariates, CBFM currently permits both parametric terms and/or smoothing terms, where the latter makes use are included in a similar manner to [mgcv::gam()]. Estimation and inference for CBFM is based on a maximum penalized quasi-likelihood (PQL) estimation approach.
#'
#' @param y A response matrix, where each row corresponds to an observational unit \eqn{i}.e, a particular space-time coordinate, and each column corresponds to a species.
#' @param formula An object of class "formula", which represents a symbolic description of the model matrix to be created (based on using this argument along with the \code{data} argument). Note there should be nothing on the left hand side of the "~". Formulas based on generalized additive models or GAMs are permitted (at least, for the basic smoothing terms we have tried so far!); please see [mgcv::formula.gam()], [mgcv::gam.models()], [mgcv::smooth.terms()], and [mgcv::s()] for more details. 
#' @param ziformula An object of class "formula", which represents a symbolic description of the model matrix to be created for the zero-inflation component (based on using this argument along with the \code{data} argument), if appropriate. Note there should be nothing on the left hand side of the "~". Formulas based on generalized additive models or GAMs are permitted (at least, for the basic smoothing terms we have tried so far!)
#' @param data A data frame containing covariate information, from which the model matrix is to be created (based on this argument along with the \code{formula} argument). 
#' @param B_space An optional matrix of spatial basis functions to be included in the CBFM. One of \code{B_space}, \code{B_time}, or \code{B_spacetime} must be supplied. The basis function matrix may be sparse or dense in form; please see the details and examples later on for illustrations of how they can constructed.
#' @param B_time An optional of matrix of temporal basis functions to be included in the CBFM. One of \code{B_space}, \code{B_time}, or \code{B_spacetime} must be supplied. The basis function matrix may be sparse or dense in form; please see the details and examples later on for illustrations of how they can constructed.
#' @param B_spacetime An optional of matrix of spatio-temporal basis functions to be included in the CBFM e.g., formed from a tensor-product of spatial and temporal basis functions. One of \code{B_space}, \code{B_time}, or \code{B_spacetime} must be supplied. The basis function matrix may be sparse or dense in form; please see the details and examples later on for illustrations of how they can constructed.
#' @param offset A matrix of offset terms to applied in association with the \code{formula} argument.
#' @param ncores To speed up fitting, parallelization can be performed, in which case this argument can be used to supply the number of cores to use in the parallelization. Defaults to \code{detectCores()-1}.
#' @param family a description of the response distribution to be used in the model, as specified by a family function. Please see details below for more information on the distributions currently permitted.
#' @param trial_size Trial sizes to use for binomial distribution. This can either equal a scalar or a matrix with the same dimension as \code{y}.
#' @param dofit Should the CBFM be fitted? If set to \code{FALSE}, then the function terminates (and return nothing) immediately after copying the C++ file to the compilation directory; please see the \code{TMB_directories} argument below.
#' @param stderrors Should standard errors of the estimates be calculated? This defaults to \code{TRUE}, but can be set of \code{FALSE} if only point estimations of the regression coefficients for the covariates and basis functions are desired. Please see details later on for more information on how standard errors are constructed. 
#' @param select For cases where \code{formula} involves smoothing terms, setting this to \code{TRUE} adds an extra penalty to each smoothing term so that it can be penalized to zero i.e., null space penalization. Please see [mgcv::gam.selection()] and [mgcv::step.gam()] for more details. Note also this argument has no effect on any parametric terms in the model i.e., it can not shrink parametric terms to zero.  
#' @param gamma For cases where \code{formula} involves smoothing terms, setting this to a value greater than one leads to smoother terms i.e., increased penalization. Note the argument can either be set to a scalar, or a vector with length equal to the number of species i.e., \code{ncol(y)}. This argument plays exactly the same role as the \code{gamma} argument in [mgcv::gam()]. Finally, note this argument has no effect on any parametric terms or the basis functions part of the CBFM.  
#' @param knots For cases where \code{formula} involves smoothing terms, this is an optional list containing user specified knot values to be used for basis construction. For most bases, the user simply supplies the knots to be used, which must match up with the k value supplied. Please see [mgcv::gam()] for more details. Note the knot values are assumed to be common across species. 
#' @param ziselect For cases where \code{ziformula} involves smoothing terms, setting this to \code{TRUE} adds an extra penalty to each smoothing term so that it can be penalized to zero i.e., null space penalization. Please see [mgcv::gam.selection()] and [mgcv::step.gam()] for more details. Note also this argument has no effect on any parametric terms in the model i.e., it can not shrink parametric terms to zero.  
#' @param zigamma For cases where \code{ziformula} involves smoothing terms, setting this to a value greater than one leads to smoother terms i.e., increased penalization. Note the argument can either be set to a scalar, or a vector with length equal to the number of species i.e., \code{ncol(y)}. This argument plays exactly the same role as the \code{gamma} argument in [mgcv::gam()]. Finally, note this argument has no effect on any parametric terms or the basis functions part of the CBFM.  
#' @param ziknots For cases where \code{ziformula} involves smoothing terms, this is an optional list containing user specified knot values to be used for basis construction. For most bases, the user simply supplies the knots to be used, which must match up with the k value supplied. Please see [mgcv::gam()] for more details. Note the knot values are assumed to be common across species. 
#' @param nonzeromean_B_space This experimental feature allows the distribution of the spatial basis function coefficients to have a non-zero mean vector. *Use at your own risk!*
#' @param nonzeromean_B_time This experimental feature allows the distribution of the temporal basis function coefficients to have a non-zero mean vector. *Use at your own risk!*
#' @param nonzeromean_B_spacetime This experimental feature allows the distribution of the spatio-temporal basis function coefficients to have a non-zero mean vector. *Use at your own risk!*
#' @param start_params Starting values for the CBFM. If desired, then a list should be supplied, which must contain at least one the following terms: 
#' \describe{
#' \item{betas: }{A matrix of starting values for the species-specific regression coefficients related to the covariates, where the number of rows is equal to the number of species.} 

#' \item{zibetas: }{A matrix of starting values for the species-specific regression coefficients related to the covariates for the zero-inflation component (if included), where the number of rows is equal to the number of species.} 

#' \item{basis_effect_mat: }{A matrix of starting values for the species-specific regression coefficients related to the combined matrix of basis functions. Again, the number of rows is equal to the number of species, while the number of columns should equal to \code{ncol(B_space, B_time, B_spacetime)} (or whatever the supplied basis functions are).}

#' \item{dispparam: }{A vector of starting values for the species-specific dispersion parameters, to be used for distributions that require one.}

#' \item{powerparam: }{A vector of starting values for the species-specific power parameters, to be used for distributions that require one.}

#' \item{custom_space_lambdas: }{Not used and can be safely ignored.}

#' \item{custom_time_lambdas: }{Not used and can be safely ignored.}

#' \item{custom_spacetime_lambdas: }{Not used and can be safely ignored.}
#' }
#' @param TMB_directories A list with two elements, identifying the directory where TMB C++ file exists (\code{cpp}), and the directory where the corresponding compiled files to be placed (\code{compile}). Unless you really want to do some real mucking around, these should be left at their default i.e., the directory where the packages were installed locally. Please note a version of the C++ file will be copied to the \code{compile} directory.
#' @param control A list of parameters for controlling the fitting process for the "outer" PQL estimation part of the CBFM. This should be a list with the following arguments:
#' \describe{
#' \item{maxit: }{The maximum number of iterations for the outer algorithm. Defaults to 100.} 

#' \item{inner_maxit: }{The maximum number of iterations for the inner (EM) algorithm. Defaults to 1, although it is recommended that this is tested out and increased if convergence presents as an issue.}  

#' \item{optim_lower/optim_upper: }{Upper and lower box constraints when updating regression coefficients related to the basis functions. Note no constraints are put in place when updating regression coefficients related to the covariates; this are controlled internally by [mgcv::gam.control()] itself.}

#' \item{convergence_type: }{The type of means by which to assess convergence. The current options are "parameters_MSE" (default), which assesses convergence based on the mean squared error of the difference between estimated parameters from successive iterations, "parameters_norm" which assesses convergence based on the sum of the squared error (i.e., the squared norm) of the difference between estimated parameters from successive iterations, "parameters_relative" which assesses convergence based on the relative change in mean squared error of the difference between estimated parameters from successive iterations, and "logLik_relative", which assess convergence based on the relative change in the PQL value between successive iterations. 
#' 
#' Although the first option is employed as a default, the second and third choices are often used in to assess convergence in other, likelihood-based optimization problems (e.g., Green, 1984).}

#' \item{tol: }{The tolerance value to use when assessing convergence.}

#' \item{final_maxit: }{The maximum number of iterations to do for the final estimation step after the PQL algorithm has converged. } 

#' \item{initial_betas_dampen: }{A dampening factor which can be used to reduce the magnitudes of the starting values obtained for the species-specific regression coefficients corresponding to the model matrix i.e., \code{betas}. This can either be set to a scalar, or a vector with length equal to the number of species i.e., \code{ncol(y)}. 
#' To elaborate, when starting values are not supplied as part of \code{start_params}, the function will attempt to obtain starting values based on fitting a stacked species distribution model. While this generally works OK, sometimes it can lead to bad starting values for the \code{betas}, due to the stacked species distribution model being severely overfitted. An *ad-hoc* fix to this is to dampen/shrink these initial values to be closer to zero, thus allowing the PQL estimation algorithm to actually "work". For instance, setting \code{initial_betas_dampen = 0.8} means the magnitudes of the staring values for the \code{betas} are reduced to 0.8 of their full values. This includes the species-specific intercepts. 
#' When \code{initial_betas_dampen} is a vector, then the dampening factor can vary with species.}

#' \item{subsequent_betas_dampen: }{A dampening factor which can be used to reduce the magnitudes of the values obtained for the species-specific regression coefficients corresponding to the model matrix i.e., \code{betas}, during the running of the PQL estimation algorithm. This can either be set to a scalar, or a vector with length equal to the number of species i.e., \code{ncol(y)}. 
#' To elaborate, during the PQL algorithm updates are made to the regression coefficients related to the combined matrix of basis functions, conditional on the regression coefficients corresponding to the model matrix. However, sometimes this updating can fails due to the latter producing non-sensible values to condition on e.g., due to severe overfitting in that component. If this occurs, then an *ad-hoc* second attempt is made, but conditioning instead on a dampened/shrunk set of the regression coefficients corresponding to the model matrix, which can often help. This amount of dampening is controlled by this argument. For instance, setting \code{subsequent_betas_dampen = 0.25} sets the magnitudes of the regression coefficients related to the model matrix to a quarter of their original size, including the intercepts. 
#' When \code{subsequent_betas_dampen} is a vector, then the dampening factor can vary with species.
#' Note that this argument *only* comes into play when the first attempt, which can be thought of as updating with \code{subsequent_betas_dampen = 1}, to update the regression coefficients associated with the combined matrix of basis functions fails. } 

#' \item{gam_method: }{When smoothing terms are included in the model, this controls the smoothing parameter estimation method. Defaults to "REML", which is maximum restricted likelihood estimation. However other options are available; please see the \code{method} argument in [mgcv::gam()] for the available options. In fact, note that [mgcv::gam()] defaults to using "GCV.Cp", which is based on generalized cross-validation. This is generally faster, but can be slightly more unstable, and hence why restricted maximum likelihood estimation is adopted as the default. }

#' \item{seed: }{The seed to use for the PQL algorithm. This is only applicable when the starting values are randomly generated, which be default should not be the case.}

#' \item{ridge: }{A additional ridge parameter that can be included to act as a ridge penalty when estimating the regression coefficients related to the covariates.}

#' \item{ziridge: }{A additional ridge parameter that can be included to act as a ridge penalty when estimating the regression coefficients related to the covariates for modeling the zero-inflation probabilities.}

#' \item{trace: }{If set to \code{TRUE} or \code{1}, then information at each iteration step of the outer algorithm will be printed. }

#' }
#' @param Sigma_control A list of parameters for controlling the fitting process for the "inner" estimation part of the CBFM pertaining to the community-level covariance matrices of the basis function regression coefficients. This should be a list with the following arguments:
#' \describe{
#' \item{rank: }{The rank of the community-level covariance matrices of the basis function regression coefficients. This either equals to a single scalar/character string equal to "full", or a vector or scalars/character strings (equal to "full") with length equal to how many of \code{B_space/B_time/B_spacetime} are supplied. If it is a single scalar or character string, then it is assumed that the same rank is used for all the community-level covariance matrices. 
#' 
#' The rank/s should be at least 2, although internal checks are also performed to assess if the rank is too large so that estimation is not feasible. If the character string "full" is used, then a full-rank covariance matrix is estimated.} 

#' \item{maxit: }{The maximum number of iterations for inner update of the community-level covariance matrices.} 

#' \item{tol: }{The tolerance value to use when assessing convergence. Convergence for the inner algorithm is assessed based on the norm of the difference between estimated parameters from successive iterations.} 

#' \item{method: }{The method by which to update the community-level covariance matrices. The current options are "REML" (default) which uses optimizing the Laplace approximated restricted maximum likelihood, "ML" which is the same but with the Laplace approximated (unrestricted) maximum likelihood, and "simple" which uses a fast large sample covariance update. *Note that the simple method is faster than the former, but is \emph{much} less accurate and we only recommend using it for pilot testing.*} 

#' \item{trace: }{If set to \code{TRUE} or \code{1}, then information at each iteration step of the inner algorithm will be printed.}

#' \item{custom_space: }{A custom, pre-specified community-level matrix for the spatial basis function regression can be supplied. If supplied, it must be a square matrix with dimension equal to the number of columns in \code{B_space}. Defaults to \code{NULL}, in which case it is estimated. Note as a side quirk, if this argument is supplied then a corresponding rank (as above) still has to be supplied, even though it is not used.}

#' \item{custom_time: }{A custom, pre-specified community-level covariance matrix for the temporal basis function regression can be supplied. If supplied, it must be a square matrix with dimension equal to the number of columns in \code{B_time}. Defaults to \code{NULL}, in which case it is estimated. Note as a side quirk, if this argument is supplied then a corresponding rank (as above) still has to be supplied, even though it is not used.}

#' \item{custom_spacetime: }{A custom, pre-specified community-level covariance matrix for the spatio-temporal basis function regression can be supplied. If supplied, it must be a square matrix with dimension equal to the number of columns in \code{B_spacetime}. Defaults to \code{NULL}, in which case it is estimated. Note as a side quirk, if this argument is supplied then a corresponding rank (as above) still has to be supplied, even though it is not used.}
#' }
#' @param G_control A list of parameters for controlling the fitting process for the "inner" estimation part of the CBFM pertaining to the so-called baseline between-species correlation matrices of the basis function regression coefficients. This should be a list with the following arguments:
#' \describe{
#' \item{rank: }{The rank of the between-species correlation matrices of the basis function regression coefficients. This either equals to a single scalar/character string equal to "full", or a vector or scalars/character strings (equal to "full") with length equal to how many of \code{B_space/B_time/B_spacetime} are supplied. If it is a scalar, then it is assumed that the same rank is used for all the correlation matrices. 
#' 
#' The rank/s should be at least 2, although internal checks are also performed to assess if the rank is too large so that estimation is not feasible. If the character string "full" is used, then a full-rank correlation matrix is estimated. 
#' 
#' Finally. note if a particular element in \code{G_control$structure} is set to "identity" or "homogeneous", then a corresponding element in \code{G_control$rank} must still be supplied, but it is ignored.} 

#' \item{structure: }{The structure to assume for the between-species correlation matrix if it is estimated. This either equals to a single character string or a vector of character strings with length equal to how many of \code{B_space/B_time/B_spacetime} are supplied. If it is a single string, then it is assumed that the same form is used for all the correlation matrices. The current options for each element are "unstructured" (default) which assumes an unstructured form subject to the corresponding element in \code{G_control$rank}, "identity" which assumes an unknown, estimated scalar multiplied by an identity matrix, and "homogeneous" which is the same as assuming an unstructured form but the entire matrix is the multiplied by an unknown, estimated scalar. 
#' *The latter two options should not be used unless you know what are doing in terms of what you want to achieve from the basis functions, especially as little checks are made in this setting to ensure parameter identifiability of the model in these settings!*} 

#' \item{nugget_profile: }{The sequence of values to try for calculating the nugget effect in each between-species correlation matrix. Please see details below for more information.} 

#' \item{maxit: }{The maximum number of iterations for inner update of the community-level covariance matrices.} 

#' \item{tol: }{The tolerance value to use when assessing convergence. Convergence for the inner algorithm is assessed based on the norm of the difference between estimated parameters from successive iterations.} 

#' \item{method: }{The method by which to update the community-level covariance matrices. The current options are "REML" (default) which uses optimizing the Laplace approximated restricted maximum likelihood, "ML" which is the same but with the Laplace approximated (unrestricted) maximum likelihood, and "simple" which uses a fast large sample covariance update. *Note that the simple method is faster than the former, but is \emph{much} less accurate and we only recommend using it for pilot testing.*} 

#' \item{trace: }{If set to \code{TRUE} or \code{1}, then information at each iteration step of the inner algorithm will be printed.}

#' \item{custom_space: }{A custom, pre-specified baseline between-species correlation matrix for the spatial basis function regression can be supplied. If supplied, it must be a square matrix with dimension equal to the number of columns in \code{B_space}. Defaults to \code{NULL}, in which case it is estimated. Note as a side quirk, if this argument is supplied then a corresponding \code{rank} and \code{structure} (as above) still has to be supplied, even though it is not used.}

#' \item{custom_time: }{A custom, pre-specified between-species correlation matrix matrix for the temporal basis function regression can be supplied. If supplied, it must be a square matrix with dimension equal to the number of columns in \code{B_time}. Defaults to \code{NULL}, in which case it is estimated. Note as a side quirk, if this argument is supplied then a corresponding \code{rank} and \code{structure} (as above) still has to be supplied, even though it is not used.}

#' \item{custom_spacetime: }{A custom, pre-specified between-species correlation matrix matrix for the spatio-temporal basis function regression can be supplied. If supplied, it must be a square matrix with dimension equal to the number of columns in \code{B_spacetime}. Defaults to \code{NULL}, in which case it is estimated. Note as a side quirk, if this argument is supplied then a corresponding rank (as above) still has to be supplied, even though it is not used.}
#' }

#' @param k_check_control A list of parameters for controlling [mgcv::k.check()] when it is applied to CBFMs involving smoothing terms for the measured covariates i.e., when smoothing terms are involved in \code{formula}. Please see [mgcv::k.check()] for more details on how this test works. This should be a list with the following two arguments:
#' \describe{
#' \item{subsample: }{If the number of observational units i.e., \code{nrow(y)} exceeds this number, then testing is done using a random sub-sample of units of this size.} 

#' \item{n.rep: }{How many re-shuffles of the residuals should be done in order to a P-value for testing. } 
#' }
#'
#'
#'
#' @details 
#' Community-level basis function models (CBFMs) are a class of joint species distribution models for spatio-temporal multivariate abundance data, which builds on the ideas of fixed rank kriging (FRK, Cressie and Johannesson, 2008; Zammit-Mangion and Cressie, 2017) and multivariate spatio-temporal mixed models (e.g., Bradley et al., 2018) and adapts them specifically for spatio-temporal multivariate abundance data in community ecology. CBFMs are a (closely connected) alternative to the popular spatio-temporal latent variable models (LVMs) approach for joint species distribution modeling, as available in a number of packages such as [Hmsc::Hmsc-package()] (Tikhonov et al., 2020), [gllvm::gllvm()] (Niku et al., 2019), and [boral::boral()] (Hui, 2016); see also Warton et al., (2015a,b), Thorson et al. (2016) and Ovaskainen and Abrego (2020) among others for general introductions to the use of LVMs in community ecology.  The key difference between LVMs and CBFMs is that rather than using a small number of latent variables (which are assumed to be random across observational units), CBFMs use a larger number of spatial and/or temporally-indexed basis functions that are specified \emph{a-priori} and remain fixed throughout the model fitting process (Hefley et al., 2017; Cressie et al., 2021). The randomness instead comes from species-specific regression coefficients related to these basis functions, which in turn induce spatio-temporal correlations within and between-species. 
#' 
#' In using a basis function approach, CBFMs can thus be framed as a type of generalized additive model (GAM, Guisan et al., 2002; Wood, 2017) . This in turn means CBFMs can leverage from the plethora of techniques that have been already developed for GAMs, with one notable benefit being that computationally, CBFMs tend to more efficient and scale better than many existing implementations of spatio-temporal LVMs (at least, at the time of writing).
#' 
#'
#' ## Some mathematics
#' Turning to some mathematical details, the CBFM in this package is characterized by the following mean regression model. For observational unit \eqn{i=1,\ldots,N} (e.g., with a corresponding space-time coordinate) and species \eqn{j=1,\ldots,m}, we have
#' 
#' \deqn{g(\mu_{ij}) = \eta_{ij} = x_i^\top\beta_j + b_i^\top a_j,}
#'
#' where \eqn{g(.)} is a known link function, \eqn{x_i} denotes a vector of predictors for unit \eqn{i} i.e., the \eqn{i}-th row from the created model matrix, \eqn{\beta_j} denotes the corresponding regression coefficients for species \eqn{j}, \eqn{b_i} denotes a vector of spatial and/or temporal basis functions for unit \eqn{i}, and \eqn{a_j} denotes the corresponding regression coefficients for species \eqn{j}. 
#' 
#' The vector of predictors \eqn{x_i} is created based on the \code{formula} and \code{data} arguments. Smoothing terms are permitted in \code{formula}, and these can be included in the same way as in [mgcv::gam.models()]; see also [mgcv::smooth.terms()]. Note smoothing terms in this context also permits the inclusion of (species-specific) random intercepts and slopes, through the use of the \code{s(..., bs = "re")}; please see [mgcv::random.effects()] and [mgcv::gam.vcomp()] for more details. These may be used, say, as a simple approach to account for nested sampling designs, multiple data sources/surveys etc..., although please note these random effects are specific to each species i.e., they are *not* random row effects as found in packages such as [boral::boral()] and [gllvm::gllvm()], and also are currently are not designed to draw the slopes from a common distribution as in [Hmsc::Hmsc-package()] or Pollock et al., (2014), say.  
#' 
#' When smoothing terms are included in the CBFM, a check of the smooth basis dimension and whether it is adequate is also automatically performed, courtesy of the [mgcv::k.check()] function; see that function's help file as well as [mgcv::choose.k()] for more general details. Furthermore, selection of smoothing terms as well as user specified knot values for basis construction are also possible; please see [mgcv::gam()], [mgcv::gam.selection()], and [mgcv::step.gam()] for more details. However, we warn the user that **some of the smoothers available as part of \code{mgcv} e.g., [mgcv::linear.functional.terms()], has not been fully tested for CBFM**, so some make may not work. If you encounter any problems, please post a Github issue on the CBFM repository!  
#' 
#' Next, the vector basis functions \eqn{b_i} is formed from the \code{B_space}, \code{B_time} and \code{B_spacetime} arguments. At least one of these arguments must be supplied. As an example, suppose we wish to fit a CBFM with spatial and temporal basis functions included in an additive manner. Then only \code{B_space} and \code{B_time} should be supplied, in which case the mean regression model for the CBFM can be written as:
#' 
#' \deqn{g(\mu_{ij}) = \eta_{ij} = x_i^\top\beta_j + b_{i,space}^\top a_{j,space} + b_{i,time}^\top a_{j,time},}
#'
#' where \eqn{b_i = (b_{i,space}, b_{i,time})} and \eqn{a_j = (a_{j,space}, a_{j,time})}. If purely spatial or temporal multivariate abundance data is recorded, then one should only supply \code{B_space} and \code{B_time}, respectively, and the above mean regression model is simplified accordingly. 
#' 
#' As another example, suppose we wish to include spatio-temporal basis functions (which are formed from a tensor-product). Then only \code{B_spacetime} should be supplied, in which case the mean regression model for the CBFM can be written as:
#' 
#' \deqn{g(\mu_{ij}) = \eta_{ij} = x_i^\top\beta_j + b_{i,spacetime}^\top a_{j,spacetime},}
#'
#' where \eqn{b_i = b_{i,spacetime}} and \eqn{a_j = a_{j,spacetime}}. More details and recommendations on how to construct this basis functions, including the tensor-product mentioned above, are provided later on. 
#' 
#' Note for zero-inflated distributions, it is the *mean of the non-zero inflated component that is modeled and not the mean of the entire distribution.*
#' 
#' 
#' ** Remark on flavors and choices of CBFMs:** 
#' For spatio-temporal multivariate abundance data, we have found that the above two examples are usually the two most appropriate "flavors" CBFM to apply, although the precise form should of course depend upon the precise interpretation and question/s of interest. As seen above, we may want the separate sets of spatial and temporal basis functions to be included in an additive manner (this is analogous to an LVM where separate spatial LVs and temporal LVMs are added together), or have a single set of spatio-temporal basis functions formed from a tensor-product say (this is analogous to an LVM with one set of spatio-temporal LVs), or have a combination of the two where (say) the basis functions included in \code{B_space} and \code{B_time} are accounting for correlations on a coarse scale while the basis functions included in \code{B_spacetime} are accounting for resolutions on a fine scale. We refer the interested reader to Thorson et al., (2016) and Thorson (2019) for examples of similar kinds of constructs and flavors within the LVM framework. 
#' 
# #' We also point out that this package only implements one possible version of a wider class of CBFMs; other potentially superior versions e.g., spatial basis functions with temporally varying corresponding regressions coefficients, are possible under the CBFM framework, but are far outside the scope of this package (sorry!).      
#'    
#' In principle, it is possible to employ a more data-driven approach such as cross-validation or information criteria to choose which type of and/or the number of basis functions to include in the arguments \code{B_space/B_time/B_spacetime}. This similar to choosing the number of latent variables in a LVM. We refer the reader to [mgcv::choose.k()] as some of the advice provided there may be applicable to CBFMs e.g., using residual analysis to informally check whether an increase the number of spatial and/or temporal basis functions is required, and echo a sentiment written there (while acknowledging things are more tricky with spatial and/or temporal basis functions, and for discrete responses!): 
#' 
#' *"So, exact choice of \eqn{k} (the number of basis functions in our situation) is not generally critical: it should be chosen to be large enough that you are reasonably sure of having enough degrees of freedom to represent the underlying 'truth' reasonably well, but small enough to maintain reasonable computational efficiency. Clearly 'large' and 'small' are dependent on the particular problem being addressed."* 
#' 
#' 
#' ## Basis function coefficients
#' In the CBFM, the basis functions \eqn{b_i} are specified \emph{a-priori} and remain fixed throughout the fitting process. Instead, it is the associated species-specific regression coefficients \eqn{a_j} which are assumed to be random. In this package, we assume follow a multivariate normal distribution as follows:
#'   
#' \deqn{(a_1,\ldots,a_m) \sim N(0, kronecker(G, \Sigma)),} 
#' 
#' where \eqn{G} is a so-called baseline between-species correlation matrix, \eqn{\Sigma} is the community-level covariance matrix for the basis function regression coefficients, and \eqn{kronecker(\cdot)} is the Kronecker product operator. When multiple sets of basis functions are included, then this carries over. For instance, in the example above involving a CBFM with separate spatial and temporal basis functions, with only \code{B_space} and \code{B_time} supplied, we have
#'   
#' \deqn{(a_{1,space},\ldots,a_{m,space}) \sim N(0, kronecker(G_{space}, \Sigma_{space})),} 
#' and
#' \deqn{(a_{1,time},\ldots,a_{m,time}) \sim N(0, kronecker(G_{time}, \Sigma_{time})).} 
#' 
#' To reduce the number of parameters needed to be estimated in both the \eqn{G}'s and \eqn{\Sigma}'s, a rank-reduced form is by default adopted in both (inspired by the rank-reduced covariance matrices characterizing LVMs). In detail, we assume \eqn{G = \Lambda_{G}\Lambda_{G}^top + \kappa_G I_m} where \eqn{\Lambda_{G}} is an \eqn{m \times d_G} loading matrix and \eqn{\kappa_G > 0} is a nugget effect, with \eqn{I_m} being an identity matrix with dimension \eqn{m}. The quantity \eqn{d_G << m} is the rank, and similar to LVMs we often choose this to be small relative to the number of species. Similarly, we have \eqn{\Sigma = \Lambda_{\Sigma}\Lambda_{\Sigma}^top + \kappa_{\Sigma} I_{q}}, where \eqn{\Lambda_{\Sigma}} is an \eqn{q \times d_{\Sigma}} loading matrix, and \eqn{q} is the number of basis functions included in the model. When multiple sets of basis functions are included e.g., both \code{B_space} and \code{B_time}, then rank-reduced structures are used accordingly. 
#'
#' The ranks \eqn{d_G} and \eqn{d_{\Sigma}} are chosen generally to be smaller than the number of species and basis functions, respectively. Generally speaking, provided the rank/s is large enough, then results should not depend much on their choice. The nugget effect is included to ensure that resulting rank-reduced forms of \eqn{G} and \eqn{\Sigma} remain positive definite and generally a bit more stable. They can also have the interpretation of adjusting for the relative strength of correlation between species, say (Shirota, 2019).
#' 
#' It is also possible to assume full-rank covariance and correlation matrices for the \eqn{G}'s and \eqn{\Sigma}'s, via the character string "full" in \code{G_control$rank} and \code{Sigma_control$rank}. Other structures may gradually be implemented in future versions of the package.
#' 
#' It is also possible for the user to specific their own custom matrices for the community-level basis function covariance matrices \eqn{\Sigma} and/or the baseline between-species correlation matrices \eqn{G}. These are supplied through the \code{custom_space/custom_time/custom_spacetime} arguments as part of \code{Sigma_control} and \code{G_control}, respectively. In \eqn{\Sigma} is supplied, then only the corresponding \eqn{G} is estimated, and note importantly that the corresponding \eqn{G} is now estimated to be a baseline between-species *covariance* matrix instead of correlation matrix (it is no longer constrained to be correlation matrix). Conversely, if \eqn{G} is supplied, then only the corresponding \eqn{\Sigma} is estimated. There is rarely any reason why one wish to supply both \eqn{\Sigma} and \eqn{G} simultaneously.
#' 
#' Using custom, pre-specified structures for \eqn{\Sigma} may be useful if the corresponding basis functions should be "equipped" with a particular structure. For example, if one or more sets of basis functions are constructed from, say, the [mgcv::smooth.terms()] package (although see later on some defaults for constructing basis functions), then they are usually equipped with a corresponding, fixed penalty (inverse \eqn{\Sigma}) matrix. Perhaps a more common example is, say, if a user wanted to include species-specific random intercepts for time. Then one would set \code{B_time} as a model-matrix formed from treating time as a factor, and set \code{Sigma_control$custom_time} to an identity matrix with dimension equal to \code{ncol(B_time)}; see Example 3b in the help file below. After fitting, the estimated diagonal elements of \eqn{G_{time}} would then be the species-specific variance components for the time random intercept, while the off-diagonal elements may be interpreted as the corresponding between-species covariation.  
#' 
#'
#' ## Distributions
#' 
#' The following response distributions are permitted: 
#' \describe{
#' \item{\code{betalogitfam()}: }{Beta distribution using a logit link. The corresponding mean-variance relationship is given by \eqn{V = \mu(1-\mu)/(1+\phi)}, where \eqn{\mu} denotes the mean and \eqn{\phi} is the dispersion parameter.}

#' \item{\code{binomial(link = "logit")}: }{Binomial distribution using a logit link. The corresponding mean-variance relationship is given by \eqn{V = N_{trial}\mu(1-\mu)}, where \eqn{\mu} denotes the mean and \eqn{N_{trial}} is the trial size.}

#' \item{\code{Gamma(link = "log")}: }{Gamma distribution using a log link. The corresponding mean-variance relationship is given by \eqn{V = \phi\mu^2}, where \eqn{\mu} denotes the mean and \eqn{\phi} is the dispersion parameter.}

#' \item{\code{gaussian(link = "identity")}: }{Gaussian or normal distribution using an identity link. The corresponding mean-variance relationship is given by \eqn{V = \phi}, where \eqn{\phi} is the dispersion parameter.}

#' \item{\code{poisson(link = "log")}: }{Poisson distribution using a log link. The corresponding mean-variance relationship is given by \eqn{V = \mu}, where \eqn{\mu} denotes the mean.}

#' \item{\code{nb2()}: }{Negative binomial distribution with the log link. The corresponding mean-variance relationship is given by \eqn{V = \mu + \phi\mu^2}, where \eqn{\mu} denotes the mean and \eqn{\phi} is the dispersion parameter.}

#' \item{\code{tweedielogfam()}: }{Tweedie distribution using log link. The corresponding mean-variance relationship is given by \eqn{V = \phi\mu^{\rho}}, where \eqn{\mu} denotes the mean, \eqn{\phi} is the dispersion parameter, and \eqn{\rho} is the power parameter.}

#' \item{\code{zipoisson()}: }{Zero-inflated Poisson distribution using a log link for the Poisson part is permitted. This partial mass function of the distribution is given by \eqn{f(y) = \pi I(y=0) + (1-\pi) f_{pois}(y)}, where \eqn{\pi} is the probability of being in the zero-inflation component, while \eqn{f_{pois}(y)} is the usual Poisson distribution. The mean of the Poisson distribution is modeled against covariates and basis functions, while the probability of zero-inflation is also modeled against covariates only via \code{ziformula}. In the case of the latter, a logit link function is used.}

#' \item{\code{zinb2()}: }{Zero-inflated negative binomial distribution using a log link for the negative binomial part is permitted. The partial mass function of the distribution is given by \eqn{f(y) = \pi I(y=0) + (1-\pi) f_{NB}(y)}, where \eqn{\pi} is the probability of being in the zero-inflation component, while \eqn{f_{NB}(y)} is the usual negative binomial distribution. The mean of the negative binomial distribution is modeled against covariates and basis functions, while the probability of zero-inflation is also modeled against covariates only via \code{ziformula}. In the case of the latter, a logit link function is used.}

#' \item{\code{ztpoisson()}: }{Zero-truncated Poisson distribution using a log link. The partial mass function of the distribution is given by \eqn{f(y) = f_{pois}(y)/(1-f_{pois}(0)}) where \eqn{f_{pois}(y)} is the usual Poisson distribution as described above. The mean of the Poisson distribution is modeled against covariates and basis functions.}
#' 
#' \item{\code{ztnb2()}: }{Zero-truncated negative binomial distribution using a log link. The partial mass function of the distribution is given by \eqn{f(y) = f_{NB}(y)/(1-f_{NB}(0)}) where \eqn{f_{NB}(y)} is the usual negative binomial distribution as described above. The mean of the negative binomial distribution is modeled against covariates and basis functions.}
#' 
#' Hurdle CBFMs are also possible; please see [makeahurdle()] for more information.
#' }
#' 
#' Missing (\code{NA}) values are permitted as part of the response matrix. These are simply passed over during the fitting process i.e., effectively the action [stats::na.omit()] is taken.
#'
#' 
#' 
#' ## Constructing basis functions
#' The CBFM approach to spatio-temporal joint species distribution modeling relies on the inclusion of the \emph{pre-specified} spatial and/or temporal basis functions to account for spatio-temporal correlations within and between species (see Hefley et al., 2017, for a general overview of using basis functions to model autocorrelation in ecological data). Currently, the package does not provide default arguments to use for this, and this is deliberately as we wish to compel the user to work and think a bit harder on designing the right basis functions for use when CBFMs to their particular analysis.
#' 
#' At the same time, it would be remiss not to provide some brief recommendations based on previous experience, and we do so below. Please also see the examples later on for some more concrete applications.
#' \describe{
#' \item{\code{B_space}: }{We have found that the resolution adaptive thin-plate spline basis functions (Tzeng and Huang, 2018), as implemented in [autoFRK::mrts()], work fairly well here as spatial basis functions. They are simple to use and require the user to only supply the number of basis functions, which itself is tied to the resolution at which the user wants to model their spatial correlations. For spatial multivariate abundance data, we have usually found that 50 or less spatial basis functions of such type suffices.
#' 
#' Another option for spatial basis functions is to use [FRK::auto_basis()], which produces basis functions that are sparse in design but consequently require many more in number compared to the resolution adaptive thin-plate splines mentioned above. This approach is more customizable however, with the choice of resolutions, basis function centers, and aperture among other choices; please see Zammit-Mangion and Cressie (2017) and Wilke et al. (2019) for more details.}

#' \item{\code{B_time}: }{Both of the approaches mentioned above for \code{B_space} can also be applied here, although with temporal basis functions we have generally found the approach implemented in [FRK::auto_basis()] to work satisfactorily in many cases, given their customizability and sparsity (local support). The resolution-adaptive thin-plate spline basis functions approach, when applied solely in the 1-D temporal dimension, can produce long-term temporal trends that may be undesirable.}

#' \item{\code{B_spacetime}: }{A general starting point for constructing spatio-temporal basis functions is to make use of a tensor-product form (analogous to [mgcv::te()]). That is, after constructing a set of spatial and a set of temporal basis functions, we can use the [tensorproduct()] function to construct the tensor-product and include them \code{B_spacetime}. 
#' 
#' It is recommended that both "ingredient" basis functions used in the tensor-product are sparse in design to facilitate computation e.g., as implemented in [FRK::auto_basis()]; see the \code{FRK} package as well as Wilke et al. (2019) for some examples. Also, we recommend you do not use these same "ingredient" basis functions in the \code{B_space} and \code{B_time} arguments, as this may lead to overfitting. Put another way, and as hinted at previously, if \code{B_spacetime} is supplied at the same time as either \code{B_space} and/or \code{B_time} is supplied, then they should generally be constructed to model different resolutions/scales of the spatio-temporal correlation.} 
#' }
#' 
#' 
#' ## A note on estimation and inference
#' Because CBFMs uses a basis function approach to model spatio-temporal correlations between and within species, then they can be thought of as a type of GAM. Similar to a common implementation of GAMs then, this package uses a maximized penalized quasi-likelihood (PQL) approach for estimation and inference (Breslow and Clayton, 1993; Wood, 2017), while the baseline between-species correlation and community-level covariance matrices are by default estimated via Laplace approximated restricted maximum likelihood estimation (Wood, 2011; Wood, 2017). Currently, CBFM makes use of both the machinery available in the [mgcv] package (Wood, 2017) as well as that of Template Model Builder (TMB, Kristensen et al., 2016) to facilitate this entire algorithm. After the PQL algorithm has converged, a final estimation step is performed purely to update regression/basis function coefficients (and dispersion parameters if appropriate).
#' 
#' If \code{start_params} is not supplied, then CBFM attempts to obtain starting values based on fitting an appropriate stacked species distribution model. This generally works alright, but can sometimes fail badly e.g., if the stacked species distribution model severely overfits for one or more species. A tell-tale sign of when it occurs is if from the returned CBFM fit, the estimates of regression coefficients corresponding to the spatial and/or temporal basis functions i.e., \code{basis_effects_mat}, are extremely close to zero for these problematic species. There are not always easy fixes for such situations (as it may reflect an underlying, perhaps intriguing feature of the proposed model for the predictors, data, or it may genuinely be that the stacked species distribution model is already fitting incredibly well!). One *ad-hoc* fix is available through \code{control$initial_betas_dampen}, but it is not guaranteed to work.  
#' 
#' Standard errors and resulting inferential tools like confidence intervals are based on the approximate large sample distribution of the regression coefficients, and use the so-called Bayesian posterior covariance matrix for the coefficients, similar to (but not as sophisticated as!) what is provided  [mgcv::summary.gam()]. Please note that **all standard errors and thus inference are currently computed without considering uncertainty in estimation of covariance \eqn{\Sigma} and correlation matrices \eqn{G}. They can lead to standard errors that are potentially too small, so please keep this in mind.** 
#' 
#' Also, the current estimation approach **does not provide uncertainty quantification of \eqn{\Sigma} and \eqn{G}**, does not provide uncertainty estimates in the smoothing parameter. This is in line with the main aims of this CBFM package, which are tailored more towards estimation and inference of regression coefficients and spatio-temporal prediction, in a relatively computationally efficient and scalable manner. Future versions of package may seek to rectify this, but for now... apologies!  
#' 
#' 
#' @return An object of class \code{CBFM} which includes the following components, not necessarily in the order below, and as appropriate:
#' \item{call: }{The matched function call.}

#' \item{family: }{The supplied response distribution i.e., family function, to be used in the model.}

#' \item{y, data, trial_size: }{The supplied response matrix, covariate information data frame, and trial size(s).}

#' \item{formula: }{The supplied symbolic description of the model matrix to be created.}

#' \item{ziformula: }{The supplied symbolic description of the model matrix to be created for modeling the probability of zero-inflation.}

#' \item{B: }{The full matrix basis functions i.e., basically the result of \code{cbind(B_space, B_time, B_spacetime)}.}

#' \item{which_B_used: }{A vector of length three, indicating which of \code{B_space, B_time, B_spacetime} was supplied. For example \code{which_B_bused = c(1,0,0)} implies only \code{B_space} was supplied.}

#' \item{num_B_space: }{The number of spatial basis functions supplied i.e., \code{ncol(B_space)}.} 

#' \item{num_B_time: }{The number of temporal basis functions supplied i.e., \code{ncol(B_time)}.} 

#' \item{num_B_spacetime: }{The number of spatio-temporal basis functions supplied i.e., \code{ncol(B_spacetime)}.} 

#' \item{num_B: }{The total number of basis functions supplied i.e., \code{ncol(cbind(B_space, B_time, B_spacetime))}.}

#' \item{converged: }{Indicates whether or not the PQL estimation algorithm converged. Note results may still be outputted (and sensible?!) even if this is \code{FALSE}.}

#' \item{logLik: }{The value of the log-likelihood (excluding the quadratic penalty term in the PQL) for the fitted model.}

#' \item{logLik_perspecies: }{The value of each species' contribution to the log-likelihood (excluding the quadratic penalty term in the PQL) for the fitted model.}

#' \item{deviance: }{The deviance for the fitted model. Note the deviance calculation here does *not* include the quadratic term of the PQL.}

#' \item{deviance_perspecies: }{The value of each species' contribution to the deviance of the fitted model. Note the deviance calculation here does *not* include the quadratic term of the PQL.}

#' \item{null_deviance: }{The null deviance i.e., deviance of a stacked model (GLM) where each species model contains only an intercept. Note the deviance calculation here does *not* include the quadratic term of the PQL}

#' \item{null_deviance_perspecies: }{The value of each species' contribution to the null deviance. Note the deviance calculation here does *not* include the quadratic term of the PQL}

#' \item{deviance_explained: }{The *percentage* of null deviance explained by the model. In community ecology this is typically not very high; please see [varpart()] for more capacity to perform variance partitioning in a CBFM. This is constructed by examining \code{deviance} to \code{null_deviance}.}

#' \item{deviance_explained_perspecies: }{The *percentage* of null deviance explained by the model, on a per-species basis. This is constructed by examining \code{deviance_perspecies} to \code{null_deviance_perspecies}.}

#' \item{pql_logLik: }{The value of the PQL i.e., the likelihood plus the quadratic penalty term, upon convergence.}

#' \item{edf/edf1: }{A matrix of estimated degrees of freedom for each model parameter in \code{formula}. The number of columns of the matrix should be equal to the number of species i.e., \code{ncol(y)}. Penalization means that many of these are less than one. \code{edf1} is an alternative estimate of EDF. Note these values are pulled straight from the GAM part of the estimation algorithm, and consequently may only be *very* approximate. }

#' \item{ziedf/ziedf1: }{A matrix of estimated degrees of freedom for each model parameter in \code{ziformula}. The number of columns of the matrix should be equal to the number of species i.e., \code{ncol(y)}. Penalization means that many of these are less than one. \code{edf1} is an alternative estimate of EDF. Note these values are pulled straight from the GAM part of the estimation algorithm, and consequently may only be *very* approximate. }

#' \item{pen_edf: }{A list with each element containing a vector of the estimated degrees of freedom associated with each smoothing term in \code{formula}. The length of the list should be equal to the number of species i.e., \code{ncol(y)}. Note these values are pulled straight from the GAM part of the estimation algorithm, and consequently may only be *very* approximate.}

#' \item{zipen_edf: }{A list with each element containing a vector of the estimated degrees of freedom associated with each smoothing term in \code{ziformula}. The length of the list should be equal to the number of species i.e., \code{ncol(y)}. Note these values are pulled straight from the GAM part of the estimation algorithm, and consequently may only be *very* approximate.}

#' \item{k_check: }{A list resulting from the application of [mgcv::k.check()], used as a diagnostic test of whether the smooth basis dimension is adequate for smoothing terms included in \code{formula}, on a per-species basis. Please see [mgcv::k.check()] for more details on the test and the output. Note that if no smoothing terms are included in \code{formula}, then this will be a list of \code{NULL} elements.}

#' \item{vcomp: }{A list with length equal to \code{ncol(y)}, where each element contains a vector of the estimated variance components (as standard deviations) associated with the smoothing terms included in \code{formula}. This output is only really useful when one or more of the smoothing terms were included in the CBFM as species-specific intercepts/slopes (see [mgcv::random.effects()] for more details), in which case the corresponding values in \code{vcomp} are the estimated variance components (estimated standard deviations to be precise) associated with these random effects; see [mgcv::random.effects()] and [mgcv::gam.vcomp()] for more details on the one-to-one relationship between smoothing parameters in GAMs and variance components in mixed models. Note that if no smoothing terms are included in \code{formula}, then this will be a list of \code{NULL} elements.}

#' \item{betas: }{The estimated matrix of species-specific regression coefficients corresponding to the model matrix created. The number of rows in \code{betas} is equal to the number of species i.e., \code{ncol(y)}.}

#' \item{zibetas: }{The estimated matrix of species-specific regression coefficients corresponding to the model matrix created for the zero-inflation component. The number of rows in \code{zibetas} is equal to the number of species i.e., \code{ncol(y)}. Recall a logit scale to model the probability of zero-inflation i.e., \eqn{log(\pi/(1-\pi))} is regressed against covariates.}

#' \item{basis_effects_mat: }{The estimated matrix of species-specific regression coefficients corresponding to the combined matrix of basis functions. The number of rows in \code{basis_effects_mat} is equal to the number of species i.e., \code{ncol(y)}.}

#' \item{dispparam: }{The estimated vector of species-specific dispersion parameters, for distributions which require one. }

#' \item{powerparam: }{The estimated vector of species-specific power parameters, for distributions which require one. }

#' \item{linear_predictors: }{The estimated matrix of linear predictors. Note that for zero-inflated distributions, the mean of the count component is modeled in CBFM, and the function returns the linear predictors corresponding to this count component in the CBFM. Similarly, for zero-truncated count distributions, the mean of the base count distribution is modeled in CBFM, and the function returns the linear predictors corresponding to this base count distribution in the CBFM (and \code{NA} values for elements corresponding to zero counts in \code{object$y}.)}

#' \item{fitted: }{The estimated matrix of fitted mean values. Note that for zero-inflated distributions, while the mean of the count component is modeled in CBFM, the fitted values are the *actual expected mean values* i.e., it returns estimated values of \eqn{(1-\pi_{ij})\mu_{ij}} where \eqn{\pi_j} is the species-specific probability of zero inflation and \eqn{\mu_{ij}} is the mean of the count component. 
#' Similarly, for zero-truncated count distributions, while the mean of the base count distribution is modeled in CBFM, the fitted values are the *actual expected mean values* i.e., it returns estimated values of \eqn{\mu_{ij}/(1-p(0,\mu_{ij}))} where \eqn{\mu_{ij}} is the mean of the base count distribution component and \eqn{p(0,\mu_{ij})} generically denotes the probability of observing a zero count for the base count distribution (and it returns \code{NA} values for elements corresponding to zero counts in \code{object$y}.)
#' }

#' \item{Sigma_space/Loading_Sigma_space/nugget_Sigma_space: }{The estimated community-level covariance matrix/loadings/nugget effect associated with the spatial basis functions, if \code{B_space} is supplied. Note if \code{Sigma_control$custom_space} was supplied, then \code{Sigma_space} would be directly the supplied matrix, while the loading and nugget arguments are set to \code{NULL}.}

#' \item{G_space/Loading_G_space/nugget_G_space: }{The estimated baseline between-species correlation matrix/loadings/nugget effect associated with the spatial basis functions, if \code{B_space} is supplied. Note if \code{Sigma_control$custom_space} was supplied, then a covariance matrix is estimated instead. Note if \code{G_control$custom_space} was supplied, then \code{G_space} would be directly the supplied matrix, while the loading and nugget arguments are set to \code{NULL}.}

#' \item{Sigma_time/Loading_Sigma_time/nugget_Sigma_time: }{The estimated community-level covariance matrix/loadings/nugget effect associated with the temporal basis functions, if \code{B_time} is supplied. Note if \code{Sigma_control$custom_time} was supplied, then \code{Sigma_time} would be directly the supplied matrix, while the loading and nugget arguments are set to \code{NULL}.}

#' \item{G_time/Loading_G_time/nugget_G_time: }{The estimated baseline between-species correlation matrix/loadings/nugget effect associated with the temporal basis functions, if \code{B_time} is supplied. Note if \code{Sigma_control$custom_tiem} was supplied, then a covariance matrix is estimated instead. Note if \code{G_control$custom_time} was supplied, then \code{G_time} would be directly the supplied matrix, while the loading and nugget arguments are set to \code{NULL}.}

#' \item{Sigma_spacetime/Loading_Sigma_spacetime/nugget_Sigma_spacetime: }{The estimated community-level covariance matrix/loadings/nugget effect associated with the spatio-temporal basis functions, if \code{B_spacetime} is supplied. Note if \code{Sigma_control$custom_spacetime} was supplied, then \code{Sigma_spacetime} would be directly the supplied matrix, while the loading and nugget arguments are set to \code{NULL}.}

#' \item{G_spacetime/Loading_G_spacetime/nugget_G_spacetime: }{The estimated baseline between-species correlation matrix/loadings/nugget effect associated with the spatio-temporal basis functions, if \code{B_spacetime} is supplied. Note if \code{Sigma_control$custom_spacetime} was supplied, then a covariance matrix is estimated instead. Note if \code{G_control$custom_spacetime} was supplied, then \code{G_spacetime} would be directly the supplied matrix, while the loading and nugget arguments are set to \code{NULL}.}

#' \item{stderrors: }{The supplied argument for \code{stderrors} i.e., whether standard errors were calculated.}

#' \item{covar_components: }{If \code{stderrors = TRUE}, then a list containing with the following components: 
#' 1) \code{topleft}, which is a matrix corresponding to the top-left block of the full Bayesian posterior covariance matrix. The top-left block specifically relates to the regression coefficients associated with the measured predictors i.e., the covariance matrix associated with \code{object$betas}, along with \code{object$zibetas} if zero-inflated distributions are used;
#' 2) \code{topright}, which is a matrix of the top-right block of the full Bayesian posterior covariance matrix. The top-right block specifically relates to the cross-covariance of the regression coefficients associated with the measured predictors (plus the coefficients associated with the probabilities of zero-inflation) and the basis functions i.e., the cross-covariance matrix between \code{object$betas} and \code{object$basis_effects_mat}; 
#' 3) \code{bottomright}, which is a matrix containing components of the bottom-right block of the full Bayesian posterior covariance matrix. The bottom-left block specifically relates to the regression coefficients associated with the basis functions i.e., the covariance matrix associated with \code{object$basis_effects_mat}.
#' 
#' Please use the [summary.CBFM()] function to obtain standard errors and confidence interval limits in a (slightly) more user-friendly form.}

#' \item{time_taken: }{The time taken to run the PQL estimation algorithm, in seconds. This is calculated simply using differences in calls of [base::proc.time()].}
#' 
#' 
#' @details # Warning
#' 1. CBFMs are designed for \emph{spatio-temporal} multivariate abundance data, such that you can sensibly construct basis functions from the space-time coordinate of each observational unit. **Please do not use them for data that are not spatially or temporally indexed**. We recommend you fit standard LVMs in those scenarios, such that made available in [gllvm::gllvm()], [boral::boral()], and [Hmsc::sampleMcmc()].
#' 
#' 2. Not for some distributions it is not the mean of the entire distribution which is modeled. For example, in zero-inflated distributions it is the mean of the non-zero-inflated component that is modeled with the regression model described above. In zero-truncated distributions, it is the mean of the base count distribution that is modeled with the regression model described above.
#' 
#' 3. Not all (in fact, not many) of the smoothing available that are available in [mgcv::gam.models()] have been fully tested out, so please be aware that some make not work well, if at all! 
#' 
#' 4. As mentioned above, all standard errors and thus inference are currently computed without considering uncertainty in estimation of covariance \eqn{\Sigma} and correlation matrices \eqn{G}, as well as the any dispersion/power parameters, analogous to default settings in [mgcv::summary.gam()]. This can lead to standard errors that are potentially too small, so please keep this in mind. Also, the current estimation approach does not provide uncertainty quantification of \eqn{\Sigma} and \eqn{G}. 
#' 
# #' Indeed, the "strength" of the CBFM approach (especially with the current approach to estimation) is its competitive predictive performance relative to computation efficiency and scalability; **estimates of \eqn{\Sigma} and \eqn{G} may not be too reliable.**
#'
#'
#' @details # CBFM is not working for my data?!
#' Once you have finished grumbling about the package and its developer, please brew some tea and make yourself comfy...debugging takes a while!
#' 
#' As with any real-life statistical modeling problem, it is almost impossible to determine what the source of the issue is without looking at the data and specific application first hand. Therefore, we can only provide some general avenues to pursue below as a first step towards making CBFM run on your data, and of course we can not guarantee that the output produced from this debugging makes any ecological sense!
#' 
#' * An initial thing to try is to bump up the number of inner iterations i.e., \code{control$inner_maxit} from the default of 1 to something like 10, 20 or even higher. This will give more opportunities for the inner estimation component of the PQL algorithm (where all regression and smoothing coefficients are updated) to try converge to a stable point, before proceeding to estimating the other model parameters. 
#' 
#' * Sometimes the starting values that CBFM constructs are not that great! A common situation where this occurs is when smoothers are employed in \code{formula} and the data are (extremely) overdispersed or show signs of complete or quasi-separation. A simple way to try and break out of bad automated starting values is to make use of the \code{control$initial_betas_dampen} argument, which as the name suggests, dampens the starting estimated coefficients from potentially extreme magnitudes, and can facilitate the underlying PQL estimation algorithm to "get going".
#' 
#' * Alternatively, supplying your own "wisely chosen" starting values is never a bad thing, plus it can often help to speed up the fitting process. In our experience, often a decent way to obtain starting values is to fit stacked GAMs using [mgcv::gam()] with the same formula as you will use in \code{formula}, plus smoothing terms to account for space and/or time. Some template code is provided as follows:
#' ```
#' manygam <- foreach::foreach(j = 1:num_spp) %dopar%
#'     gam(response ~ s(temp) + s(depth) + s(chla) + s(O2) + s(x,y), data = data.frame(response = simy_train[,j], dat_train), family = xxx)
#' start_params = list(betas = t(sapply(manygam, coef)[1:37,])) # Or as appropriate the number of coefficients excluding the spatial-temporal smoothing terms
#' ```
#' 
#' * Analogously, the argument \code{control$subsequent_betas_dampen} can be used to dampen the values obtained for the species-specific regression coefficients during subsequent running of the PQL estimation algorithm. Basically, it is an *ad-hoc* solution to when updates potentially fail due to e.g., overfitting causing coefficients to become extremely large in magnitude.  
#' 
#' * If, after multiple attempts at debugging and you CBF'd anymore (pun-intended), then you can post the issue up on [CBFM Github page](https://github.com/fhui28/CBFM) *if* you think there is a genuine bug in the package. Otherwise, you can email the authors of this package on potential general statistical modeling questions, although we may not be able to get to them soon let along have the time to get to them at all (we are paid to be statistical consultants...apologies in advance!). Please do not post general statistical modeling issues on Github issues, as they will likely be ignored or deleted without prior consent.    
#'  
#'
#' @author Francis K.C. Hui <fhui28@gmail.com>, Chris Haak
#'
#'
#' @references
#' Bradley, J. R., Holan, S. H., and Wikle, C. K. (2018). Computationally efficient multivariate spatio-temporal models for high-dimensional count-valued data (with discussion). Bayesian Analysis, 13, 253-310.
#' 
#' Breslow, N. E., and Clayton, D. G. (1993). Approximate inference in generalized linear mixed models. Journal of the American statistical Association, 88, 9-25.
#' 
#' Brooks, M. E., Kristensen, K., Van Benthem, K. J., Magnusson, A., Berg, C. W., Nielsen, A., and Bolker, B. M. (2017). glmmTMB balances speed and flexibility among packages for zero-inflated generalized linear mixed modeling. The R journal, 9, 378-400.
#' 
#' Cressie, N., and Johannesson, G. (2008). Fixed rank kriging for very large spatial data sets. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 70, 209-226.
#' 
#' Green, P. J. (1984). Iteratively reweighted least squares for maximum likelihood estimation, and some robust and resistant alternatives. Journal of the Royal Statistical Society: Series B (Methodological), 46, 149-170.
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
#' Pollock, L. J., Tingley, R., Morris, W. K., Golding, N., O'Hara, R. B., Parris, K. M., Vesk, P. A., and McCarthy, M. A. (2014). Understanding co‐occurrence by modelling species simultaneously with a Joint Species Distribution Model (JSDM). Methods in Ecology and Evolution, 5, 397-406.
#' 
#' Shirota, S., Gelfand, A. E., & Banerjee, S. (2019). Spatial joint species distribution modeling using Dirichlet processes. Statistica Sinica, 29, 1127-1154.
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
#' Warton, D. I., Blanchet, F. G., O'Hara, R., Ovaskainen, O., Taskinen, S., Walker, S. C., and Hui, F. K. C. (2016). Extending joint models in community ecology: A response to Beissinger et al. Trends in ecology & evolution, 31, 737-738.
#'
#' Wood, S. N. (2011). Fast stable restricted maximum likelihood and marginal likelihood estimation of semiparametric generalized linear models. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 73, 3-36.
#' 
#' Wood, S. N. (2017). Generalized additive models: An introduction with R. CRC press.
#' 
#' Zammit-Mangion, A., and Cressie, N. (2017). FRK: An R package for spatial and spatio-temporal prediction with large datasets. arXiv preprint arXiv:1705.08105.
#'
#' Cressie, N., Sainsbury-Dale, M., and Zammit-Mangion, A. (2021). Basis-Function Models in Spatial Statistics. Annual Review of Statistics and Its Application, 9.
#'
#' @seealso [corX()] for calculating between-species (cross-)correlations due to measured covariates, [corB()] for calculating residual between-species (cross-)correlations due to the basis functions, [fitted.CBFM()] for extracting the fitted values from a CBFM fit, [gratia_effects()] for calculating smooth estimates and parametric effects suitable for plotting/visualization; [influence.CBFM()] for calculating some basic influence measures from a CBFM fit, [ordinate.CBFM()] for an *ad-hoc* approach to constructing spatio-temporal ordinations from a CBFM fit, [plot.CBFM()] for basic residual diagnostics from a CBFM fit, [predict.CBFM()] for constructing predictions from a CBFM fit, [residuals.CBFM()] for calculating residuals from a CBFM fit, [simulate.CBFM()] for simulating spatio-temporal multivariate abundance data from a CBFM fit, [summary.CBFM()] for summaries including standard errors and confidence intervals, and [varpart()] for variance partitioning of a CBFM fit.
#' 
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
#' rm(X, spp_loadings, true_lvs, xy, simy, dat)
#' 
#' 
#' # Fit stacked GLM as a baseline
#' fitstacked <- manyglm(simy_train ~ temp + depth + chla + O2, family = binomial(), data = dat_train)
#' 
#' 
#' # Set up spatial basis functions for CBFM -- Most users will start here! 
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
#' 
#' # Fit CBFM 
#' tic <- proc.time()
#' useformula <- ~ temp + depth + chla + O2
#' fitcbfm <- CBFM(y = simy_train, formula = useformula, data = dat_train,,
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
#' # Set up spatial basis functions for CBFM -- Most users will start here! 
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
#' useformula <- ~ temp + s(depth) + chla + s(O2)
#' fitcbfm_gam <- CBFM(y = simy_train, formula = useformula, 
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
#' # Fit Poisson CBFM 
#' tic <- proc.time()
#' useformula <- ~ temp + depth + chla + O2
#' fitcbfm <- CBFM(y = simy_train, formula = useformula, data = dat_train, 
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
#' # Fit negative binomial CBFM
#' # Please make sure you set up spatial basis functions for CBFM, we do so in Example 1a
#' # Please note the negative binomial distribution is not necessarily required for 
#' # overdispersed count data, since the latent variables can to some degree account for
#' # overdispersion. 
#' tic <- proc.time()
#' useformula <- ~ temp + depth + chla + O2
#' fitcbfm <- CBFM(y = simy_train, formula = useformula, data = dat_train, 
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
#' library(gamlss)
#' library(pscl)
#' 
#' # Probability of zero-inflation 
#' spp_zislopes <- matrix(runif(num_spp * num_X, -1, 1), nrow = num_spp)
#' spp_ziintercepts <- runif(num_spp, -0.5, 0)
#' 
#' # Simulate spatial multivariate abundance data
#' # Note the deliberate "+2" on the linear predictor: This creates data that is a bit more 
#' # clearly overdispersed and zero-inflated...
#' zieta <- tcrossprod(cbind(1,mm), cbind(spp_ziintercepts,spp_zislopes))
#' component_ind <- matrix(rbinom(num_sites * num_spp, size = 1, 
#' prob = matrix(plogis(zieta), num_sites, num_spp, byrow = TRUE)), num_sites,num_spp)
#' simy <- matrix(rpois(num_sites * num_spp, lambda = exp(eta+2) * (1-component_ind)), 
#' num_sites, num_spp)
#' rm(component_ind, zieta)
#' 
#' # Form training and test sets
#' simy_train <- simy[1:500,]
#' simy_test <- simy[501:1000,]
#' 
#' 
#' # Fit stacked zero-inflated Poisson regression models as a baseline
#' fitstacked <- NULL 
#' for(j in 1:num_spp) {
#' fitstacked[[j]] <- zeroinfl(resp ~ temp + depth + chla + O2, 
#' data = data.frame(resp = simy_train[,j], dat_train))
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
#' 
#' # Fit zero-inflated Poisson CBFM
#' tic <- proc.time()
#' useformula <- ~ temp + depth + chla + O2
#' fitcbfm <- CBFM(y = simy_train, formula = useformula, ziformula = useformula, data = dat_train, 
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
#' predictions_cbfm <- exp(predict(fitcbfm, newdata = dat_test, type = "response", 
#' new_B_space = test_basisfunctions))
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
#' countpredictions_stacked <- sapply(1:num_spp, function(j) predict(fitstacked[[j]], 
#' newdata = dat_test, type = "count"))
#' zipredictions_stacked <- sapply(1:num_spp, function(j) predict(fitstacked[[j]], 
#' newdata = dat_test, type = "zero"))
#' countpredictions_cbfm <- exp(predict(fitcbfm, newdata = dat_test, type = "link", 
#' new_B_space = test_basisfunctions))
#' zipredictions_cbfm <- plogis(tcrossprod(
#' predict(fitcbfm, newdata = dat_test, type = "zilpmatrix"), fitcbfm$zibetas))
#' 
#' preddeviance <- data.frame(
#' stacked = sapply(1:num_spp, function(j) { 
#' -2*sum(dZIP(simy_test[,j], mu = countpredictions_stacked[,j], 
#' sigma = zipredictions_stacked[,j], log = TRUE))
#' }),
#' cbfm = sapply(1:num_spp, function(j) { 
#' -2*sum(dZIP(simy_test[,j], mu = countpredictions_cbfm[,j], 
#' sigma = zipredictions_cbfm[,j], log = TRUE))
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
#' rm(countpredictions_stacked, zipredictions_stacked, 
#' countpredictions_cbfm, zipredictions_cbfm)
#' 
#' 
#' ##------------------------------
#' ## **Example 1f: Repeat Example 1a but illustrate applications to ZINB count data**
#' ## **This time, we use constant species-specific probabilities of zero inflation**
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
#' # Fit stacked zero-inflated NB regression models as a baseline
#' fitstacked <- NULL 
#' for(j in 1:num_spp) {
#' fitstacked[[j]] <- zeroinfl(resp ~ temp + depth + chla + O2 | 1, 
#' dist = "negbin", data = data.frame(resp = simy_train[,j], dat_train))
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
#' # Fit zero-inflated negative binomial CBFM
#' tic <- proc.time()
#' useformula <- ~ temp + depth + chla + O2
#' fitcbfm <- CBFM(y = simy_train, formula = useformula, ziformula = ~ 1, data = dat_train, 
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
#' countpredictions_stacked <- sapply(1:num_spp, function(j) predict(fitstacked[[j]], 
#' newdata = dat_test, type = "count"))
#' zipredictions_stacked <- sapply(1:num_spp, function(j) predict(fitstacked[[j]], 
#' newdata = dat_test, type = "zero"))
#' countpredictions_cbfm <- exp(predict(fitcbfm, newdata = dat_test, type = "link", 
#' new_B_space = test_basisfunctions))
#' zipredictions_cbfm <- plogis(tcrossprod(
#' predict(fitcbfm, newdata = dat_test, type = "zilpmatrix"), fitcbfm$zibetas))
#' 
#' preddeviance <- data.frame(
#' stacked = sapply(1:num_spp, function(j) { 
#' -2*sum(dZINBI(simy_test[,j], mu = countpredictions_stacked[,j], 
#' nu = zipredictions_stacked[,j], sigma = 1/fitstacked[[j]]$theta), log = TRUE)
#' }),
#' cbfm = sapply(1:num_spp, function(j) { 
#' -2*sum(dZINBI(simy_test[,j], mu = countpredictions_cbfm[,j], 
#' nu = zipredictions_cbfm[,j], sigma = fitcbfm$dispparam[j]), log = TRUE)
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
#' # Fit Tweedie CBFM
#' tic <- proc.time()
#' useformula <- ~ temp + depth + chla + O2
#' fitcbfm <- CBFM(y = simy_train, formula = useformula, data = dat_train, 
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
#' num_X <- 4 # Number of regression slope
#' spp_slopes <- matrix(runif(num_spp * num_X, -1, 1), nrow = num_spp)
#' spp_intercepts <- runif(num_spp, -2, 0)
#' 
#' # Simulate spatial coordinates and environmental covariate components
#' # We will use this information in later examples as well
#' xy <- data.frame(x = runif(num_sites, 0, 5), y = runif(num_sites, 0, 5))
#' X <- mvtnorm::rmvnorm(num_sites, mean = rep(0,4))
#' colnames(X) <- c("temp", "depth", "chla", "O2")
#' dat <- data.frame(xy, X)
#' 
#' # Simulate latent variable component
#' # We will use this information in later examples as well
#' true_lvs <- grf(grid = cbind(xy$x, xy$y), nsim = 2, cov.model = "exponential", 
#' cov.pars = c(1, 2))$data %>% 
#'      as.matrix
#' spp_loadings <- matrix(runif(num_spp * 2, -1, 1), nrow = num_spp) 
#' set.seed(NULL)
#' 
#' # Simulate spatial multivariate abundance data (presence-absence)
#' # We will use this information in later examples as well
#' eta <- tcrossprod(cbind(1,mm), cbind(spp_intercepts,spp_slopes)) + 
#' tcrossprod(true_lvs, spp_loadings)
#' 
#' spp_dispersion <- runif(num_spp)
#' # Simulate spatial multivariate abundance data
#' simy <- matrix(rNBItr(num_sites * num_spp, mu = exp(eta), sigma = 
#' matrix(spp_dispersion, num_sites, num_spp, byrow = TRUE)), num_sites, num_spp)
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
#' # Set up spatial basis functions for CBFM -- Most users will start here! 
#' num_basisfunctions <- 20 # Number of spatial basis functions to use
#' basisfunctions <- mrts(all.data$xy, num_basisfunctions) %>% 
#' as.matrix %>%
#' {.[,-(1)]} # Remove the first intercept column
#' 
#' # Fit CBFM 
#' tic <- proc.time()
#' fitcbfm <- CBFM(y = all.data$Y, formula = all.data$X.formula, data = all.data$X.data, 
#' B_space = basisfunctions, family = gaussian(), control = list(trace = 1))
#' toc <- proc.time()
#' toc - tic
#' 
#' summary(fitcbfm) %>% 
#' str
#' 
#' 
#' # Evaluate in-sample performance (similar to what was done in the Hmsc vignette)
#' # with the evaluateModelFit() function 
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
#' # Set up spatial basis functions for CBFM -- Most users will start here! 
#' num_basisfunctions <- 20 # Number of spatial basis functions to use
#' basisfunctions <- mrts(all.data$xy, num_basisfunctions) %>% 
#' as.matrix %>%
#' {.[,-(1)]} # Remove the first intercept column
#' 
#' # Fit CBFM 
#' # Note also that Hmsc generates and fits models assuming a probit link, 
#' # but CBFM uses a logit link
#' tic <- proc.time()
#' fitcbfm <- CBFM(y = all.data$Y, formula = all.data$X.formula, data = all.data$X.data, 
#' B_space = basisfunctions, family = binomial(), control = list(trace = 1))
#' toc <- proc.time()
#' toc - tic
#' 
#' 
#' # Evaluate in-sample performance (similar to what was done in the Hmsc vignette)
#' # with the evaluateModelFit() function 
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
#' ## Simulated from a spatio-temporal latent variable model
#' ## Please note the data generation process (thus) differs from CBFM.
#' ## The additive CBFM might take a while to fit...grab a cuppa!
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
#' X <- mvtnorm::rmvnorm(num_sites, mean = rep(0,4))
#' colnames(X) <- c("temp", "depth", "chla", "O2")
#' dat <- data.frame(xy, time = sort(runif(1000, 0, 10)) , X)
#' mm <- model.matrix(~ temp + depth + chla + O2 - 1, data = dat) %>% 
#' scale %>% 
#' as.matrix
#' 
#' # Simulate latent variable component
#' # We will also use this information in examples below
#' true_space_lvs <- grf(grid = cbind(xy$x, xy$y), nsim = 2, cov.model = "exponential", 
#' cov.pars = c(1, 2))$data %>% 
#'      as.matrix
#' true_time_lvs <- grf(grid = cbind(dat$time, 0), nsim = 2, cov.model = "gaussian", 
#' cov.pars = c(1, 1))$data %>% 
#'      as.matrix
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
#' rm(X, eta, xy, simy, true_time_lvs, spp_time_loadings)
#' 
#' 
#' # Fit stacked GLM as a baseline
#' fitstacked <- manyglm(simy_train ~ temp + depth + chla + O2, family = binomial(), data = dat_train)
#' 
#' 
#' # Set up spatial basis functions for CBFM -- Most users will start here! 
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
#' time_basisfunctions <- eval_basis(time_basisfunctions, s = as.matrix(dat$time)) %>%
#' as.matrix
#' train_time_basisfunctions <- time_basisfunctions[1:500,] 
#' test_time_basisfunctions <- time_basisfunctions[501:1000,] 
#' rm(time_basisfunctions, time_knots)
#' 
#' # Fit CBFM with additive spatial and temporal basis functions 
#' tic <- proc.time()
#' useformula <- ~ temp + depth + chla + O2
#' fitcbfm_additive <- CBFM(y = simy_train, formula = useformula, data = dat_train, 
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
#' ## **Example 3b: Repeat Example 3a but use a species-specific random intercept for time**
#' ## Simulated from a spatio-temporal latent variable model with random intercept for time
#' ## Please note the data generation process (thus) differs from CBFM.
#' ## The additive CBFM might take a while to fit...grab a cuppa!
#' ##------------------------------
#' # Consider 10 time points, 100 spatial locations sampled per time point
#' dat$time <- rep(1:10, each = 100) %>% factor
#' spp_time_variances <- runif(num_spp, 0.5, 2)
#' spp_time_randint <- sapply(spp_time_variances, function(x) rnorm(10, 0, sd = sqrt(x))) %>%
#' t
#' mm_time <- model.matrix(~ time - 1, data = dat)
#' 
#' # Simulate spatial multivariate abundance data (presence-absence)
#' # We will also use this information in examples below
#' eta <- tcrossprod(cbind(1,mm), cbind(spp_intercepts,spp_slopes)) + 
#' tcrossprod(true_space_lvs, spp_space_loadings) +
#' tcrossprod(mm_time, spp_time_randint)
#' simy <- matrix(rbinom(num_sites * num_spp, size = 1, 
#' prob = plogis(eta)), nrow = num_sites)
#' 
#' # Form training and test sets
#' simy_train <- simy[1:500,]
#' simy_test <- simy[501:1000,]
#' dat_train <- dat[1:500,]
#' dat_test <- dat[501:1000,]
#' rm(eta, mm, mm_time, simy, true_space_lvs, spp_space_loadings)
#' 
#' 
#' # Fit stacked GLM as a baseline
#' fitstacked <- manyglm(simy_train ~ temp + depth + chla + O2, family = binomial(), data = dat_train)
#' 
#' 
#' # Set up spatial basis functions for CBFM -- Most users will start here! 
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
#' # Training temporal basis functions. For species-specific random intercepts, 
#' # set up the model matrix reflecting a random intercept for time. Then set up a custom Sigma,
#' # which in this case will be equal to an identity matrix 
#' # (since within each species, the random intercepts are drawn independently; see above)
#' dat_train$time <- factor(dat_train$time) # Re-factor so that it now only has five time points
#' train_time_basisfunctions <- model.matrix(~ time - 1, data = dat_train)
#' custom_Sigma_time <- diag(nrow = ncol(train_time_basisfunctions))
#' 
#' 
#' # Fit CBFM with additive spatial and temporal basis functions 
#' tic <- proc.time()
#' useformula <- ~ temp + depth + chla + O2
#' fitcbfm_additive <- CBFM(y = simy_train, formula = useformula, data = dat_train, 
#' B_space = train_space_basisfunctions, B_time = train_time_basisfunctions, family = binomial(), 
#' G_control = list(rank = c(5,5)), 
#' Sigma_control = list(rank = c(5,"full"), custom_time = custom_Sigma_time), 
#' control = list(trace = 1))
#' toc <- proc.time()
#' toc - tic
#' 
#' 
#' # Calculate predictions onto test dataset
#' predictions_stacked <- predict(fitstacked, newdata = dat_test, type = "response")
#' # Note the test data contains completely different time points to the training data, and
#' # so the test set of temporal basis functions is then just a matrix of zeros. This is 
#' # analogous to, say, how mgcv handles random effects; see ?smooth.construct.re.smooth.spec
#' test_time_basisfunctions <- matrix(0, nrow = 500, ncol = ncol(train_time_basisfunctions))
#' predictions_cbfm_additive <- predict(fitcbfm_additive, newdata = dat_test, type = "response", 
#' new_B_space = test_space_basisfunctions, new_B_time = test_time_basisfunctions)
#' 
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
#' ## **Example 3c: Repeat Example 3a but with tensor product basis functions**
#' ## Please note the data generation process (thus) differs from CBFM.
#' ## To save some time, and for illustrative purposes, we will use the fast method estimating
#' ## the covariance matrices
#' # Nevertheless, please note this might take quite a while...grab a big cuppa!
#' ##------------------------------
#' tic <- proc.time()
#' useformula <- ~ temp + depth + chla + O2
#' train_st_basisfunctions <- tensorproduct(train_space_basisfunctions, train_time_basisfunctions)
#' dim(train_st_basisfunctions)
#' 
#' fitcbfm_tensor <- CBFM(y = simy_train, formula = useformula, data = dat_train, 
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
#' @importFrom foreach foreach %dopar% 
#' @import Matrix 
#' @importFrom compiler cmpfun
#' @importFrom doParallel registerDoParallel
#' @importFrom MASS theta.mm
#' @importFrom methods as
#' @importFrom mgcv betar gam gam.vcomp k.check ldTweedie logLik.gam model.matrix.gam notExp notLog pen.edf predict.gam nb rTweedie Tweedie tw ziplss
#' @importFrom numDeriv grad
#' @importFrom parallel detectCores
#' @importFrom stats fitted dnorm pnorm qnorm rnorm dbinom pbinom rbinom dnbinom pnbinom rnbinom dbeta pbeta rbeta dexp pexp rexp dgamma pgamma rgamma dlogis plogis qlogis dpois ppois rpois runif dchisq pchisq qchisq qqnorm as.formula binomial formula Gamma logLik model.matrix na.omit optim nlminb residuals 
#' @importFrom TMB MakeADFun
#' @importFrom utils capture.output
#' @md

CBFM <- function(y, formula, ziformula = NULL, data, B_space = NULL, B_time = NULL, B_spacetime = NULL, 
     offset = NULL, ncores = NULL, family = stats::gaussian(), trial_size = 1, dofit = TRUE, stderrors = TRUE, 
     select = FALSE, gamma = 1, knots = NULL, 
     ziselect = FALSE, zigamma = 1, ziknots = NULL,
     nonzeromean_B_space = FALSE, nonzeromean_B_time = FALSE, nonzeromean_B_spacetime = FALSE,
     start_params = list(betas = NULL, zibetas = NULL, basis_effects_mat = NULL, dispparam = NULL, powerparam = NULL, 
                         custom_space_lambdas = NULL, custom_time_lambdas = NULL, custom_spacetime_lambdas = NULL),
     TMB_directories = list(cpp = system.file("executables", package = "CBFM"), compile = system.file("executables", package = "CBFM")),
     control = list(maxit = 100, inner_maxit = 1, optim_lower = -50, optim_upper = 50, convergence_type = "parameters_MSE", tol = 1e-4, final_maxit = 100,   
                    initial_beta_dampen = 1, subsequent_betas_dampen = 0.25, 
                    gam_method = "REML", seed = NULL, ridge = 0, ziridge = 0, trace = 0),
     Sigma_control = list(rank = 5, maxit = 100, tol = 1e-4, method = "REML", trace = 0, 
                          custom_space = NULL, custom_time = NULL, custom_spactime = NULL), 
     G_control = list(rank = 5, structure = "unstructured", nugget_profile = seq(0.05, 0.95, by = 0.05), maxit = 100, 
                      tol = 1e-4, method = "REML", trace = 0, 
                      custom_space = NULL, custom_time = NULL, custom_spactime = NULL),
     k_check_control = list(subsample = 5000, n.rep = 400)) { 
     
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
          colnames(y) <- paste0("species", 1:ncol(y))
     if(is.null(rownames(y)))
          rownames(y) <- paste0("units", 1:nrow(y))
          
     full_gamma <- gamma
     if(length(full_gamma) == 1)
          full_gamma <- rep(full_gamma, ncol(y))
     if(!(length(full_gamma) %in% c(1, ncol(y))))
            stop("gamma should either be a scalar or a vector equal to the number of species i.e., ncol(y).")
     
     full_zigamma <- zigamma
     if(length(full_zigamma) == 1)
          full_zigamma <- rep(full_zigamma, ncol(y))
     if(!(length(full_zigamma) %in% c(1, ncol(y))))
          stop("zigamma should either be a scalar or a vector equal to the number of species i.e., ncol(y).")
     
     ## Form full basis function matrix B
     #.check_B_forms(B_space = B_space, B_time = B_time, B_spacetime = B_spacetime)
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
     underlying_B_colnames <- colnames(B)
     colnames(B) <- paste0("B_", 1:ncol(B))
     
     
     control <- .fill_control(control = control, num_spp = ncol(y), which_B_used = which_B_used, nonzeromean_B_space = nonzeromean_B_space, nonzeromean_B_time = nonzeromean_B_time, nonzeromean_B_spacetime = nonzeromean_B_spacetime)
     G_control <- .fill_G_control(control = G_control, which_B_used = which_B_used, num_spp = ncol(y), Sigma_control = Sigma_control)
     Sigma_control <- .fill_Sigma_control(control = Sigma_control, which_B_used = which_B_used, 
                                          num_spacebasisfns = num_spacebasisfns, 
                                          num_timebasisfns = num_timebasisfns, 
                                          num_spacetimebasisfns = num_spacetimebasisfns, 
                                          G_control = G_control) # Sigma_control must be done after G_control, as the latter needs the former supplied to check G_control$structure
     
     which_nonzeromean_B <- 1*c(nonzeromean_B_space, nonzeromean_B_time, nonzeromean_B_spacetime)

     ## Form covariate model matrix
     formula <- .check_X_formula(formula = formula, data = as.data.frame(data))          
     tmp_formula <- as.formula(paste("response", paste(as.character(formula),collapse = " ") ) )
     nullfit <- gam(tmp_formula, data = data.frame(data, response = runif(nrow(y))), knots = knots, fit = TRUE, control = list(maxit = 1))
     X <- model.matrix(nullfit)
     rm(tmp_formula, nullfit)
     rownames(X) <- rownames(y)

     .check_BX(B = B, X = X) 
     .check_offset(offset = offset, y = y) 
     .check_nonzeromeans(nonzeromean_B_space = nonzeromean_B_space, nonzeromean_B_time = nonzeromean_B_time, nonzeromean_B_spacetime = nonzeromean_B_spacetime) 
     .check_customSigma_Gstructure(Sigma_control = Sigma_control, G_control = G_control, which_B_used = which_B_used)
     
     
     ## Some other checks
     if(Sigma_control$method == "simple" || G_control$method == "simple")
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
     
     
     ## If required, form covariate model matrix for zero-inflation component
     if(is.null(ziformula)) {
          ziX <- NULL          
          zioffset <- NULL
          }
     if(!is.null(ziformula)) {
          ziformula <- .check_X_formula(formula = ziformula, data = as.data.frame(data))          
          tmp_formula <- as.formula(paste("response", paste(as.character(ziformula),collapse = " ") ) )
          nullfit <- gam(tmp_formula, data = data.frame(data, response = runif(nrow(y))), knots = ziknots, fit = TRUE, control = list(maxit = 1))

          ziX <- model.matrix(nullfit)
          zioffset <- model.offset(model.frame(nullfit))
          rm(tmp_formula, nullfit)
          rownames(ziX) <- rownames(y)
          }
     
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
     TMB::compile(paste0(getDLL, ".cpp")) #, flags = "-Wno-ignored-attributes -O2 -mfpmath=sse -msse2 -mstackrealign" 
     dyn.load(paste0(TMB_directories$compile, "/", TMB::dynlib(getDLL)))
     setwd(origwd)
     
     
     ##----------------
     ## Obtain starting values -- No selection is attempted here  
     ##----------------
     .check_start_params(start_params = start_params, num_spp = num_spp, num_basisfns = num_basisfns, num_X = num_X)
     initfit_fn <- function(j, formula) {
          tmp_formula <- as.formula(paste("response", paste(as.character(formula),collapse = " ") ) )

                    
          if(family$family %in% c("binomial","gaussian","poisson","Gamma","negative.binomial","Beta","tweedie")) {
               cw_family <- family
               if(family$family %in% c("negative.binomial"))
                    cw_family <- nb()
               if(family$family %in% c("Beta"))
                    cw_family <- betar(link = "logit")
               if(family$family %in% c("tweedie"))
                    cw_family <- Tweedie(p = 1.6, link = "log")
               
               use_size <- 0
               if(family$family %in% c("binomial")) {
                    tmp_formula <- as.formula(paste("cbind(response, size - response)", paste(as.character(formula),collapse = " ") ) )
                    use_size <- .ifelse_size(trial_size = trial_size, trial_size_length = trial_size_length, j = j, num_units = num_units)
                    }
                    
               fit0 <- try(gam(tmp_formula, 
                               data = as.data.frame(cbind(response = y[,j], data, size = use_size)), 
                               offset = offset[,j], 
                               knots = knots, 
                               method = control$gam_method, 
                               family = cw_family, 
                               gamma = full_gamma[j]), 
                           silent = TRUE)
               fit0$logLik <-  try(logLik(fit0), silent = TRUE)
               rm(cw_family)
               }
          if(family$family == "zipoisson") {
               tmp_ziformula <- as.formula(paste("taus", paste(as.character(ziformula), collapse = " ") ) )
               cw_offset <- offset[,j]
               if(is.null(cw_offset))
                    cw_offset <- numeric(num_units)
               
               init_pi <- rep(mean(y[,j] == 0, na.rm = TRUE), num_units) # Initial weights/posterior probabilities of being in zero-inflation component
               init_lambda <- mean(y[,j], na.rm = TRUE)
               w <- ifelse(y[,j] == 0, init_pi / (init_pi + (1-init_pi) * dpois(0, init_lambda)), 0)
               fit0 <- gam(tmp_formula, data = data.frame(response = y[,j], data), weights = 1-w, offset = cw_offset, knots = knots, method = control$gam_method, family = "poisson", gamma = full_gamma[j])
               fitzi <- suppressWarnings(gam(tmp_ziformula, data = data.frame(taus = w, data), knots = ziknots, method = control$gam_method, family = "binomial", gamma = full_zigamma[j]))
               MM <- model.matrix(fit0)
               cw_eta <- MM %*% fit0$coefficients
               if(!is.null(model.offset(model.frame(fit0))))
                    cw_eta <- cw_eta + model.offset(model.frame(fit0))
               cw_inner_logL <- .dzipoisson_log(y = na.omit(y[,j]), eta = cw_eta, zeroinfl_prob = fitted(fitzi)) 
               cw_inner_logL <- sum(cw_inner_logL[is.finite(cw_inner_logL)])

               inner_err <- Inf
               inner_inner_counter <- 0
               while(inner_err > 1e-4) {
                    if(inner_inner_counter > 20)
                         break;
                    
                    cw_eta <- MM %*% fit0$coefficients
                    if(!is.null(model.offset(model.frame(fit0))))
                         cw_eta <- cw_eta + model.offset(model.frame(fit0))
                    w <- ifelse(y[,j] == 0, fitted(fitzi) / (fitted(fitzi) + (1-fitted(fitzi)) * dpois(0, lambda = exp(cw_eta))), 0) # Posterior probabilities of being in zero-inflation component
                    fit0 <- gam(tmp_formula, data = data.frame(response = y[,j], data), weights = 1-w, offset = cw_offset, knots = knots, method = control$gam_method, family = "poisson", gamma = full_gamma[j])
                    fitzi <- suppressWarnings(gam(tmp_ziformula, data = data.frame(taus = w, data), knots = ziknots, method = control$gam_method, family = "binomial", gamma = full_zigamma[j]))
                    cw_eta <- MM %*% fit0$coefficients
                    if(!is.null(model.offset(model.frame(fit0))))
                         cw_eta <- cw_eta + model.offset(model.frame(fit0))
                    new_inner_logL <- .dzipoisson_log(y = na.omit(y[,j]), eta = cw_eta, zeroinfl_prob = fitted(fitzi))
                    new_inner_logL <- sum(new_inner_logL[is.finite(new_inner_logL)])
                    inner_err <- abs(new_inner_logL/cw_inner_logL-1)
                    cw_inner_logL <- new_inner_logL
                    inner_inner_counter <- inner_inner_counter + 1
                    #print(new_inner_logL)
                    }
               
               fit0$logLik <- new_inner_logL
               fit0$zicoefficients <- fitzi$coefficients
               }
          if(family$family[1] == "zinegative.binomial") {
               tmp_ziformula <- as.formula(paste("taus", paste(as.character(ziformula), collapse = " ") ) )
               cw_offset <- offset[,j]
               if(is.null(cw_offset))
                    cw_offset <- numeric(num_units)
               
               init_pi <- rep(mean(y[,j] == 0, na.rm = TRUE), num_units) # Initial weights/posterior probabilities of being in zero-inflation component
               init_lambda <- mean(y[,j], na.rm = TRUE)
               w <- ifelse(y[,j] == 0, init_pi / (init_pi + (1-init_pi) * dnbinom(0, mu = init_lambda, size = 1/0.2)), 0)
               fit0 <- gam(tmp_formula, data = data.frame(response = y[,j], data), weights = 1-w, offset = cw_offset, knots = knots, method = control$gam_method, family = nb(), gamma = full_gamma[j])
               fitzi <- suppressWarnings(gam(tmp_ziformula, data = data.frame(taus = w, data), knots = ziknots, method = control$gam_method, family = "binomial", gamma = full_zigamma[j]))
               MM <- model.matrix(fit0)
               cw_eta <- MM %*% fit0$coefficients
               if(!is.null(model.offset(model.frame(fit0))))
                    cw_eta <- cw_eta + model.offset(model.frame(fit0))
               cw_inner_logL <- .dzinegativebinomial_log(y = na.omit(y[,j]), eta = cw_eta, zeroinfl_prob = fitted(fitzi), phi = 1/fit0$family$getTheta(TRUE))
               cw_inner_logL <- sum(cw_inner_logL[is.finite(cw_inner_logL)])
               
               inner_err <- Inf
               inner_inner_counter <- 0
               while(inner_err > 1e-4) {
                    if(inner_inner_counter > 20)
                         break;

                    cw_eta <- MM %*% fit0$coefficients
                    if(!is.null(model.offset(model.frame(fit0))))
                         cw_eta <- cw_eta + model.offset(model.frame(fit0))
                    w <- ifelse(y[,j] == 0, fitted(fitzi) / (fitted(fitzi) + (1-fitted(fitzi)) * dnbinom(0, mu = exp(cw_eta), size = fit0$family$getTheta(TRUE))), 0) # Posterior probabilities of being in zero-inflation component
                    fit0 <- gam(tmp_formula, data = data.frame(response = y[,j], data), weights = 1-w, offset = cw_offset, knots = knots, method = control$gam_method, family = nb(), gamma = full_gamma[j])
                    fitzi <- suppressWarnings(gam(tmp_ziformula, data = data.frame(taus = w, data), knots = ziknots, method = control$gam_method, family = "binomial", gamma = full_zigamma[j]))
                    cw_eta <- MM %*% fit0$coefficients
                    if(!is.null(model.offset(model.frame(fit0))))
                         cw_eta <- cw_eta + model.offset(model.frame(fit0))
                    new_inner_logL <- .dzinegativebinomial_log(y = na.omit(y[,j]), eta = cw_eta, zeroinfl_prob = fitted(fitzi), phi = 1/fit0$family$getTheta(TRUE))
                    new_inner_logL <- sum(new_inner_logL[is.finite(new_inner_logL)])
                    inner_err <- abs(new_inner_logL/cw_inner_logL-1)
                    cw_inner_logL <- new_inner_logL
                    inner_inner_counter <- inner_inner_counter + 1
                    #print(new_inner_logL)
                    }
               
               fit0$logLik <- new_inner_logL
               fit0$zicoefficients <- fitzi$coefficients
               }
          if(family$family %in% c("ztpoisson")) {
               tmp_formula <- as.formula(paste("response", paste(as.character(formula), collapse = " "), "+ offset(off)" ) )
               cw_offset <- offset[,j]
               if(is.null(cw_offset))
                    cw_offset <- numeric(num_units)
               tmp_dat <- data.frame(response = c(y[,j], numeric(10)), data[c(1:nrow(data),1:10),], off = c(cw_offset, numeric(10))) # Append some zeros speeds ziplss up a heck of a lot!
               
               fit0 <- gam(list(tmp_formula, ~1), data = tmp_dat, method = control$gam_method, family = ziplss(), knots = knots, gamma = full_gamma[j])
               if(!inherits(fit0, "try-error")) {
                    MM <- model.matrix(fit0)
                    MM <- MM[1:(nrow(MM)-10),,drop=FALSE]
                    MM <- MM[, -ncol(MM), drop=FALSE]
                    fit0$coefficients <- fit0$coefficients[1:ncol(MM)] 
                    cw_eta <- MM %*% fit0$coefficients
                    if(!is.null(model.offset(model.frame(fit0))))
                         cw_eta <- cw_eta + model.offset(model.frame(fit0))
                    fit0$logLik <- .dztpois(y[,j], lambda = exp(cw_eta), log = TRUE) # .dztpois y = 0 values to -Inf, and handles NA values; because offset is in formula, then offset is contained in linear.predictors
                    fit0$logLik <- sum(fit0$logLik[is.finite(fit0$logLik)])
                    rm(MM)
                    }
               }
          if(family$family %in% c("ztnegative.binomial")) {
               cw_offset <- offset[,j]
               if(is.null(cw_offset))
                    cw_offset <- numeric(num_units)

               find_nonzeros <- which(y[,j] > 0) # This automatically excludes NA values
               init_lambda <- mean(y[,j], na.rm = TRUE)
               initw <- dnbinom(0, mu = init_lambda, size = 1/0.2) / (1-dnbinom(0, mu = init_lambda, size = 1/0.2))
               w1 <- c(rep(1,num_units))
               w1[which(y[,j]==0)] <- 0
               w2 <- rep(initw, num_units)
               w2[which(y[,j]==0)] <- 0
               w <- c(w1, w2)
               rm(w2) ## This way of constructing the weights and GAM set up is the *only* one I have constructed so far that ensures MM = X
               fit0 <- gam(tmp_formula, data = data.frame(response = c(y[,j], numeric(num_units)), data[c(1:num_units,1:num_units),]), 
                               weights = w, offset = cw_offset[c(1:num_units, 1:num_units)], knots = knots, method = control$gam_method, family = nb(), gamma = full_gamma[j])
               MM <- predict.gam(fit0, newdata = data, type = "lpmatrix")
               cw_eta <- MM[find_nonzeros,,drop=FALSE] %*% fit0$coefficients
               if(!is.null(model.offset(model.frame(fit0))))
                    cw_eta <- cw_eta + model.offset(model.frame(fit0))[find_nonzeros]
               cw_inner_logL <- .dztnbinom(y[find_nonzeros,j], mu = exp(cw_eta), size = fit0$family$getTheta(TRUE), log = TRUE)
               cw_inner_logL <- sum(cw_inner_logL[is.finite(cw_inner_logL)])
               
               inner_err <- Inf
               inner_inner_counter <- 0
               while(inner_err > 1e-4) { 
                    if(inner_inner_counter > 20)
                         break;
                    
                    cw_eta <- MM[find_nonzeros,,drop=FALSE] %*% fit0$coefficients
                    if(!is.null(model.offset(model.frame(fit0))))
                         cw_eta <- cw_eta + model.offset(model.frame(fit0))[find_nonzeros]
                    initw <- dnbinom(0, mu = exp(cw_eta), size = fit0$family$getTheta(TRUE)) / (1-dnbinom(0, mu = exp(cw_eta), size = fit0$family$getTheta(TRUE)) + 1e-4)
                    initw[initw > 0.1/1e-4] <- 0.1/1e-4
                    w2 <- numeric(num_units)
                    w2[find_nonzeros] <- initw
                    w <- c(w1, w2)
                    rm(w2) 
                    fit0 <- try(gam(tmp_formula, data = data.frame(response = c(y[,j],numeric(num_units)), data[c(1:num_units,1:num_units),]), 
                                    weights = w, offset = cw_offset[c(1:num_units, 1:num_units)], knots = knots, method = control$gam_method, family = nb(), gamma = full_gamma[j]), silent = TRUE)
                    if(inherits(fit0, "try-error"))
                         break;
                    cw_eta <- MM[find_nonzeros,,drop=FALSE] %*% fit0$coefficients
                    if(!is.null(model.offset(model.frame(fit0))))
                         cw_eta <- cw_eta + model.offset(model.frame(fit0))[find_nonzeros]
                    new_inner_logL <- .dztnbinom(y[find_nonzeros,j], mu = exp(cw_eta), size = fit0$family$getTheta(TRUE), log = TRUE)
                    new_inner_logL <- sum(new_inner_logL[is.finite(new_inner_logL)])
                    inner_err <- abs(new_inner_logL/cw_inner_logL-1)
                    cw_inner_logL <- new_inner_logL
                    inner_inner_counter <- inner_inner_counter + 1
                    #print(new_inner_logL)
                    }
               
               fit0$logLik <- new_inner_logL
               }
                   
          if(inherits(fit0, "try-error")) {
               fit0 <- list(coefficients = runif(num_X), dispparam = 1)
               if(family$family[1] %in% c("zipoisson","zinegative.binomial"))
                    fit0$zicoefficients <- runif(ncol(ziX))
               }
          
          return(fit0)
          }

     if(is.null(start_params$betas)) { 
          if(control$trace > 0)
               message("Calculating starting values...")
          
          all_start_fits <- foreach(j = 1:num_spp) %dopar% initfit_fn(j = j, formula = formula)              
          start_params$betas <- do.call(rbind, lapply(all_start_fits, function(x) x$coefficients))
          start_params$betas <- start_params$betas * control$initial_betas_dampen # Should be OK even if control$initial_betas_dampen is vector equal to number of species
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
          if(family$family[1] %in% c("ztnegative.binomial", "zinegative.binomial"))
               start_params$dispparam <- rep(0.2, num_spp)
          if(family$family[1] %in% c("Beta")) 
               start_params$dispparam <- rep(0.5, num_spp)          
          }
     if(is.null(start_params$powerparam))
          start_params$powerparam <- rep(ifelse(family$family[1] == "tweedie", 1.6, 0), num_spp)
     if(is.null(start_params$zibetas)) {
          if(family$family[1] %in% c("zipoisson","zinegative.binomial")) {
               start_params$zibetas <- do.call(rbind, lapply(all_start_fits, function(x) x$zicoefficients))
               start_params$zibetas <- start_params$zibetas * control$initial_betas_dampen 
               }
          if(!(family$family[1] %in% c("zipoisson","zinegative.binomial")))
               start_params$zibetas <- NULL
          }
     
     
     if(which_B_used[1]) {
          if(is.null(Sigma_control[["custom_space"]])) {
               start_params$Sigma_space <- diag(1, nrow = num_spacebasisfns)
               new_LoadingnuggetSigma_space <- list(invcov = start_params$Sigma_space)
               }
          if(!is.null(Sigma_control[["custom_space"]])) { 
               if(is.matrix(Sigma_control[["custom_space"]])) {
                    start_params$Sigma_space <- NULL
                    new_LoadingnuggetSigma_space <- list(invcov = .pinv(V = Sigma_control[["custom_space"]]), cov = Sigma_control[["custom_space"]])
                    }
               if(is.list(Sigma_control[["custom_space"]])) {
                    start_params$Sigma_space <- NULL
                    if(!is.null(start_params$custom_space_lambdas))
                         new_LoadingnuggetSigma_space <- list(lambdas = start_params$custom_space_lambdas) 
                    if(is.null(start_params$custom_space_lambdas))
                         new_LoadingnuggetSigma_space <- list(lambdas = rep(10, length(Sigma_control[["custom_space"]]))) #' This is an arbitrary starting value!
                    new_LoadingnuggetSigma_space$invcov <- Reduce("+", lapply(1:length(Sigma_control[["custom_space"]]), function(x) .pinv(Sigma_control[["custom_space"]][[x]]) / new_LoadingnuggetSigma_space$lambdas[x])) 
                    # Note this object will not contain the individual constituent Sigma_space matrices. 
                    new_LoadingnuggetSigma_space$cov <- .pinv(new_LoadingnuggetSigma_space$invcov)
                    }
               }
          if(is.null(G_control[["custom_space"]])) {
               start_params$G_space <- diag(1, nrow = num_spp)
               new_LoadingnuggetG_space <- list(invcov = start_params$G_space)
               }
          if(!is.null(G_control[["custom_space"]])) {
               start_params$G_space <- NULL
               new_LoadingnuggetG_space <- list(invcov = .pinv(V = G_control[["custom_space"]]), cov = G_control[["custom_space"]])
               }
          start_params[["mean_B_space"]] <- numeric(num_spacebasisfns)
          }
     if(which_B_used[2]) {
          if(is.null(Sigma_control[["custom_time"]])) {
               start_params$Sigma_time <- diag(1, nrow = num_timebasisfns)
               new_LoadingnuggetSigma_time <- list(invcov = start_params$Sigma_time)
               }
          if(!is.null(Sigma_control[["custom_time"]])) { 
               if(is.matrix(Sigma_control[["custom_time"]])) {
                    start_params$Sigma_time <- NULL
                    new_LoadingnuggetSigma_time <- list(invcov = .pinv(V = Sigma_control[["custom_time"]]), cov = Sigma_control[["custom_time"]])
                    }
               if(is.list(Sigma_control[["custom_time"]])) {
                    start_params$Sigma_time <- NULL
                    if(!is.null(start_params$custom_time_lambdas))
                         new_LoadingnuggetSigma_time <- list(lambdas = start_params$custom_time_lambdas) 
                    if(is.null(start_params$custom_time_lambdas))
                         new_LoadingnuggetSigma_time <- list(lambdas = rep(10, length(Sigma_control[["custom_time"]]))) #' This is an arbitrary starting value!
                    new_LoadingnuggetSigma_time$invcov <- Reduce("+", lapply(1:length(Sigma_control[["custom_time"]]), function(x) .pinv(Sigma_control[["custom_time"]][[x]]) / new_LoadingnuggetSigma_time$lambdas[x])) 
                    # Note this object will not contain the individual constituent Sigma_time matrices. 
                    new_LoadingnuggetSigma_time$cov <- .pinv(new_LoadingnuggetSigma_time$invcov)
                    }
               }
          if(is.null(G_control[["custom_time"]])) {
               start_params$G_time <- diag(1, nrow = num_spp)
               new_LoadingnuggetG_time <- list(invcov = start_params$G_time)
               }
          if(!is.null(G_control[["custom_time"]])) {
               start_params$G_time <- NULL
               new_LoadingnuggetG_time <- list(invcov = .pinv(V = G_control[["custom_time"]]), cov = G_control[["custom_time"]])
               }
          start_params[["mean_B_time"]] <- numeric(num_timebasisfns)
          }
     if(which_B_used[3]) {
          if(is.null(Sigma_control[["custom_spacetime"]])) {
               start_params$Sigma_spacetime <- diag(1, nrow = num_spacetimebasisfns)
               new_LoadingnuggetSigma_spacetime <- list(invcov = start_params$Sigma_spacetime)
               }
          if(!is.null(Sigma_control[["custom_spacetime"]])) { 
               if(is.matrix(Sigma_control[["custom_spacetime"]])) {
                    start_params$Sigma_spacetime <- NULL
                    new_LoadingnuggetSigma_spacetime <- list(invcov = .pinv(Sigma_control[["custom_spacetime"]]), cov = Sigma_control[["custom_spacetime"]]) 
                    }
               if(is.list(Sigma_control[["custom_spacetime"]])) {
                    start_params$Sigma_spacetime <- NULL
                    if(!is.null(start_params$custom_spacetime_lambdas))
                         new_LoadingnuggetSigma_spacetime <- list(lambdas = start_params$custom_spacetime_lambdas) 
                    if(is.null(start_params$custom_spacetime_lambdas))
                         new_LoadingnuggetSigma_spacetime <- list(lambdas = rep(0.1, length(Sigma_control[["custom_spacetime"]]))) #' This is an arbitrary starting value!
                    new_LoadingnuggetSigma_spacetime$invcov <- Reduce("+", lapply(1:length(Sigma_control[["custom_spacetime"]]), function(x) .pinv(Sigma_control[["custom_spacetime"]][[x]]) / new_LoadingnuggetSigma_spacetime$lambdas[x])) 
                    # Note this object will not contain the individual constituent Sigma_spacetime matrices. 
                    new_LoadingnuggetSigma_spacetime$cov <- .pinv(new_LoadingnuggetSigma_spacetime$invcov)
                    }
               }
          if(is.null(G_control[["custom_spacetime"]])) {
               start_params$G_spacetime <- diag(1, nrow = num_spp)
               new_LoadingnuggetG_spacetime <- list(invcov = start_params$G_spacetime)
               }
          if(!is.null(G_control[["custom_spacetime"]])) {
               start_params$G_spacetime <- NULL
               new_LoadingnuggetG_spacetime <- list(invcov = .pinv(V = G_control[["custom_spacetime"]]), cov = G_control[["custom_spacetime"]])
               }
          start_params[["mean_B_spacetime"]] <- numeric(num_spacetimebasisfns)
          }
     start_params$logLik <- -Inf
     
     
     ##----------------
     ## Run PQL algorithm
     ##----------------
     tic <- proc.time()
     counter <- 0
     diff <- 10
     converged <- FALSE
     
     ## Only make offset here, after initial fitting done
     if(is.null(offset)) { 
          offset <- Matrix(0, nrow = nrow(y), ncol = ncol(y), sparse = TRUE)
          }
     formula_offset <- numeric(nrow(y))
     if(!is.null(model.offset(model.frame(all_start_fits[[1]]))))
          formula_offset <- model.offset(model.frame(all_start_fits[[1]]))
     rm(all_start_fits)
     
     if(control$trace > 0)
          message("Commencing model fitting...")
     
     while(diff > control$tol & counter <= control$maxit) {
          make_tidibits_data <- function() {
               out <- list(B = B)
               if(identical(which_B_used, c(1,0,0))) {
                    out$Sigmainv <- as.matrix(new_LoadingnuggetSigma_space$invcov)
                    out$Ginv <- as.matrix(new_LoadingnuggetG_space$invcov)
                    }
               if(identical(which_B_used, c(0,1,0))) {
                    out$Sigmainv <- as.matrix(new_LoadingnuggetSigma_time$invcov)
                    out$Ginv <- as.matrix(new_LoadingnuggetG_time$invcov)
                    }
               if(identical(which_B_used, c(0,0,1))) {
                    out$Sigmainv <- as.matrix(new_LoadingnuggetSigma_spacetime$invcov)
                    out$Ginv <- as.matrix(new_LoadingnuggetG_spacetime$invcov)
                    }
                    
               if(identical(which_B_used, c(1,1,0))) {
                    out$Sigmainv_B1 <- as.matrix(new_LoadingnuggetSigma_space$invcov)
                    out$Ginv_B1 <- as.matrix(new_LoadingnuggetG_space$invcov)
                    out$Sigmainv_B2 <- as.matrix(new_LoadingnuggetSigma_time$invcov)
                    out$Ginv_B2 <- as.matrix(new_LoadingnuggetG_time$invcov)
                    }
               if(identical(which_B_used, c(1,0,1))) {
                    out$Sigmainv_B1 <- as.matrix(new_LoadingnuggetSigma_space$invcov)
                    out$Ginv_B1 <- as.matrix(new_LoadingnuggetG_space$invcov)
                    out$Sigmainv_B2 <- as.matrix(new_LoadingnuggetSigma_spacetime$invcov)
                    out$Ginv_B2 <- as.matrix(new_LoadingnuggetG_spacetime$invcov)
                    }
               if(identical(which_B_used, c(0,1,1))) {
                    out$Sigmainv_B1 <- as.matrix(new_LoadingnuggetSigma_time$invcov)
                    out$Ginv_B1 <- as.matrix(new_LoadingnuggetG_time$invcov)
                    out$Sigmainv_B2 <- as.matrix(new_LoadingnuggetSigma_spacetime$invcov)
                    out$Ginv_B2 <- as.matrix(new_LoadingnuggetG_spacetime$invcov)
                    }
               if(identical(which_B_used, c(1,1,1))) {
                    out$Sigmainv_B1 <- as.matrix(new_LoadingnuggetSigma_space$invcov)
                    out$Ginv_B1 <- as.matrix(new_LoadingnuggetG_space$invcov)
                    out$Sigmainv_B2 <- as.matrix(new_LoadingnuggetSigma_time$invcov)
                    out$Ginv_B2 <- as.matrix(new_LoadingnuggetG_time$invcov)
                    out$Sigmainv_B3 <- as.matrix(new_LoadingnuggetSigma_spacetime$invcov)
                    out$Ginv_B3 <- as.matrix(new_LoadingnuggetG_spacetime$invcov)
                    }                    

               return(out)
               }

          if(counter == 0) {
               tidbits_data <- make_tidibits_data()
                                   
               new_fit_CBFM_ptest <- list(
                    betas = start_params$betas,
                    basis_effects_mat = start_params$basis_effects_mat,
                    zibetas = start_params$zibetas,
                    mean_B_space = start_params[["mean_B_space"]], # NULL if which_B_used[1] == 0
                    mean_B_time = start_params[["mean_B_time"]], # NULL if which_B_used[2] == 0
                    mean_B_spacetime = start_params[["mean_B_spacetime"]], # ]NULL if which_B_used[3] == 0
                    dispparam = start_params$dispparam,
                    powerparam = start_params$powerparam,
                    logLik = start_params$logLik
                    )
                    
               cw_params <- as.vector(unlist(new_fit_CBFM_ptest))
               cw_params <- cw_params[-length(cw_params)] # Get rid of the logLik term
               cw_inner_params <- cw_params
               cw_logLik <- new_fit_CBFM_ptest$logLik
               rm(start_params)
               }
          if(counter > 0) {          
               tidbits_data <- make_tidibits_data()
               cw_inner_params <- cw_params
               }

          cw_inner_logL <- -Inf
          inner_err <- 100
          inner_counter <- 0
          if(control$trace > 0)
               message("Updating all coefficients and dispersion/power parameters (this includes running an inner EM algorithm if appropriate).")         
          while(inner_err > 1e-3 & inner_counter <= control$inner_maxit) {

               ##-------------------------
               ## For zero-inflated distributions, E-step + updating zero-inflation probabilities for distributions that require it. 
               ## Otherwise, effectively do nothing
               ##-------------------------
               getweights <- .estep_fn(family = family,
                                       cwfit = new_fit_CBFM_ptest,
                                       y = y,
                                       X = X,
                                       offset = offset,
                                       formula_offset = formula_offset,
                                       B = B,
                                       ziX = ziX,
                                       zioffset = zioffset) # Posterior probabilities of zero-inflation

               ##-------------------------
               ## Update smoothing coefficients for all basis functions, one species at a time
               ##-------------------------
               update_basiscoefsspp_fn <- function(j) {
                    dyn.load(paste0(TMB_directories$compile, "/", TMB::dynlib(getDLL)))
                    
                    tidbits_parameters <- list(basis_effects = new_fit_CBFM_ptest$basis_effects_mat[j,])
                    
                    tidbits_data$y <- y[,j]
                    tidbits_data$Xbeta <- as.vector(X %*% new_fit_CBFM_ptest$betas[j,])
                    tidbits_data$dispparam <- new_fit_CBFM_ptest$dispparam[j]
                    tidbits_data$powerparam <- new_fit_CBFM_ptest$powerparam[j]
                    tidbits_data$estep_weights <- as.vector(getweights[,j])
                    tidbits_data$offset <- offset[,j] + formula_offset
                    if(sum(which_B_used) == 1) {
                         tidbits_data$mean_basis_effects <- c(new_fit_CBFM_ptest[["mean_B_space"]], new_fit_CBFM_ptest[["mean_B_time"]], new_fit_CBFM_ptest[["mean_B_spacetime"]])
                         tidbits_data$other_centered_basis_effects_mat <- new_fit_CBFM_ptest$basis_effects_mat - matrix(tidbits_data$mean_basis_effects, nrow = num_spp, ncol = num_basisfns, byrow = TRUE)
                         }
                    if(sum(which_B_used) == 2) {
                         if(identical(which_B_used, c(1,1,0))) { 
                              tidbits_data$mean_basis_effects_B1 <- new_fit_CBFM_ptest[["mean_B_space"]]
                              tidbits_data$mean_basis_effects_B2 <- new_fit_CBFM_ptest[["mean_B_time"]]
                              tidbits_data$other_centered_basis_effects_mat_B1 <- new_fit_CBFM_ptest$basis_effects_mat[,1:num_spacebasisfns,drop=FALSE] - matrix(tidbits_data$mean_basis_effects_B1, nrow = num_spp, ncol = num_spacebasisfns, byrow = TRUE)
                              tidbits_data$other_centered_basis_effects_mat_B2 <- new_fit_CBFM_ptest$basis_effects_mat[,num_spacebasisfns+(1:num_timebasisfns),drop=FALSE] - matrix(tidbits_data$mean_basis_effects_B2, nrow = num_spp, ncol = num_timebasisfns, byrow = TRUE)
                              }
                         if(identical(which_B_used, c(1,0,1))) { 
                              tidbits_data$mean_basis_effects_B1 <- new_fit_CBFM_ptest[["mean_B_space"]]
                              tidbits_data$mean_basis_effects_B2 <- new_fit_CBFM_ptest[["mean_B_spacetime"]]
                              tidbits_data$other_centered_basis_effects_mat_B1 <- new_fit_CBFM_ptest$basis_effects_mat[,1:num_spacebasisfns,drop=FALSE] - matrix(tidbits_data$mean_basis_effects_B1, nrow = num_spp, ncol = num_spacebasisfns, byrow = TRUE)
                              tidbits_data$other_centered_basis_effects_mat_B2 <- new_fit_CBFM_ptest$basis_effects_mat[,num_spacebasisfns+num_timebasisfns+(1:num_spacetimebasisfns),drop=FALSE] - matrix(tidbits_data$mean_basis_effects_B2, nrow = num_spp, ncol = num_spacetimebasisfns, byrow = TRUE)
                              }
                         if(identical(which_B_used, c(0,1,1))) { 
                              tidbits_data$mean_basis_effects_B1 <- new_fit_CBFM_ptest[["mean_B_time"]]
                              tidbits_data$mean_basis_effects_B2 <- new_fit_CBFM_ptest[["mean_B_spacetime"]]
                              tidbits_data$other_centered_basis_effects_mat_B1 <- new_fit_CBFM_ptest$basis_effects_mat[,num_spacebasisfns+(1:num_timebasisfns),drop=FALSE] - matrix(tidbits_data$mean_basis_effects_B1, nrow = num_spp, ncol = num_timebasisfns, byrow = TRUE)
                              tidbits_data$other_centered_basis_effects_mat_B2 <- new_fit_CBFM_ptest$basis_effects_mat[,num_spacebasisfns+num_timebasisfns+(1:num_spacetimebasisfns),drop=FALSE] - matrix(tidbits_data$mean_basis_effects_B2, nrow = num_spp, ncol = num_spacetimebasisfns, byrow = TRUE)
                              }
                         }
                    if(sum(which_B_used) == 3) {
                         tidbits_data$mean_basis_effects_B1 <- new_fit_CBFM_ptest[["mean_B_space"]]
                         tidbits_data$mean_basis_effects_B2 <- new_fit_CBFM_ptest[["mean_B_time"]]
                         tidbits_data$mean_basis_effects_B3 <- new_fit_CBFM_ptest[["mean_B_spacetime"]]
                         tidbits_data$other_centered_basis_effects_mat_B1 <- new_fit_CBFM_ptest$basis_effects_mat[,1:num_spacebasisfns,drop=FALSE] - matrix(tidbits_data$mean_basis_effects_B1, nrow = num_spp, ncol = num_spacebasisfns, byrow = TRUE)
                         tidbits_data$other_centered_basis_effects_mat_B2 <- new_fit_CBFM_ptest$basis_effects_mat[,num_spacebasisfns+(1:num_timebasisfns),drop=FALSE] - matrix(tidbits_data$mean_basis_effects_B1, nrow = num_spp, ncol = num_timebasisfns, byrow = TRUE)
                         tidbits_data$other_centered_basis_effects_mat_B3 <- new_fit_CBFM_ptest$basis_effects_mat[,num_spacebasisfns+num_timebasisfns+(1:num_spacetimebasisfns),drop=FALSE] - matrix(tidbits_data$mean_basis_effects_B2, nrow = num_spp, ncol = num_spacetimebasisfns, byrow = TRUE)
                         }
                    
                    tidbits_data$spp_ind <- j
                    tidbits_data$family <- .family_to_counter(family = family)
                    tidbits_data$trial_size <- .ifelse_size(trial_size = trial_size, trial_size_length = trial_size_length, j = j, num_units = num_units, family = family)
     
                    tidbits_constraints <- list(
                         lower = rep(control$optim_lower, sum(sapply(tidbits_parameters, length))),
                         upper = rep(control$optim_upper, sum(sapply(tidbits_parameters, length)))
                         )
                         
                    CBFM_objs <- TMB::MakeADFun(data = tidbits_data, parameters = tidbits_parameters, DLL = getDLL, hessian = FALSE, silent = TRUE)
                    
                    new_fit_CBFM <- try(nlminb(start = CBFM_objs$par, objective = CBFM_objs$fn, gradient = CBFM_objs$gr, 
                                               lower = tidbits_constraints$lower, upper = tidbits_constraints$upper, 
                                               control = list(iter.max = 500, eval.max = 1000)), silent = TRUE) #' Note choices in control here are pretty arbitrary!
                    
                    # Dampen Xbeta component and run it again...it is kind of ad-hoc but has been shown to be helpful especially with GAM fits to extremely overdispersed counts 
                    if(inherits(new_fit_CBFM, "try-error")) {
                         if(length(control$subsequent_betas_dampen) == 1)
                              tidbits_data$Xbeta <- tidbits_data$Xbeta * control$subsequent_betas_dampen 
                         if(length(control$subsequent_betas_dampen) == num_spp)
                              tidbits_data$Xbeta <- tidbits_data$Xbeta * matrix(control$subsequent_betas_dampen, nrow = num_units, ncol = num_spp, byrow = TRUE)
                         
                         CBFM_objs <- TMB::MakeADFun(data = tidbits_data, parameters = tidbits_parameters, DLL = getDLL, hessian = FALSE, silent = TRUE)
                         new_fit_CBFM <- try(nlminb(start = CBFM_objs$par, objective = CBFM_objs$fn, gradient = CBFM_objs$gr,
                              lower = tidbits_constraints$lower, upper = tidbits_constraints$upper, 
                              control = list(iter.max = 500, eval.max = 1000)), silent = TRUE)
                         }
                    if(inherits(new_fit_CBFM, "try-error")) {
                         new_fit_CBFM <- list(par = tidbits_parameters$basis_effects)
                    }
                    #new_fit_CBFM <- optim(par = CBFM_objs$par, fn = CBFM_objs$fn, gr = CBFM_objs$gr, method = "BFGS")
                    
                    return(new_fit_CBFM)
                    }
               update_basiscoefsspp_cmpfn <- compiler::cmpfun(update_basiscoefsspp_fn)
                    
               all_update_coefs <- foreach(j = 1:num_spp, .export = c("tidbits_data")) %dopar% update_basiscoefsspp_cmpfn(j = j)
               for(j in 1:num_spp) {
                    new_fit_CBFM_ptest$basis_effects_mat[j,] <- all_update_coefs[[j]]$par 
                    }
               rm(all_update_coefs, update_basiscoefsspp_fn)

     
               ##-------------------------
               ## Update coefficients related to covariate model matrix X, and other nuisance parameters, one response at a time
               #' new_offset contains the explicit offset argument, while formula_offset contains the offset from the formula, and this has to be included manually to the linear predictor 
               #' The point is that with both new_offset and formula_offset, there will *always* be an offset included, meaning fit0$offset will never be NULL
               ##-------------------------
               update_Xcoefsspp_fn <- function(j) {
                    tmp_formula <- as.formula(paste("response", paste(as.character(formula), collapse = " ") ) )
                    new_offset <- offset[,j] + as.vector(B %*% new_fit_CBFM_ptest$basis_effects_mat[j,])
                    Hmat <- diag(control$ridge+1e-15, nrow = num_X)
                    
                    if(family$family %in% c("binomial","gaussian","poisson","Gamma","negative.binomial","Beta","tweedie")) {
                         cw_family <- family
                         if(family$family == "negative.binomial")
                              cw_family <- nb()
                         if(family$family == "Beta")
                              cw_family <- betar(link = "logit")
                         if(family$family == "tweedie")
                              cw_family <- tw(link = "log")
                         
                         use_size <- 0
                         if(family$family %in% c("binomial")) {
                              tmp_formula <- as.formula(paste("cbind(response, size - response)", paste(as.character(formula), collapse = " ") ) )
                              use_size <- .ifelse_size(trial_size = trial_size, trial_size_length = trial_size_length, j = j, num_units = num_units)
                              }
                              
                         if(control$ridge > 0)
                              fit0 <- gam(tmp_formula, 
                                          data = as.data.frame(cbind(response = y[,j], data, size = use_size)), 
                                          offset = new_offset, 
                                          knots = knots, 
                                          method = control$gam_method, 
                                          H = Hmat, 
                                          family = cw_family, 
                                          select = select, 
                                          gamma = full_gamma[j])
                         if(control$ridge == 0)
                              fit0 <- gam(tmp_formula, 
                                          data = as.data.frame(cbind(response = y[,j], data, size = use_size)), 
                                          offset = new_offset, 
                                          knots = knots, 
                                          method = control$gam_method, 
                                          family = cw_family, 
                                          select = select, 
                                          gamma = full_gamma[j])
                         
                         fit0$logLik <- as.vector(logLik(fit0))
                         fit0$linear.predictors <- X %*% fit0$coefficients + new_offset + formula_offset
                         rm(cw_family)
                         }
     
                    if(family$family %in% c("zipoisson")) { # M-step
                         tmp_ziformula <- as.formula(paste("taus", paste(as.character(ziformula), collapse = " ") ) )
                         ziHmat <- diag(control$ziridge+1e-15, nrow = ncol(ziX))
                         
                         if(control$ridge > 0)
                              fit0 <- gam(tmp_formula, data = data.frame(response = y[,j], data), offset = new_offset, knots = knots, 
                                          method = control$gam_method, weights = 1-getweights[,j], H = Hmat, family = "poisson", 
                                          select = select, gamma = full_gamma[j])
                         if(control$ridge == 0)
                              fit0 <- gam(tmp_formula, data = data.frame(response = y[,j], data), offset = new_offset, knots = knots, 
                                          method = control$gam_method, weights = 1-getweights[,j], family = "poisson", 
                                          select = select, gamma = full_gamma[j])
                         
                         if(control$ziridge > 0)
                              fitzi <- suppressWarnings(gam(tmp_ziformula, data = data.frame(taus = getweights[,j], data), knots = ziknots, 
                                                            method = control$gam_method, H = ziHmat, family = "binomial", 
                                                            select = ziselect, gamma = full_zigamma[j]))
                         if(control$ziridge == 0)
                              fitzi <- suppressWarnings(gam(tmp_ziformula, data = data.frame(taus = getweights[,j], data), knots = ziknots, 
                                                            method = control$gam_method, family = "binomial", select = ziselect, gamma = full_zigamma[j]))

                         
                         fit0$logLik <- sum(.dzipoisson_log(y = na.omit(y[,j]), eta = model.matrix(fit0) %*% fit0$coefficients + fit0$offset, zeroinfl_prob = fitted(fitzi)))
                         fit0$linear.predictors <- X %*% fit0$coefficients + new_offset + formula_offset
                         }
                    if(family$family %in% c("zinegative.binomial")) { # M-step
                         tmp_ziformula <- as.formula(paste("taus", paste(as.character(ziformula), collapse = " ") ) )
                         ziHmat <- diag(control$ziridge+1e-15, nrow = ncol(ziX))
                         
                         if(control$ridge > 0)
                              fit0 <- gam(tmp_formula, data = data.frame(response = y[,j], data), offset = new_offset, knots = knots, 
                                          method = control$gam_method, weights = 1-getweights[,j], H = Hmat, family = nb(), 
                                          select = select, gamma = full_gamma[j])
                         if(control$ridge == 0)
                              fit0 <- gam(tmp_formula, data = data.frame(response = y[,j], data), offset = new_offset, knots = knots, 
                                          method = control$gam_method,weights = 1-getweights[,j], family = nb(), 
                                          select = select, gamma = full_gamma[j])
                         if(control$ziridge > 0)
                              fitzi <- suppressWarnings(gam(tmp_ziformula, data = data.frame(taus = getweights[,j], data), knots = ziknots, 
                                                            method = control$gam_method, H = ziHmat, family = "binomial", 
                                                            select = ziselect, gamma = full_zigamma[j]))
                         if(control$ziridge == 0)
                              fitzi <- suppressWarnings(gam(tmp_ziformula, data = data.frame(taus = getweights[,j], data), knots = ziknots, 
                                                            method = control$gam_method, family = "binomial", select = ziselect, 
                                                            gamma = full_zigamma[j]))

                         fit0$logLik <- sum(.dzinegativebinomial_log(y = na.omit(y[,j]), eta = model.matrix(fit0) %*% fit0$coefficients + fit0$offset, zeroinfl_prob = fitted(fitzi), phi = 1/fit0$family$getTheta(TRUE)))
                         fit0$linear.predictors <- X %*% fit0$coefficients + new_offset + formula_offset
                         }
                    if(family$family %in% c("ztpoisson")) {
                         Hmat <- diag(control$ridge+1e-15, nrow = num_X+1)
                         tmp_dat <- data.frame(response = c(y[,j], numeric(10)), data[c(1:nrow(data), 1:10),], off = c(new_offset, numeric(10))) # Append some zeros avoids a non-convergence problem in mgcv
                         tmp_formula <- as.formula(paste("response", paste(as.character(formula), collapse = " "), "+ offset(off)" ) )
                         if(control$ridge > 0)
                              fit0 <- gam(list(tmp_formula, ~1), data = tmp_dat, knots = knots, method = control$gam_method, 
                                          H = Hmat, family = ziplss(), select = select, gamma = full_gamma[j])
                         if(control$ridge == 0)
                              fit0 <- gam(list(tmp_formula, ~1), data = tmp_dat, knots = knots, method = control$gam_method, 
                                          family = ziplss(), select = select, gamma = full_gamma[j])
                         fit0$coefficients <- fit0$coefficients[1:num_X] 
                         fit0$linear.predictors <- X %*% fit0$coefficients + new_offset + formula_offset
                         fit0$logLik <- .dztpois(y[,j], lambda = exp(fit0$linear.predictors), log = TRUE) # .dztpois sets y = 0 values to -Inf, and handles NA values
                         fit0$logLik <- sum(fit0$logLik[is.finite(fit0$logLik)])                               
                         }
                    if(family$family %in% c("ztnegative.binomial")) { # Do EM here to update betas (since E-step is also specifically performed within here)
                         find_nonzeros <- which(y[,j] > 0)
                         inner_err <- 100
                         inner_inner_counter <- 0
                         cw_inner_logL <- -Inf
                         
                         while(inner_err > 1e-3) {
                              if(inner_inner_counter > 20)
                                   break;

                              initw <- dnbinom(0, 
                                               mu = as.vector(exp(X %*% new_fit_CBFM_ptest$betas[j,] + new_offset + formula_offset)), 
                                               size = 1/new_fit_CBFM_ptest$dispparam[j])
                              initw <- initw / (1 - initw + 1e-4)
                              initw[initw > 0.1/1e-4] <- 0.1/1e-4
                              initw <- initw[find_nonzeros]
                              w1 <- rep(1, num_units)
                              w1[which(y[,j]==0)] <- 0
                              w2 <- numeric(num_units)
                              w2[find_nonzeros] <- initw
                              w <- c(w1, w2)
                              rm(w2)
                              if(control$ridge > 0) {
                                   fit0 <- try(gam(tmp_formula, data = data.frame(response = c(y[,j],numeric(num_units)), data[c(1:num_units,1:num_units),]), 
                                                   offset = new_offset[c(1:num_units,1:num_units)], knots = knots, method = control$gam_method, 
                                                   weights = w, H = Hmat, family = nb(), select = select, gamma = full_gamma[j]), silent = TRUE)
                                   }
                              if(control$ridge == 0) {
                                   fit0 <- try(gam(tmp_formula, data = data.frame(response = c(y[,j],numeric(num_units)), data[c(1:num_units,1:num_units),]), 
                                                   offset = new_offset[c(1:num_units,1:num_units)], knots = knots, method = control$gam_method,
                                                   weights = w, family = nb(), select = select, gamma = full_gamma[j]), silent = TRUE)
                                   }
                              if(inherits(fit0, "try-error"))
                                   break;
                              
                              new_inner_logL <- .dztnbinom(y[find_nonzeros,j], 
                                                           mu = exp(X[find_nonzeros,,drop=FALSE] %*% fit0$coefficients + new_offset[find_nonzeros] + formula_offset[find_nonzeros]), 
                                                           size = fit0$family$getTheta(TRUE), log = TRUE) 
                              new_inner_logL <- sum(new_inner_logL[is.finite(new_inner_logL)]) 
                              inner_err <- abs(new_inner_logL/cw_inner_logL - 1)
                              cw_inner_logL <- new_inner_logL
                              new_fit_CBFM_ptest$betas[j,] <- fit0$coefficients
                              new_fit_CBFM_ptest$dispparam[j] <- 1/fit0$family$getTheta(TRUE)
                              inner_inner_counter <- inner_inner_counter + 1
                              }
                         
                         fit0$logLik <- new_inner_logL
                         fit0$linear.predictors <- X %*% fit0$coefficients + new_offset + formula_offset
                         }

                    
                    out <- list(coefficients = fit0$coefficients, 
                                linear.predictors = fit0$linear.predictors, 
                                logLik = fit0$logLik, 
                                fit = fit0, 
                                S = .get_bigS(fit_gam = fit0, num_X = num_X))
                    if(family$family %in% c("gaussian","Gamma"))                        
                         out$dispparam <- fit0$sig2
                    if(family$family %in% c("negative.binomial","zinegative.binomial","ztnegative.binomial"))                        
                         out$dispparam <- 1/fit0$family$getTheta(TRUE)
                    if(family$family == "Beta")                        
                         out$dispparam <- exp(fit0$family$getTheta(TRUE))
                    if(family$family == "tweedie") {                        
                         out$dispparam <- fit0$sig2
                         out$powerparam <- fit0$family$getTheta(TRUE)
                         }
                    if(family$family %in% c("zipoisson","zinegative.binomial")) {
                         out$zicoefficients <- fitzi$coefficients
                         out$fitzi <- fitzi
                         out$ziS <- .get_bigS(fit_gam = fitzi, num_X = ncol(ziX))
                         }
                    
                    return(out)
                    }               
               update_Xcoefsspp_cmpfn <- compiler::cmpfun(update_Xcoefsspp_fn)
               
               all_update_coefs <- foreach(j = 1:num_spp, .export = c("new_fit_CBFM_ptest")) %dopar% update_Xcoefsspp_cmpfn(j = j)
               new_fit_CBFM_ptest$betas <- do.call(rbind, lapply(all_update_coefs, function(x) x$coefficients))
               new_fit_CBFM_ptest$linear_predictors <- sapply(all_update_coefs, function(x) x$linear.predictors)         
               for(j in 1:num_spp) {
                    if(family$family[1] %in% c("gaussian","Gamma","negative.binomial","tweedie","Beta","zinegative.binomial","ztnegative.binomial"))
                         new_fit_CBFM_ptest$dispparam[j] <- all_update_coefs[[j]]$dispparam
                    if(family$family[1] == "tweedie") {                        
                         new_fit_CBFM_ptest$powerparam[j] <- all_update_coefs[[j]]$powerparam
                         }
                    }
               if(family$family[1] %in% c("zipoisson","zinegative.binomial")) {                        
                    new_fit_CBFM_ptest$zibetas <- do.call(rbind, lapply(all_update_coefs, function(x) x$zicoefficients))
                    }
               new_fit_CBFM_ptest$logLik <- sum(sapply(all_update_coefs, function(x) x$logLik))          

                              
               ##-------------------------
               ## Check whether to finish inner EM algorithm
               ##-------------------------
               new_inner_params <- c(c(new_fit_CBFM_ptest$betas), c(new_fit_CBFM_ptest$basis_effects_mat), c(new_fit_CBFM_ptest$zibetas), 
                               new_fit_CBFM_ptest[["mean_B_space"]], new_fit_CBFM_ptest[["mean_B_time"]], new_fit_CBFM_ptest[["mean_B_spacetime"]]) 
               inner_params_diff <- suppressWarnings(new_inner_params - cw_inner_params) # Need to suppress warnings because in first iteration, cw_params is longer than new_params 
               if(control$convergence_type == "parameters_MSE")
                    inner_err <- mean((inner_params_diff)^2) 
               if(control$convergence_type == "parameters_norm")
                    inner_err <- sum((inner_params_diff)^2) 
               if(control$convergence_type == "parameters_relative")
                    inner_err <- mean((inner_params_diff)^2)/mean(cw_inner_params^2) 
               if(control$convergence_type == "logLik_relative")
                    inner_err <- abs(new_fit_CBFM_ptest$logLik/cw_inner_logL-1)  
               
               cw_inner_logL <- new_fit_CBFM_ptest$logLik
               cw_inner_params <- new_inner_params
               inner_counter <- inner_counter + 1
               rm(all_update_coefs, update_Xcoefsspp_fn, inner_params_diff)
               }          
          
          
          ##-------------------------
          ## Update mean vector for normal distribution of basis functions coefficients, if required 
          ##-------------------------
          if(sum(which_nonzeromean_B) > 0)
               message("Updating all non-zero mean vectors in the distribution of the basis effect coefficients.")
          
          three_nonzeromean_options <- c(nonzeromean_B_space, nonzeromean_B_time, nonzeromean_B_spacetime)
          three_mean_options <- c("mean_B_space", "mean_B_time", "mean_B_spacetime")
          three_LoadingnuggetSigma_options <- c("new_LoadingnuggetSigma_space", "new_LoadingnuggetSigma_time", "new_LoadingnuggetSigma_spacetime")
          three_LoadingnuggetG_options <- c("new_LoadingnuggetG_space", "new_LoadingnuggetG_time", "new_LoadingnuggetG_spacetime")
          three_num_B <- c(num_spacebasisfns, num_timebasisfns, num_spacetimebasisfns)
          
          for(j1 in 1:3) {
               if(three_nonzeromean_options[j1]) {
                    onemIq <- kronecker(matrix(1, nrow = num_spp, ncol = 1), Diagonal(n = three_num_B[j1]))
                    MM <- crossprod(onemIq, kronecker(get(three_LoadingnuggetG_options[j1])$invcov, get(three_LoadingnuggetSigma_options[j1])$invcov))       
                    rhs <- as.vector(t(new_fit_CBFM_ptest$basis_effects_mat[, three_num_B[1]*(j1>1) + three_num_B[2]*(j1>2) + 1:three_num_B[j1], drop = FALSE])) 
                    
                    new_nonzeromean_current <- as.vector(Matrix::solve(a = MM %*% onemIq, b = MM %*% rhs))
                    assign(new_fit_CBFM_ptest[[ three_mean_options[j1] ]], new_nonzeromean_current)
                    rm(rhs, MM, onemIq, new_nonzeromean_current)
                    }
               }
          
          rm(three_nonzeromean_options, three_mean_options, three_LoadingnuggetSigma_options, three_LoadingnuggetG_options, three_num_B)
               
          ##-------------------------
          ## Update between spp correlation matrix G. First assume unstructured, then cov2cor, then update loading and nugget. 
          ##-------------------------
          if(control$trace == 1)
               message("Updating between response correlation (covariance) matrices, G.")
          
          three_custom_options <- c("custom_space", "custom_time", "custom_spacetime")
          three_mean_options <- c("mean_B_space", "mean_B_time", "mean_B_spacetime")
          three_B_options <- c("B_space", "B_time", "B_spacetime")
          three_LoadingnuggetSigma_options <- c("new_LoadingnuggetSigma_space", "new_LoadingnuggetSigma_time", "new_LoadingnuggetSigma_spacetime")
          three_LoadingnuggetG_options <- c("new_LoadingnuggetG_space", "new_LoadingnuggetG_time", "new_LoadingnuggetG_spacetime")
          three_num_B <- c(num_spacebasisfns, num_timebasisfns, num_spacetimebasisfns)
          
          for(j1 in 1:3) {
               if(which_B_used[j1]) {
                    if(is.null(G_control[[ three_custom_options[j1] ]])) {
                         
                         centered_BF_mat <- new_fit_CBFM_ptest$basis_effects_mat[, three_num_B[1]*(j1>1) + three_num_B[2]*(j1>2) + 1:three_num_B[j1], drop = FALSE] + .Machine$double.eps
                         if(c(nonzeromean_B_space, nonzeromean_B_time, nonzeromean_B_spacetime)[j1])
                              centered_BF_mat <- centered_BF_mat - matrix(new_fit_CBFM_ptest[[ three_mean_options[j1] ]], nrow = num_spp, ncol = three_num_B[j1], byrow = TRUE)
                         
                         estimate_G_as_correlation <- .check_G_correlation(custom_Sigma = Sigma_control[[ three_custom_options[j1] ]], 
                                                                           G_structure = G_control$structure[sum(which_B_used[1:j1])])
                         
                         
                         new_G_current <- .update_G_fn(Ginv = get(three_LoadingnuggetG_options[j1])$invcov, 
                                                       basis_effects_mat = centered_BF_mat, 
                                                       Sigmainv = get(three_LoadingnuggetSigma_options[j1])$invcov, 
                                                       B = get(three_B_options[j1]), 
                                                       X = X, 
                                                       ziX = ziX, 
                                                       zioffset = zioffset,
                                                       y_vec = as.vector(y), 
                                                       linpred_vec = c(new_fit_CBFM_ptest$linear_predictors), 
                                                       dispparam = new_fit_CBFM_ptest$dispparam, 
                                                       powerparam = new_fit_CBFM_ptest$powerparam, 
                                                       zibetas = new_fit_CBFM_ptest$zibetas, 
                                                       trial_size = trial_size, 
                                                       family = family, 
                                                       G_control = G_control,
                                                       use_rank_element = sum(which_B_used[1:j1]),
                                                       return_correlation = estimate_G_as_correlation)
                         
                         new_LoadingnuggetG_current <- .update_LoadingG_fn(G = new_G_current, 
                                                                           G_control = G_control, 
                                                                           use_rank_element = sum(which_B_used[1:j1]), 
                                                                           correlation = estimate_G_as_correlation)
                         
                         assign(three_LoadingnuggetG_options[j1], new_LoadingnuggetG_current)
                         rm(new_G_current, new_LoadingnuggetG_current, centered_BF_mat, estimate_G_as_correlation)
                         }
                    }     
               }
          rm(three_custom_options, three_mean_options, three_B_options, three_LoadingnuggetSigma_options, three_LoadingnuggetG_options, three_num_B)
          
          
          ##-------------------------
          ## Update covariance function for basis functions Sigma. First assuming unstructured, then update loading and nugget 
          #' For G_control$structure = "identity" and "homogeneous", the estimation of lambda in G = lambda * R is also performed here
          ##-------------------------
          if(control$trace == 1)
               message("Updating covariance matrices for basis functions, Sigma, if required.")
          
          three_custom_options <- c("custom_space", "custom_time", "custom_spacetime")
          three_mean_options <- c("mean_B_space", "mean_B_time", "mean_B_spacetime")
          three_B_options <- c("B_space", "B_time", "B_spacetime")
          three_LoadingnuggetSigma_options <- c("new_LoadingnuggetSigma_space", "new_LoadingnuggetSigma_time", "new_LoadingnuggetSigma_spacetime")
          three_LoadingnuggetG_options <- c("new_LoadingnuggetG_space", "new_LoadingnuggetG_time", "new_LoadingnuggetG_spacetime")
          three_num_B <- c(num_spacebasisfns, num_timebasisfns, num_spacetimebasisfns)
          
          for(j1 in 1:3) {
               if(which_B_used[j1]) {
                    estimate_lambda_not_Sigma <- as.numeric(G_control$structure[sum(which_B_used[1:j1])] %in% c("identity", "homogeneous"))
                    
                    if(is.null(Sigma_control[[ three_custom_options[j1] ]]) | estimate_lambda_not_Sigma == 1) {
                         
                         centered_BF_mat <- new_fit_CBFM_ptest$basis_effects_mat[, three_num_B[1]*(j1>1) + three_num_B[2]*(j1>2) + 1:three_num_B[j1], drop = FALSE] + .Machine$double.eps
                         if(c(nonzeromean_B_space, nonzeromean_B_time, nonzeromean_B_spacetime)[j1])
                              centered_BF_mat <- centered_BF_mat - matrix(new_fit_CBFM_ptest[[ three_mean_options[j1] ]], nrow = num_spp, ncol = three_num_B[j1], byrow = TRUE)
                         
                         new_Sigma_current <- .update_Sigma_fn(Sigmainv = get(three_LoadingnuggetSigma_options[j1])$invcov, 
                                                               lambdas = get(three_LoadingnuggetSigma_options[j1])$lambdas, 
                                                               basis_effects_mat = centered_BF_mat, 
                                                               Ginv = get(three_LoadingnuggetG_options[j1])$invcov, 
                                                               B = get(three_B_options[j1]), 
                                                               X = X, 
                                                               ziX = ziX, 
                                                               zioffset = zioffset,
                                                               y_vec = as.vector(y), 
                                                               linpred_vec = c(new_fit_CBFM_ptest$linear_predictors), 
                                                               dispparam = new_fit_CBFM_ptest$dispparam, 
                                                               powerparam = new_fit_CBFM_ptest$powerparam, 
                                                               zibetas = new_fit_CBFM_ptest$zibetas, 
                                                               trial_size = trial_size, 
                                                               family = family, 
                                                               Sigma_control = Sigma_control, 
                                                               estimate_lambda = estimate_lambda_not_Sigma, 
                                                               which_B = j1)
                         
                         new_LoadingnuggetSigma_current <- .update_LoadingSigma_fn(Sigma = new_Sigma_current, 
                                                                                   Sigma_control = Sigma_control, 
                                                                                   use_rank_element = sum(which_B_used[1:j1]), 
                                                                                   estimate_lambda_not_Sigma = estimate_lambda_not_Sigma, 
                                                                                   which_B = j1)
                         
                         assign(three_LoadingnuggetSigma_options[j1], new_LoadingnuggetSigma_current)
                         rm(new_Sigma_current, new_LoadingnuggetSigma_current, centered_BF_mat, estimate_lambda_not_Sigma)
                         }
                    }     
               }
          rm(three_custom_options, three_mean_options, three_B_options, three_LoadingnuggetSigma_options, three_LoadingnuggetG_options, three_num_B)
          

          ##-------------------------
          ## Finish iteration 
          ##-------------------------
          new_params <- c(c(new_fit_CBFM_ptest$betas), c(new_fit_CBFM_ptest$basis_effects_mat), c(new_fit_CBFM_ptest$zibetas), 
                          new_fit_CBFM_ptest[["mean_B_space"]], new_fit_CBFM_ptest[["mean_B_time"]], new_fit_CBFM_ptest[["mean_B_spacetime"]]) # Stop checking dispersion and power parameters
          
          new_logLik <- new_fit_CBFM_ptest$logLik
          if(which_B_used[1]) {
               centered_BF_mat <- new_fit_CBFM_ptest$basis_effects_mat[,1:num_spacebasisfns,drop=FALSE]
               if(nonzeromean_B_space)
                    centered_BF_mat <- centered_BF_mat - matrix(new_fit_CBFM_ptest[["mean_B_space"]], nrow = num_spp, ncol = num_spacebasisfns, byrow = TRUE)
               new_logLik <- new_logLik + .calc_pqlquadraticterm_basiseffects(basis_effects_mat = centered_BF_mat, 
                                                                              Ginv = new_LoadingnuggetG_space$invcov, Sigmainv = new_LoadingnuggetSigma_space$invcov)
               }
          if(which_B_used[2]) {
               centered_BF_mat <- new_fit_CBFM_ptest$basis_effects_mat[,num_spacebasisfns + (1:num_timebasisfns),drop=FALSE]
               if(nonzeromean_B_time)
                    centered_BF_mat <- centered_BF_mat - matrix(new_fit_CBFM_ptest[["mean_B_time"]], nrow = num_spp, ncol = num_timebasisfns, byrow = TRUE)
               new_logLik <- new_logLik + .calc_pqlquadraticterm_basiseffects(basis_effects_mat = centered_BF_mat, 
                                                                              Ginv = new_LoadingnuggetG_time$invcov, Sigmainv = new_LoadingnuggetSigma_time$invcov)
               }
          if(which_B_used[3]) {
               centered_BF_mat <- new_fit_CBFM_ptest$basis_effects_mat[,num_spacebasisfns + num_timebasisfns + (1:num_spacetimebasisfns),drop=FALSE]
               if(nonzeromean_B_spacetime)
                    centered_BF_mat <- centered_BF_mat - matrix(new_fit_CBFM_ptest[["mean_B_spacetime"]], nrow = num_spp, ncol = num_spacetimebasisfns, byrow = TRUE)
               new_logLik <- new_logLik + .calc_pqlquadraticterm_basiseffects(basis_effects_mat = centered_BF_mat,
                                                                              Ginv = new_LoadingnuggetG_spacetime$cov, Sigmainv = new_LoadingnuggetSigma_spacetime$invcov)
               }
          rm(centered_BF_mat)

          
          params_diff <- suppressWarnings(new_params - cw_params) # Need to suppress warnings because in first iteration, cw_params is longer than new_params 
          if(control$convergence_type == "parameters_MSE")
               diff <- mean((params_diff)^2) 
          if(control$convergence_type == "parameters_norm")
               diff <- sum((params_diff)^2) 
          if(control$convergence_type == "parameters_relative")
               diff <- mean((params_diff)^2)/mean(cw_params^2) 
          if(control$convergence_type == "logLik_relative")
               diff <- abs(new_logLik/cw_logLik-1)  
          #if(control$convergence_type == "linear_predictors")
          #     diff <- mean((new_fit_CBFM_ptest$linear_predictors - cw_linpred)^2) 
          if(control$trace > 0) {
               if(control$convergence_type == "parameters_MSE")
                    message("Iteration: ", counter, "\t Difference in parameter estimates (mean squared error): ", round(diff,5))
               if(control$convergence_type == "parameters_norm")
                    message("Iteration: ", counter, "\t Difference in parameter estimates (squared norm): ", round(diff,5))
               if(control$convergence_type == "parameters_relative")
                    message("Iteration: ", counter, "\t Relative difference in parameter estimates: ", round(diff,5))
               if(control$convergence_type == "logLik_relative")
                    message("Iteration: ", counter, "\t Relative difference in PQL values: ", round(diff,5))
               #if(control$convergence_type == "linear_predictors")
               #     message("Iteration: ", counter, "\t MSE of linear predictors: ", round(mean((new_fit_CBFM_ptest$linear_predictors - cw_linpred)^2),5))
               }
               
          cw_params <- new_params
          cw_logLik <- new_logLik
          counter <- counter + 1
          }
     toc <- proc.time()
     gc()

     
     ##----------------
     ## Do final fit -- just for coefficients and dispersion/power parameters only. 
     ##----------------
     if(diff < control$tol)
          converged <- TRUE
     tidbits_data <- make_tidibits_data()
     inner_err <- 100
     final_counter <- 0
     cw_inner_logL <- -Inf
     cw_inner_params <- cw_params
     
     while(inner_err > control$tol & final_counter <= control$final_maxit) {
          getweights <- .estep_fn(family = family, cwfit = new_fit_CBFM_ptest, y = y, X = X, B = B, ziX = ziX, zioffset = zioffset)

          all_update_coefs <- foreach(j = 1:num_spp, .export = c("tidbits_data")) %dopar% update_basiscoefsspp_cmpfn(j = j)
          for(j in 1:num_spp) {
               new_fit_CBFM_ptest$basis_effects_mat[j,] <- all_update_coefs[[j]]$par
               }
          rm(all_update_coefs)
          
          all_update_coefs <- foreach(j = 1:num_spp, .export = c("new_fit_CBFM_ptest")) %dopar% update_Xcoefsspp_cmpfn(j = j)
          new_fit_CBFM_ptest$betas <- do.call(rbind, lapply(all_update_coefs, function(x) x$coefficients))
          new_fit_CBFM_ptest$linear_predictors <- sapply(all_update_coefs, function(x) x$linear.predictors)          
          for(j in 1:num_spp) {
               if(family$family %in% c("gaussian","Gamma","negative.binomial","tweedie","Beta", "zinegative.binomial","ztnegative.binomial"))
                    new_fit_CBFM_ptest$dispparam[j] <- all_update_coefs[[j]]$dispparam
               if(family$family == "tweedie") {                        
                    new_fit_CBFM_ptest$powerparam[j] <- all_update_coefs[[j]]$powerparam
                    }
               }
          if(family$family[1] %in% c("zipoisson","zinegative.binomial")) {                        
               new_fit_CBFM_ptest$zibetas <- do.call(rbind, lapply(all_update_coefs, function(x) x$zicoefficients))
          }
          
          # Degrees of freedom calculations below may be shoddy for zero-inflated and zero-truncated distributions 
          new_fit_CBFM_ptest$edf <- sapply(all_update_coefs, function(x) x$fit$edf) 
          new_fit_CBFM_ptest$edf1 <- sapply(all_update_coefs, function(x) x$fit$edf1)    
          if(family$family[1] %in% c("ztpoisson")) { # Need to make this adjustment due to the way ziplss works in mgcv
               new_fit_CBFM_ptest$edf <- new_fit_CBFM_ptest$edf[1:num_X,]     
               new_fit_CBFM_ptest$edf1 <- new_fit_CBFM_ptest$edf1[1:num_X,]     
               }
          new_fit_CBFM_ptest$pen_edf <- lapply(all_update_coefs, function(x) pen.edf(x$fit)) 
          if(family$family[1] %in% c("zipoisson","zinegative.binomial")) {                        
               new_fit_CBFM_ptest$ziedf <- sapply(all_update_coefs, function(x) x$fitzi$edf)
               new_fit_CBFM_ptest$ziedf1 <- sapply(all_update_coefs, function(x) x$fitzi$edf1)    
               new_fit_CBFM_ptest$zipen_edf <- lapply(all_update_coefs, function(x) pen.edf(x$fitzi)) 
               }
          
          new_fit_CBFM_ptest$logLik <- sum(sapply(all_update_coefs, function(x) x$logLik))          
          new_fit_CBFM_ptest$logLik_perspp <- sapply(all_update_coefs, function(x) x$logLik)         
          
          new_inner_params <- c(c(new_fit_CBFM_ptest$betas), c(new_fit_CBFM_ptest$basis_effects_mat), c(new_fit_CBFM_ptest$zibetas), 
                          new_fit_CBFM_ptest[["mean_B_space"]], new_fit_CBFM_ptest[["mean_B_time"]], new_fit_CBFM_ptest[["mean_B_spacetime"]]) 
          inner_params_diff <- suppressWarnings(new_inner_params - cw_inner_params) # Need to suppress warnings because in first iteration, cw_params is longer than new_params 
          if(control$convergence_type == "parameters_MSE")
               inner_err <- mean((inner_params_diff)^2) 
          if(control$convergence_type == "parameters_norm")
               inner_err <- sum((inner_params_diff)^2) 
          if(control$convergence_type == "parameters_relative")
               inner_err <- mean((inner_params_diff)^2)/mean(cw_inner_params^2) 
          if(control$convergence_type == "logLik_relative")
               inner_err <- abs(new_fit_CBFM_ptest$logLik/cw_inner_logL-1)  
          
          cw_inner_logL <- new_fit_CBFM_ptest$logLik
          cw_inner_params <- new_inner_params
          final_counter <- final_counter + 1
          
          #if(!(family$family[1] %in% c("zipoisson","zinegative.binomial")))
          #     break;
          }

     all_S <- sapply(all_update_coefs, function(x) x$S)
     all_ziS <- sapply(all_update_coefs, function(x) x$ziS)
     all_k_check <- foreach(j = 1:num_spp) %dopar% k.check(all_update_coefs[[j]]$fit, subsample = k_check_control$subsample, n.rep = k_check_control$n.rep)
     names(all_k_check) <- colnames(y)
     if(family$family[1] %in% c("zipoisson", "zinegative.binomiak", "ztpoisson","ztnegative.binomial"))
          warning("k_check may not be terrible or not available for zero-inflated and zero-truncated distributions. Please take any results given here with a big grain of salt!")
     invisible(capture.output( all_vcomp <- lapply(1:num_spp, function(j) {
          if(class(all_update_coefs[[j]]$fit)[1] == "gam")
               out <- gam.vcomp(all_update_coefs[[j]]$fit)
          if(is.matrix(out))
               return(out[,1])
          if(is.list(out))
               return(out$all)
          if(is.numeric(out))
               return(out)
          }
          )))
     names(all_vcomp) <- colnames(y)
     rm(all_update_coefs, tidbits_data, inner_err, inner_params_diff, cw_inner_logL, cw_logLik, 
        cw_params, cw_inner_params, new_params, new_inner_params, diff, inner_counter, counter, final_counter)
     gc()
     
     
     # Calculate deviance, null deviance etc...note deviance calculation **excludes** the quadratic term in the PQL
     nulldeviance <- foreach(j = 1:num_spp) %dopar% initfit_fn(j = j, formula = ~ 1)
     nulldeviance_perspp <- sapply(nulldeviance, function(x) -2*x$logLik)
     nulldeviance <- sum(sapply(nulldeviance, function(x) -2*x$logLik))
     rm(initfit_fn)
     
     new_logLik <- new_fit_CBFM_ptest$logLik
     new_logLik_perspp <- new_fit_CBFM_ptest$logLik_perspp
     if(which_B_used[1]) {
          centered_BF_mat <- new_fit_CBFM_ptest$basis_effects_mat[,1:num_spacebasisfns,drop=FALSE]
          if(nonzeromean_B_space)
               centered_BF_mat <- centered_BF_mat - matrix(new_fit_CBFM_ptest[["mean_B_space"]], nrow = num_spp, ncol = num_spacebasisfns, byrow = TRUE)
          new_logLik_pql <- new_logLik + .calc_pqlquadraticterm_basiseffects(basis_effects_mat = centered_BF_mat, 
                                                                             Ginv = new_LoadingnuggetG_space$invcov, Sigmainv = new_LoadingnuggetSigma_space$invcov)
          }
     if(which_B_used[2]) {
          centered_BF_mat <- new_fit_CBFM_ptest$basis_effects_mat[,num_spacebasisfns + (1:num_timebasisfns),drop=FALSE]
          if(nonzeromean_B_time)
               centered_BF_mat <- centered_BF_mat - matrix(new_fit_CBFM_ptest[["mean_B_time"]], nrow = num_spp, ncol = num_timebasisfns, byrow = TRUE)
          new_logLik_pql <- new_logLik + .calc_pqlquadraticterm_basiseffects(basis_effects_mat = centered_BF_mat, 
                                                                             Ginv = new_LoadingnuggetG_time$invcov, Sigmainv = new_LoadingnuggetSigma_time$invcov)
          }
     if(which_B_used[3]) {
          centered_BF_mat <- new_fit_CBFM_ptest$basis_effects_mat[,num_spacebasisfns + num_timebasisfns + (1:num_spacetimebasisfns),drop=FALSE]
          if(nonzeromean_B_time)
               centered_BF_mat <- centered_BF_mat - matrix(new_fit_CBFM_ptest[["mean_B_spacetime"]], nrow = num_spp, ncol = num_spacetimebasisfns, byrow = TRUE)
          new_logLik_pql <- new_logLik + .calc_pqlquadraticterm_basiseffects(basis_effects_mat = centered_BF_mat,
                                                                             Ginv = new_LoadingnuggetG_spacetime$cov, Sigmainv = new_LoadingnuggetSigma_spacetime$invcov)
          }
     rm(centered_BF_mat)


     ##-----------------
     ## Format output
     ##-----------------
     colnames(B) <- underlying_B_colnames
     out_CBFM <- list(call = match.call())
     out_CBFM$family <- family
     out_CBFM$y <- y
     out_CBFM$data <- data
     out_CBFM$trial_size <- trial_size
     out_CBFM$formula <- formula
     out_CBFM$ziformula <- ziformula
     out_CBFM$select <- select
     out_CBFM$gamma <- full_gamma
     out_CBFM$knots <- knots
     out_CBFM$ziselect <- ziselect
     out_CBFM$zigamma <- full_zigamma
     out_CBFM$ziknots <- ziknots
     out_CBFM$B <- B
     out_CBFM$which_B_used <- which_B_used
     out_CBFM$which_custom_Sigma_used <- Sigma_control$which_custom_Sigma_used
     out_CBFM$which_custom_G_used <- G_control$which_custom_G_used
     out_CBFM$which_nonzeromean_B <- which_nonzeromean_B
     out_CBFM$num_B_space <- num_spacebasisfns
     out_CBFM$num_B_time <- num_timebasisfns
     out_CBFM$num_B_spacetime <- num_spacetimebasisfns
     out_CBFM$num_B <- num_basisfns
     out_CBFM$converged <- converged
     out_CBFM$logLik <- new_logLik
     out_CBFM$logLik_perspecies <- new_logLik_perspp
     out_CBFM$deviance <- -2*out_CBFM$logLik
     out_CBFM$deviance_perspecies <- -2*out_CBFM$logLik_perspecies
     out_CBFM$null_deviance <- nulldeviance
     out_CBFM$null_deviance_perspecies <- nulldeviance_perspp
     out_CBFM$deviance_explained <- 100*(out_CBFM$null_deviance - out_CBFM$deviance)/out_CBFM$null_deviance
     out_CBFM$deviance_explained_perspecies <- 100*(out_CBFM$null_deviance_perspecies - out_CBFM$deviance_perspecies)/out_CBFM$null_deviance_perspecies
     out_CBFM$pql_logLik <- new_logLik_pql
     out_CBFM$edf <- matrix(new_fit_CBFM_ptest$edf, ncol = num_spp)
     out_CBFM$edf1 <- matrix(new_fit_CBFM_ptest$edf1, ncol = num_spp)
     out_CBFM$pen_edf <- new_fit_CBFM_ptest$pen_edf
     if(family$family[1] %in% c("zipoisson","zinegative.binomial")) {
          out_CBFM$ziedf <- matrix(new_fit_CBFM_ptest$ziedf, ncol = num_spp)
          out_CBFM$ziedf1 <- matrix(new_fit_CBFM_ptest$ziedf1, ncol = num_spp)
          out_CBFM$zipen_edf <- new_fit_CBFM_ptest$zipen_edf
          }
     out_CBFM$k_check <- all_k_check
     out_CBFM$vcomp <- all_vcomp

     out_CBFM$betas <- new_fit_CBFM_ptest$betas
     out_CBFM$zibetas <- new_fit_CBFM_ptest$zibetas
     out_CBFM$basis_effects_mat <- new_fit_CBFM_ptest$basis_effects_mat
     out_CBFM$dispparam <- new_fit_CBFM_ptest$dispparam
     out_CBFM$powerparam <- new_fit_CBFM_ptest$powerparam
     out_CBFM$mean_B_space <- new_fit_CBFM_ptest[["mean_B_space"]]
     out_CBFM$mean_B_time <- new_fit_CBFM_ptest[["mean_B_time"]]
     out_CBFM$mean_B_spacetime <- new_fit_CBFM_ptest[["mean_B_spacetime"]]
     out_CBFM$linear_predictors <- new_fit_CBFM_ptest$linear_predictors
     rm(new_fit_CBFM_ptest, getweights, converged, all_k_check)

     # ## Get restricted CBFM estimates of the coefficients (really more for exploration at this point in time)
     #OLSmatrix_transpose <- X %*% solve(crossprod(X))
     #out_CBFM$partitioned_beta <- out_CBFM$betas + tcrossprod(out_CBFM$basis_effects_mat, B) %*% OLSmatrix_transpose
     #rm(OLSmatrix_transpose)
     out_CBFM$linear_predictors[is.na(out_CBFM$y)] <- NA

     if(!(family$family %in% c("zipoisson","zinegative.binomial","ztpoisson","ztnegative.binomial"))) 
          out_CBFM$fitted <- family$linkinv(out_CBFM$linear_predictors)
     if(family$family %in% c("zipoisson","zinegative.binomial")) {
          zieta <- tcrossprod(ziX, out_CBFM$zibetas)
          if(!is.null(zieta))
               zieta <- zieta + zioffset
          zeroinfl_prob <- plogis(zieta)
          rm(zieta)
          out_CBFM$fitted <- family$linkinv(out_CBFM$linear_predictors) * (1-zeroinfl_prob)
          rm(zeroinfl_prob)
          }
     if(family$family == "ztpoisson") {
          out_CBFM$fitted <- exp(out_CBFM$linear_predictors) / (1 - dpois(0, lambda = exp(out_CBFM$linear_predictors)))
          out_CBFM$fitted[out_CBFM$fitted < 1 | out_CBFM$fitted == Inf] <- 1 # Predictions less than 1 or at infinity will be almost certainly due to the linear predictor being so close to zero that you get underflow issues...unless the linear predictor is stupidly large
          out_CBFM$linear_predictors[which(out_CBFM$y == 0)] <- NA
          out_CBFM$fitted[which(out_CBFM$y == 0)] <- NA
          }
     if(family$family == "ztnegative.binomial") {
          out_CBFM$fitted <- exp(out_CBFM$linear_predictors) / (1 - dnbinom(0, mu = exp(out_CBFM$linear_predictors), size = matrix(1/out_CBFM$dispparam, nrow = num_units, ncol = num_spp, byrow = TRUE)))
          out_CBFM$fitted[out_CBFM$fitted < 1 | out_CBFM$fitted == Inf] <- 1 # Predictions less than 1 will be almost certainly due to the linear predictor being so close to zero that you get underflow issues...unless the linear predictor is stupidly large
          out_CBFM$linear_predictors[which(out_CBFM$y == 0)] <- NA
          out_CBFM$fitted[which(out_CBFM$y == 0)] <- NA
          }
     
     names(out_CBFM$dispparam) <- names(out_CBFM$powerparam) <- colnames(out_CBFM$edf) <- colnames(out_CBFM$edf1) <- names(out_CBFM$logLik_perspecies) <- 
          names(out_CBFM$deviance_perspecies) <- names(out_CBFM$null_deviance_perspecies) <- names(out_CBFM$deviance_explained_perspecies) <- 
          names(out_CBFM$pen_edf) <- names(out_CBFM$k_check) <- names(out_CBFM$vcomp) <- colnames(y)
     if(family$family %in% c("zipoisson","zinegative.binomial")) {
          colnames(out_CBFM$ziedf) <- colnames(out_CBFM$ziedf1) <- names(out_CBFM$zipen_edf) <- colnames(y)
          }
     
     if(which_B_used[1]) {
          if(is.null(Sigma_control[["custom_space"]])) {
               out_CBFM$Sigma_space <- new_LoadingnuggetSigma_space$cov
               rownames(out_CBFM$Sigma_space) <- colnames(out_CBFM$Sigma_space) <- colnames(B_space)
               if(Sigma_control$rank[1] != "full") {
                    out_CBFM$Loading_Sigma_space <- as.matrix(new_LoadingnuggetSigma_space$Loading)
                    out_CBFM$nugget_Sigma_space <- new_LoadingnuggetSigma_space$nugget
                    rownames(out_CBFM$Loading_Sigma_space) <- colnames(B_space)
                    if(num_spacebasisfns > 2)
                         colnames(out_CBFM$Loading_Sigma_space) <- paste0("Loading", 1:Sigma_control$rank[1])
                    }
               rm(new_LoadingnuggetSigma_space)
               }
          if(!is.null(Sigma_control[["custom_space"]])) {
               out_CBFM$Sigma_space <- new_LoadingnuggetSigma_space$cov
               out_CBFM$Loading_Sigma_space <- out_CBFM$nugget_Sigma_space <- NULL
               }
          
          if(is.null(G_control[["custom_space"]])) {
               out_CBFM$G_space <- new_LoadingnuggetG_space$cov
               rownames(out_CBFM$G_space) <- colnames(out_CBFM$G_space) <- colnames(y)
               if(G_control$rank[1] != "full" & G_control$structure[1] != "identity") {
                    if(num_spp > 1) {
                         out_CBFM$Loading_G_space <- as.matrix(new_LoadingnuggetG_space$Loading)
                         out_CBFM$nugget_G_space <- new_LoadingnuggetG_space$nugget
                         rownames(out_CBFM$Loading_G_space) <- colnames(y)
                         }
                    if(num_spp > 2)
                         colnames(out_CBFM$Loading_G_space) <- paste0("Loading", 1:G_control$rank[1])
                    }
               rm(new_LoadingnuggetG_space)
               }
          if(!is.null(G_control[["custom_space"]])) {
               out_CBFM$G_space <- new_LoadingnuggetG_space$cov
               out_CBFM$Loading_G_space <- out_CBFM$nugget_G_space <- NULL
               }
          
          if(nonzeromean_B_space)
               names(out_CBFM[["mean_B_space"]]) <- colnames(B_space)
          }          
     if(which_B_used[2]) {
          if(is.null(Sigma_control[["custom_time"]])) {
               out_CBFM$Sigma_time <- new_LoadingnuggetSigma_time$cov
               rownames(out_CBFM$Sigma_time) <- colnames(out_CBFM$Sigma_time) <- colnames(B_time)
               if(Sigma_control$rank[sum(which_B_used[1:2])] != "full") {
                    out_CBFM$Loading_Sigma_time <- as.matrix(new_LoadingnuggetSigma_time$Loading)
                    out_CBFM$nugget_Sigma_time <- new_LoadingnuggetSigma_time$nugget
                    rownames(out_CBFM$Loading_Sigma_time) <- colnames(B_time)
                    if(num_timebasisfns > 2)
                         colnames(out_CBFM$Loading_Sigma_time) <- paste0("Loading", 1:Sigma_control$rank[sum(which_B_used[1:2])])
                    }
               rm(new_LoadingnuggetSigma_time)
               }
          if(!is.null(Sigma_control[["custom_time"]])) {
               out_CBFM$Sigma_time <- new_LoadingnuggetSigma_time$cov
               out_CBFM$Loading_Sigma_time <- out_CBFM$nugget_Sigma_time <- NULL
               }

          if(is.null(G_control[["custom_time"]])) {
               out_CBFM$G_time <- new_LoadingnuggetG_time$cov
               rownames(out_CBFM$G_time) <- colnames(out_CBFM$G_time) <- colnames(y)
               if(G_control$rank[sum(which_B_used[1:2])] != "full" & G_control$structure[sum(which_B_used[1:2])] != "identity") {
                    if(num_spp > 1) {
                         out_CBFM$Loading_G_time <- as.matrix(new_LoadingnuggetG_time$Loading)
                         out_CBFM$nugget_G_time <- new_LoadingnuggetG_time$nugget
                         rownames(out_CBFM$Loading_G_time) <- colnames(y)
                         }
                    if(num_spp > 2)
                         colnames(out_CBFM$Loading_G_time) <- paste0("Loading", 1:G_control$rank[sum(which_B_used[1:2])])
                    }
               rm(new_LoadingnuggetG_time)
               }
          if(!is.null(G_control[["custom_time"]])) {
               out_CBFM$G_time <- new_LoadingnuggetG_time$cov
               out_CBFM$Loading_G_time <- out_CBFM$nugget_G_time <- NULL
               }
          
          if(nonzeromean_B_time)
               names(out_CBFM[["mean_B_time"]]) <- colnames(B_time)
          }
     if(which_B_used[3]) {
          if(is.null(Sigma_control[["custom_spacetime"]])) {
               out_CBFM$Sigma_spacetime <- new_LoadingnuggetSigma_spacetime$cov
               rownames(out_CBFM$Sigma_spacetime) <- colnames(out_CBFM$Sigma_spacetime) <- colnames(B_spacetime)          
               if(Sigma_control$rank[sum(which_B_used[1:3])] != "full") {
                    out_CBFM$Loading_Sigma_spacetime <- as.matrix(new_LoadingnuggetSigma_spacetime$Loading)
                    out_CBFM$nugget_Sigma_spacetime <- new_LoadingnuggetSigma_spacetime$nugget
                    rownames(out_CBFM$Loading_Sigma_spacetime) <- colnames(B_spacetime)          
                    if(num_spacetimebasisfns > 2)
                         colnames(out_CBFM$Loading_Sigma_spacetime) <- paste0("Loading", 1:Sigma_control$rank[sum(which_B_used[1:3])])
                    }
               rm(new_LoadingnuggetSigma_spacetime)
               }
          if(!is.null(Sigma_control[["custom_spacetime"]])) {
               out_CBFM$Sigma_spacetime <- new_LoadingnuggetSigma_spacetime$cov
               out_CBFM$Loading_Sigma_spacetime <- out_CBFM$nugget_Sigma_spacetime <- NULL
               }
          
          if(is.null(G_control[["custom_spacetime"]])) {
               out_CBFM$G_spacetime <- new_LoadingnuggetG_spacetime$cov
               rownames(out_CBFM$G_spacetime) <- colnames(out_CBFM$G_spacetime) <- colnames(y)
               if(G_control$rank[sum(which_B_used[1:3])] != "full" & G_control$structure[sum(which_B_used[1:3])] != "identity") {
                    if(num_spp > 1) {
                         out_CBFM$Loading_G_spacetime <- as.matrix(new_LoadingnuggetG_spacetime$Loading)
                         out_CBFM$nugget_G_spacetime <- new_LoadingnuggetG_spacetime$nugget
                         rownames(out_CBFM$Loading_G_spacetime) <- colnames(y)
                         }
                    if(num_spp > 2)
                         colnames(out_CBFM$Loading_G_spacetime) <- paste0("Loading", 1:G_control$rank[sum(which_B_used[1:3])])
                    }
               rm(new_LoadingnuggetG_spacetime)
               }
          if(!is.null(G_control[["custom_spacetime"]])) {
               out_CBFM$G_spacetime <- new_LoadingnuggetG_spacetime$cov
               out_CBFM$Loading_G_spacetime <- out_CBFM$nugget_G_spacetime <- NULL
               }
          
          if(nonzeromean_B_spacetime)
               names(out_CBFM[["mean_B_spacetime"]]) <- colnames(B_spacetime)
          }

     rownames(out_CBFM$betas) <- rownames(out_CBFM$basis_effects_mat) <- colnames(out_CBFM$linear_predictors) <- colnames(out_CBFM$fitted) <- colnames(y)
     rownames(out_CBFM$linear_predictors) <- rownames(out_CBFM$fitted) <- rownames(X)
     colnames(out_CBFM$betas) <- colnames(X)
     colnames(out_CBFM$basis_effects_mat) <- colnames(B)     
     if(family$family[1] %in% c("zipoisson","zinegative.binomial")) {
          rownames(out_CBFM$zibetas) <- colnames(y)
          colnames(out_CBFM$zibetas) <- colnames(ziX)
          }

     
     ##-----------------
     ## Calculate structures needed for producing standard errors for coefficients
     ## The Bayesian posterior covariance matrix used, as opposed to the frequentist sandwich form. This is consistent with the default available in mgcv     
     ## Similar to the default with summary.gam in mgcv, the uncertainty in the nuisance parameters or the covariance matrix is not accounted for! 
     ## Make use of blockwise inversion 
     ##-----------------
     out_CBFM$stderrors <- stderrors
     if(stderrors) {          
          if(control$trace)
               message("Calculating (components of) the covariance (standard error) matrix...")
          
          zieta <- NULL
          if(family$family[1] %in% c("zipoisson","zinegative.binomial")) {
               zieta <- tcrossprod(ziX, out_CBFM$zibetas)
               if(!is.null(zioffset))
                    zieta <- zieta + zioffset
               }
          
          weights_mat <- .neghessfamily(family = family, eta = out_CBFM$linear_predictors, y = y, 
                                        phi = matrix(out_CBFM$dispparam, num_units, num_spp, byrow = TRUE), 
                                        powerparam = matrix(out_CBFM$powerparam, num_units, num_spp, byrow = TRUE),
                                        zieta = zieta, trial_size = trial_size, domore = TRUE)
          if(family$family[1] %in% c("ztpoisson", "ztnegative.binomial"))
               weights_mat$out[is.na(weights_mat$out)] <- 0
          if(!(family$family[1] %in% c("zipoisson","zinegative.binomial"))) {
               weights_mat <- matrix(weights_mat$out, nrow = num_units, ncol = num_spp) # Overwrite weights_mat since only one quantity needed
               weights_mat[is.na(out_CBFM$y)] <- 0
               }
          if(family$family[1] %in% c("zipoisson","zinegative.binomial")) {
               weights_mat_betabeta <- matrix(weights_mat$out, nrow = num_units, ncol = num_spp)
               weights_mat_betabeta[is.na(out_CBFM$y)] <- 0
               weights_mat$out_zeroinflzeroinfl[is.na(out_CBFM$y)] <- 0
               weights_mat$out_zeroinflbetas[is.na(out_CBFM$y)] <- 0
               }
          
          # Bottom right of covariance matrix
          D1minusCAinvB_fn <- function(j) {                
               if(!(family$family[1] %in% c("zipoisson","zinegative.binomial"))) {
                    XTWX_inv <- crossprod(X*sqrt(weights_mat[,j])) + Diagonal(x = control$ridge+1e-15, n = num_X) + all_S[[j]]
                    XTWX_inv <- .pinv( 0.5*(XTWX_inv + t(XTWX_inv)) )
                    BTWX <- crossprod(B, X*weights_mat[,j])               
                    return(crossprod(B*sqrt(weights_mat[,j])) - BTWX %*% tcrossprod(XTWX_inv, BTWX))
                    }

               if(family$family[1] %in% c("zipoisson","zinegative.binomial")) {
                    Xi <- bdiag(ziX, X)
                    bigW <- cbind(Diagonal(x = weights_mat$out_zeroinflzeroinfl[,j]),  Diagonal(x = weights_mat$out_zeroinflbetas[,j]))
                    bigW <- rbind(bigW, cbind(Diagonal(x = weights_mat$out_zeroinflbetas[,j]),  Diagonal(x = weights_mat_betabeta[,j])))
                    XTWX_inv <- crossprod(Xi, bigW) %*% Xi + Diagonal(x = rep(c(control$ridge+1e-15,control$ziridge+1e-15), c(num_X,ncol(ziX)))) + bdiag(all_ziS[[j]], all_S[[j]])
                    XTWX_inv <- .pinv( 0.5*(XTWX_inv + t(XTWX_inv)) )
                    BTWX <- crossprod(B, cbind(Diagonal(x = weights_mat$out_zeroinflbetas[,j]), Diagonal(x = weights_mat_betabeta[,j]))) %*% Xi               
                    return(crossprod(B*sqrt(weights_mat_betabeta[,j])) - BTWX %*% tcrossprod(XTWX_inv, BTWX))
                    }
               }
          all_D1minusCAinvB <- foreach(j = 1:num_spp) %dopar% D1minusCAinvB_fn(j = j)
          all_D1minusCAinvB <- bdiag(all_D1minusCAinvB)
          if(identical(which_B_used, c(1,0,0)))
               DminusCAinvB_inv <- forceSymmetric(all_D1minusCAinvB + kronecker(.pinv(out_CBFM$G_space), .pinv(out_CBFM$Sigma_space)))
          if(identical(which_B_used, c(0,1,0)))
               DminusCAinvB_inv <- forceSymmetric(all_D1minusCAinvB + kronecker(.pinv(out_CBFM$G_time), .pinv(out_CBFM$Sigma_time)))
          if(identical(which_B_used, c(0,0,1)))
               DminusCAinvB_inv <- forceSymmetric(all_D1minusCAinvB + kronecker(.pinv(out_CBFM$G_spacetime), .pinv(out_CBFM$Sigma_spacetime)))
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
          DminusCAinvB_inv <- .pinv(DminusCAinvB_inv) ## Bottleneck! 

          # Top right of covariance matrix
          AinvandB_fn <- function(j) {
               if(!(family$family[1] %in% c("zipoisson","zinegative.binomial"))) {
                    XTWX_inv <- crossprod(X*sqrt(weights_mat[,j])) + Diagonal(x=control$ridge+1e-15, n = num_X) + all_S[[j]]
                    XTWX_inv <- .pinv( 0.5*(XTWX_inv + t(XTWX_inv)) )
                    XTWB <- crossprod(X*weights_mat[,j], B)
                    }
               if(family$family[1] %in% c("zipoisson","zinegative.binomial")) {
                    Xi <- bdiag(ziX, X)
                    bigW <- cbind(Diagonal(x = weights_mat$out_zeroinflzeroinfl[,j]),  Diagonal(x = weights_mat$out_zeroinflbetas[,j]))
                    bigW <- rbind(bigW, cbind(Diagonal(x = weights_mat$out_zeroinflbetas[,j]),  Diagonal(x = weights_mat_betabeta[,j])))
                    XTWX_inv <- crossprod(Xi, bigW) %*% Xi + Diagonal(x = rep(c(control$ridge+1e-15,control$ziridge+1e-15), c(num_X,ncol(ziX)))) + bdiag(all_ziS[[j]], all_S[[j]])
                    XTWX_inv <- .pinv( 0.5*(XTWX_inv + t(XTWX_inv)) )
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
                    
                    
          if(!(family$family[1] %in% c("zipoisson","zinegative.binomial")))
               rownames(out_CBFM$covar_components$topleft) <- colnames(out_CBFM$covar_components$topleft) <- rownames(out_CBFM$covar_components$topright) <-
                    apply(as.data.frame.table(t(out_CBFM$betas))[,1:2],1,function(x) paste(x, collapse = ":"))
          if(family$family[1] %in% c("zipoisson","zinegative.binomial")) {
               make_tab <- cbind(out_CBFM$zibetas, out_CBFM$betas)
               colnames(make_tab)[1:ncol(ziX)] <- paste0("ZeroInfl", colnames(make_tab)[1:ncol(ziX)])
               rownames(out_CBFM$covar_components$topleft) <- colnames(out_CBFM$covar_components$topleft) <- rownames(out_CBFM$covar_components$topright) <-
                    apply(as.data.frame.table(t(make_tab))[,1:2],1,function(x) paste(x, collapse = ":"))
               rm(make_tab) 
               }
          
          colnames(out_CBFM$covar_components$topright) <- rownames(out_CBFM$covar_components$bottomright) <- colnames(out_CBFM$covar_components$bottomright) <-
               apply(as.data.frame.table(t(out_CBFM$basis_effects_mat))[,1:2],1,function(x) paste(x, collapse = ":"))
          } 
 
 
     
     ##-----------------
     ## Final touches!
     ##-----------------
     if(!(family$family[1] %in% c("Beta","gaussian","Gamma","negative.binomial","tweedie","zinegative.binomial","ztnegative.binomial")))
          out_CBFM$dispparam <- NULL
     if(!(family$family %in% c("tweedie")))                        
          out_CBFM$powerparam <- NULL
     if(!(family$family %in% c("zipoisson","zinegative.binomial")))                        
          out_CBFM$zibetas <- NULL
     if(!nonzeromean_B_space)
          out_CBFM[["mean_B_space"]] <- NULL
     if(!nonzeromean_B_time)
          out_CBFM[["mean_B_time"]] <- NULL
     if(!nonzeromean_B_spacetime)
          out_CBFM[["mean_B_spacetime"]] <- NULL
          
     
     out_CBFM$time_taken <- toc[3] - tic[3] 
     class(out_CBFM) <- "CBFM"
     return(out_CBFM)
     }
     
