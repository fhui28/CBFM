#' @title Create a hurdle CBFM 
#' 
#' @description 
#' `r lifecycle::badge("experimental")`
#' 
#' Given one CBFM fitted to presence-absence data, and a second CBFM fitted to zero-truncated count data, this function simply links the two fitted CBFM to create a hurdle CBFM object. 
#'
#' @param pa_object An object of class \code{CBFM}. A check is made to see if the CBFM was fitted using the binary/binomial family.
#' @param count_object An object of class \code{CBFM}. A check is made to see if the CBFM was fitted using a zero-truncated count family.
#' @param ... Not used.
#'
#'  
#' @details 
#' For modeling count data, a hurdle regression model comprises two independent, component models: 
#' 
#' 1. A presence-absence model fitted to the presence-absence form of the data i.e., the data is converted to ones and zeros.   
#' 2. A count model fitted to the zero-truncated count form the data i.e., a subset of the data including only non-zero counts. 
#' 
#' As the name suggests, this class of models therefore separately models the probability of overcoming the "hurdle" i.e., being present, and then distribution of the counts conditional on overcoming the hurdle i.e., the observed count conditional on the counts being positive; we refer to Cragg (1971) for the original hurdle regression model idea, and Welsh et al., (1996); Barry and Welsh (2002); Potts and Elith (2006); Cantoni et al., (2017); Sadykova et al., (2017); Smith et al., (2019) among many others for their use in (potentially spatial and/or temporal) species distribution modeling in ecology. Two well established packages used to two hurdle regressions include [pscl::hurdle()], and zero-altered distributions available [gamlss::gamlss()]. Also, motivated strongly by fisheries data, the excellent VAST package by Thorson (2019) has hurdle count families available for spatio-temporal joint species distribution modeling, using an latent variable modeling (LVM) framework.
#' 
#' Because a hurdle model consists of two separate models, and assuming the spatio-temporal fields for the presence-absence and count components are independent of each other, then a hurdle CBFM can be constructed by "simply" fitting a presence-absence CBFM, then a zero-truncated count CBFM, and then linking the two model together. This function's sole aim is to perform the linkage: it creates an object of class "CBFM_hurdle" containing the two CBFMs themselves. This object can be passed downstream for use in selected functions such as [AIC.CBFM_hurdle()], [fitted.CBFM_hurdle()], and [predict.CBFM_hurdle()]. Note the object still carries the two component CBFMs, so the user can operate on each of the two components using all standard CBFM functions.
#' 
#' Mathematically, and for completeness, the hurdle CBFM is characterized by the following two regression models: for observational unit \eqn{i=1,\ldots,N} and species \eqn{j=1,\ldots,m}, the presence-absence component CBFM is characterized by
#' 
#' \deqn{logit(\mu_{ij,pa}) = x_{i,pa}^\top\beta_{j,pa} + b_{i,pa}^\top a_{j,pa},}
#'
#' where \eqn{logit(.)} is the logit link function, \eqn{x_{i,pa}} denotes a vector of predictors for unit \eqn{i} i.e., the \eqn{i}-th row from the created model matrix, \eqn{\beta_{j,pa}} denotes the corresponding regression coefficients for species \eqn{j}, \eqn{b_{i,pa}} denotes a vector of spatial and/or temporal basis functions for unit \eqn{i}, and \eqn{a_{j,pa}} denotes the corresponding regression coefficients for species \eqn{j}. Analogously, for the count component CBFM is characterized by
#' 
#' \deqn{log(\mu_{ij,count}) = x_{i,count}^\top\beta_{j,count} + b_{i,count}^\top a_{j,count},}
#'
#' where a log link function is used and all the covariates, basis functions, and coefficients are defined similarly to the above. For observational unit \eqn{i} and species \eqn{j}, the mean probability of presence is given by \eqn{\mu_{ij,pa}}, while the mean count conditional on the species being presence is given by \eqn{\mu_{ij,count}}. 
#' 
#' Note because the user is "forced" into fitting two components CBFMs separately (yes we purposefully do this!), then different covariates and basis functions are allowed for the two components.
#'     
#' 
#' @return An object of class "CBFM_hurdle" which includes the following two components:
#' \item{pa_fit: }{Equivalent to \code{pa_object} i.e., a CBFM fitted using the binary/binomial family.}
#' 
#' \item{count_fit: }{Equivalent to \code{count_object} i.e., a CBFM fitted using a zero-truncated count family.}
#' 
#' @details # Warning
#' Just because a hurdle CBFM is available does not necessarily mean you should always use it! It is a relatively sophisticated model with a large number of parameters, and when appropriate it may be better to adopt a simpler model such as the negative binomial CBFM instead; please see [CBFM()] for other families for handling spatio-temporal multivariate count data. 
#' 
#' 
#' @author Francis K.C. Hui <fhui28@gmail.com>, Chris Haak
#' 
#' 
#' @references
#' Barry, S. C., and Welsh, A. H. (2002). Generalized additive modelling and zero inflated count data. Ecological Modelling, 157, 179-188.
#' 
#' Cantoni, E., Flemming, J. M., and Welsh, A. H. (2017). A random-effects hurdle model for predicting bycatch of endangered marine species. The Annals of Applied Statistics, 11, 2178-2199.
#' 
#' Cragg, J. G. (1971). Some statistical models for limited dependent variables with application to the demand for durable goods. Econometrica: Journal of the Econometric Society, 39, 829-844.
#' 
#' Potts, J. M., and Elith, J. (2006). Comparing species abundance models. Ecological modelling, 199, 153-163.
#' 
#' Sadykova, D., Scott, B. E., De Dominicis, M., Wakelin, S. L., Sadykov, A., and Wolf, J. (2017). Bayesian joint models with INLA exploring marine mobile predator–prey and competitor species habitat overlap. Ecology and Evolution, 7, 5212-5226.
#' 
#' Smith, A., Hofner, B., Lamb, J. S., Osenkowski, J., Allison, T., Sadoti, G., and Paton, P. (2019). Modeling spatio-temporal abundance of mobile wildlife in highly variable environments using boosted GAMLSS hurdle models. Ecology and evolution, 9, 2346-2364.
#' 
#' Thorson, J. T. (2019). Guidance for decisions using the Vector Autoregressive Spatio-Temporal (VAST) package in stock, ecosystem, habitat and climate assessments. Fisheries Research, 210, 143-161.
#' 
#' Welsh, A. H., Cunningham, R. B., Donnelly, C. F., and Lindenmayer, D. B. (1996). Modelling the abundance of rare species: statistical models for counts with extra zeros. Ecological Modelling, 88, 297-308.
#' 
#' 
#' @seealso [AIC.CBFM_hurdle()] and [AICc.CBFM_hurdle()] for calculation various information criteria from a hurdle CBFM fit, [fitted.CBFM_hurdle()] for extracting the fitted values from a hurdle CBFM fit, [logLik.CBFM_hurdle()] for extracting the log-likelihood from a hurdle CBFM fit, [plot.CBFM_hurdle()] for basic residual diagnostics from a hurdle CBFM fit, [predict.CBFM_hurdle()] for constructing predictions from a hurdle CBFM fit, [residuals.CBFM_hurdle()] for constructing residuals from a hurdle CBFM fit, [simulate.CBFM_hurdle()] for simulating spatio-temporal multivariate abundance data from a hurdle CBFM fit. Note the component CBFMs making up the object can be operated on using standard CBFM functions; see the [CBFM()] help file for more information.
#' 
#' 
#' @examples
#' \dontrun{
#' library(autoFRK)
#' library(gamlss.tr)
#' library(FRK)
#' library(MASS)
#' library(mvtnorm)
#' library(sp)
#' library(geoR)
#' library(tidyverse)
#' 
#' 
#' ##------------------------------
#' ## **Example 1: Fitting a CBFM to spatial multivariate hurdle Poisson data**
#' ## simulated from a spatial hurdle latent variable model
#' ## Please note the data generation process (thus) differs from CBFM.
#' ##------------------------------
#' set.seed(2021)
#' num_sites <- 1000 # 500 (units) sites for training set + 500 sites for testing.
#' num_spp <- 50 # Number of species
#' num_X <- 4 # Number of regression slopes
#' 
#' # Generate and combine latent variables model for presence-absence and zero truncated Poisson X
#' spp_slopes_pa <- matrix(runif(num_spp * num_X, -1, 1), nrow = num_spp)
#' spp_slopes_ztp <- matrix(runif(num_spp * num_X, -1, 1), nrow = num_spp)
#' spp_intercepts_pa <- runif(num_spp, -2, 0)
#' spp_intercepts_ztp <- runif(num_spp, -2, 0)
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
#' # Simulate latent variable components
#' true_lvs_pa <- grf(grid = cbind(xy$x, xy$y), nsim = 2, cov.model = "exponential",
#' cov.pars = c(1, 2))$data %>%
#'      as.matrix
#' spp_loadings_pa <- matrix(runif(num_spp * 2, -1, 1), nrow = num_spp)
#' true_lvs_ztp <- grf(grid = cbind(xy$x, xy$y), nsim = 2, cov.model = "exponential",
#' cov.pars = c(1, 2.5))$data %>%
#'      as.matrix
#' spp_loadings_ztp <- matrix(runif(num_spp * 2, -1, 1), nrow = num_spp) 
#' set.seed(NULL)
#' 
#' # Simulate spatial multivariate abundance data
#' eta <- tcrossprod(cbind(1,mm), cbind(spp_intercepts_pa,spp_slopes_pa)) + 
#' tcrossprod(true_lvs_pa, spp_loadings_pa)
#' simy_pa <- matrix(rbinom(num_sites * num_spp, size = 1, 
#' prob = plogis(eta)), nrow = num_sites)
#' 
#' # Now simulate spatial count data from a truncated Poisson distribution
#' # Note the use of sapply as ztpR behaves oddly when using vectorized arugments
#' eta <- tcrossprod(cbind(1,mm), cbind(spp_intercepts_ztp,spp_slopes_ztp)) + 
#' tcrossprod(true_lvs_ztp, spp_loadings_ztp)
#' ztpR <- trun.r(par = 0, family = "PO", type = "left") 
#' simy_ztp <- matrix(sapply(1:(num_sites*num_spp), 
#' function(x) ztpR(1, mu = exp(eta[x]))), nrow = num_sites)
#' 
#' # Spatial multivariate count data from a hurdle Poisson model is then the product of the two
#' simy <- simy_pa *  simy_ztp
#' 
#' # Form training and test sets
#' dat_train <- dat[1:500,]
#' dat_test <- dat[501:1000,]
#' simy_train <- simy[1:500,]
#' simy_test <- simy[501:1000,]
#' 
#' # Delete the "component" responses and present you only observe the final response
#' rm(eta, simy_pa, simy_ztp, X, mm, spp_loadings_pa, spp_loadings_ztp, true_lvs_pa, true_lvs_ztp, 
#' xy, simy, dat, ztpR)
#' 
#' 
#' 
#' # Set up spatial basis functions for CBFM -- Most users will start here! 
#' # We will use this basis functions in both CBFM fits
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
#' # Fit presence-absence CBFM 
#' simy_pa_train <- (simy_train > 0)*1 
#' useformula <- ~ temp + depth + chla + O2
#' fitcbfm_pa <- CBFM(y = simy_pa_train, formula = useformula, data = dat_train, 
#' B_space = train_basisfunctions, family = binomial(), control = list(trace = 1))
#' rm(simy_pa_train)
#' 
#' # Now fit zero-truncated Poisson CBFM
#' # This could take a while...
#' fitcbfm_ztp <- CBFM(y = simy_train, formula = useformula, data = dat_train, 
#' B_space = train_basisfunctions, family = ztpoisson(), control = list(trace = 1))
#' 
#' 
#' # Finally, form the hurdle CBFM, and remove the component CBFMs to save space!
#' fitcbfm_hurdle <- makeahurdle(pa_object = fitcbfm_pa, count_object = fitcbfm_ztp)
#' rm(fitcbfm_pa, fitcbfm_ztp)
#' 
#' 
#' fitcbfm_hurdle
#' 
#' fitted(fitcbfm_hurdle)
#' 
#' residuals(fitcbfm_hurdle, type = "dunnsmyth")
#' 
#' plot(fitcbfm_hurdle)
#' 
#' simulate(fitcbfm_hurdle)
#' 
#' logLik(fitcbfm_hurdle, use_edf = TRUE) 
#' ## Notice the (effective) degrees of freedom are quite large!
#' 
#' AIC(fitcbfm_hurdle, use_edf = TRUE) 
#' AICc(fitcbfm_hurdle, use_edf = TRUE)
#' 
#' 
#' # Examples of in-sample prediction first
#' predict(fitcbfm_hurdle) # One can confirm this is (basically) the same as fitted(fitcbfm_hurdle)
#' 
#' predict(fitcbfm_hurdle, se_fit = TRUE) 
#' # Can take a while, as it uses simulation...
#' # Also, the results will differ from predict(fitcbfm_hurdle) due to simulation error
#' 
#' 
#' # Now compare prediction on test set
#' library(pscl)
#' library(mvabund)
#' 
#' # Fit stacked hurdle Poisson regression models as a baseline
#' fitstacked <- NULL 
#' for(j in 1:num_spp) {
#' fitstacked[[j]] <- hurdle(resp ~ temp + depth + chla + O2, 
#' data = data.frame(resp = simy_train[,j], dat_train), dist = "poisson")
#' }
#' 
#' # Fit stacked NB GLM as another baseline
#' fitstacked2 <- manyglm(simy_train ~ temp + depth + chla + O2, family = "negative.binomial",
#' data = dat_train)
#' 
#' # Fit a negative binomial CBFM as another baseline
#' useformula <- ~ temp + depth + chla + O2
#' fitcbfm_nb <- CBFM(y = simy_train, formula = useformula, data = dat_train, 
#' B_space = train_basisfunctions, family = nb2(), control = list(trace = 1))
#' 
#' 
#' # Calculate predictions onto test dataset
#' predictions_stacked <- sapply(1:num_spp, function(j) predict(fitstacked[[j]], 
#' newdata = dat_test, type = "response"))
#' predictions_stacked2 <- predict(fitstacked2, newdata = dat_test, type = "response")
#' predictions_cbfm_nb <- predict(fitcbfm_nb, newdata = dat_test, type = "response", 
#' new_B_space = test_basisfunctions)
#' predictions_cbfm_hurdle <- predict(fitcbfm_hurdle, newdata_pa = dat_test,  
#' new_B_space_pa = test_basisfunctions, newdata_count = dat_test,
#' new_B_space_count = test_basisfunctions)
#' 
#' # Evaluation predictions
#' # Pseudo R-squared across species
#' pseudoR2 <- data.frame(
#' stacked_hurdle = sapply(1:num_spp, function(j) { 
#' out <- cor(predictions_stacked[,j], simy_test[,j], method = "spearman")
#' out^2 * sign(out)     
#' }),
#' stacked_nb = sapply(1:num_spp, function(j) { 
#' out <- cor(predictions_stacked2[,j], simy_test[,j], method = "spearman")
#' out^2 * sign(out)     
#' }),
#' cbfm_nb = sapply(1:num_spp, function(j) { 
#' out <- cor(predictions_cbfm_nb[,j], simy_test[,j], method = "spearman")
#' out^2 * sign(out)     
#' }),
#' cbfm_hurdle = sapply(1:num_spp, function(j) { 
#' out <- cor(predictions_cbfm_hurdle[,j], simy_test[,j], method = "spearman")
#' out^2 * sign(out)     
#' })
#' )
#' 
#' boxplot(pseudoR2, main = "Pseudo-R2", 
#' names = c("Stacked hurdle Poisson", "Stacked NB", "NB CBFM", "Hurdle Poisson CBFM"))
#' 
#' 
#' 
#' ##------------------------------
#' ## **Example 2: Fitting a CBFM to spatial multivariate hurdle negative binomial data**
#' ## simulated from a spatial hurdle latent variable model
#' ## Please note the data generation process (thus) differs from CBFM.
#' ##------------------------------
#' set.seed(2021)
#' num_sites <- 1000 # 500 (units) sites for training set + 500 sites for testing.
#' num_spp <- 50 # Number of species
#' num_X <- 4 # Number of regression slopes
#' 
#' # Generate and combine latent variables model for presence-absence and zero truncated Poisson X
#' spp_slopes_pa <- matrix(runif(num_spp * num_X, -1, 1), nrow = num_spp)
#' spp_slopes_ztnb <- matrix(runif(num_spp * num_X, -1, 1), nrow = num_spp)
#' spp_intercepts_pa <- runif(num_spp, -2, 0)
#' spp_intercepts_ztnb <- runif(num_spp, -2, 0)
#' spp_dispersion_ztnb <- runif(num_spp, 0, 5)
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
#' # Simulate latent variable components
#' true_lvs_pa <- grf(grid = cbind(xy$x, xy$y), nsim = 2, cov.model = "exponential",
#' cov.pars = c(1, 2))$data %>%
#'      as.matrix
#' spp_loadings_pa <- matrix(runif(num_spp * 2, -1, 1), nrow = num_spp)
#' true_lvs_ztnb <- grf(grid = cbind(xy$x, xy$y), nsim = 2, cov.model = "exponential",
#' cov.pars = c(1, 2.5))$data %>%
#'      as.matrix
#' spp_loadings_ztnb <- matrix(runif(num_spp * 2, -1, 1), nrow = num_spp) 
#' set.seed(NULL)
#' 
#' # Simulate spatial multivariate abundance data
#' eta <- tcrossprod(cbind(1,mm), cbind(spp_intercepts_pa,spp_slopes_pa)) + 
#' tcrossprod(true_lvs_pa, spp_loadings_pa)
#' simy_pa <- matrix(rbinom(num_sites * num_spp, size = 1, 
#' prob = plogis(eta)), nrow = num_sites)
#' 
#' # Now simulate spatial count data from a truncated NB distribution
#' eta <- tcrossprod(cbind(1,mm), cbind(spp_intercepts_ztnb,spp_slopes_ztnb)) + 
#' tcrossprod(true_lvs_ztnb, spp_loadings_ztnb)
#' ztNBR <- trun.r(par = 0, family = "NBI", type = "left") 
#' simy_ztnb <- matrix(ztNBR(num_sites * num_spp, mu = exp(eta)), nrow = num_sites)
#' 
#' # Spatial multivariate count data from a hurdle NB model is then the product of the two
#' simy <- simy_pa *  simy_ztnb
#' 
#' # Form training and test sets
#' dat_train <- dat[1:500,]
#' dat_test <- dat[501:1000,]
#' simy_train <- simy[1:500,]
#' simy_test <- simy[501:1000,]
#' 
#' # Delete the "component" responses and present you only observe the final response
#' rm(eta, simy_pa, simy_ztnb, X, mm, spp_loadings_pa, spp_loadings_ztnb, true_lvs_pa, true_lvs_ztnb, 
#' xy, simy, dat, ztNBR)
#' 
#' 
#' 
#' # Set up spatial basis functions for CBFM -- Most users will start here! 
#' # We will use this basis functions in both CBFM fits
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
#' # Fit presence-absence CBFM 
#' simy_pa_train <- (simy_train > 0)*1 
#' useformula <- ~ temp + depth + chla + O2
#' fitcbfm_pa <- CBFM(y = simy_pa_train, formula = useformula, data = dat_train, 
#' B_space = train_basisfunctions, family = binomial(), control = list(trace = 1))
#' rm(simy_pa_train)
#' 
#' # Now fit zero-truncated negative binomial CBFM
#' # This could take a while...
#' fitcbfm_ztnb <- CBFM(y = simy_train, formula = useformula, data = dat_train, 
#' B_space = train_basisfunctions, family = ztnb2(), control = list(trace = 1))
#' 
#' 
#' # Finally, form the hurdle CBFM, and remove the component CBFMs to save space!
#' fitcbfm_hurdle <- makeahurdle(pa_object = fitcbfm_pa, count_object = fitcbfm_ztnb)
#' rm(fitcbfm_pa, fitcbfm_ztnb)
#' 
#' 
#' fitcbfm_hurdle
#' 
#' fitted(fitcbfm_hurdle)
#' 
#' residuals(fitcbfm_hurdle, type = "dunnsmyth")
#' 
#' plot(fitcbfm_hurdle)
#' 
#' simulate(fitcbfm_hurdle)
#' 
#' logLik(fitcbfm_hurdle, use_edf = TRUE) 
#' ## Notice the (effective) degrees of freedom are quite large!
#' 
#' AIC(fitcbfm_hurdle, use_edf = TRUE) 
#' AICc(fitcbfm_hurdle, use_edf = TRUE)
#' 
#' 
#' # Examples of in-sample prediction first
#' predict(fitcbfm_hurdle) # One can confirm this is (basically) the same as fitted(fitcbfm_hurdle)
#' 
#' predict(fitcbfm_hurdle, se_fit = TRUE) 
#' # Can take a while, as it uses simulation...
#' # Also, the results will differ from predict(fitcbfm_hurdle) due to simulation error
#' 
#' 
#' # Now compare prediction on test set
#' library(pscl)
#' library(mvabund)
#' 
#' # Fit stacked hurdle negative binomial regression models as a baseline
#' fitstacked <- NULL 
#' for(j in 1:num_spp) {
#' fitstacked[[j]] <- hurdle(resp ~ temp + depth + chla + O2, 
#' data = data.frame(resp = simy_train[,j], dat_train), dist = "negbin")
#' }
#' 
#' # Fit stacked negative binomial GLM as another baseline
#' fitstacked2 <- manyglm(simy_train ~ temp + depth + chla + O2, family = "negative.binomial",
#' data = dat_train)
#' 
#' # Fit a negative binomial CBFM as another baseline
#' useformula <- ~ temp + depth + chla + O2
#' fitcbfm_nb <- CBFM(y = simy_train, formula = useformula, data = dat_train, 
#' B_space = train_basisfunctions, family = nb2(), control = list(trace = 1))
#' 
#' 
#' # Calculate predictions onto test dataset
#' predictions_stacked <- sapply(1:num_spp, function(j) predict(fitstacked[[j]], 
#' newdata = dat_test, type = "response"))
#' predictions_stacked2 <- predict(fitstacked2, newdata = dat_test, type = "response")
#' predictions_cbfm_nb <- predict(fitcbfm_nb, newdata = dat_test, type = "response", 
#' new_B_space = test_basisfunctions)
#' predictions_cbfm_hurdle <- predict(fitcbfm_hurdle, newdata_pa = dat_test,  
#' new_B_space_pa = test_basisfunctions, newdata_count = dat_test,
#' new_B_space_count = test_basisfunctions)
#' 
#' # Evaluation predictions
#' # Pseudo R-squared across species
#' pseudoR2 <- data.frame(
#' stacked_hurdle = sapply(1:num_spp, function(j) { 
#' out <- cor(predictions_stacked[,j], simy_test[,j], method = "spearman")
#' out^2 * sign(out)     
#' }),
#' stacked_nb = sapply(1:num_spp, function(j) { 
#' out <- cor(predictions_stacked2[,j], simy_test[,j], method = "spearman")
#' out^2 * sign(out)     
#' }),
#' cbfm_nb = sapply(1:num_spp, function(j) { 
#' out <- cor(predictions_cbfm_nb[,j], simy_test[,j], method = "spearman")
#' out^2 * sign(out)     
#' }),
#' cbfm_hurdle = sapply(1:num_spp, function(j) { 
#' out <- cor(predictions_cbfm_hurdle[,j], simy_test[,j], method = "spearman")
#' out^2 * sign(out)     
#' })
#' )
#' 
#' boxplot(pseudoR2, main = "Pseudo-R2", 
#' names = c("Stacked hurdle Poisson", "Stacked NB", "NB CBFM", "Hurdle NB CBFM"))
#' }
#' 
#' 
#' @export 
#' @md

makeahurdle <- function(pa_object, count_object, ...) {
     if(!inherits(pa_object, "CBFM")) 
        stop("`pa_object' is not of class \"CBFM\"")
     if(!inherits(count_object, "CBFM")) 
        stop("`count_object' is not of class \"CBFM\"")
     
     if(pa_object$family$family[1] != "binomial" | pa_object$trial_size != 1)
          stop("pa_object must corresponding to a binary CBFM fitted to presence absence data.")
     if(!(count_object$family$family[1] %in% c("ztpoisson", "ztnegative.binomial")))
          stop("count_object must corresponding to a zero-truncated CBFM.")

     
     out <- list(pa_fit = pa_object, count_fit = count_object)
     class(out) <- "CBFM_hurdle"
     message("Hurdle CBFM formed.")
     return(out)
     }



# Adapted from and acknowledgement goes to the authors of gamlss.dist
# Note zeroprob is the probability of a zero
.phurdlepoisson <- function(q, lambda = 5, zeroprob, lower.tail = TRUE, log.p = FALSE) {
        if(any(lambda <= 0)) 
                stop(paste("lambda must be greater than 0", "\n", ""))
        if(any(zeroprob <= 0) | any(zeroprob >= 1)) 
                stop(paste("zeroprob must be between 0 and 1", "\n", ""))
        if(any(q < 0)) 
                stop(paste("y must be 0 or greater than 0", "\n", ""))
    
        ly <- max(length(q), length(lambda), length(zeroprob))
        q <- rep(q, length = ly)
        lambda <- rep(lambda, length = ly)
        zeroprob <- rep(zeroprob, length = ly)
        
        cdf <- rep(0, ly)
        cdf1 <- ppois(q, lambda = lambda, lower.tail = TRUE, log.p = FALSE)
        cdf2 <- ppois(0, lambda = lambda, lower.tail = TRUE, log.p = FALSE)
        cdf3 <- zeroprob + ((1 - zeroprob) * (cdf1 - cdf2)/(1 - cdf2))
        cdf <- ifelse((q == 0), zeroprob, cdf3)
        
        if(lower.tail == TRUE) 
                cdf <- cdf
        else 
                cdf <- 1 - cdf
        
        if(log.p == FALSE) 
                cdf <- cdf
        else
                cdf <- log(cdf)
        
        return(cdf)
        }




# Adapted from and acknowledgement goes to the authors of gamlss.dist
# Note zeroprob is the probability of a zero
.phurdlenb2 <- function(q, mu, phi, zeroprob, lower.tail = TRUE, log.p = FALSE) {
        if(any(mu <= 0)) 
                stop(paste("mu must be greater than 0 ", "\n", ""))
        if(any(phi <= 0)) 
                stop(paste("phi must be greater than 0 ", "\n", ""))
        if(any(zeroprob <= 0) | any(zeroprob >= 1)) 
                stop(paste("zeroprob must be between 0 and 1 ", "\n", ""))
        if(any(q < 0))
                stop(paste("y must be >=0", "\n", ""))
    
        
        ly <- max(length(q), length(mu), length(zeroprob), length(phi))
        q <- rep(q, length = ly)
        phi <- rep(phi, length = ly)
        mu <- rep(mu, length = ly)
        zeroprob <- rep(zeroprob, length = ly)
    
        cdf0 <- pnbinom(0, mu = mu, size = 1/phi)
        cdf1 <- pnbinom(q, mu = mu, size = 1/phi)
        cdf3 <- zeroprob + ((1 - zeroprob) * (cdf1 - cdf0)/(1 - cdf0))
        cdf <- ifelse((q == 0), zeroprob, cdf3)
        
        if(lower.tail == TRUE) 
                cdf <- cdf
        else 
                cdf <- 1 - cdf
        if(log.p == FALSE) 
                cdf <- cdf
        else 
                cdf <- log(cdf)
        return(cdf)
        }
