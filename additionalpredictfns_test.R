#' ---
#' title: Testing out some additional prediction functions purely for Chris Haak's application
#' author: FKCH
#' date: Code started Dec 2022
#' ---
#' 
#' 
##------------------------------
# Setup
##------------------------------
rm(list = ls())
library(tidyverse)
library(autoFRK)
library(FRK)
library(MASS)
library(mvabund)
library(mvtnorm)
library(ROCR)
library(sp)
library(RandomFields)
library(gamlss)
source("additionalpredictfns.R")
library(CBFM)

##------------------------------
## Pulled from help file for makeahurdle function 
## **Example 2: Fitting a CBFM to spatial multivariate hurdle negative binomial data**
## simulated from a spatial hurdle latent variable model
## Please note the data generation process (thus) differs from CBFM.
##------------------------------
set.seed(2021)
num_sites <- 1000 # 500 (units) sites for training set + 500 sites for testing.
num_spp <- 50 # Number of species
num_X <- 4 # Number of regression slopes

# Generate and combine latent variables model for presence-absence and zero truncated Poisson X
spp_slopes_pa <- matrix(runif(num_spp * num_X, -1, 1), nrow = num_spp)
spp_slopes_ztnb <- matrix(runif(num_spp * num_X, -1, 1), nrow = num_spp)
spp_intercepts_pa <- runif(num_spp, -2, 0)
spp_intercepts_ztnb <- runif(num_spp, -2, 0)
spp_dispersion_ztnb <- runif(num_spp, 0, 5)

# Simulate spatial coordinates and environmental covariate components
xy <- data.frame(x = runif(num_sites, 0, 5), y = runif(num_sites, 0, 5))
X <- rmvnorm(num_sites, mean = rep(0,4))
colnames(X) <- c("temp", "depth", "chla", "O2")
dat <- data.frame(xy, X)
mm <- model.matrix(~ temp + depth + chla + O2 - 1, data = dat) %>%
scale %>%
as.matrix

# Simulate latent variable components
true_lvs_pa <- RFsimulate(model = RMexp(var=1, scale=2),
x = xy$x, y = xy$y, n = 2)@data %>%
as.matrix
spp_loadings_pa <- matrix(runif(num_spp * 2, -1, 1), nrow = num_spp)
true_lvs_ztnb <- RFsimulate(model = RMexp(var=1, scale=2.5),
x = xy$x, y = xy$y, n = 2)@data %>%
as.matrix
spp_loadings_ztnb <- matrix(runif(num_spp * 2, -1, 1), nrow = num_spp)
set.seed(NULL)

# Simulate spatial multivariate abundance data
eta <- tcrossprod(cbind(1,mm), cbind(spp_intercepts_pa,spp_slopes_pa)) +
tcrossprod(true_lvs_pa, spp_loadings_pa)
simy_pa <- matrix(rbinom(num_sites * num_spp, size = 1,
prob = plogis(eta)), nrow = num_sites)

# Now simulate spatial count data from a truncated NB distribution
eta <- tcrossprod(cbind(1,mm), cbind(spp_intercepts_ztnb,spp_slopes_ztnb)) +
tcrossprod(true_lvs_ztnb, spp_loadings_ztnb)
ztNBR <- trun.r(par = 0, family = "NBI", type = "left")
simy_ztnb <- matrix(ztNBR(num_sites * num_spp, mu = exp(eta)), nrow = num_sites)

# Spatial multivariate count data from a hurdle NB model is then the product of the two
simy <- simy_pa *  simy_ztnb

# Form training and test sets
dat_train <- dat[1:500,]
dat_test <- dat[501:1000,]
simy_train <- simy[1:500,]
simy_test <- simy[501:1000,]

# Delete the "component" responses and present you only observe the final response
rm(eta, simy_pa, simy_ztnb, X, mm, spp_loadings_pa, spp_loadings_ztnb, true_lvs_pa, true_lvs_ztnb,
xy, simy, dat, ztNBR)



# Set up spatial basis functions for CBFM -- Most users will start here!
# We will use this basis functions in both CBFM fits
num_basisfunctions <- 25 # Number of spatial basis functions to use
# Training set basis functions
train_basisfunctions <- mrts(dat_train[,c("x","y")], num_basisfunctions) %>%
as.matrix %>%
{.[,-(1)]} # Remove the first intercept column
# Testing set basis functions
test_basisfunctions <- mrts(dat_train[,c("x","y")], num_basisfunctions) %>%
predict(newx = dat_test[,c("x","y")]) %>%
as.matrix %>%
{.[,-c(1)]}

# Fit presence-absence CBFM
simy_pa_train <- (simy_train > 0)*1
useformula <- ~ temp + depth + chla + O2
fitcbfm_pa <- CBFM(y = simy_pa_train, formula = useformula, data = dat_train,
B_space = train_basisfunctions, family = binomial(), control = list(trace = 1))
rm(simy_pa_train)


# Now fit zero-truncated negative binomial CBFM
# This could take a while...
fitcbfm_ztnb <- CBFM(y = simy_train, formula = useformula, data = dat_train,
                     B_space = train_basisfunctions, family = ztnb2(), control = list(trace = 1))


# Finally, form the hurdle CBFM, and remove the component CBFMs to save space!
fitcbfm_hurdle <- makeahurdle(pa_object = fitcbfm_pa, count_object = fitcbfm_ztnb)
rm(fitcbfm_pa, fitcbfm_ztnb)

fitcbfm_hurdle



##------------------------------
# Comparing different forms of predictions
##------------------------------
# oospred_cbfm_timespacetime_hurdle_mean <- predict(fitcbfm_hurdle, newdata_pa = dat_test, newdata_count = dat_test,
#                                    new_B_space_pa = test_basisfunctions, new_B_space_count = test_basisfunctions,
#                                    type = "response")
oospred_cbfm_timespacetime_hurdle_mean <- predict_hurdle_CBFM(fitcbfm_hurdle, prediction_form = "mean",
                                                                newdata_pa = dat_test, new_B_space_pa = test_basisfunctions,
                                                                newdata_count = dat_test, new_B_space_count = test_basisfunctions,
                                                                parameter_uncertainty = FALSE, num_sims = 10, ncores = 10)

oospred_cbfm_timespacetime_hurdle_median <- predict_hurdle_CBFM(fitcbfm_hurdle, prediction_form = "median",
                                                                newdata_pa = dat_test, new_B_space_pa = test_basisfunctions,
                                                                newdata_count = dat_test, new_B_space_count = test_basisfunctions,
                                                                parameter_uncertainty = FALSE, num_sims = 10, ncores = 10)

oospred_cbfm_timespacetime_hurdle_threshold <- predict_hurdle_CBFM(fitcbfm_hurdle, prediction_form = "threshold", quantile_threshold = 0.9,
                                                                newdata_pa = dat_test, new_B_space_pa = test_basisfunctions,
                                                                newdata_count = dat_test, new_B_space_count = test_basisfunctions,
                                                                parameter_uncertainty = FALSE, num_sims = 10, ncores = 10)

oospred_cbfm_timespacetime_hurdle_linex <- predict_hurdle_CBFM(fitcbfm_hurdle, prediction_form = "linex", 
                                                                   newdata_pa = dat_test, new_B_space_pa = test_basisfunctions,
                                                                   newdata_count = dat_test, new_B_space_count = test_basisfunctions,
                                                                   parameter_uncertainty = FALSE, num_sims = 10, ncores = 10)



##-----------------------
sessionInfo()
##-----------------------
