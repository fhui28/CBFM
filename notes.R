## See github issues for outstanding things to do!!!

## GENERAL APPROACH TO MODIFYING PACKAGE
## 1. Do all changes
## 1. Consider moving calls of other packages from :: to using importFrom
## 2. Delete all Rd files in man (but not figures!), the NAMESPACE file, and any compiled C++ files
## 3. Check package (Ctrl+Shift+E)
##4. Build (optional)
## 5. Push to github

#Removed from DESCRIPTION since we now turn off compilation of the TMB C++ files [learning from how James Thorson does it with VAST!]
LinkingTo: 
    TMB,
    RcppEigen


##--------------------------------------
## Random testing jazz
##-------------------------------------
library(autoFRK)
library(FRK)
library(MASS)
library(mvabund)
library(mvtnorm)
library(ROCR)
library(sp)
library(RandomFields)
library(tidyverse)

##------------------------------
## **Example 1a: Fitting a CBFM to spatial multivariate presence-absence data**
## simulated from a spatial latent variable model
## Please note the data generation process (thus) differs from CBFM.
##------------------------------
set.seed(2021)
num_sites <- 1000 # 500 (units) sites for training set + 500 sites for testing.
num_spp <- 50 # Number of species
num_X <- 4 # Number of regression slopes

#spp_slopes <- matrix(runif(num_spp * num_X, -1, 1), nrow = num_spp)
spp_slopes <- cbind(rnorm(num_spp, -0.5, sd = 0.25), rnorm(num_spp, 0.5, sd = 0.25), rnorm(num_spp, -0.25, sd = 0.1), rnorm(num_spp, 0.25, sd = 0.1))
spp_intercepts <- runif(num_spp, -2, 0)

# Simulate spatial coordinates and environmental covariate components
# We will use this information in later examples as well
xy <- data.frame(x = runif(num_sites, 0, 5), y = runif(num_sites, 0, 5))
X <- rmvnorm(num_sites, mean = rep(0,4))
colnames(X) <- c("temp", "depth", "chla", "O2")
dat <- data.frame(xy, X)
mm <- model.matrix(~ temp + depth + chla + O2 - 1, data = dat) %>%
scale %>%
as.matrix

# Simulate latent variable component
# We will use this information in later examples as well
true_lvs <- RFsimulate(model = RMexp(var=1, scale=2),
x = xy$x, y = xy$y, n = 2)@data %>%
as.matrix
spp_loadings <- matrix(runif(num_spp * 2, -1, 1), nrow = num_spp)
set.seed(NULL)

# Simulate spatial multivariate abundance data (presence-absence)
# We will use this information in later examples as well
eta <- tcrossprod(cbind(1,mm), cbind(spp_intercepts,spp_slopes)) +
tcrossprod(true_lvs, spp_loadings)
simy <- matrix(rbinom(num_sites * num_spp, size = 1,
prob = plogis(eta)), nrow = num_sites)

# Form training and test sets
dat_train <- dat[1:500,]
dat_test <- dat[501:1000,]
simy_train <- simy[1:500,]
simy_test <- simy[501:1000,]
rm(X, mm, spp_loadings, true_lvs, xy, simy, dat)


# Fit stacked GLM as a baseline
fitstacked <- manyglm(simy_train ~ temp + depth + chla + O2, family = binomial(), data = dat_train)


# Set up spatial basis functions for CBFM -- Most users will start here!
# We will also use this basis functions in some later examples
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



# Fit CBFM
mm <- model.matrix(~ temp + depth + chla + O2, data = dat_train)[,-1,drop=FALSE]
useformula <- ~ 1
fitcbfm <- CBFM(y = simy_train, formula_X = useformula, data = dat_train,
B_space = train_basisfunctions, B_time = mm,
family = binomial(), 
control = list(trace = 1, nonzeromean_B_time = TRUE),
Sigma_control = list(rank = c(5,1))
)


fitcbfm_pure <- CBFM(y = simy_train, formula_X = ~ temp + depth + chla + O2, 
                     data = dat_train,
                     B_space = train_basisfunctions, 
                     family = binomial(), control = list(trace = 1))



ggmatplot(spp_slopes, fitcbfm_pure$betas[,-1]) + geom_abline(intercept = 0, slope = 1)
ggmatplot(spp_slopes, fitcbfm$basis_effects_mat[,fitcbfm$num_B_space+(1:fitcbfm$num_B_time)]) + geom_abline(intercept = 0, slope = 1)




y = simy_train
useformula <- ~ 1
mm <- model.matrix(~ temp + depth + chla + O2, data = dat_train)[,-1,drop=FALSE]
formula_X = useformula
data = dat_train
family =  binomial()
B_space = train_basisfunctions
B_time = mm
B_spacetime = NULL
offset = NULL
ncores = NULL
gamma = 1
trial_size = 1
dofit = TRUE
stderrors = TRUE
select = FALSE
start_params = list(betas = NULL, basis_effects_mat = NULL, dispparam = NULL, powerparam = NULL, zeroinfl_prob = NULL)
TMB_directories = list(cpp = system.file("executables", package = "CBFM"), compile = system.file("executables", package = "CBFM"))
control = list(maxit = 100, convergence_type = "parameters", tol = 1e-4, seed = NULL, trace = 1, ridge = 0, nonzeromean_B_time = TRUE)
Sigma_control = list(rank = c(5,1), maxit = 100, tol = 1e-4, method = "LA", trace = 0)
G_control = list(rank = c(5,5), nugget_profile = seq(0.05, 0.95, by = 0.05), maxit = 100, tol = 1e-4, method = "LA", trace = 0)
k_check_control = list(subsample = 5000, n.rep = 400)



Ginv = new_LoadingnuggetG_space$covinv
basis_effects_mat = new_fit_CBFM_ptest$basis_effects_mat[,1:num_spacebasisfns,drop=FALSE]+G_control$tol
Sigmainv = new_LoadingnuggetSigma_space$covinv
B = B_space
y_vec = as.vector(y)
linpred_vec = c(new_fit_CBFM_ptest$linear_predictor)
dispparam = new_fit_CBFM_ptest$dispparam
powerparam = new_fit_CBFM_ptest$powerparam
zeroinfl_prob_intercept = new_fit_CBFM_ptest$zeroinfl_prob_intercept
return_correlation = TRUE


#-----------------------------
#-----------------------------
library(autoFRK)
library(FRK)
library(MASS)
library(mvtnorm)
library(sp)
library(RandomFields)
library(tidyverse)
library(gamlss.tr)
library(gamlss.add)
 
set.seed(2021)
num_sites <- 500 # 500 (units) sites 
num_spp <- 50 # Number of species
num_X <- 4 # Number of regression slopes
 
# Simulate spatial coordinates and environmental covariate components
xy <- data.frame(x = runif(num_sites, 0, 5), y = runif(num_sites, 0, 5))
X <- rmvnorm(num_sites, mean = rep(0,4)) 
colnames(X) <- c("temp", "depth", "chla", "O2")
dat <- data.frame(xy, X)
useformula <- ~ temp + depth + chla + O2
 
# Set up spatial basis functions for CBFM 
num_basisfunctions <- 25 # Number of spatial basis functions to use
basisfunctions <- mrts(dat[,c("x","y")], num_basisfunctions) %>% 
as.matrix %>%
{.[,-(1)]} # Remove the first intercept column


spp_slopes_ztnb <- matrix(runif(num_spp * num_X, -0.5, 0.5), nrow = num_spp)
#spp_basis_effects_mat_ztnb <- matrix(runif(num_spp * (num_basisfunctions-1), -0.5, 0.5), nrow = num_spp)
#spp_basis_effects_mat_ztnb[,15:19] <- 0
spp_intercepts_ztnb <- runif(num_spp, -2, 0)
 
true_Sigma_space_ztnb <- rWishart(1, num_basisfunctions+1, diag(x = 0.1, nrow = num_basisfunctions-1))[,,1]/10
true_G_space_ztnb <- rWishart(1, num_spp+1, diag(x = 0.1, nrow = num_spp))[,,1] %>% cov2cor
spp_dispersion <- runif(num_spp, 0, 5)

 
simy_count <- create_CBFM_life(family = nb2(), formula_X = useformula, data = dat,
                              B_space = basisfunctions, betas = cbind(spp_intercepts_ztnb, spp_slopes_ztnb),
                              G = list(space = true_G_space_ztnb), Sigma = list(space = true_Sigma_space_ztnb),
                              dispparam = spp_dispersion,  max_resp = 20000) 



# Now fit a zero-truncated count distribution, and CBFM has to ignore the zeros in the data
y = simy_count$y
useformula <- ~ temp + depth + chla + O2
formula_X = useformula
data = dat
B_space = basisfunctions
family =  ztnb2()
B_time = NULL
B_spacetime = NULL
offset = NULL
ncores = NULL
gamma = 1
trial_size = 1
dofit = TRUE
stderrors = TRUE
select = FALSE
start_params = list(betas = NULL, basis_effects_mat = NULL, dispparam = NULL, powerparam = NULL, zeroinfl_prob = NULL)
TMB_directories = list(cpp = system.file("executables", package = "CBFM"), compile = system.file("executables", package = "CBFM"))
control = list(maxit = 100, convergence_type = "parameters", tol = 1e-4, initial_betas_dampen = 0.05, seed = NULL, trace = 1, ridge = 0) 
Sigma_control = list(rank = 5, maxit = 100, tol = 1e-4, method = "LA", trace = 0)
G_control = list(rank = 5, nugget_profile = seq(0.05, 0.95, by = 0.05), maxit = 100, tol = 1e-4, method = "LA", trace = 0)
k_check_control = list(subsample = 5000, n.rep = 400)




##---------------------
