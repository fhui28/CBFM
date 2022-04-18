## See github issues for outstanding things to do!!!

## GENERAL APPROACH TO MODIFYING PACKAGE
## 1. Do all changes
## 1. Consider moving calls of other packages from :: to using importFrom
## 2. Delete all Rd files in man (but not figures!), the NAMESPACE file, and any compiled C++ files
## 3. Check package (Ctrl+Shift+E)
## 4. Build (optional)
## 5. Push to github

#Removed from DESCRIPTION since we now turn off compilation of the TMB C++ files [learning from how James Thorson does it with VAST!]
LinkingTo: 
    TMB,
    RcppEigen


##--------------------------------------
## Random testing jazz
##-------------------------------------

y = simy_train
useformula <- ~ 1
formula_X = useformula
data = dat_train
B_space = train_basisfunctions
family =  binomial()
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
control = list(maxit = 100, convergence_type = "parameters", tol = 1e-4, seed = NULL, trace = 1, ridge = 0)
Sigma_control = list(rank = c(2,1), maxit = 100, tol = 1e-4, method = "LA", trace = 0)
G_control = list(rank = 5, nugget_profile = seq(0.05, 0.95, by = 0.05), maxit = 100, tol = 1e-4, method = "LA", trace = 0)
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

set.seed(2022)
num_sites <- 500 # 500 (units) sites 
num_spp <- 50 # Number of species
num_X <- 4 # Number of regression slopes
 
spp_slopes <- matrix(runif(num_spp * num_X, -1, 1), nrow = num_spp)
spp_intercepts <- runif(num_spp, -2, 0)
spp_gear <- rnorm(num_spp, mean = 1.5, sd = 0.2)

# Simulate spatial coordinates and environmental covariate components
xy <- data.frame(x = runif(num_sites, 0, 5), y = runif(num_sites, 0, 5))
X <- rmvnorm(num_sites, mean = rep(0,4))
colnames(X) <- c("temp", "depth", "chla", "O2")
dat <- data.frame(xy, X)
mm <- model.matrix(~ temp + depth + chla + O2 - 1, data = dat) %>%
scale %>%
as.matrix
X_year <- matrix(rep(c(0,1), c(400,100)), ncol = 1)

 
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
    tcrossprod(X_year, spp_gear) +
    tcrossprod(true_lvs, spp_loadings)
simy <- matrix(rbinom(num_sites * num_spp, size = 1, prob = plogis(eta)), nrow = num_sites)

# Set up spatial basis functions for CBFM 
num_basisfunctions <- 25 # Number of spatial basis functions to use
basisfunctions <- mrts(dat[,c("x","y")], num_basisfunctions) %>% 
as.matrix %>%
{.[,-(1)]} # Remove the first intercept column



# Fit CBFM
useformula <- ~ temp + depth + chla + O2
fitcbfm <- CBFM(y = simy, formula_X = useformula, data = dat,
B_space = basisfunctions, B_time = X_year, 
family = binomial(), control = list(trace = 1),
Sigma_control = list(rank = c(5,1)),
G_control = list(rank = c(5,5))
)

 
y = simy
useformula <- ~ temp + depth + chla + O2
formula_X = useformula
data = dat
family =  binomial()
B_space = basisfunctions
B_time = X_year
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
Sigma_control = list(rank = c(5,1), maxit = 100, tol = 1e-4, method = "LA", trace = 0)
G_control = list(rank = c(5,5), nugget_profile = seq(0.05, 0.95, by = 0.05), maxit = 100, tol = 1e-4, method = "LA", trace = 0)
k_check_control = list(subsample = 5000, n.rep = 400)




library(ggmatplot)
ggmatplot(spp_slopes, fitcbfm$betas[,-1]) + geom_abline(intercept = 0, slope = 1)
qplot(spp_intercepts, fitcbfm$betas[,1]) + geom_abline(intercept = 0, slope = 1)
qplot(spp_gear, fitcbfm$basis_effects_mat[,25]) + geom_abline(intercept = 0, slope = 1)
fitcbfm$basis_effects_mat[,25] %>% summary
fitcbfm$basis_effects_mat[,25] %>% sd
fitcbfm$G_time
fitcbfm$Sigma_time

