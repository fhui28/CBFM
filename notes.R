## See github issues for outstanding things to do -- Sept 18 2021

## GENERAL APPROACH TO MODIFYING PACKAGE
## 1. Do all changes
## 1. Consider moving calls of other packages from :: to using importFrom
## 2. Delete all Rd files in man (but not figures!), the NAMESPACE file, and any compiled C++ files
## 3. Check package (Ctrl+Shift+E)
## 4. Build (optional)
## 5. Push to github

# Could use remotes to do installaton for countreg, but then likely can not officially put it up on CRAN
Remotes: svn:://svn.r-forge.r-project.org/svnroot/countreg/

#Removed from DESCRIPTION since we now turn off compilation of the TMB C++ files [learning from how James Thorson does it with VAST!]
LinkingTo: 
    TMB,
    RcppEigen


##--------------------------------------
## Random testing jazz
##-------------------------------------

# 
# dat2 <- data.frame(response = y[,2], dat)
# fit_trun <- gamlss(formula = response ~ temp + depth + chla + O2, 
#                    data = dat2, 
#                    family = NBI,
#                    #family = trun(0, family = "NBI"),
#                    control = gamlss.control(n.cyc = 30))
# 
# 
# fit_trun2 <- gamlss(formula = response ~ temp + depth + chla + O2, 
#                    data = dat2, 
#                    family = trun(0, family = "NBI"),
#                    control = gamlss.control(n.cyc = 30))
# 
# 
# 


y = simy_train
useformula <- ~ temp + depth + chla + O2
formula_X = useformula
data = dat_train
B_space = train_basisfunctions
family =  nb2()
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
control = list(maxit = 100, convergence_type = "parameters", tol = 1e-4, seed = NULL, trace = 1, ridge = 0)
Sigma_control = list(rank = 5, maxit = 100, tol = 1e-4, method = "LA", trace = 0)
G_control = list(rank = 5, nugget_profile = seq(0.05, 0.95, by = 0.05), maxit = 100, tol = 1e-4, method = "LA", trace = 0)
k_check_control = list(subsample = 5000, n.rep = 400)



Ginv = new_LoadingnuggetG_space$covinv
basis_effects_mat = new_fit_CBFM_ptest$basis_effects_mat[,1:num_spacebasisfns,drop=FALSE]
Sigmainv = new_LoadingnuggetSigma_space$covinv
B = B_space
y_vec = as.vector(y)
linpred_vec = c(new_fit_CBFM_ptest$linear_predictor)
dispparam = new_fit_CBFM_ptest$dispparam
powerparam = new_fit_CBFM_ptest$powerparam
zeroinfl_prob_intercept = new_fit_CBFM_ptest$zeroinfl_prob_intercept
return_correlation = TRUE




#------------------------
y = Y[sel_training_units,]
data = X[sel_training_units,]
formula_X = ~ TOW_EFFECT + s(DEPTH) + s(SURFTEMP) + s(BOTTEMP) + s(STRESS_Q95_YR)
B_space =  train_sp_basisfunctions
B_time = train_time_basisfunctions
B_spacetime = NULL
offset = NULL
ncores = NULL
gamma = 1
trial_size = 1
dofit = TRUE
stderrors = TRUE
select = FALSE
family = nb2()
control = list(trace = 1, initial_betas_dampen = 1)
start_params = list(betas = t(sapply(stackedgams_spacetime, coef)[1:38,]))
G_control = list(rank = c(5,5), method = "LA")
Sigma_control = list(rank = c(5,2), method = "LA")
#TMB_directories = list(cpp = system.file("executables", package = "CBFM"), compile = system.file("executables", package = "CBFM"))
TMB_directories = list(cpp = "/home/fh/Dropbox/private/Maths/ANU/Rpackage_CBFM/inst/executables", 
                       compile = "/home/fh/Dropbox/private/Maths/ANU/Rpackage_CBFM/inst/executables")






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
 
# Now generate data from a count distribution 
spp_dispersion <- runif(num_spp, 0, 5)
simy_ztpois <- create_CBFM_life(family = poisson(), formula_X = useformula, data = dat,
                                B_space = basisfunctions, betas = cbind(spp_intercepts_ztnb, spp_slopes_ztnb),
                                G = list(space = true_G_space_ztnb), Sigma = list(space = true_Sigma_space_ztnb),
                                #basis_effects_mat = spp_basis_effects_mat_ztnb,
                                max_resp = 20000) 


# manygam <- foreach(j = 1:num_spp) %dopar%
#     gam(list(response ~ s(temp) + s(depth) + chla + s(O2) + s(x,y), ~1), 
#         data = data.frame(response = c(simy_ztpois$y[,j], numeric(20)), dat[c(1:nrow(dat),1:20),]),
#         family = ziplss())
# 
# start_params = list(betas = t(sapply(manygam, coef)[1:29,]))
# rm(manygam)


# Now fit a zero-truncated count distribution, and CBFM has to ignore the zeros in the data
y = simy_ztpois$y
useformula <- ~ temp + depth + chla + O2
formula_X = useformula
data = dat
B_space = basisfunctions
family =  ztpoisson()
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

useformula <- ~ temp + depth + chla + O2
fitcbfm <- CBFM(y = simy_ztpois$y, formula_X = useformula, data = dat, B_space = basisfunctions, 
                #start_params = start_params, 
                family = ztpoisson(), control = list(trace = 1, initial_betas_dampen = 0.05))

plot(spp_slopes_ztnb, fitcbfm$betas[,2:5])
abline(0,1)

plot(fitcbfm)
summary(fitcbfm)



plot(eta, fitcbfm$linear_predictor)
abline(0,1)
