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
y = simy_train
useformula <- ~ s(temp) + s(depth) + chla + s(O2)
formula_X = useformula
data = dat_train
B_space = train_basisfunctions
family =  binomial()
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
control = list(maxit = 1000, optim_lower = -5, optim_upper = 5, convergence_type = "parameters", tol = 1e-4, seed = NULL, trace = 1, ridge = 0)
Sigma_control = list(rank = 5, maxit = 1000, tol = 1e-4, method = "LA", trace = 0)
G_control = list(rank = 5, nugget_profile = seq(0.05, 0.95, by = 0.05), maxit = 1000, tol = 1e-4, method = "LA", trace = 0)
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
library(tidyverse)
library(mgcv)
#library(countreg)
library(gamlss)
library(gamlss.tr)
library(gamlss.add)

data("CrabSatellites", package = "countreg")
cs2 <- CrabSatellites
cs2$color <- as.numeric(cs2$color)
cs2 <- subset(cs2, subset = satellites > 0)

fit_tnb1 <- zerotrunc(satellites ~ width + color, data = cs2, dist = "negbin")

ZTNBI <- trun(par = 0, family = "NBI", type = "left", local = FALSE)
fit_tnb2 <- gamlss(satellites~ga(~s(width) + color, method = "ML"), data = cs2, family = ZTNBI) # OK

fit_tnb3 <- gamlss(satellites ~ pb(width) + color, data = cs2, family = ZTNBI)

getSmo(fit_tnb2) %>% plot

s <- getSmo(fit_tnb3)
s$coef



simy <- create_CBFM_life(family = ztpoisson(), formula_X = useformula, data = dat,
                         B_space = basisfunctions, betas = cbind(spp_intercepts, spp_slopes),
                         Sigma = list(space = true_Sigma_space), G = list(space = true_G_space))


spp_dispersion <- runif(num_spp)
simy <- create_CBFM_life(family = ztnb2, formula_X = useformula, data = dat, B_space = basisfunctions, betas = cbind(spp_intercepts, spp_slopes), dispparam = spp_dispersion, max_resp = 20000, Sigma = list(space = true_Sigma_space), G = list(space = true_G_space))



