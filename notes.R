## See github issues for outstanding things to do!!!

## GENERAL APPROACH TO MODIFYING PACKAGE
## 1. Do all changes
## 2. Delete all Rd files in man (but not figures!), the NAMESPACE file, and any compiled C++ files
## 3. Check package
##4. Build (optional)
## 5. Push to github

#Removed from DESCRIPTION since we now turn off compilation of the TMB C++ files [learning from how James Thorson does it with VAST!]
LinkingTo: 
    TMB,
    RcppEigen


##--------------------------------------
#' # Random testing jazz
##-------------------------------------
rm(list = ls())
library(autoFRK)
library(FRK)
library(MASS)
library(mvabund)
library(mvtnorm)
library(ROCR)
library(sp)
library(geoR)
library(tidyverse)
#library(ggmatplot)
library(doParallel)
library(foreach)
registerDoParallel(cores = detectCores()-2)


##------------------------------
#' # **Example 1: Fitting a CBFM to data from a spatial CBFM**
#' This is a just a general testing function
##------------------------------
set.seed(022025)
num_sites <- 1000 # Number of sites
num_spp <- 10 # Number of species
num_X <- 4 # Number of regression slopes
 
spp_slopes <- matrix(runif(num_spp * num_X, -1, 1), nrow = num_spp)
spp_intercepts <- runif(num_spp, -2, 0)

# Simulate spatial coordinates and environmental covariate components
xy <- data.frame(x = runif(num_sites, 0, 5), y = runif(num_sites, 0, 5))
X <- cbind(rmvnorm(num_sites, mean = rep(0,4)), rlnorm(n = num_sites, meanlog = 2, sdlog = 2))
colnames(X) <- c("temp", "depth", "chla", "O2", "areaswept")
dat <- data.frame(xy, X)
useformula <- ~ temp + depth + chla + O2 

# Set up spatial basis functions for CBFM 
num_basisfunctions <- 25 # Number of spatial basis functions to use
basisfunctions <- mrts(dat[,c("x","y")], num_basisfunctions) %>%
as.matrix %>%
{.[,-(1)]} # Remove the first intercept column

true_Sigma_space <- rWishart(1, num_basisfunctions+1, diag(x = 0.1, nrow = num_basisfunctions-1))[,,1]/10
true_G_space <- rWishart(1, num_spp+1, diag(x = 0.1, nrow = num_spp))[,,1] %>%
cov2cor

simy <- create_CBFM_life(family = nb2(), 
                         formula = useformula, 
                         data = dat,
                         B_space = basisfunctions, 
                         betas = cbind(spp_intercepts, spp_slopes),
                         Sigma = list(space = true_Sigma_space), 
                         G = list(space = true_G_space))



#' ## Fit different flavors of CBFMs
fitcbfm_offset1 <- CBFM(y = simy$y[1:500,], 
                      formula = ~ temp + depth + chla + O2, 
                      data = dat[1:500,],
                      offset = matrix(log(dat$areaswept)[1:500], nrow = nrow(simy$y[1:500,]), ncol = ncol(simy$y), byrow = FALSE),
                      B_space = basisfunctions[1:500,], 
                      family = binomial(), 
                      control = list(trace = 1), 
                      G_control = list(rank = 2),
                      Sigma_control = list(rank = 2))

fitcbfm_offset2 <- CBFM(y = simy$y[1:500,], 
                        formula = ~ temp + depth + chla + O2 + offset(log(areaswept)), 
                        data = dat[1:500,],
                        B_space = basisfunctions[1:500,], 
                        family = binomial(), 
                        control = list(trace = 1), 
                        G_control = list(rank = 2),
                        Sigma_control = list(rank = 2))

fitcbfm_nooffset <- CBFM(y = simy$y[1:500,], 
                        formula = ~ temp + depth + chla + O2, 
                        data = dat[1:500,],
                        B_space = basisfunctions[1:500,], 
                        family = binomial(), 
                        control = list(trace = 1), 
                        G_control = list(rank = 2),
                        Sigma_control = list(rank = 2))


fitcbfm_offset1
fitcbfm_offset2

fitcbfm_offset1$betas - fitcbfm_offset2$betas
fitcbfm_nooffset$betas - fitcbfm_offset2$betas

#' Should be the same; for the first model predict function does not know about the offset, but the second model does since offset is included in formula
plot(plogis(predict(fitcbfm_offset1, newdata = dat[501:1000,]) + matrix(log(dat$areaswept)[501:1000], nrow = 500, ncol = 10, byrow = FALSE)),
     predict(fitcbfm_offset2, newdata = dat[501:1000,], type = "response")) 



##------------------------------
#' # **Example 2 modified: Fitting a CBFM to spatial multivariate presence-absence data simulated from a spatial latent variable model**
#' Please note the data generation process (thus) differs from CBFM.
#' This example was constructed to test the G_control$structure = "identity" argument, which is like fitting a hierarchical GAM
##------------------------------
set.seed(082022)
num_sites <- 1000 # 500 (units) sites for training set + 500 sites for testing.
num_spp <- 20 # Number of species. Need to keep this lower than usual as factor smooths do not scale very well with the number of factor levels
num_X <- 4 # Number of regression slopes

# Simulate spatial coordinates and environmental covariate components
# We will use this information in later examples as well
xy <- data.frame(x = runif(num_sites, 0, 5), y = runif(num_sites, 0, 5))
X <- rmvnorm(num_sites, mean = rep(0,4))
colnames(X) <- c("temp", "depth", "chla", "O2")
dat <- data.frame(xy, X)
f1 <- function(x, a, b) a * sin(pi * b * x)

spp_slopes <- cbind(a = rnorm(num_spp)+2, b = rnorm(num_spp, sd = 0.5))
spp_intercepts <- rnorm(num_spp, 0, sd = 0.5)


# Simulate spatial multivariate abundance data
eta <- matrix(spp_intercepts, nrow = num_sites, ncol = num_spp, byrow = TRUE) + sapply(1:num_spp, function(j) f1(dat$temp, a = spp_slopes[j,1], b = spp_slopes[j,2]))
simy <- matrix(rbinom(num_sites * num_spp, size = 1, prob = plogis(eta)), num_sites, num_spp)
colnames(simy) <- paste0("spp", 1:num_spp)

# Form training and test sets
dat_train <- dat[1:500,]
dat_test <- dat[501:1000,]
simy_train <- simy[1:500,]
simy_test <- simy[501:1000,]
rm(X, mm, spp_loadings, true_lvs, xy, simy, dat)


# Fit a hierarchical GAM
dat_long <- data.frame(simy_train, dat_train) %>% 
     pivot_longer(spp1:spp20, names_to = "species") %>% 
     mutate(species = fct_inorder(species))

hgam <- gam(value ~ s(temp, by = species, bs = "gp", id = 1), data = dat_long, family = "binomial", method = "REML")
#hgam_select <- mgcv::gam(value ~ s(temp, by = species, id = 1), data = dat_long, select = TRUE, family = "binomial", method = "REML")
#hgam_fs <- mgcv::gam(value ~ s(temp, species, bs = "fs"), data = dat_long, family = "binomial", method = "REML")

hgam$sp
data.frame(temp = dat_train$temp, hgam = matrix(hgam$linear.predictors, ncol = num_spp, byrow = TRUE)) %>% 
     pivot_longer(-temp) %>% 
     mutate(species = rep(paste0("spp",1:num_spp),500)) %>% 
     ggplot(., aes(x = temp, y = value)) +
     geom_line() +
     facet_wrap(. ~ species, nrow = 5)


sm <- smoothCon(s(temp, bs = "gp"), data = dat_train, knots = NULL, absorb.cons = TRUE)[[1]]

hgam$smooth[[1]]$S[[1]] - sm$S[[1]] ## Basically matches
hgam_mm <- model.matrix(hgam)[,2:12] # Extract only the first species (first column is an intercept), since a factor smooth basically Kronecker products (up to a permutation of the levels of the factor) this matrix num_spp times
hgam_mm <- hgam_mm[seq(1,nrow(hgam_mm),by=20),]
norm(hgam_mm - sm$X) ## matches

## So hgam matches construction of a smoother for each species, where the identifiability constraints are absorbed into the basis, and the penalty is scaled. You should be able to get away without the latter (see scale.penalty argument in ?smoothCon). The former I am less sure about, but would be useful to do without in order to work with nicer penalty matrices?
rm(sm)



# Fit CBFMs
sm_train_useincbfm <- smoothCon(s(temp, bs = "gp"), data = dat_train, knots = NULL, absorb.cons = TRUE)[[1]] 
sm_train_useincbfm$X # Compare summary(sm$X) and summary(sm_useincbfm$X) and you can see the identifiability not being absorbed as the matrix is no longer mean centered. But note that by not absorbing the constraint, the intercept term of the smooth becomes clear

# For CBFM, drop intercept as it is already contained in the formula argument -- needed when absorb.cons = FALSE
# mm_train_useincbfm <- sm_train_useincbfm$X[,!apply(sm_train_useincbfm$X==1,2,all)]
# penmat_useincbfm <- sm_train_useincbfm$S[[1]][!apply(sm_train_useincbfm$X==1,2,all),!apply(sm_train_useincbfm$X==1,2,all)]
# mm_test_useincbfm <- PredictMat(sm_train_useincbfm, data = dat_test)
# mm_test_useincbfm <- mm_test_useincbfm[,!apply(sm_train_useincbfm$X==1,2,all)]
# rm(sm_train_useincbfm)
     
mm_train_useincbfm <- sm_train_useincbfm$X
penmat_useincbfm <- sm_train_useincbfm$S[[1]]
mm_test_useincbfm <- PredictMat(sm_train_useincbfm, data = dat_test)
rm(sm_train_useincbfm)


fitcbfm <- CBFM(y = simy_train, 
                    formula = ~ 1, 
                    data = dat_train,
                    B_space = mm_train_useincbfm, 
                    family = binomial(), 
                    control = list(trace = 1, optim_lower = -5000, optim_upper = 5000),
                    Sigma_control = list(rank = "full", custom_space = .pinv(penmat_useincbfm)), 
                    G_control = list(rank = "full", structure = "identity")) #custom_space = diag(gam.vcomp(hgam, rescale = FALSE)[[1]][1]^2, nrow = num_spp)


fitcbfm$G_space
fitcbfm$basis_effects_mat

data.frame(temp = dat_train$temp, cbfm = fitcbfm$linear_predictors, hgam = matrix(hgam$linear.predictors, ncol = num_spp, byrow = TRUE)) %>% 
     pivot_longer(-temp) %>% 
     mutate(species = rep(rep(paste0("spp",1:num_spp),2),500)) %>% 
     mutate(model = rep(rep(c("cbfm","gam"),each=num_spp),500)) %>% 
     ggplot(., aes(x = temp, y = value, color = model)) +
          geom_line() +
     facet_wrap(. ~ species, nrow = 5)
     
## Not perfect but it's not too bad either? Setting optim_lower/upper is pretty critical for this




##----------------------------------
#' Custom testing
##----------------------------------
function() {
     y = Y[sel_training_units,c(1,5), drop=FALSE]
     formula = 
          ~ offset(log(AREA_SWEPT_EST)) + 
          SOURCE + 
          s(BOTTEMP_DOPGLO_monthly_MEAN, bs="tp", m=2, k=5) +
          s(BOTSALIN_DOPGLO_monthly_MEAN, bs="tp", m=2, k=5) +
          s(logDEPTH_1km, bs="tp", m=2, k=5) +
          s(logTideVel_max, bs="tp", m=2, k=5) +
          s(logBPI_1km, bs="tp", m=2, k=5) +
          s(logCOMPLEX_1km, bs="tp", m=2, k=5) +
          s(GRAIN_1km, bs="tp", m=2, k=5)
     ziformula <- NULL
     data = X[sel_training_units,]
     family =  CBFM::nb2() 
     B_space = NULL
     B_time = NULL
     B_spacetime = MM_train_year_gp
     knots = NULL
     ziknots = NULL
     offset = NULL
     ncores = NULL
     gamma = 1
     zigamma = 1
     trial_size = 1
     nonzeromean_B_space = FALSE
     nonzeromean_B_time = FALSE
     nonzeromean_B_spacetime = FALSE
     dofit = TRUE
     stderrors = TRUE
     select = FALSE
     ziselect = FALSE
     start_params = list(betas = (sapply(GAMs_simple_yeargp_ident[c(1,5)], coef)[1:33,]*0.1) %>% t,
                         dispparam = 1/sapply(GAMs_simple_yeargp_ident[c(1,5)], function(x) x$family$getTheta(TRUE)))
     #start_params = list(betas = NULL, zibetas = NULL, basis_effects_mat = NULL, dispparam = NULL, powerparam = NULL)
     TMB_directories = list(cpp = system.file("executables", package = "CBFM"), compile = system.file("executables", package = "CBFM"))
     control = list(trace = 1, maxit=2, final_maxit=100, initial_betas_dampen = c(1, 1)) 
     G_control = list(rank = c("full"), structure = c("identity"))
     Sigma_control = list(rank = c("full"), 
                          custom_spacetime = Sigma_year_gp, 
                          control = list(trace = 1))     
     k_check_control = list(subsample = 5000, n.rep = 400)
     
     
     }


function() {
     y = simy_train
     formula = useformula
     data = dat_train
     ziformula <- NULL
     B_space = train_basisfunctions
     B_time = NULL
     B_spacetime = NULL
     family = nb2() 
     ncores = detectCores() - 2
     control = list(trace = 1, initial_ridge = 0.5)
     offset = NULL
     gamma = 1
     zigamma = 1
     knots = NULL
     trial_size = 1
     nonzeromean_B_space = FALSE
     nonzeromean_B_time = FALSE
     nonzeromean_B_spacetime = FALSE
     dofit = TRUE
     stderrors = TRUE
     select = FALSE
     ziselect = FALSE
     start_params = list(betas = NULL, zibetas = NULL, basis_effects_mat = NULL, dispparam = NULL, powerparam = NULL)
     TMB_directories = list(cpp = system.file("executables", package = "CBFM"), compile = system.file("executables", package = "CBFM"))
     G_control = list()
     Sigma_control = list()
     k_check_control = list(subsample = 5000, n.rep = 400)
     
     
     }



function() {
     Sigmainv = new_LoadingnuggetSigma_spacetime$invcov
     basis_effects_mat = centered_BF_mat 
     Ginv = new_LoadingnuggetG_spacetime$invcov
     B = B_spacetime
     X = X
     ziX = ziX
     y_vec = as.vector(y)
     linpred_vec = c(new_fit_CBFM_ptest$linear_predictors)
     dispparam = new_fit_CBFM_ptest$dispparam
     powerparam = new_fit_CBFM_ptest$powerparam
     zibetas = new_fit_CBFM_ptest$zibetas 
     estimate_lambda = estimate_lambda_not_Sigma
     which_B = 3
     }
