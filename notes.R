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
library(autoFRK)
library(FRK)
library(MASS)
library(mvabund)
library(mvtnorm)
library(ROCR)
library(sp)
library(RandomFields)
library(tidyverse)
library(ggmatplot)
library(doParallel)
library(foreach)
registerDoParallel(cores = detectCores()-2)


##------------------------------
#' # **Example 1: Fitting a CBFM to data from a spatial CBFM**
#' This is a just a general testing function
##------------------------------
set.seed(2023)
num_sites <- 500 # 500 (units) sites 
num_spp <- 10 # Number of species
num_X <- 4 # Number of regression slopes
 
spp_slopes <- matrix(runif(num_spp * num_X, -1, 1), nrow = num_spp)
spp_intercepts <- runif(num_spp, -2, 0)
spp_gear <- rnorm(num_spp, mean = 1.5, sd = 1)

# Simulate spatial coordinates and environmental covariate components
xy <- data.frame(x = runif(num_sites, 0, 5), y = runif(num_sites, 0, 5))
X <- cbind(rmvnorm(num_sites, mean = rep(0,4)), rep(c(0,1), c(450,50)))
colnames(X) <- c("temp", "depth", "chla", "O2", "gear")
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

simy <- create_CBFM_life(family = binomial(), formula = useformula, data = dat,
                         B_space = basisfunctions, betas = cbind(spp_intercepts, spp_slopes),
                         Sigma = list(space = true_Sigma_space), G = list(space = true_G_space))


#' ## Fit different flavors of CBFMs
fitcbfm_fixed <- CBFM(y = simy$y, formula = ~ s(temp) + depth + chla + O2, data = dat,
                   B_space = basisfunctions, family = binomial(), 
                   control = list(trace = 1), G_control = list(rank = 5), Sigma_control = list(rank = 5))



fake_gam <- gam(simy$y[,1] ~ s(temp), data = dat, family = binomial(), method = "REML", fit = FALSE)
MM_temp <- fake_gam$X[,-1]
Sigma_temp <- .pinv(fake_gam$smooth[[1]]$S[[1]])
fitcbfm_random <- CBFM(y = simy$y, formula = ~ depth + chla + O2, data = dat,
                   B_space = basisfunctions, B_time = MM_temp, family = binomial(), 
                   control = list(trace = 1, inner_maxit = 10),
                   G_control = list(rank = c(5,"full"), structure = c("unstructured", "homogeneous")), 
                   Sigma_control = list(rank = c(5,1), custom_time = Sigma_temp))


Sigma_temp
fitcbfm_random$Sigma_time
fitcbfm_random$Sigma_time/Sigma_temp
fitcbfm_random$G_time

matplot(simy$linear_predictors, fitcbfm_fixed$linear_predictors, pch = 19); abline(0,1)
matplot(simy$linear_predictors, fitcbfm_random$linear_predictors, pch = 19); abline(0,1)



fake_gam <- gam(simy$y[,1] ~ te(temp, depth), data = dat, family = binomial(), method = "REML", fit = FALSE)
MM_tensorprod <- fake_gam$X[,-1]
Sigma_tensorprod <- lapply(1:length(fake_gam$smooth[[1]]$S), function(x) .pinv(fake_gam$smooth[[1]]$S[[x]]))
fitcbfm_tensorprod <- CBFM(y = simy$y, formula = ~ chla + O2, data = dat,
                       B_space = basisfunctions, B_spacetime = MM_tensorprod, family = binomial(), 
                       control = list(trace = 1),
                       G_control = list(rank = c(5,2), structure = c("unstructured", "homogeneous")), 
                       Sigma_control = list(rank = c(5,1), custom_spacetime = Sigma_tensorprod))


fitcbfm_random$basis_effects_mat
Sigma_tensorprod
fitcbfm_tensorprod$Sigma_spacetime
fitcbfm_tensorprod$G_spacetime

dev.off()
matplot(simy$linear_predictors, fitcbfm_fixed$linear_predictors, pch = 19); abline(0,1)
matplot(simy$linear_predictors, fitcbfm_random$linear_predictors, pch = 19); abline(0,1)
matplot(simy$linear_predictors, fitcbfm_tensorprod$linear_predictors, pch = 19); abline(0,1)






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
     y = simy_train
     formula <- useformula
     ziformula <- NULL
     data = dat_train
     family =  stats::binomial() 
     B_space = train_basisfunctions
     B_time = NULL
     knots = NULL
     ziknots = NULL
     B_spacetime = NULL
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
     start_params = list(betas = NULL, zibetas = NULL, basis_effects_mat = NULL, dispparam = NULL, powerparam = NULL)
     TMB_directories = list(cpp = system.file("executables", package = "CBFM"), compile = system.file("executables", package = "CBFM"))
     control = list(trace = 1)
     G_control = list(rank = c(5), structure = c("unstructured"))
     Sigma_control = list(rank = c(5))
     k_check_control = list(subsample = 5000, n.rep = 400)
     }


function() {
     y = simy$y[,1,drop=FALSE]
     formula = ~ chla + O2
     data = dat
     ziformula <- NULL
     B_space = basisfunctions
     B_time = NULL
     B_spacetime = MM_tensorprod
     family = binomial() 
     ncores = 8
     control = list(trace = 1)
     offset = NULL
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
     start_params = list(betas = NULL, zibetas = NULL, basis_effects_mat = NULL, dispparam = NULL, powerparam = NULL)
     TMB_directories = list(cpp = system.file("executables", package = "CBFM"), compile = system.file("executables", package = "CBFM"))
     G_control = list(rank = c(5,"full"), structure = c("unstructured", "homogeneous"))
     Sigma_control = list(rank = c(5,1), custom_spacetime = Sigma_tensorprod)
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
