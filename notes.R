## See github issues for outstanding things to do!!!

## GENERAL APPROACH TO MODIFYING PACKAGE
## 1. Do all changes
## 1. Consider moving calls of other packages from :: to using importFrom
## 2. Delete all Rd files in man (but not figures!), the NAMESPACE file, and any compiled C++ files
## 3. Check package
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
library(ggmatplot)
library(doParallel)
library(foreach)
registerDoParallel(cores = detectCores()-2)


##------------------------------
## **Example 0: Fitting a CBFM to data from a spatial CBFM
## Estimate betas and basis functions alright, but as expected (although it's worse than I thought) estimation of G and Sigma is pretty poor. It maybe something to do with estimating the scale issue when they are estimated separately
##------------------------------
set.seed(2022)
num_sites <- 500 # 500 (units) sites 
num_spp <- 50 # Number of species
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


# Fit CBFMs
fitcbfm_fixed <- CBFM(y = simy$y, formula = ~ s(temp) + depth + chla + O2, data = dat,
                   B_space = basisfunctions, family = binomial(), 
                   control = list(trace = 1), G_control = list(rank = 5), Sigma_control = list(rank = 5))

fake_gam <- gam(simy$y[,1] ~ s(temp), data = dat)
MM_temp <- model.matrix(fake_gam)[,-1]
Sigma_temp <- .pinv(fake_gam$smooth[[1]]$S[[1]])

fitcbfm_random <- CBFM(y = simy$y, formula = ~ depth + chla + O2, data = dat,
                   B_space = basisfunctions, B_time = MM_temp, family = binomial(), 
                   control = list(trace = 1, initial_betas_dampen = 1), 
                   G_control = list(rank = c(5,"full")), 
                   Sigma_control = list(rank = c(5,1), custom_time = Sigma_temp))


ggmatplot(fitcbfm_fixed$betas[,-c(1:4)], fitcbfm_random$basis_effects_mat[,-c(1:24)], shape = 19) + 
     geom_abline(intercept = 0, slope = 1) +
     scale_color_viridis_d()



ggmatplot(cbind(spp_slopes,spp_gear), fitcbfm_fixed$betas[,-1]) + geom_abline(intercept = 0, slope = 1)
ggmatplot(cbind(spp_slopes,spp_gear), cbind(fitcbfm_random$betas[,-1], fitcbfm_random$basis_effects_mat[,1])) + geom_abline(intercept = 0, slope = 1)


ggmatplot(simy$basis_effects_mat, fitcbfm_sp$basis_effects_mat, shape = 1) + geom_abline(intercept = 0, slope = 1)
ggmatplot(true_Sigma_space[lower.tri(true_Sigma_space)], fitcbfm_sp$Sigma_space[lower.tri(fitcbfm_sp$Sigma_space)]) + geom_abline(intercept = 0, slope = 1)
ggmatplot(true_G_space[lower.tri(true_G_space)], fitcbfm_sp$G_space[lower.tri(fitcbfm_sp$G_space)]) + geom_abline(intercept = 0, slope = 1)



##------------------------------
## **Example 1 modified: Fitting a CBFM to spatial multivariate data**
## simulated from a spatial latent variable model
## Please note the data generation process (thus) differs from CBFM.
##------------------------------
set.seed(2022)
num_sites <- 1000 # 500 (units) sites for training set + 500 sites for testing.
num_spp <- 50 # Number of species
num_X <- 4 # Number of regression slopes

spp_slopes <- matrix(runif(num_spp * num_X, -1, 1), nrow = num_spp)
spp_intercepts <- rnorm(num_spp, -3, sd = 0.5)
spp_zislopes <- matrix(runif(num_spp * num_X, -1, 1), nrow = num_spp)
spp_ziintercepts <- runif(num_spp, -0.5, 0)

# Simulate spatial coordinates and environmental covariate components
# We will use this information in later examples as well
xy <- data.frame(x = runif(num_sites, 0, 5), y = runif(num_sites, 0, 5))
X <- rmvnorm(num_sites, mean = rep(0,4))
colnames(X) <- c("temp", "depth", "chla", "O2")
# X <- matrix(rep(c(0,1), num_sites*c(0.4,0.6)), ncol = 1)
# colnames(X) <- "gear"
dat <- data.frame(xy, X)
mm <- model.matrix(~ temp + depth + chla + O2 - 1, data = dat) %>%
 scale %>%
 as.matrix
# mm <- model.matrix(~ gear - 1, data = dat) %>%
#       as.matrix

# Simulate latent variable component
# We will use this information in later examples as well
true_lvs <- RFsimulate(model = RMexp(var=1, scale=2),
x = xy$x, y = xy$y, n = 2)@data %>%
as.matrix
spp_loadings <- matrix(runif(num_spp * 2, -1, 1), nrow = num_spp)
set.seed(NULL)

# Simulate spatial multivariate abundance data (presence-absence)
eta <- tcrossprod(cbind(1,mm), cbind(spp_intercepts,spp_slopes)) + tcrossprod(true_lvs, spp_loadings)
zieta <- tcrossprod(cbind(1,mm), cbind(spp_ziintercepts,spp_zislopes))
component_ind <- matrix(rbinom(num_sites * num_spp, size = 1, prob = matrix(plogis(zieta), num_sites, num_spp, byrow = TRUE)), num_sites,num_spp)
simy <- matrix(rpois(num_sites * num_spp, lambda = exp(eta+2) * (1-component_ind)), num_sites, num_spp)

# Form training and test sets
dat_train <- dat[1:500,]
dat_test <- dat[501:1000,]
simy_train <- simy[1:500,]
simy_test <- simy[501:1000,]
mm_train <- mm[1:500,,drop=FALSE]
mm_test <- mm[501:1000,,drop=FALSE]
rm(X, mm, spp_loadings, true_lvs, xy, simy, dat)



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


# Fit CBFMs
fitcbfm_zip <- CBFM(y = simy_train, 
                     formula = ~ s(temp) + depth + s(chla) + O2, 
                     ziformula = ~ s(temp) + depth + s(chla) + O2, 
                     data = dat_train,
                     B_space = train_basisfunctions, 
                     family = zipoisson(), 
                     control = list(trace = 1),
                     Sigma_control = list(rank = 5), 
                     G_control = list(rank = 2))



# Calculate predictions onto test dataset
predictions_cbfm_gold <- predict(fitcbfm_zip, type = "response")
predictions_cbfm_test <- predict(fitcbfm_zip, type = "response", se_fit = TRUE)
matplot(predictions_cbfm_gold, predictions_cbfm_test$fit, log = "|xy"); abline(0,1)


predictions_cbfm_pure <- predict(fitcbfm_zip, newdata = dat_test, type = "response", new_B_space = test_basisfunctions)


# Evaluation predictions
# Tjur R-squared across species
tjurR2 <- data.frame(
cbfm_pure = sapply(1:num_spp, function(j) {
m1 <- predictions_cbfm_pure[which(simy_test[,j] > 0),j] %>%
mean(na.rm = TRUE)
m0 <- predictions_cbfm_pure[which(simy_test[,j] == 0),j] %>%
mean(na.rm = TRUE)
m1 - m0
}),
cbfm = sapply(1:num_spp, function(j) {
m1 <- predictions_cbfm[which(simy_test[,j] > 0),j] %>%
mean(na.rm = TRUE)
m0 <- predictions_cbfm[which(simy_test[,j] == 0),j] %>%
mean(na.rm = TRUE)
m1 - m0
})
)

boxplot(tjurR2, main = "Tjur-R2", names = c("CBFM_pure", "CBFM"))

ggplot(tjurR2, aes(x = cbfm_pure, y = cbfm)) +
geom_point() +
geom_abline(intercept = 0, slope = 1, linetype = 2) +
labs(x = "Stacked SDM", y = "CBFM", main = "Tjur-R2") +
theme_bw()

# AUC across species
aucs <- data.frame(
cbfm_pure = sapply(1:num_spp, function(j) {
pred <- prediction(predictions_cbfm_pure[,j], labels = simy_test[,j]) %>%
performance(measure = "auc")
pred@y.values[[1]]
}),
cbfm = sapply(1:num_spp, function(j) {
pred <- prediction(predictions_cbfm[,j], labels = simy_test[,j]) %>%
performance(measure = "auc")
pred@y.values[[1]]
})
)

boxplot(aucs, main = "AUC", names = c("CBFM_pure", "CBFM"))

ggplot(aucs, aes(x = cbfm_pure, y = cbfm)) +
geom_point() +
geom_abline(intercept = 0, slope = 1, linetype = 2) +
labs(x = "Stacked SDM", y = "CBFM", main = "AUC") +
theme_bw()





##------------------------------
## **Example 2 modified: Fitting a CBFM to spatial multivariate presence-absence data**
## simulated from a spatial latent variable model
## Please note the data generation process (thus) differs from CBFM.
## This example was constructed to test the G_control$structure = "identity" argument, which is like fitting a hierarchical GAM
##------------------------------
set.seed(2022)
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

hgam <- mgcv::gam(value ~ s(temp, by = species, bs = "gp", id = 1), data = dat_long, family = "binomial", method = "REML")
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
# Set up spatial basis functions for CBFM -- Most users will start here!
# We will also use this basis functions in some later examples
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
     



function() {
     y = simy_train
     useformula <- ~ 1
     formula <- useformula
     ziformula <- NULL
     data = dat_train
     family =  binomial() 
     B_space = mm_train_useincbfm
     B_time = NULL
     B_spacetime = NULL
     offset = NULL
     ncores = NULL
     gamma = 1
     zigamma = 1
     trial_size = 1
     dofit = TRUE
     stderrors = TRUE
     select = FALSE
     ziselect = FALSE
     start_params = list(betas = NULL, zibetas = NULL, basis_effects_mat = NULL, dispparam = NULL, powerparam = NULL)
     TMB_directories = list(cpp = system.file("executables", package = "CBFM"), compile = system.file("executables", package = "CBFM"))
     control = list(maxit = 100, convergence_type = "parameters_MSE", tol = 1e-4, seed = NULL, trace = 1, ridge = 0)
     G_control = list(rank = c("full"), structure = "identity")
     Sigma_control = list(rank = c("full"), custom_space = penmat_useincbfm)
     k_check_control = list(subsample = 5000, n.rep = 400)
}


# Calculate predictions onto test dataset
predictions_cbfm_gold <- predict(fitcbfm_zip, type = "response")
predictions_cbfm_test <- predict(fitcbfm_zip, type = "response", se_fit = TRUE)
matplot(predictions_cbfm_gold, predictions_cbfm_test$fit, log = "|xy"); abline(0,1)


predictions_cbfm_pure <- predict(fitcbfm_zip, newdata = dat_test, type = "response", new_B_space = test_basisfunctions)



##----------------------------------
## Custom testing
##----------------------------------
y = simy_train
useformula <- ~ temp + depth + chla + O2
formula <- useformula
ziformula <- NULL
data = dat_train
family =  ztnb2() 
B_space = train_basisfunctions
B_time = NULL
B_spacetime = NULL
offset = NULL
ncores = NULL
gamma = 1
zigamma = 1
trial_size = 1
dofit = TRUE
stderrors = TRUE
select = FALSE
ziselect = FALSE
start_params = list(betas = NULL, zibetas = NULL, basis_effects_mat = NULL, dispparam = NULL, powerparam = NULL)
TMB_directories = list(cpp = system.file("executables", package = "CBFM"), compile = system.file("executables", package = "CBFM"))
control = list(maxit = 100, convergence_type = "parameters_MSE", tol = 1e-4, seed = NULL, trace = 1, ridge = 0)
G_control = list(rank = c(5), structure = "unstructured")
Sigma_control = list(rank = c(5))
k_check_control = list(subsample = 5000, n.rep = 400)





Ginv = new_LoadingnuggetG_space$covinv
basis_effects_mat = centered_BF_mat
Sigmainv = new_LoadingnuggetSigma_space$covinv
B = B_space
y_vec = as.vector(y)
linpred_vec = c(new_fit_CBFM_ptest$linear_predictors)
dispparam = new_fit_CBFM_ptest$dispparam
powerparam = new_fit_CBFM_ptest$powerparam
zibetas = new_fit_CBFM_ptest$zibetas
return_correlation = is.null(Sigma_control$custom_space)