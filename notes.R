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
## **Example 0: Fitting a CBFM to data from a spatial CBFM[]
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
                         B_space = basisfunctions*0, betas = cbind(spp_intercepts, spp_slopes),
                         Sigma = list(space = true_Sigma_space), G = list(space = true_G_space))


##------------------------------
## Example 0.5: Try to address the above problem by exploring a single species GAM versus CBFM fits
##------------------------------
stackedgams_fn <- function(j, y, formula_X, data) {
    tmp_formula <- as.formula(paste("response", paste(as.character(formula_X),collapse="") ) )
    fit0 <- gam(tmp_formula, data = data.frame(response = y[,j], data), method = "REML", family = binomial())
    return(fit0)
    }

stackedgams <- foreach(j = 1:ncol(simy$y)) %dopar% stackedgams_fn(j = j, y = simy$y, formula_X = ~ depth + chla + O2 + s(temp), data = dat)
stackedgams_coef <- sapply(stackedgams, coef)[-(1:4),] %>% t



MM_temp <- model.matrix(stackedgams[[1]])[,-c(1:4),drop=FALSE]
Sigma_temp <- .pinv(stackedgams[[1]]$smooth[[1]]$S[[1]])
G_temp <- sapply(stackedgams, function(x) gam.vcomp(x, rescale = FALSE)[[1]][1]^2)

fitcbfm <- CBFM(y = simy$y[,1:5,drop=FALSE], formula = ~ 1, data = dat,
                B_space = MM_temp, 
                family = binomial(), 
                control = list(trace = 1, convergence_type = "parameters_norm", tol = 1e-6),
                #G_control = list(rank = "full", custom_space = diag(x = G_temp[1:10])),
                G_control = list(rank = "full", trace = 1),
                Sigma_control = list(rank = "full", custom_space = Sigma_temp))


fitcbfm_singlespp <- CBFM(y = simy$y[,5,drop=FALSE], formula = ~ depth + chla + O2, data = dat,
                B_space = MM_temp, 
                family = binomial(), 
                control = list(trace = 1, convergence_type = "parameters_norm", tol = 1e-6),
                #G_control = list(rank = "full", custom_space = diag(x = G_temp[1:10])),
                G_control = list(rank = "full", trace = 1),
                Sigma_control = list(rank = 1, custom_space = Sigma_temp))

fitcbfm_singlespp$basis_effects_mat

fitcbfm$basis_effects_mat
fitcbfm$G_space
fitcbfm$Sigma_space


ggmatplot(t(stackedgams_coef[1:5,]), t(fitcbfm$basis_effects_mat), shape = 19) + 
     geom_abline(intercept = 0, slope = 1) +
     scale_color_viridis_d()



y = simy$y
useformula <- ~ depth + chla + O2
formula <- useformula
ziformula <- NULL
data = dat
family =  binomial() 
B_space = MM_temp
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
control = list(maxit = 100, convergence_type = "parameters_norm", tol = 1e-6, seed = NULL, trace = 1, ridge = 0)
G_control = list(rank = c("full"), trace = 1, method = "ML", tol = 1e-6)
Sigma_control = list(rank = c(1), custom_space = Sigma_temp)
k_check_control = list(subsample = 5000, n.rep = 400)




##------------------------------
## END Example 0.5. Return to CBFM and stacked fits
##------------------------------
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
## **Example 1 modified: Fitting a CBFM to spatial multivariate presence-absence data**
## simulated from a spatial latent variable model
## Please note the data generation process (thus) differs from CBFM.
##------------------------------
set.seed(072022)
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




##----------------------------------
## Custom testing
##----------------------------------
y = simy$y
useformula <- ~ depth + chla + O2
formula <- useformula
ziformula <- NULL
data = dat
family =  binomial() 
B_space = basisfunctions
B_time = MM_temp
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
control = list(maxit = 100, convergence_type = "parameters", tol = 1e-4, seed = NULL, trace = 1, ridge = 0)
G_control = list(rank = c(5,"full"))
Sigma_control = list(rank = c(5,1), custom_time = Sigma_temp)
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

