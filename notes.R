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
## **Example 0: Fitting a CBFM to data from a spatial CBFM
## Estimate betas and basis functions alright, but as expected (although it's worse than I thought) estimation of G and Sigma is pretty poor. It maybe something to do with estimating the scale issue when they are estimated separately
##------------------------------
set.seed(2022)
num_sites <- 500 # 500 (units) sites 
num_spp <- 50 # Number of species
num_X <- 4 # Number of regression slopes
 
spp_slopes <- matrix(runif(num_spp * num_X, -1, 1), nrow = num_spp)
spp_intercepts <- runif(num_spp, -2, 0)
#spp_gear <- rnorm(num_spp, mean = 1.5, sd = 0.2)

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

true_Sigma_space <- rWishart(1, num_basisfunctions+1, diag(x = 0.1, nrow = num_basisfunctions-1))[,,1]/10
true_G_space <- rWishart(1, num_spp+1, diag(x = 0.1, nrow = num_spp))[,,1] %>%
cov2cor

simy <- create_CBFM_life(family = binomial(), formula_X = useformula, data = dat,
                         B_space = basisfunctions, betas = cbind(spp_intercepts, spp_slopes),
                         Sigma = list(space = true_Sigma_space), G = list(space = true_G_space))


fitcbfm_sp <- CBFM(y = simy$y, formula_X = ~ temp + depth + chla + O2, data = dat,
                   B_space = basisfunctions, 
                   family = binomial(), control = list(trace = 1),
                   G_control = list(rank = 5), Sigma_control = list(rank = 5, custom_space = true_Sigma_space))

# Fit CBFM
useformula <- ~ temp + depth + chla + O2
fitcbfm <- CBFM(y = simy, formula_X = useformula, data = dat,
B_space = basisfunctions, B_time = X_year, 
family = binomial(), control = list(trace = 1),
Sigma_control = list(rank = c(5,1)),
G_control = list(rank = c(5,5))
)

library(ggmatplot)
ggmatplot(spp_slopes, fitcbfm_sp$betas[,-1]) + geom_abline(intercept = 0, slope = 1)
ggmatplot(simy$basis_effects_mat, fitcbfm_sp$basis_effects_mat, shape = 1) + geom_abline(intercept = 0, slope = 1)

ggmatplot(true_Sigma_space[lower.tri(true_Sigma_space)], fitcbfm_sp$Sigma_space[lower.tri(fitcbfm_sp$Sigma_space)]) + geom_abline(intercept = 0, slope = 1)
ggmatplot(true_G_space[lower.tri(true_G_space)], fitcbfm_sp$G_space[lower.tri(fitcbfm_sp$G_space)]) + geom_abline(intercept = 0, slope = 1)




##------------------------------
## **Example 1a: Fitting a CBFM to spatial multivariate presence-absence data**
## simulated from a spatial latent variable model
## Please note the data generation process (thus) differs from CBFM.
##------------------------------
set.seed(2022)
num_sites <- 1000 # 500 (units) sites for training set + 500 sites for testing.
num_spp <- 50 # Number of species
num_X <- 4 # Number of regression slopes

#spp_slopes <- matrix(runif(num_spp * num_X, -1, 1), nrow = num_spp)
spp_slopes <- cbind(rnorm(num_spp, -1, sd = 0.25), rnorm(num_spp, 1, sd = 0.25), rnorm(num_spp, -0.25, sd = 0.1), rnorm(num_spp, 0.25, sd = 0.1))
#spp_slopes <- cbind(rnorm(num_spp, 1, sd = 0.25))
spp_intercepts <- rnorm(num_spp, -3, sd = 0.5)

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
#      as.matrix

# Simulate latent variable component
# We will use this information in later examples as well
true_lvs <- RFsimulate(model = RMexp(var=1, scale=2),
x = xy$x, y = xy$y, n = 2)@data %>%
as.matrix
spp_loadings <- matrix(runif(num_spp * 2, -1, 1), nrow = num_spp)
set.seed(NULL)

# Simulate spatial multivariate abundance data (presence-absence)
# We will use this information in later examples as well
eta <- tcrossprod(cbind(1,mm), cbind(spp_intercepts,spp_slopes)) + tcrossprod(true_lvs, spp_loadings)
simy <- matrix(rbinom(num_sites * num_spp, size = 1, prob = plogis(eta)), nrow = num_sites)

# Form training and test sets
dat_train <- dat[1:500,]
dat_test <- dat[501:1000,]
simy_train <- simy[1:500,]
simy_test <- simy[501:1000,]
mm_train <- mm[1:500,,drop=FALSE]
mm_test <- mm[501:1000,,drop=FALSE]
rm(X, mm, spp_loadings, true_lvs, xy, simy, dat)


#-----------------------------------------
# Fit stacked GLM as a baseline
fitstacked <- manyglm(simy_train ~ temp + depth + chla + O2, family = binomial(), data = dat_train)
#fitstacked <- manyglm(simy_train ~ gear, family = binomial(), data = dat_train)


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
fitcbfm_pure <- CBFM(y = simy_train, 
                     formula_X = ~ s(temp) + depth + s(chla) + O2, 
                     data = dat_train,
                     B_space = train_basisfunctions, 
                     family = binomial(), control = list(trace = 1),
                     Sigma_control = list(rank = 5), 
                     G_control = list(rank = 2))


fitcbfm_pure$all_parametric_effects$response <- fitcbfm_pure$all_parametric_effects$response %>% 
     fct_inorder()
ggplot(fitcbfm_pure$all_parametric_effects, aes(x = value, y = partial, group = response, color = response)) +
     geom_line() +
     facet_wrap(. ~ term, nrow = 2) +
     geom_rug(data = fitcbfm_pure$all_parametric_effects %>% dplyr::filter(response == "response1"), 
              aes(x = value, y = partial), sides = "b", alpha = 0.2, color = "black") +
     theme_bw() +
     theme(legend.position = "bottom")


fitcbfm_pure$all_smooth_estimates$response <- fitcbfm_pure$all_smooth_estimates$response %>% 
     fct_inorder()
ggplot(fitcbfm_pure$all_smooth_estimates %>% dplyr::filter(smooth == "s(temp)"), 
       aes(x = temp, y = est, group = response, color = response)) +
     geom_line() +
     theme_bw() +
     theme(legend.position = "bottom")
ggplot(fitcbfm_pure$all_smooth_estimates %>% dplyr::filter(smooth == "s(chla)"), 
       aes(x = chla, y = est, group = response, color = response)) +
     geom_line() +
     theme_bw() +
     theme(legend.position = "bottom")



fitcbfm <- CBFM(y = simy_train, formula_X = ~ depth + chla + O2, data = dat_train,
                B_space = train_basisfunctions, B_time = mm_train,
                family = binomial(),
                control = list(trace = 1),
                #control = list(trace = 1, nonzeromean_B_time = TRUE),
                Sigma_control = list(rank = c(5,"full"), custom_time = Sinv),
                G_control = list(rank = c(2,"full"))
               )


library(ggmatplot)
ggmatplot(spp_slopes, fitcbfm_pure$betas[,-1]) + geom_abline(intercept = 0, slope = 1)
ggmatplot(spp_slopes, fitcbfm$basis_effects_mat[,fitcbfm$num_B_space+(1:fitcbfm$num_B_time)]) + geom_abline(intercept = 0, slope = 1)

(fitcbfm_pure$betas[,-c(1:4)]-fitcbfm$basis_effects_mat[,fitcbfm$num_B_space+(1:fitcbfm$num_B_time)]) %>% summary


(spp_slopes - fitcbfm_pure$betas[,-1]) %>% norm("F")
(spp_slopes - fitcbfm$basis_effects_mat[,fitcbfm$num_B_space+(1:fitcbfm$num_B_time)]) %>% norm("F")


colMeans(fitcbfm_pure$betas[,-1])
fitcbfm$mean_B_time


apply(fitcbfm_pure$betas[,-1], 2, sd)
fitcbfm$Sigma_time %>% diag %>% sqrt




# Calculate predictions onto test dataset
predictions_cbfm_pure <- predict(fitcbfm_pure, newdata = dat_test, type = "response", new_B_space = test_basisfunctions)
predictions_cbfm <- predict(fitcbfm, newdata = dat_test, type = "response", new_B_space = test_basisfunctions, new_B_time = mm_test, se_fit = FALSE)

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
y = simy_train
useformula <- ~ 1
formula_X = ~ temp + depth + chla + O2
data = dat_train
family =  binomial()
B_space = train_basisfunctions
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
Sigma_control = list(rank = c(5), trace = 0)
G_control = list(rank = c(2), trace = 0)
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

library(ggmatplot)
ggmatplot(spp_slopes, fitcbfm$betas[,-1]) + geom_abline(intercept = 0, slope = 1)
qplot(spp_intercepts, fitcbfm$betas[,1]) + geom_abline(intercept = 0, slope = 1)
qplot(spp_gear, fitcbfm$basis_effects_mat[,25]) + geom_abline(intercept = 0, slope = 1)
fitcbfm$basis_effects_mat[,25] %>% summary
fitcbfm$basis_effects_mat[,25] %>% sd
fitcbfm$G_time
fitcbfm$Sigma_time

