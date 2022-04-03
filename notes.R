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
set.seed(2021)
num_sites <- 500 # 500 (units) sites
num_spp <- 20 # Number of species
num_X <- 4 # Number of regression slopes

spp_slopes <- matrix(runif(num_spp * num_X, -1, 1), nrow = num_spp)
spp_intercepts <- runif(num_spp, -2, 0)

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
set.seed(2021)
num_sites <- 1000 # 500 (units) sites for training set + 500 sites for testing.
num_spp <- 50 # Number of species
num_X <- 4 # Number of regression slopes

#spp_slopes <- matrix(runif(num_spp * num_X, -1, 1), nrow = num_spp)
spp_slopes <- cbind(rnorm(num_spp, -1, sd = 0.25), rnorm(num_spp, 1, sd = 0.25), rnorm(num_spp, -0.25, sd = 0.1), rnorm(num_spp, 0.25, sd = 0.1))
# true_G_betas <- rWishart(1, num_spp, diag(num_spp))[,,1] %>% cov2cor 
# true_Sigma_betas <- diag(x = c(0.5,0.5,0.1,0.1))
# true_mean_betas <- c(-1,1,-0.25,0.25)
# true_mean_betas <- rep(true_mean_betas, num_spp)
# spp_slopes <- rmvnorm(n = 1, mean = true_mean_betas, sigma = kronecker(true_G_betas, true_Sigma_betas)) %>% 
#     matrix(nrow = num_spp, byrow = TRUE)
spp_intercepts <- rnorm(num_spp, -3, sd = 0.5)

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
eta <- tcrossprod(cbind(1,mm), cbind(spp_intercepts,spp_slopes)) + tcrossprod(true_lvs, spp_loadings)
simy <- matrix(rbinom(num_sites * num_spp, size = 1, prob = plogis(eta)), nrow = num_sites)

# Form training and test sets
dat_train <- dat[1:500,]
dat_test <- dat[501:1000,]
simy_train <- simy[1:500,]
simy_test <- simy[501:1000,]
mm_train <- mm[1:500,]
mm_test <- mm[501:1000,]
rm(X, mm, spp_loadings, true_lvs, xy, simy, dat)


#-----------------------------------------
# Fit stacked GLM as a baseline
fitstacked <- manyglm(simy_train ~ temp + depth + chla + O2, family = binomial(), data = dat_train)


# # Fit HMSC (gold standard) -- Apply to probit data only
# library(Hmsc)
# studyDesign = data.frame(sample = factor(rownames(dat_train)))
# rL.nngp = HmscRandomLevel(sData = dat_train %>% dplyr::select(x:y), sMethod = 'NNGP', nNeighbours = 20)
# rL.nngp = setPriors(rL.nngp, nfMin = 2, nfMax = 2) # Note Hmsc scales really poorly with increasing number of latent variables!
# m.nngp = Hmsc(Y = simy_train, XData = dat_train, XScale = FALSE, XFormula = ~ temp + depth + chla + O2, 
#               studyDesign = studyDesign, ranLevels = list("sample" = rL.nngp), distr = "probit")
# 
# m.nngp <- sampleMcmc(m.nngp, thin = 10, samples = 100, transient = 100, nChains = 3, updater = list(GammaEta=FALSE), nParallel = 3)  
# 
# est_coefs <- getPostEstimate(m.nngp, parName = "Beta")$mean %>% 
#     t %>% 
#     {.[,-1]}
# 
# 
# ggmatplot(spp_slopes, est_coefs) + geom_abline(intercept = 0, slope = 1)
# 
# (spp_slopes - est_coefs) %>% norm


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
fitcbfm_pure <- CBFM(y = simy_train, formula_X = ~ temp + depth + chla + O2, 
                     data = dat_train,
                     B_space = train_basisfunctions, 
                     family = binomial(), control = list(trace = 1),
                     G_control = list(rank = 2), Sigma_control = list(rank = 5))


fitcbfm <- CBFM(y = simy_train, formula_X = ~ 1, data = dat_train,
                B_space = train_basisfunctions, B_time = mm_train,
                family = binomial(),
                control = list(trace = 1, nonzeromean_B_time = TRUE),
                Sigma_control = list(rank = c(5,"full")),
                G_control = list(rank = c(2,"full"))
                )


library(ggmatplot)
ggmatplot(spp_slopes, fitcbfm_pure$betas[,-1]) + geom_abline(intercept = 0, slope = 1)
ggmatplot(spp_slopes, fitcbfm$basis_effects_mat[,fitcbfm$num_B_space+(1:fitcbfm$num_B_time)]) + geom_abline(intercept = 0, slope = 1)
ggmatplot(spp_slopes, out_CBFM$basis_effects_mat[,fitcbfm$num_B_space+(1:fitcbfm$num_B_time)]) + geom_abline(intercept = 0, slope = 1)

(spp_slopes - fitcbfm_pure$betas[,-1]) %>% norm("F")
(spp_slopes - fitcbfm$basis_effects_mat[,fitcbfm$num_B_space+(1:fitcbfm$num_B_time)]) %>% norm("F")
(spp_slopes - out_CBFM$basis_effects_mat[,out_CBFM$num_B_space+(1:out_CBFM$num_B_time)]) %>% norm("F")


colMeans(fitcbfm_pure$betas[,-1])
fitcbfm$mean_B_time
out_CBFM$mean_B_time

apply(fitcbfm_pure$betas[,-1], 2, sd)
fitcbfm$Sigma_time %>% diag %>% sqrt
new_LoadingnuggetSigma_time$cov %>% diag %>% sqrt




# Calculate predictions onto test dataset
predictions_cbfm_pure <- predict(fitcbfm_pure, newdata = dat_test, type = "response", new_B_space = test_basisfunctions)
predictions_cbfm <- predict(fitcbfm2, newdata = dat_test, type = "response", new_B_space = test_basisfunctions, new_B_time = mm_test, se_fit = FALSE)

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





y = simy_train
useformula <- ~ 1
formula_X = useformula
data = dat_train
family =  binomial()
B_space = train_basisfunctions
B_time = mm_train
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
Sigma_control = list(rank = c(5,"full"), maxit = 100, tol = 1e-4, method = "LA", trace = 0)
G_control = list(rank = c(2,"full"), nugget_profile = seq(0.05, 0.95, by = 0.05), maxit = 100, tol = 1e-4, method = "LA", trace = 0)
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


