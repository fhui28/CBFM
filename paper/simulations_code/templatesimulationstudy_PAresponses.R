##--------------------
## Template code for simulation study -- presence-absence responses. The script can be straightforwardly modified to work for Poisson responses, different sample sizes, 
## Spatial multivariate abundance data are generated from a stacked species distribution model
##--------------------
rm(list = ls())
library(autoFRK)
library(here)
library(doParallel)
library(foreach)
library(mvabund)
library(tidyverse)
library(abind)
library(ROCR)
library(sjSDM)
library(Hmsc)
library(CBFM) # devtools::install_github("fhui28/CBFM")     
source("simulatefunctions.R")


# Number of simulated datasets 
nsims <- 200

# Number of spatial locations
nlocs <- 500/0.8


# Parameters in true model
set.seed(2021)
num_spp <- 30
candidate_coefs <- data.frame(
     big_coefs = runif(num_spp, -1.5, 1.5), 
     small_coefs = runif(num_spp, -0.5, 0.5)
     )
spp_slopes <- data.frame(
     temp = sample(candidate_coefs$small_coefs), 
     depth = sample(candidate_coefs$big_coefs), 
     chla = sample(candidate_coefs$small_coefs),
     O2 = sample(candidate_coefs$big_coefs)
     ) 
rownames(spp_slopes) <- paste0("spp", 1:num_spp)
set.seed(NULL)


# Get coefficients corresponding to the covariates that are actually included in the models to be fitted
spp_slopes_long <- spp_slopes %>%
     as.data.frame %>%
     rownames_to_column(var = "species") %>%
     dplyr::select(c(species, temp, depth)) %>%
     pivot_longer(-species)
spp_slopes <- as.matrix(spp_slopes)



# Commence simuulation -- Note depending on the computing resources available, the user may wish to parallelize the simulations across multiple machines. 
for (i in 1:nsims) {
     
     #------------------
     # Generate data
     # Assume that only the first two covariates (temp, depth) are observed, while the second two (chla, O2) are assumed to be missing. In the below example, the two missing covariates are generated from stationary spatial fields (and set to "nonstationary" in the other simulation setting...
     # Note new spatial locations and covariate values are generated for each simulated dataset
     #------------------
     num_blocks <- 25
     env_dat <- sim_covariates(n = nlocs, missing_type = "stationary", num_blocks = num_blocks) 
     env <- env_dat$data
     coordinates(env) <- ~ x+y
     env_dat$data <- NULL

     # Simulate spatial multivariate presence-absence data, with a canonical logit link function
     myfam <- binomial(link = "logit")
     simdat <- create_life(env  = env, formula_X = ~ temp + depth + chla + O2, family = myfam, spp_slopes = spp_slopes)     


     #------------------
     # Split into training and test
     # Survey region is split into 25 squares in a checkerboard type pattern. Then 80% i.e., 20 of the 25 blocks are randomly chosen as training blocks, with the remaining five chosen as test blocks. Note that because the spatial locations are simulated uniformly in the survey region, then on average we expect 0.8*n of the observations to be in the training dataset, with some variation across simulated datasets
     # The splitting is also done to ensure there are "sufficient" presences in both training and test blocks
     #------------------
     training_blocks <- sample(1:num_blocks, ceiling(0.8*num_blocks)) %>%
          sort
     test_blocks <- (1:num_blocks)[-training_blocks]
     simdat_train <- list(resp = simdat$resp[env@data$block_number %in% training_blocks,], eta = simdat$eta[env@data$block_number %in% training_blocks,])
     simdat_test <- list(resp = simdat$resp[env@data$block_number %in% test_blocks,], eta = simdat$eta[env@data$block_number %in% test_blocks,])

     while(!all(colSums(simdat_train$resp>0) > 20) || !all(colSums(simdat_test$resp>0) > 10)) {
          training_blocks <- sample(1:num_blocks, ceiling(0.8*num_blocks)) %>%
               sort
          test_blocks <- (1:num_blocks)[-training_blocks]
          simdat_train <- list(resp = simdat$resp[env@data$block_number %in% training_blocks,], eta = simdat$eta[env@data$block_number %in% training_blocks,])
          simdat_test <- list(resp = simdat$resp[env@data$block_number %in% test_blocks,], eta = simdat$eta[env@data$block_number %in% test_blocks,])
          }

     rm(simdat)     
     n_train <- nrow(simdat_train$resp)
     n_test <- nrow(simdat_test$resp)
     
     env_train <- cbind(env@coords, env@data)[env@data$block_number %in% training_blocks,] 
     rownames(env_train) <- 1:nrow(env_train)
     coordinates(env_train) <- ~ x+y 

     env_test <- cbind(env@coords, env@data)[env@data$block_number %in% test_blocks,] 
     rownames(env_test) <- nrow(env_train) + (1:nrow(env_test))
     coordinates(env_test) <- ~ x+y 
     rm(env)

     # Form some quantities needed for assessing performance later on.
     rownames(simdat_test$eta) <- paste0("unit", 1:nrow(simdat_test$eta))
     test_eta <-  simdat_test$eta %>%
          as.data.frame %>%
          rownames_to_column(var = "unit") %>%
          pivot_longer(spp1:spp30) 

     # Construct linear predictor corresponding to missing component
     missing_linear_predictors <- list(true = tcrossprod(model.matrix(~ chla + O2 - 1, data = env_train@data), spp_slopes[,-c(1:2)]) %>% cov, 
          betabetaT = tcrossprod(spp_slopes[,-c(1:2)]))


     #------------------
     # Method 1: Stacked SDMs 
     #------------------
     tic <- proc.time()
     fit_stackedsdm <- lapply(1:num_spp, function(j) glm(resp ~ temp + depth, data = data.frame(resp = as.vector(simdat_train$resp[,j]), env_train@data), 
          family = myfam))
     toc <- proc.time()
     computing_time <- data.frame(stackedsdm = (toc-tic)[3])
     
     # Comparing true and estimated coefficients
     spp_slopes_long$stackedsdm <- as.vector(sapply(fit_stackedsdm, function(x) x$coefficients[-1]))      
     
     # Compute confidence interval coverage and width
     spp_slopes_ciwidth <- spp_slopes_coverage <- spp_slopes_long
     cw_cis_stackedsdm <- lapply(fit_stackedsdm, function(x) confint.default(x)[1+(1:(length(fit_stackedsdm[[1]]$coefficients)-1)),,drop=FALSE]) %>%
          abind(along = 1) %>%
          as.data.frame %>%
          rownames_to_column(var = "covariate") %>%
          mutate(species = rep(paste0("spp",1:num_spp), each = length(fit_stackedsdm[[1]]$coefficients)-1))
     colnames(cw_cis_stackedsdm)[2:3] <- c("lower","upper")
     spp_slopes_ciwidth$stackedsdm<- cw_cis_stackedsdm$upper - cw_cis_stackedsdm$lower
     spp_slopes_coverage$stackedsdm <- as.numeric(cw_cis_stackedsdm$lower < spp_slopes_long$value & cw_cis_stackedsdm$upper > spp_slopes_long$value)     

     # Out of sample performance measures
     getpreds <- sapply(fit_stackedsdm, function(x) predict(x, newdata = env_test@data, type = "response"))
     brier_score <- data.frame(stackedsdm = colMeans((simdat_test$resp - getpreds)^2))
     if(myfam$family == "binomial") {
          tjur_r2 <- data.frame(stackedsdm = sapply(1:num_spp, function(j) { 
               m1 <- getpreds[which(simdat_test$resp[,j] == 1),j] %>%
                    mean(na.rm = TRUE)
               m0 <- getpreds[which(simdat_test$resp[,j] == 0),j] %>%
                    mean(na.rm = TRUE)
               m1 - m0     
               }))
          log_score <- data.frame(stackedsdm = sapply(1:num_spp, function(j) { 
                    sum(dbinom(simdat_test$resp[,j], size = 1 , prob = getpreds[,j], log = TRUE))
                    })
               )
          aucs <- data.frame(stackedsdm = sapply(1:num_spp, function(j) { 
                    pred <- ROCR::prediction(getpreds[,j], labels = simdat_test$resp[,j]) %>%
                         ROCR::performance(measure = "auc")
                    pred@y.values[[1]]
                    })
               )
          all_predictions <- list(stackedsdm = getpreds)
          }
#      if(myfam$family == "poisson") {
#           log_score <- data.frame(stackedsdm = sapply(1:num_spp, function(j) { 
#                     sum(dpois(simdat_test$resp[,j], lambda = getpreds[,j], log = TRUE))
#                     })
#                )
#           pseudoR2 <- data.frame(stackedsdm = sapply(1:num_spp, function(j) {
#                     out <- cor(simdat_test$resp[,j], getpreds[,j], method = "spearman")
#                     out^2*sign(out)
#                     })
#                )
#           }
     rm(getpreds)
     

     #------------------
     # Method 2: HMSC: Hierarchical Modelling of Species Communities with Nearest Neighbour Gaussian Process. 
     # The code below is adapted from <https://cran.r-project.org/web/packages/Hmsc/vignettes/vignette_4_spatial.pdf>, with additional help from personal communication from Gleb Tikhonov
     #------------------
     studyDesign = data.frame(sample = factor(rownames(env_train@coords)))
     rL.nngp = HmscRandomLevel(sData = env_train@coords, sMethod = 'NNGP', nNeighbours = 20)
     rL.nngp = setPriors(rL.nngp, nfMin = 2, nfMax = 2) # Note Hmsc scales really poorly with increasing number of latent variables!
     env_trainX <- model.matrix(~ temp + depth, data = env_train@data)
     m.nngp = Hmsc(Y = simdat_train$resp, X = env_trainX, XScale = FALSE, studyDesign = studyDesign, ranLevels = list("sample" = rL.nngp), distr = "probit")

     tic <- proc.time()
     m.nngp = sampleMcmc(m.nngp, thin = 10, samples = 1000, transient = 1000, nChains = 3, updater = list(GammaEta=FALSE), nParallel = 3) 
     toc <- proc.time()
     computing_time$hmsc_nngp <- (toc-tic)[3]

     # Comparing true and estimated coefficients -- Note due to misspecification of the link function, then assessment of performance on coefficients is not performed.
     est_coefs <- getPostEstimate(m.nngp, parName = "Beta")$mean %>% 
          t
     spp_slopes_long$hmsc_nngp<- est_coefs[,-1] %>%
          t %>%
          as.vector

     # Compute confidence interval coverage and width -- Note due to misspecification of the link function, then assessment of performance on coefficients is not performed.
     cw_cis_hmscnngp <- getPostEstimate(m.nngp, parName = "Beta", q = c(0.025,0.975))$q[,2:3,,drop=FALSE] 
     dimnames(cw_cis_hmscnngp) <- list(quantile = c("lower","upper"), covariate = spp_slopes_long$name[1:2], species = paste0("spp",1:num_spp))
     cw_cis_hmscnngp <- as.data.frame.table(cw_cis_hmscnngp) %>%
          pivot_wider(names_from = quantile, values_from = Freq)
     spp_slopes_ciwidth$hmsc_nngp <- cw_cis_hmscnngp$upper - cw_cis_hmscnngp$lower
     spp_slopes_coverage$hmsc_nngp <- as.numeric(cw_cis_hmscnngp$lower < spp_slopes_long$value & cw_cis_hmscnngp$upper > spp_slopes_long$value)     

     # Out of sample performance measures -- based on posterior mean
     rL.nngp_test = HmscRandomLevel(sData = rbind(env_train@coords, env_test@coords), sMethod = 'NNGP', nNeighbours = 20)
     test_studyDesign = data.frame(sample = factor(rownames(env_test@coords)))
     env_testX <- model.matrix(~ temp + depth, data = env_test@data)
     getpreds <- predict(m.nngp, X = env_testX, studyDesign = test_studyDesign, ranLevels = list("sample" = rL.nngp_test), expected = TRUE) %>%
          abind(along = 3) %>%
          apply(., c(1,2), mean, na.rm = TRUE)
     brier_score$hmsc_nngp <- colMeans((simdat_test$resp - getpreds)^2)
     if(myfam$family == "binomial") {     
          tjur_r2$hmsc_nngp <- sapply(1:num_spp, function(j) { 
               m1 <- getpreds[which(simdat_test$resp[,j] == 1),j] %>%
                    mean(na.rm = TRUE)
               m0 <- getpreds[which(simdat_test$resp[,j] == 0),j] %>%
                    mean(na.rm = TRUE)
               m1 - m0     
               })
          log_score$hmsc_nngp <- sapply(1:num_spp, function(j) { 
               sum(dbinom(simdat_test$resp[,j], size = 1 , prob = getpreds[,j], log = TRUE))
               })
          aucs$hmsc_nngp <- sapply(1:num_spp, function(j) { 
               pred <- ROCR::prediction(getpreds[,j], labels = simdat_test$resp[,j]) %>%
                    ROCR::performance(measure = "auc")
               pred@y.values[[1]]
               })
          all_predictions$hmsc_nngp <- getpreds
          }
#      if(myfam$family == "poisson") {
#           log_score$hmsc_nngp <- sapply(1:num_spp, function(j) { 
#                sum(dbinom(simdat_test$resp[,j], size = 1 , prob = getpreds[,j], log = TRUE))
#                })
#           pseudoR2$hmsc_nngp <- sapply(1:num_spp, function(j) {
#                out <- cor(simdat_test$resp[,j], getpreds[,j], method = "spearman")
#                out^2*sign(out)
#                })
#           }
     rm(getpreds, est_coefs)

     missing_linear_predictors$hmsc_nngp_Associations <- computeAssociations(m.nngp)[[1]]$meanp
     gc()
     

     #------------------
     # Method 3: HMSC: Hierarchical Modelling of Species Communities with Gaussian predictive process
     # The code below is adapted from <https://cran.r-project.org/web/packages/Hmsc/vignettes/vignette_4_spatial.pdf>, with additional help from personal communication from Gleb Tikhonov
     #------------------
     studyDesign = data.frame(sample = factor(rownames(env_train@coords)))
     Knots = constructKnots(env_train@coords, knotDist = 0.5, minKnotDist = 1) #nKnots = nrow(env_train@coords)*0.5
     rL.gpp = HmscRandomLevel(sData = env_train@coords, sMethod = 'GPP', sKnot = Knots)
     rL.gpp = setPriors(rL.gpp, nfMin = 2, nfMax = 2) # Note Hmsc scales really poorly with increasing number of latent variables!
     env_trainX <- model.matrix(~ temp + depth, data = env_train@data)
     m.gpp = Hmsc(Y = simdat_train$resp, X= env_trainX, XScale = FALSE, studyDesign = studyDesign, ranLevels = list("sample" = rL.gpp), distr = "probit")
     
     tic <- proc.time()
     m.gpp = sampleMcmc(m.gpp, thin = 10, samples = 1000, transient = 1000, nChains = 3, updater = list(GammaEta=FALSE), nParallel = 3) 
     toc <- proc.time()
     computing_time$hmsc_gpp <- (toc-tic)[3]

     # Comparing true and estimated coefficients -- Note due to misspecification of the link function, then assessment of performance on coefficients is not performed.
     est_coefs <- getPostEstimate(m.gpp, parName = "Beta")$mean %>% 
          t
     spp_slopes_long$hmsc_gpp <- est_coefs[,-1] %>%
          t %>%
          as.vector

     # Compute confidence interval coverage and width -- Note due to misspecification of the link function, then assessment of performance on coefficients is not performed.
     cw_cis_hmscgpp <- getPostEstimate(m.gpp, parName = "Beta", q = c(0.025,0.975))$q[,2:3,,drop=FALSE] 
     dimnames(cw_cis_hmscgpp) <- list(quantile = c("lower","upper"), covariate = spp_slopes_long$name[1:2], species = paste0("spp",1:num_spp))
     cw_cis_hmscgpp <- as.data.frame.table(cw_cis_hmscgpp) %>%
          pivot_wider(names_from = quantile, values_from = Freq)
     spp_slopes_ciwidth$hmsc_gpp <- cw_cis_hmscgpp$upper - cw_cis_hmscgpp$lower
     spp_slopes_coverage$hmsc_gpp <- as.numeric(cw_cis_hmscgpp$lower < spp_slopes_long$value & cw_cis_hmscgpp$upper > spp_slopes_long$value)     

     # Out of sample performance measures
     Knots = constructKnots(rbind(env_train@coords, env_test@coords), knotDist = 0.5, minKnotDist = 1) #nKnots = nrow(env_train@coords)*0.5
     rL.gpp_test = HmscRandomLevel(sData = rbind(env_train@coords, env_test@coords), sMethod = 'GPP', sKnot = Knots)
     test_studyDesign = data.frame(sample = factor(rownames(env_test@coords)))
     env_testX <- model.matrix(~ temp + depth, data = env_test@data)
     getpreds <- predict(m.gpp, X = env_testX, studyDesign = test_studyDesign, ranLevels = list("sample" = rL.nngp_test), expected = TRUE) %>%
          abind(along = 3) %>%
          apply(., c(1,2), mean, na.rm = TRUE)
     brier_score$hmsc_gpp <- colMeans((simdat_test$resp - getpreds)^2)
     if(myfam$family == "binomial") {     
          tjur_r2$hmsc_gpp <- sapply(1:num_spp, function(j) { 
               m1 <- getpreds[which(simdat_test$resp[,j] == 1),j] %>%
                    mean(na.rm = TRUE)
               m0 <- getpreds[which(simdat_test$resp[,j] == 0),j] %>%
                    mean(na.rm = TRUE)
               m1 - m0     
               })
          log_score$hmsc_gpp <- sapply(1:num_spp, function(j) { 
               sum(dbinom(simdat_test$resp[,j], size = 1 , prob = getpreds[,j], log = TRUE))
               })
          aucs$hmsc_gpp <- sapply(1:num_spp, function(j) { 
               pred <- ROCR::prediction(getpreds[,j], labels = simdat_test$resp[,j]) %>%
                    ROCR::performance(measure = "auc")
               pred@y.values[[1]]
               })
          all_predictions$hmsc_gpp <- getpreds
          }
#      if(myfam$family == "poisson") {
#           log_score$hmsc_gpp <- sapply(1:num_spp, function(j) { 
#                sum(dbinom(simdat_test$resp[,j], size = 1 , prob = getpreds[,j], log = TRUE))
#                })
#           pseudoR2$hmsc_nngp <- sapply(1:num_spp, function(j) {
#                out <- cor(simdat_test$resp[,j], getpreds[,j], method = "spearman")
#                out^2*sign(out)
#                })
#           }
     rm(getpreds, est_coefs)

     missing_linear_predictors$hmsc_gpp_Associations <- computeAssociations(m.gpp)[[1]]$meanp
     gc()

     
     #------------------
     # Method 4: sjSDM
     # This is the method of Pichler and Hartig (2020), <https://arxiv.org/abs/2003.05331>, including a spatial component that is modelled via neural network. 
     # devtools::install_github("https://github.com/TheoreticalEcology/s-jSDM", subdir = "sjSDM")
     # sjSDM::install_sjSDM(version = "cpu")
     #------------------
     tic <- proc.time()
     fit_sjsdm <- sjSDM(Y = simdat_train$resp, env = linear(data = env_train@data, formula = ~ temp + depth), 
          spatial = DNN(env_train@coords, hidden = c(10L, 10L, 10L), ~0+.), se = TRUE, family = myfam, parallel = detectCores()-2)
     toc <- proc.time()
     computing_time$sjsdm <- (toc-tic)[3]

     # Comparing true and estimated coefficients
     spp_slopes_long$sjsdm <- coef(fit_sjsdm)$env[[1]][,-1] %>%
          t %>%
          as.vector
          
     # Compute confidence interval coverage and width
     s <- summary(fit_sjsdm)$coefmat %>%
          as.data.frame %>%
          mutate(covariate = rep(c("Intercept",spp_slopes_long$name[1:2]), num_spp)) %>%
          dplyr::filter(covariate != "Intercept") %>%
          mutate(species = rep(paste0("spp",1:num_spp), each = 2))
     cw_cis_sjsdm <- s %>%
          mutate(lower = Estimate - qnorm(0.975)*Std.Err, upper = Estimate + qnorm(0.975)*Std.Err) %>%
          dplyr::select(covariate:upper)
     rm(s)
     spp_slopes_ciwidth$sjsdm <- cw_cis_sjsdm$upper - cw_cis_sjsdm$lower
     spp_slopes_coverage$sjsdm <- as.numeric(cw_cis_sjsdm$lower < spp_slopes_long$value & cw_cis_sjsdm$upper > spp_slopes_long$value)     

     # Out of sample performance measures
     getpreds <- predict(fit_sjsdm, newdata = env_test@data, SP = env_test@coords)
     brier_score$sjsdm <- colMeans((simdat_test$resp - getpreds)^2)
     if(myfam$family == "binomial") {
          tjur_r2$sjsdm <- sapply(1:num_spp, function(j) { 
                    m1 <- getpreds[which(simdat_test$resp[,j] > 0),j] %>%
                         mean(na.rm = TRUE)
                    m0 <- getpreds[which(simdat_test$resp[,j] == 0),j] %>%
                         mean(na.rm = TRUE)
                    m1 - m0     
                    })
          log_score$sjsdm <- sapply(1:num_spp, function(j) { 
                    sum(dbinom(simdat_test$resp[,j], size = 1, prob = getpreds[,j], log = TRUE))
                    })
          aucs$sjsdm <- sapply(1:num_spp, function(j) { 
                    pred <- ROCR::prediction(getpreds[,j], labels = simdat_test$resp[,j]) %>%
                         ROCR::performance(measure = "auc")
                    pred@y.values[[1]]
                    })
          all_predictions$sjsdm <- getpreds
          }
#      if(myfam$family == "poisson") {
#           log_score$sjsdm <- sapply(1:num_spp, function(j) { 
#                     sum(dpois(simdat_test$resp[,j], lambda = getpreds[,j], log = TRUE))
#                     })
#           pseudoR2$sjsdm <- sapply(1:num_spp, function(j) {
#                out <- cor(simdat_test$resp[,j], getpreds[,j], method = "spearman")
#                out^2*sign(out)
#                })
#           }
     rm(getpreds)

     missing_linear_predictors$sjsdm <- getCov(fit_sjsdm)

     
     #------------------
     # Method 5: Community Basis Function model (CBFM) 
     #------------------
     # 25 basis functions, rank of G = 5, rank of Sigma = 5
     num_spbasisfunctions <- 25 
     sp_basisfunctions <- mrts(env_train@coords, num_spbasisfunctions) %>% 
          as.matrix %>%
          as(., "sparseMatrix") %>%
          {.[,-(1)]} # Remove intercept 
     test_basisfunctions <- mrts(env_train@coords, num_spbasisfunctions) %>% 
          predict(newx = env_test@coords) %>% 
          as.matrix %>%
          as(., "sparseMatrix") %>%
          {.[,-(1)]} # Remove intercept 
     
     
     tic <- proc.time()     
     CBFMmod <- CBFM(y = simdat_train$resp, formula = ~ temp + depth, data = env_train@data, B_space = as.matrix(sp_basisfunctions), 
                     ncores = detectCores()-2, family = myfam, G_control = list(rank = 5), Sigma_control = list(rank = 5))
     toc <- proc.time()     
     computing_time$cbfm255 <- (toc-tic)[3]
     
     # Comparing true and estimated coefficients
     spp_slopes_long$cbfm255 <- CBFMmod$beta[,-1] %>%
          t %>%
          as.vector
     
     # Compute confidence interval coverage and width
     cw_cis_cbfm <- lapply(summary(CBFMmod)$summary_tables, function(x) x$parametric_coefs[-1,]) %>% 
          do.call(rbind, .) %>%
          dplyr::mutate(Species = rep(1:num_spp, each = 2), Predictor = rep(c("Temp", "Depth"), num_spp)) 
     spp_slopes_ciwidth$cbfm255 <- cw_cis_cbfm$Upper - cw_cis_cbfm$Lower
     spp_slopes_coverage$cbfm255 <- as.numeric(cw_cis_cbfm$Lower < spp_slopes_long$value & cw_cis_cbfm$Upper > spp_slopes_long$value)     

     # Out of sample performance measures
     getpreds <- predict(CBFMmod, newdata = env_test@data, new_B_space = test_basisfunctions, type = "response", se_fit = FALSE)
     brier_score$cbfm255 <- colMeans((simdat_test$resp - getpreds)^2)
     if(myfam$family == "binomial") {
          tjur_r2$cbfm255 <- sapply(1:num_spp, function(j) { 
                    m1 <- getpreds[,j][which(simdat_test$resp[,j] == 1)] %>%
                         mean(na.rm = TRUE)
                    m0 <- getpreds[,j][which(simdat_test$resp[,j] == 0)] %>%
                         mean(na.rm = TRUE)
                    m1 - m0     
                    })
          log_score$cbfm255 <- sapply(1:num_spp, function(j) { 
                    sum(dbinom(simdat_test$resp[,j], size = 1, prob = getpreds[,j], log = TRUE))
                    })
          aucs$cbfm255 <- sapply(1:num_spp, function(j) { 
                    pred <- ROCR::prediction(getpreds[,j], labels = simdat_test$resp[,j]) %>%
                         ROCR::performance(measure = "auc")
                    pred@y.values[[1]]
                    })
          all_predictions$cbfm255 <- getpreds
          }
#      if(myfam$family == "poisson") {
#           log_score$cbfm255 <- sapply(1:num_spp, function(j) { 
#                     sum(dpois(simdat_test$resp[,j], lambda = getpreds[,j]), log = TRUE)
#                     })
#           pseudoR2$cbfm255 <- sapply(1:num_spp, function(j) {
#                out <- cor(simdat_test$resp[,j], getpreds[,j], method = "spearman")
#                out^2*sign(out)
#                })
#           }
     rm(getpreds)     

     missing_linear_predictors$cbfm255_G <- CBFMmod$G_space
     missing_linear_predictors$cbfm255_corB <- corB(CBFMmod)

     
     # 25 basis functions, rank of G = 5, rank of Sigma = 10
     tic <- proc.time()     
     CBFMmod <- CBFM(y = simdat_train$resp, formula = ~ temp + depth, data = env_train@data, B_space = as.matrix(sp_basisfunctions), 
                     ncores = detectCores()-2, family = myfam, G_control = list(rank = 5), Sigma_control = list(rank = 10))
     toc <- proc.time()     
     computing_time$cbfm2510 <- (toc-tic)[3]
     
     # Comparing true and estimated coefficients
     spp_slopes_long$cbfm2510 <- CBFMmod$beta[,-1] %>%
          t %>%
          as.vector
     
     # Compute confidence interval coverage and width
     cw_cis_cbfm <- lapply(summary(CBFMmod)$summary_tables, function(x) x$parametric_coefs[-1,]) %>% 
          do.call(rbind, .) %>%
          dplyr::mutate(Species = rep(1:num_spp, each = 2), Predictor = rep(c("Temp", "Depth"), num_spp)) 
     spp_slopes_ciwidth$cbfm2510 <- cw_cis_cbfm$Upper - cw_cis_cbfm$Lower
     spp_slopes_coverage$cbfm2510 <- as.numeric(cw_cis_cbfm$Lower < spp_slopes_long$value & cw_cis_cbfm$Upper > spp_slopes_long$value)     

     # Out of sample performance measures
     getpreds <- predict(CBFMmod, newdata = env_test@data, new_B_space = test_basisfunctions, type = "response", se_fit = FALSE)
     brier_score$cbfm2510 <- colMeans((simdat_test$resp - getpreds)^2)
     if(myfam$family == "binomial") {
          tjur_r2$cbfm2510 <- sapply(1:num_spp, function(j) { 
                    m1 <- getpreds[,j][which(simdat_test$resp[,j] == 1)] %>%
                         mean(na.rm = TRUE)
                    m0 <- getpreds[,j][which(simdat_test$resp[,j] == 0)] %>%
                         mean(na.rm = TRUE)
                    m1 - m0     
                    })
          log_score$cbfm2510 <- sapply(1:num_spp, function(j) { 
                    sum(dbinom(simdat_test$resp[,j], size = 1, prob = getpreds[,j], log = TRUE))
                    })
          aucs$cbfm2510 <- sapply(1:num_spp, function(j) { 
                    pred <- ROCR::prediction(getpreds[,j], labels = simdat_test$resp[,j]) %>%
                         ROCR::performance(measure = "auc")
                    pred@y.values[[1]]
                    })
          all_predictions$cbfm2510 <- getpreds
          }
#      if(myfam$family == "poisson") {
#           log_score$cbfm2510 <- sapply(1:num_spp, function(j) { 
#                     sum(dpois(simdat_test$resp[,j], lambda = getpreds[,j]), log = TRUE)
#                     })
#           pseudoR2$cbfm2510 <- sapply(1:num_spp, function(j) {
#                out <- cor(simdat_test$resp[,j], getpreds[,j], method = "spearman")
#                out^2*sign(out)
#                })
#           }
     rm(getpreds)     

     missing_linear_predictors$cbfm2510_G <- CBFMmod$G_space
     missing_linear_predictors$cbfm2510_corB <- corB(CBFMmod)

     
     #-----------------------
     # Finish and save results
     #-----------------------          
     results <- list(
          slopes_estimates = spp_slopes_long, 
          slopes_CIwidth = spp_slopes_ciwidth, 
          slopes_CIcoverage = spp_slopes_coverage, 
          tjur_r2 = tjur_r2, 
          aucs = aucs, 
          log_score = log_score, 
          brier_score = brier_score, 
          time = computing_time,
          all_predictions = all_predictions, 
          test_set_responses = simdat_test$resp,
          missing_linear_predictors = missing_linear_predictors,
          m_nngp = m.nngp, 
          m_gpp = m.gpp
          )
     
     
     # Save as RDS file (this can be changed depending on user's desire)
     filename <- paste0("simstudy_PA500training_stationary_dataset",i,".rds")
     saveRDS(results, file = filename)
     rm(results, m.nngp, m.gpp, 
        cw_cis_cbfm, cw_cis_hmscgpp, cw_cis_hmscnngp, cw_cis_sjsdm, cw_cis_sjsdm, cw_cis_stackedsdm,
        missing_linear_predictors, spp_slopes_ciwidth, spp_slopes_coverage, 
        tjur_r2, aucs, log_score, brier_score, all_predictions, 
        simdat_test, simdat_train, env_test, env_train, CBFMmod, fit_stackedsdm, test_eta, 
        training_blocks, test_blocks)
     }




#-----------------------
# Versions of packages used
#-----------------------          
sessionInfo()
## Results of sessionInfo() from machine running CBFM simulation study, executed on December 13 2022 ##

# R version 4.2.2 Patched (2022-11-10 r83330)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 20.04.5 LTS
# 
# Matrix products: default
# BLAS/LAPACK: /opt/intel/compilers_and_libraries_2020.4.304/linux/mkl/lib/intel64_lin/libmkl_rt.so
# 
# locale:
#  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
#  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
#  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
#  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
#  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
# [1] parallel  stats     graphics  grDevices utils     datasets  methods  
# [8] base     
# 
# other attached packages:
#  [1] RandomFields_3.3.14     RandomFieldsUtils_1.2.5 sp_1.5-1               
#  [4] scales_1.2.1            mvtnorm_1.1-3           CBFM_0.1               
#  [7] TMB_1.9.1               Hmsc_3.0-13             coda_0.19-4            
# [10] sjSDM_1.0.3             ROCR_1.0-11             abind_1.4-5            
# [13] forcats_0.5.2           stringr_1.4.1           dplyr_1.0.10           
# [16] purrr_0.3.5             readr_2.1.3             tidyr_1.2.1            
# [19] tibble_3.1.8            ggplot2_3.4.0           tidyverse_1.3.2        
# [22] mvabund_4.2.1           doParallel_1.0.17       iterators_1.0.14       
# [25] foreach_1.5.2           here_1.0.1              autoFRK_1.4.3          
# [28] spam_2.9-1             
# 
# loaded via a namespace (and not attached):
#   [1] filematrix_1.3       googledrive_2.0.0    colorspace_2.0-3    
#   [4] gamlss.tr_5.1-7      ellipsis_0.3.2       LatticeKrig_8.4     
#   [7] tdigest_0.4.1        rprojroot_2.0.3      fs_1.5.2            
#  [10] rstudioapi_0.14      MatrixModels_0.5-1   bit64_4.0.5         
#  [13] fansi_1.0.3          filehashSQLite_0.2-6 lubridate_1.9.0     
#  [16] mathjaxr_1.6-0       xml2_1.3.3           codetools_0.2-18    
#  [19] splines_4.2.2        cachem_1.0.6         jsonlite_1.8.3      
#  [22] pROC_1.18.0          mcmc_0.9-7           broom_1.0.1         
#  [25] dbplyr_2.2.1         png_0.1-7            gamlss.dist_6.0-5   
#  [28] compiler_4.2.2       httr_1.4.4           backports_1.4.1     
#  [31] assertthat_0.2.1     Matrix_1.5-1         fastmap_1.1.0       
#  [34] gargle_1.2.1         cli_3.4.1            quantreg_5.94       
#  [37] tools_4.2.2          dotCall64_1.0-2      gtable_0.3.1        
#  [40] glue_1.6.2           maps_3.4.1           rappdirs_0.3.3      
#  [43] Rcpp_1.0.9           cellranger_1.1.0     vctrs_0.5.1         
#  [46] filehash_2.4-3       ape_5.6-2            nlme_3.1-160        
#  [49] rvest_1.0.3          timechange_0.1.1     lifecycle_1.0.3     
#  [52] statmod_1.4.37       googlesheets4_1.0.1  MASS_7.3-58.1       
#  [55] hms_1.1.2            gamlss.data_6.0-2    SparseM_1.81        
#  [58] fields_14.1          mvnfast_0.2.7        memoise_2.0.1       
#  [61] reticulate_1.26      gridExtra_2.3        stringi_1.7.8       
#  [64] RSQLite_2.2.18       gratia_0.7.3         checkmate_2.1.0     
#  [67] truncnorm_1.0-8      rlang_1.0.6          pkgconfig_2.0.3     
#  [70] matrixStats_0.63.0   pracma_2.4.2         lattice_0.20-45     
#  [73] patchwork_1.1.2      bit_4.0.5            tidyselect_1.2.0    
#  [76] plyr_1.8.8           magrittr_2.0.3       R6_2.5.1            
#  [79] generics_0.1.3       BayesLogit_2.1       DBI_1.1.3           
#  [82] pillar_1.8.1         haven_2.5.1          withr_2.5.0         
#  [85] mgcv_1.8-41          survival_3.4-0       nnet_7.3-18         
#  [88] modelr_0.1.10        crayon_1.5.2         gamlss_5.4-10       
#  [91] utf8_1.2.2           tzdb_0.3.0           viridis_0.6.2       
#  [94] grid_4.2.2           readxl_1.4.1         blob_1.2.3          
#  [97] FNN_1.1.3.1          reprex_2.0.2         numDeriv_2016.8-1.1 
# [100] MCMCpack_1.6-3       munsell_0.5.0        viridisLite_0.4.1   
# [103] tweedie_2.3.5 
