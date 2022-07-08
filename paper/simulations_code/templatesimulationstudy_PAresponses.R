##--------------------
## Template code for simulation study -- presence-absence responses. The script can be straightforwardly modified to work for Poisson responses, different sample sizes, 
## Spatial multivariate abundance data are generated from a stacked species distribution model
##--------------------
here::i_am("desired_directory/simulationstudy_presenceabsenceresponses.R") # install.packages("here")

rm(list = ls())
library(autoFRK)
library(here)
library(doParallel)
library(foreach)
library(mvabund)
library(tidyverse)
library(mvtnorm)
library(scales)
library(RandomFields)
library(sp)
library(abind)
library(ROCR)
library(sjSDM)
library(Hmsc)
library(CBFM)
source(here("simulatefunctions.R"))


# Number of simulated datasets 
nsims <- 200

# Number of spatial locations
nlocs <- 625


# Create list to store results
result_list <- result_CIs <- vector(mode = "list", length = nsims)

# Parameters in true model
set.seed(2021)
num_spp <- 30
candidate_coefs <- data.frame(
     big_coefs = runif(num_spp, -1.5, 1.5), 
     small_coefs = runif(num_spp, -0.5, 0.5), 
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
     # Assume that only the first two covariates (temp, depth) are observed, while the second two (chla, O2) are assumed to be missing. In the below example, the two missing covariates are generated from stationary spatial fields
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

     #' Form some quantities needed for assessing performance later on.
     rownames(simdat_test$eta) <- paste0("unit", 1:nrow(simdat_test$eta))
     test_eta <-  simdat_test$eta %>%
          as.data.frame %>%
          rownames_to_column(var = "unit") %>%
          pivot_longer(spp1:spp30) 



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
     model_error <- data.frame(stackedsdm = sqrt(colMeans((simdat_test$resp - getpreds)^2)))
     if(myfam$family == "binomial") {
          tjur_r2 <- data.frame(stackedsdm = sapply(1:num_spp, function(j) { 
               m1 <- getpreds[which(simdat_test$resp[,j] == 1),j] %>%
                    mean(na.rm = TRUE)
               m0 <- getpreds[which(simdat_test$resp[,j] == 0),j] %>%
                    mean(na.rm = TRUE)
               m1 - m0     
               }))
          pred_logL <- data.frame(stackedsdm = sapply(1:num_spp, function(j) { 
                    sum(dbinom(simdat_test$resp[,j], size = 1 , prob = getpreds[,j], log = TRUE))
                    })
               )
          aucs <- data.frame(stackedsdm = sapply(1:num_spp, function(j) { 
                    pred <- ROCR::prediction(getpreds[,j], labels = simdat_test$resp[,j]) %>%
                         ROCR::performance(measure = "auc")
                    pred@y.values[[1]]
                    })
               )
          }
#      if(myfam$family == "poisson") {
#           pred_logL <- data.frame(stackedsdm = sapply(1:num_spp, function(j) { 
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
     rL.nngp = setPriors(rL.nngp, nfMin = 2, nfMax = 2) 
     m.nngp = Hmsc(Y = simdat_train$resp, XData = env_train@data, XScale = FALSE, XFormula = ~ temp + depth, 
          studyDesign = studyDesign, ranLevels = list("sample" = rL.nngp), distr = "probit")

     tic <- proc.time()
     m.nngp = sampleMcmc(m.nngp, thin = 10, samples = 1000, transient = 1000, nChains = 3, nParallel = 3)  ## Need to bump this up!
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
     getpreds <- predict(m.nngp, XData = env_test@data, studyDesign = test_studyDesign, ranLevels = list("sample" = rL.nngp_test), expected = TRUE) %>%
          abind(along = 3) %>%
          apply(., c(1,2), mean, na.rm = TRUE)
     model_error$hmsc_nngp <- sqrt(colMeans((simdat_test$resp - getpreds)^2))
     if(myfam$family == "binomial") {     
          tjur_r2$hmsc_nngp <- sapply(1:num_spp, function(j) { 
               m1 <- getpreds[which(simdat_test$resp[,j] == 1),j] %>%
                    mean(na.rm = TRUE)
               m0 <- getpreds[which(simdat_test$resp[,j] == 0),j] %>%
                    mean(na.rm = TRUE)
               m1 - m0     
               })
          pred_logL$hmsc_nngp <- sapply(1:num_spp, function(j) { 
               sum(dbinom(simdat_test$resp[,j], size = 1 , prob = getpreds[,j], log = TRUE))
               })
          aucs$hmsc_nngp <- sapply(1:num_spp, function(j) { 
               pred <- ROCR::prediction(getpreds[,j], labels = simdat_test$resp[,j]) %>%
                    ROCR::performance(measure = "auc")
               pred@y.values[[1]]
               })
          }
#      if(myfam$family == "poisson") {
#           pred_logL$hmsc_nngp <- sapply(1:num_spp, function(j) { 
#                sum(dbinom(simdat_test$resp[,j], size = 1 , prob = getpreds[,j], log = TRUE))
#                })
#           pseudoR2$hmsc_nngp <- sapply(1:num_spp, function(j) {
#                out <- cor(simdat_test$resp[,j], getpreds[,j], method = "spearman")
#                out^2*sign(out)
#                })
#           }
     rm(getpreds, est_coefs)

     

     #------------------
     # Method 3: HMSC: Hierarchical Modelling of Species Communities with Gaussian predictive process
     # The code below is adapted from <https://cran.r-project.org/web/packages/Hmsc/vignettes/vignette_4_spatial.pdf>, with additional help from personal communication from Gleb Tikhonov
     #------------------
     studyDesign = data.frame(sample = factor(rownames(env_train@coords)))
     Knots = constructKnots(env_train@coords, knotDist = 0.5, minKnotDist = 1) 
     rL.gpp = HmscRandomLevel(sData = env_train@coords, sMethod = 'GPP', sKnot = Knots)
     rL.gpp = setPriors(rL.gpp, nfMin = 2, nfMax = 2) 
     m.gpp = Hmsc(Y = simdat_train$resp, XData = env_train@data, XScale = FALSE, XFormula = ~ temp + depth, 
          studyDesign = studyDesign, ranLevels = list("sample" = rL.gpp), distr = "probit")
     
     tic <- proc.time()
     m.gpp = sampleMcmc(m.gpp, thin = 10, samples = 1000, transient = 1000, nChains = 3, nParallel = 3) ## Should bump this up
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
     getpreds <- predict(m.gpp, XData = env_test@data, studyDesign = test_studyDesign, ranLevels = list("sample" = rL.gpp_test), expected = TRUE) %>%
          abind(along = 3) %>%
          apply(., c(1,2), mean, na.rm = TRUE)
     model_error$hmsc_gpp <- sqrt(colMeans((simdat_test$resp - getpreds)^2))
     if(myfam$family == "binomial") {     
          tjur_r2$hmsc_gpp <- sapply(1:num_spp, function(j) { 
               m1 <- getpreds[which(simdat_test$resp[,j] == 1),j] %>%
                    mean(na.rm = TRUE)
               m0 <- getpreds[which(simdat_test$resp[,j] == 0),j] %>%
                    mean(na.rm = TRUE)
               m1 - m0     
               })
          pred_logL$hmsc_gpp <- sapply(1:num_spp, function(j) { 
               sum(dbinom(simdat_test$resp[,j], size = 1 , prob = getpreds[,j], log = TRUE))
               })
          aucs$hmsc_gpp <- sapply(1:num_spp, function(j) { 
               pred <- ROCR::prediction(getpreds[,j], labels = simdat_test$resp[,j]) %>%
                    ROCR::performance(measure = "auc")
               pred@y.values[[1]]
               })
          }
#      if(myfam$family == "poisson") {
#           pred_logL$hmsc_gpp <- sapply(1:num_spp, function(j) { 
#                sum(dbinom(simdat_test$resp[,j], size = 1 , prob = getpreds[,j], log = TRUE))
#                })
#           pseudoR2$hmsc_nngp <- sapply(1:num_spp, function(j) {
#                out <- cor(simdat_test$resp[,j], getpreds[,j], method = "spearman")
#                out^2*sign(out)
#                })
#           }
     rm(getpreds, est_coefs)


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
     model_error$sjsdm <- sqrt(colMeans((simdat_test$resp - getpreds)^2))
     if(myfam$family == "binomial") {
          tjur_r2$sjsdm <- sapply(1:num_spp, function(j) { 
                    m1 <- getpreds[which(simdat_test$resp[,j] > 0),j] %>%
                         mean(na.rm = TRUE)
                    m0 <- getpreds[which(simdat_test$resp[,j] == 0),j] %>%
                         mean(na.rm = TRUE)
                    m1 - m0     
                    })
          pred_logL$sjsdm <- sapply(1:num_spp, function(j) { 
                    sum(dbinom(simdat_test$resp[,j], size = 1, prob = getpreds[,j], log = TRUE))
                    })
          aucs$sjsdm <- sapply(1:num_spp, function(j) { 
                    pred <- ROCR::prediction(getpreds[,j], labels = simdat_test$resp[,j]) %>%
                         ROCR::performance(measure = "auc")
                    pred@y.values[[1]]
                    })
          }
#      if(myfam$family == "poisson") {
#           pred_logL$sjsdm <- sapply(1:num_spp, function(j) { 
#                     sum(dpois(simdat_test$resp[,j], lambda = getpreds[,j], log = TRUE))
#                     })
#           pseudoR2$sjsdm <- sapply(1:num_spp, function(j) {
#                out <- cor(simdat_test$resp[,j], getpreds[,j], method = "spearman")
#                out^2*sign(out)
#                })
#           }
     rm(getpreds)


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
     CBFMmod <- CBFM(y = simdat_train$resp, formula_X = ~ temp + depth, data = env_train@data, B_space = sp_basisfunctions, 
          ncores = detectCores()-2, family = myfam, G_control = list(rank = 5), Sigma_control = list(rank = 5))
     toc <- proc.time()     
     computing_time$cbfm255 <- (toc-tic)[3]
     

     # Comparing true and estimated coefficients
     spp_slopes_long$cbfm255 <- CBFMmod$beta[,-1] %>%
          t %>%
          as.vector
     
     # Compute confidence interval coverage and width
     cw_cis_cbfm <- summary(CBFMmod)$betas_results %>%
          as.data.frame %>%
          dplyr::filter(Predictor != "(Intercept)") %>%
          arrange(Response, Predictor)
     spp_slopes_ciwidth$cbfm255 <- cw_cis_cbfm$upper - cw_cis_cbfm$lower
     spp_slopes_coverage$cbfm255 <- as.numeric(cw_cis_cbfm$lower < spp_slopes_long$value & cw_cis_cbfm$upper > spp_slopes_long$value)     
     spp_slopes_ciwidth$cbfm255_plusinterceptci <- cw_cis_cbfm$upper_plusintercept - cw_cis_cbfm$lower_plusintercept
     spp_slopes_coverage$cbfm255_plusinterceptci <- as.numeric(cw_cis_cbfm$lower_plusintercept < spp_slopes_long$value & cw_cis_cbfm$upper_plusintercept > spp_slopes_long$value)     

     #' Out of sample performance measures
     getpreds <- predict(CBFMmod, newdata = env_test@data, newB_space = test_basisfunctions, type = "response", se_fit = FALSE)
     model_error$cbfm255 <- sqrt(colMeans((simdat_test$resp - getpreds)^2))
     if(myfam$family == "binomial") {
          tjur_r2$cbfm255 <- sapply(1:num_spp, function(j) { 
                    m1 <- getpreds[,j][which(simdat_test$resp[,j] == 1)] %>%
                         mean(na.rm = TRUE)
                    m0 <- getpreds[,j][which(simdat_test$resp[,j] == 0)] %>%
                         mean(na.rm = TRUE)
                    m1 - m0     
                    })
          pred_logL$cbfm255 <- sapply(1:num_spp, function(j) { 
                    sum(dbinom(simdat_test$resp[,j], size = 1, prob = getpreds[,j], log = TRUE))
                    })
          aucs$cbfm255 <- sapply(1:num_spp, function(j) { 
                    pred <- ROCR::prediction(getpreds[,j], labels = simdat_test$resp[,j]) %>%
                         ROCR::performance(measure = "auc")
                    pred@y.values[[1]]
                    })
          }
#      if(myfam$family == "poisson") {
#           pred_logL$cbfm255 <- sapply(1:num_spp, function(j) { 
#                     sum(dpois(simdat_test$resp[,j], lambda = getpreds[,j]), log = TRUE)
#                     })
#           pseudoR2$cbfm255 <- sapply(1:num_spp, function(j) {
#                out <- cor(simdat_test$resp[,j], getpreds[,j], method = "spearman")
#                out^2*sign(out)
#                })
#           }
     rm(getpreds)     
            

     # 25 basis functions, rank of G = 5, rank of Sigma = 10
     tic <- proc.time()     
     CBFMmod <- fit_CBFM(y = simdat_train$resp, formula_X = ~ temp + depth, data = env_train@data, B_space = sp_basisfunctions, 
          ncores = detectCores()-2, family = myfam, G_control = list(rank = 5), Sigma_control = list(rank = 10))
     toc <- proc.time()     
     computing_time$cbfm2510 <- (toc-tic)[3]
     

     #' Comparing true and estimated coefficients
     spp_slopes_long$cbfm2510 <- CBFMmod$beta[,-1] %>%
          t %>%
          as.vector
     
     # Compute confidence interval coverage and width
     cw_cis_cbfm <- summary(CBFMmod)$betas_results %>%
          as.data.frame %>%
          dplyr::filter(Predictor != "(Intercept)") %>%
          arrange(Response, Predictor)
     spp_slopes_ciwidth$cbfm2510 <- cw_cis_cbfm$upper - cw_cis_cbfm$lower
     spp_slopes_coverage$cbfm2510 <- as.numeric(cw_cis_cbfm$lower < spp_slopes_long$value & cw_cis_cbfm$upper > spp_slopes_long$value)     
     spp_slopes_ciwidth$cbfm2510_plusinterceptci <- cw_cis_cbfm$upper_plusintercept - cw_cis_cbfm$lower_plusintercept
     spp_slopes_coverage$cbfm2510_plusinterceptci <- as.numeric(cw_cis_cbfm$lower_plusintercept < spp_slopes_long$value & cw_cis_cbfm$upper_plusintercept > spp_slopes_long$value)     

     # Out of sample performance measures
     getpreds <- predict(CBFMmod, newdata = env_test@data, newB_space = test_basisfunctions, type = "response", se_fit = FALSE)
     model_error$cbfm2510 <- sqrt(colMeans((simdat_test$resp - getpreds)^2))
     if(myfam$family == "binomial") {
          tjur_r2$cbfm2510 <- sapply(1:num_spp, function(j) { 
                    m1 <- getpreds[,j][which(simdat_test$resp[,j] == 1)] %>%
                         mean(na.rm = TRUE)
                    m0 <- getpreds[,j][which(simdat_test$resp[,j] == 0)] %>%
                         mean(na.rm = TRUE)
                    m1 - m0     
                    })
          pred_logL$cbfm2510 <- sapply(1:num_spp, function(j) { 
                    sum(dbinom(simdat_test$resp[,j], size = 1, prob = getpreds[,j], log = TRUE))
                    })
          aucs$cbfm2510 <- sapply(1:num_spp, function(j) { 
                    pred <- ROCR::prediction(getpreds[,j], labels = simdat_test$resp[,j]) %>%
                         ROCR::performance(measure = "auc")
                    pred@y.values[[1]]
                    })
          }
#      if(myfam$family == "poisson") {
#           pred_logL$cbfm2510 <- sapply(1:num_spp, function(j) { 
#                     sum(dpois(simdat_test$resp[,j], lambda = getpreds[,j]), log = TRUE)
#                     })
#           pseudoR2$cbfm2510 <- sapply(1:num_spp, function(j) {
#                out <- cor(simdat_test$resp[,j], getpreds[,j], method = "spearman")
#                out^2*sign(out)
#                })
#           }
     rm(getpreds)     

               
     #-----------------------
     # Finish and save results
     #-----------------------          
     result_list[[i]] <- list(slopes_estimates = spp_slopes_long, slopes_CIwidth = spp_slopes_ciwidth, slopes_CIcoverage = spp_slopes_coverage, 
          tjur_r2 = tjur_r2, aucs = aucs, predlogLik = pred_logL, model_error = model_error, time = computing_time)
     }


# Save as RDS file (this can be changed depending on user's desire)
saveRDS(result_list, file = paste0("simstudy_PAresponses_",nlocs,".rds"))

