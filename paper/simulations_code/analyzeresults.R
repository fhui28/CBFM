##--------------------
## Template code for processing results from simulation study -- presence-absence responses. 
## It is assumed that the user has already run the script templatesimulationstudy_PAresponses.R several times using different numbers of spatial locations (N_total = 625, 1250, 2500), and consequently obtained results stored in simstudy_PA500training_stationary_datasetxxx.rds, simstudy_PA1000training_stationary_datasetxxx.rds, and simstudy_PA2000training_stationary_datasetxxx.rds
## All results are assuming the two missing covariates are spatially stationary. The script can be straightforwardly altered below to accomodate results from missing spatially non-stationary covariates, Poisson responses, and so on.
##--------------------
rm(list = ls())
library(autoFRK)
library(here)
library(doParallel)
library(foreach)
library(mvabund)
library(tidyverse)
library(patchwork)
library(abind)
library(ROCR)
library(sjSDM)
library(Hmsc)
library(CBFM) # devtools::install_github("fhui28/CBFM")     
# Note not all of the packages above are actually needed, but they are included to be consistent with templatesimulationstudy_PAresponses.R

num_spp <- 30
num_dataset <- 200

bind0 <- function(...) {
   abind::abind(..., along = 0)
   }
 

#--------------------------
# Collect all results
#--------------------------
# Note results associated with estimation and inference or coefficients are not calculated for HMSC and presence-absence responses, since the link is misspecified by construction i.e., the spatial multivariate presence-absence data are generated assumed a logit link, and all methods use the logit link except for HMSC which can only use the probit link.

results_nobs <- results_nobs_ci <- results_nobs_ciwidth <- results_nobs_tjuR2 <- 
   results_nobs_log_score <- results_nobs_brier_score <- results_nobs_time <- 
   results_nobs_residualcor <- NULL

for(k0 in 1:num_dataset) {
   raw_results <- readRDS(file = paste0("simstudy_PA500training_stationary_dataset", k0, ".rds"))
   cw_results500 <- raw_results$slopes_estimates %>%
      mutate(diff_stackedsdm = stackedsdm - value,
      diff_sjsdm = sjsdm - value,
      diff_cbfm255 = cbfm255 - value,
      diff_cbfm2510 = cbfm2510 - value,
      ) %>%
      dplyr::select(species, name, value, starts_with("diff"))
   cw_results500_ci <- cbind(raw_results$slopes_CIcoverage)
   cw_results500_ciwidth <- cbind(raw_results$slopes_CIwidth)
   cw_results500_tjurR2 <- raw_results$tjur_r2 %>% rownames_to_column(var = "species")
   cw_results500_log_score <- (raw_results$log_score/500) %>% rownames_to_column(var = "species")
   cw_results500_brier_score <- raw_results$brier_score %>% rownames_to_column(var = "species")
   cw_results500_time <- raw_results$time %>% rownames_to_column(var = "species")
   
   true_residual_correlation_missing <- cov2cor(raw_results$missing_linear_predictors$true)
   raw_results$missing_linear_predictors$nngp_omega <- computeAssociations(raw_results$m_nngp)[[1]]$mean
   raw_results$missing_linear_predictors$gpp_omega <- computeAssociations(raw_results$m_gpp)[[1]]$mean
   allEta <- lapply(poolMcmcChains(raw_results$m_nngp$postList), function(a) a[["Eta"]][[1]]) %>% do.call(bind0, .)
   allLambda <- lapply(poolMcmcChains(raw_results$m_nngp$postList), function(a) a[["Lambda"]][[1]]) %>% do.call(bind0, .)
   raw_results$meanresidualcor_nngp <- lapply(1:dim(allEta)[1], function(k) cov(allEta[k,,] %*% allLambda[k,,]) %>% cov2cor) %>% 
      abind::abind(., along = 3) %>% 
      apply(., c(1,2), mean)
   allEta <- lapply(poolMcmcChains(raw_results$m_gpp$postList), function(a) a[["Eta"]][[1]]) %>% do.call(bind0, .)
   allLambda <- lapply(poolMcmcChains(raw_results$m_gpp$postList), function(a) a[["Lambda"]][[1]]) %>% do.call(bind0, .)
   raw_results$meanresidualcor_gpp <- lapply(1:dim(allEta)[1], function(k) cov(allEta[k,,] %*% allLambda[k,,]) %>% cov2cor) %>% 
      abind::abind(., along = 3) %>% 
      apply(., c(1,2), mean)
   rm(allEta, allLambda)
   cw_results500_residualcor <- data.frame(
      hmsc_nngp_omega = norm(raw_results$missing_linear_predictors$nngp_omega - true_residual_correlation_missing, "F"),
      hmsc_nngp_corB = norm(raw_results$meanresidualcor_nngp - true_residual_correlation_missing, "F"),
      hmsc_gpp_omega = norm(raw_results$missing_linear_predictors$gpp_omega - true_residual_correlation_missing, "F"),
      hmsc_gpp_corB = norm(raw_results$meanresidualcor_gpp - true_residual_correlation_missing, "F"),
      sdjsm = norm(cov2cor(raw_results$missing_linear_predictors$sjsdm) - true_residual_correlation_missing, "F"),
      cbfm_25_5_G = norm(raw_results$missing_linear_predictors$cbfm255_G - true_residual_correlation_missing, "F"),
      cbfm_25_5_corB = norm(raw_results$missing_linear_predictors$cbfm255_corB$correlation - true_residual_correlation_missing, "F"),
      cbfm_25_10_G = norm(raw_results$missing_linear_predictors$cbfm2510_G - true_residual_correlation_missing, "F"),
      cbfm_25_10_corB = norm(raw_results$missing_linear_predictors$cbfm2510_corB$correlation - true_residual_correlation_missing, "F")
      )
    
    
   raw_results <- readRDS(file = paste0("simstudy_PA1000training_stationary_dataset", k0, ".rds"))
   cw_results1000 <- raw_results$slopes_estimates %>%
      mutate(diff_stackedsdm = stackedsdm - value,
      diff_sjsdm = sjsdm - value,
      diff_cbfm255 = cbfm255 - value,
      diff_cbfm2510 = cbfm2510 - value,
      ) %>%
      dplyr::select(species, name, value, starts_with("diff"))
   cw_results1000_ci <- cbind(raw_results$slopes_CIcoverage)
   cw_results1000_ciwidth <- cbind(raw_results$slopes_CIwidth)
   cw_results1000_tjurR2 <- raw_results$tjur_r2 %>% rownames_to_column(var = "species")
   cw_results1000_log_score <- (raw_results$log_score/1000) %>% rownames_to_column(var = "species")
   cw_results1000_brier_score <- raw_results$brier_score %>% rownames_to_column(var = "species")
   cw_results1000_time <- raw_results$time %>% rownames_to_column(var = "species")
   
   true_residual_correlation_missing <- cov2cor(raw_results$missing_linear_predictors$true)
   raw_results$missing_linear_predictors$nngp_omega <- computeAssociations(raw_results$m_nngp)[[1]]$mean
   raw_results$missing_linear_predictors$gpp_omega <- computeAssociations(raw_results$m_gpp)[[1]]$mean
   allEta <- lapply(poolMcmcChains(raw_results$m_nngp$postList), function(a) a[["Eta"]][[1]]) %>% do.call(bind0, .)
   allLambda <- lapply(poolMcmcChains(raw_results$m_nngp$postList), function(a) a[["Lambda"]][[1]]) %>% do.call(bind0, .)
   raw_results$meanresidualcor_nngp <- lapply(1:dim(allEta)[1], function(k) cov(allEta[k,,] %*% allLambda[k,,]) %>% cov2cor) %>% 
      abind::abind(., along = 3) %>% 
      apply(., c(1,2), mean)
   allEta <- lapply(poolMcmcChains(raw_results$m_gpp$postList), function(a) a[["Eta"]][[1]]) %>% do.call(bind0, .)
   allLambda <- lapply(poolMcmcChains(raw_results$m_gpp$postList), function(a) a[["Lambda"]][[1]]) %>% do.call(bind0, .)
   raw_results$meanresidualcor_gpp <- lapply(1:dim(allEta)[1], function(k) cov(allEta[k,,] %*% allLambda[k,,]) %>% cov2cor) %>% 
      abind::abind(., along = 3) %>% 
      apply(., c(1,2), mean)
   rm(allEta, allLambda)
   cw_results1000_residualcor <- data.frame(
      hmsc_nngp_omega = norm(raw_results$missing_linear_predictors$nngp_omega - true_residual_correlation_missing, "F"),
      hmsc_nngp_corB = norm(raw_results$meanresidualcor_nngp - true_residual_correlation_missing, "F"),
      hmsc_gpp_omega = norm(raw_results$missing_linear_predictors$gpp_omega - true_residual_correlation_missing, "F"),
      hmsc_gpp_corB = norm(raw_results$meanresidualcor_gpp - true_residual_correlation_missing, "F"),
      sdjsm = norm(cov2cor(raw_results$missing_linear_predictors$sjsdm) - true_residual_correlation_missing, "F"),
      cbfm_25_5_G = norm(raw_results$missing_linear_predictors$cbfm255_G - true_residual_correlation_missing, "F"),
      cbfm_25_5_corB = norm(raw_results$missing_linear_predictors$cbfm255_corB$correlation - true_residual_correlation_missing, "F"),
      cbfm_25_10_G = norm(raw_results$missing_linear_predictors$cbfm2510_G - true_residual_correlation_missing, "F"),
      cbfm_25_10_corB = norm(raw_results$missing_linear_predictors$cbfm2510_corB$correlation - true_residual_correlation_missing, "F")
      )

      
      
   raw_results <- readRDS(file = paste0("simstudy_PA2000training_stationary_dataset", k0, ".rds"))
   cw_results2000 <- raw_results$slopes_estimates %>%
      mutate(diff_stackedsdm = stackedsdm - value,
      diff_sjsdm = sjsdm - value,
      diff_cbfm255 = cbfm255 - value,
      diff_cbfm2510 = cbfm2510 - value,
      ) %>%
      dplyr::select(species, name, value, starts_with("diff"))
   cw_results2000_ci <- cbind(raw_results$slopes_CIcoverage)
   cw_results2000_ciwidth <- cbind(raw_results$slopes_CIwidth)
   cw_results2000_tjurR2 <- raw_results$tjur_r2 %>% rownames_to_column(var = "species")
   cw_results2000_log_score <- (raw_results$log_score/2000) %>% rownames_to_column(var = "species")
   cw_results2000_brier_score <- raw_results$brier_score %>% rownames_to_column(var = "species")
   cw_results2000_time <- raw_results$time %>% rownames_to_column(var = "species")
   
   true_residual_correlation_missing <- cov2cor(raw_results$missing_linear_predictors$true)
   raw_results$missing_linear_predictors$nngp_omega <- computeAssociations(raw_results$m_nngp)[[1]]$mean
   raw_results$missing_linear_predictors$gpp_omega <- computeAssociations(raw_results$m_gpp)[[1]]$mean
   allEta <- lapply(poolMcmcChains(raw_results$m_nngp$postList), function(a) a[["Eta"]][[1]]) %>% do.call(bind0, .)
   allLambda <- lapply(poolMcmcChains(raw_results$m_nngp$postList), function(a) a[["Lambda"]][[1]]) %>% do.call(bind0, .)
   raw_results$meanresidualcor_nngp <- lapply(1:dim(allEta)[1], function(k) cov(allEta[k,,] %*% allLambda[k,,]) %>% cov2cor) %>% 
      abind::abind(., along = 3) %>% 
      apply(., c(1,2), mean)
   allEta <- lapply(poolMcmcChains(raw_results$m_gpp$postList), function(a) a[["Eta"]][[1]]) %>% do.call(bind0, .)
   allLambda <- lapply(poolMcmcChains(raw_results$m_gpp$postList), function(a) a[["Lambda"]][[1]]) %>% do.call(bind0, .)
   raw_results$meanresidualcor_gpp <- lapply(1:dim(allEta)[1], function(k) cov(allEta[k,,] %*% allLambda[k,,]) %>% cov2cor) %>% 
      abind::abind(., along = 3) %>% 
      apply(., c(1,2), mean)
   rm(allEta, allLambda)
   cw_results2000_residualcor <- data.frame(
      hmsc_nngp_omega = norm(raw_results$missing_linear_predictors$nngp_omega - true_residual_correlation_missing, "F"),
      hmsc_nngp_corB = norm(raw_results$meanresidualcor_nngp - true_residual_correlation_missing, "F"),
      hmsc_gpp_omega = norm(raw_results$missing_linear_predictors$gpp_omega - true_residual_correlation_missing, "F"),
      hmsc_gpp_corB = norm(raw_results$meanresidualcor_gpp - true_residual_correlation_missing, "F"),
      sdjsm = norm(cov2cor(raw_results$missing_linear_predictors$sjsdm) - true_residual_correlation_missing, "F"),
      cbfm_25_5_G = norm(raw_results$missing_linear_predictors$cbfm255_G - true_residual_correlation_missing, "F"),
      cbfm_25_5_corB = norm(raw_results$missing_linear_predictors$cbfm255_corB$correlation - true_residual_correlation_missing, "F"),
      cbfm_25_10_G = norm(raw_results$missing_linear_predictors$cbfm2510_G - true_residual_correlation_missing, "F"),
      cbfm_25_10_corB = norm(raw_results$missing_linear_predictors$cbfm2510_corB$correlation - true_residual_correlation_missing, "F")
      )
  
  
   results_nobs <- rbind(results_nobs, 
                         data.frame(dataset = k0, model = "Stationary covariates", N = "N = 500", cw_results500),
                         data.frame(dataset = k0, model = "Stationary covariates", N = "N = 1000", cw_results1000),
                         data.frame(dataset = k0, model = "Stationary covariates", N = "N = 2000", cw_results2000)
                         )
   results_nobs_ci <- rbind(results_nobs_ci, 
                            data.frame(dataset = k0, N = "N = 500", model = "Stationary covariates", cw_results500_ci),
                            data.frame(dataset = k0, N = "N = 1000", model = "Stationary covariates", cw_results1000_ci),
                            data.frame(dataset = k0, N = "N = 2000", model = "Stationary covariates", cw_results2000_ci)
                            )
   results_nobs_ciwidth <- rbind(results_nobs_ciwidth, 
                            data.frame(dataset = k0, N = "N = 500", model = "Stationary covariates", cw_results500_ciwidth),
                            data.frame(dataset = k0, N = "N = 1000", model = "Stationary covariates", cw_results1000_ciwidth),
                            data.frame(dataset = k0, N = "N = 2000", model = "Stationary covariates", cw_results2000_ciwidth)
                            )
   results_nobs_tjuR2 <- rbind(results_nobs_tjuR2, 
                         data.frame(dataset = k0, N = "N = 500", model = "Stationary covariates", cw_results500_tjurR2),
                         data.frame(dataset = k0, N = "N = 1000", model = "Stationary covariates", cw_results1000_tjurR2),
                         data.frame(dataset = k0, N = "N = 2000", model = "Stationary covariates", cw_results2000_tjurR2)
                         )
   results_nobs_log_score <- rbind(results_nobs_log_score, 
                               data.frame(dataset = k0, N = "N = 500", model = "Stationary covariates", cw_results500_log_score),
                               data.frame(dataset = k0, N = "N = 1000", model = "Stationary covariates", cw_results1000_log_score),
                               data.frame(dataset = k0, N = "N = 2000", model = "Stationary covariates", cw_results2000_log_score)
                               )
   results_nobs_brier_score <- rbind(results_nobs_brier_score, 
                                   data.frame(dataset = k0, N = "N = 500", model = "Stationary covariates", cw_results500_brier_score),
                                   data.frame(dataset = k0, N = "N = 1000", model = "Stationary covariates", cw_results1000_brier_score),
                                   data.frame(dataset = k0, N = "N = 2000", model = "Stationary covariates", cw_results2000_brier_score)
                                   )
   results_nobs_time <- rbind(results_nobs_time, 
                              data.frame(dataset = k0, N = "N = 500", model = "Stationary covariates", cw_results500_time),
                              data.frame(dataset = k0, N = "N = 1000", model = "Stationary covariates", cw_results1000_time),
                              data.frame(dataset = k0, N = "N = 2000", model = "Stationary covariates", cw_results2000_time)
                              )
   results_nobs_residualcor <- rbind(results_nobs_residualcor, 
                              data.frame(dataset = k0, N = "N = 500", model = "Stationary covariates", cw_results500_residualcor),
                              data.frame(dataset = k0, N = "N = 1000", model = "Stationary covariates", cw_results1000_residualcor),
                              data.frame(dataset = k0, N = "N = 2000", model = "Stationary covariates", cw_results2000_residualcor)
                              )
  
  }
    

rm(list = ls(pattern = "cw_results"))
rm(list = ls(pattern = "meanresidualcor"))
rm(raw_results)
gc()


#--------------------------
# Bias of species-specific coefficients 
#--------------------------
method_names <- c("Stacked SDM", "sjSDM", "CBFM_25_5", "CBFM_25_10", "CBFM_50_5", "CBFM_50_10")

biasresults <- results_nobs %>% 
    mutate(species = sapply(species, function(x) strsplit(x, "spp")[[1]][2] %>% as.numeric) %>% factor %>% fct_inorder()) %>% 
    rename(covariate = name) %>% 
    mutate(covariate = factor(covariate) %>% fct_inorder()) %>% 
    group_by(N, model, species, covariate)
biasresults <- biasresults %>% 
    summarise(stackedsdm = mean(diff_stackedsdm),
              sjsdm = mean(diff_sjsdm),
              cbfm255 = mean(diff_cbfm255),
              cbfm2510 = mean(diff_cbfm2510)
              ) %>% 
    arrange(N, model, species) %>% 
    pivot_longer(stackedsdm:last_col(), values_to = "bias", names_to = "method") %>% 
    mutate(method = factor(method) %>% fct_inorder())
levels(biasresults$method) <- method_names


bias_final <- ggplot(biasresults, aes(x = method, y = bias)) +
    geom_boxplot() +
    geom_hline(yintercept = 0, linetype = 2, color = "darkgrey") +
    facet_grid(model ~ N, scales = "free_y") +
    scale_x_discrete(labels = c("Stacked SDM", "sjSDM", "CBFM (25/5)", "CBFM (25/10)")) +
    labs(x = "Method", y = "Bias", fill = "Method", title = "Empirical bias") +
    #scale_color_brewer(palette = "Set2") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank())

bias_final


#--------------------------
# RMSE of species-specific coefficients 
#--------------------------
rmseresults <- results_nobs %>% 
    mutate(species = sapply(species, function(x) strsplit(x, "spp")[[1]][2] %>% as.numeric) %>% factor() %>% fct_inorder()) %>% 
    rename(covariate = name) %>% 
    mutate(covariate = factor(covariate) %>% fct_inorder()) %>% 
    group_by(N, model, species, covariate)
rmseresults <- rmseresults %>% 
   summarise(stackedsdm = sqrt(mean(diff_stackedsdm^2)),
             sjsdm = sqrt(mean(diff_sjsdm^2)),
             cbfm255 = sqrt(mean(diff_cbfm255^2)),
             cbfm2510 = sqrt(mean(diff_cbfm2510^2))
             ) %>% 
    arrange(N, model, species) %>% 
    pivot_longer(stackedsdm:last_col(), values_to = "RMSE", names_to = "method") %>% 
    mutate(method = factor(method) %>% fct_inorder())
levels(rmseresults$method) <- method_names


rmse_final <- ggplot(rmseresults, aes(x = method, y = RMSE)) +
    geom_boxplot() +
    facet_grid(model ~ N, scales = "free_y") +
    labs(x = "Method", y = "RMSE", fill = "Method", title = "Empirical RMSE") +
    scale_x_discrete(labels = c("Stacked SDM", "sjSDM", "CBFM (25/5)", "CBFM (25/10)")) +
    #scale_color_brewer(palette = "Set2") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank())

rmse_final


p <- bias_final | rmse_final
p


#--------------------------
# Coverage probability of species-specific coefficients 
#--------------------------
method_names_ci <- c("Stacked SDM", "sjSDM", "CBFM_25_5", "CBFM_25_10")

coverage_results <- results_nobs_ci %>% 
    mutate(species = sapply(species, function(x) strsplit(x, "spp")[[1]][2] %>% as.numeric) %>% factor() %>% fct_inorder()) %>% 
    rename(covariate = name) %>% 
    mutate(covariate = factor(covariate) %>% fct_inorder()) %>% 
    group_by(N, model, species, covariate)
coverage_results <- coverage_results %>% 
    summarise(stackedsdm = mean(stackedsdm),
              sjsdm = mean(sjsdm),
              cbfm255 = mean(cbfm255),
              cbfm2510 = mean(cbfm2510)
              ) %>% 
    arrange(N, model, species) %>% 
    pivot_longer(stackedsdm:last_col(), values_to = "coverage", names_to = "method") %>% 
    mutate(method = factor(method) %>% fct_inorder())
levels(coverage_results$method) <- method_names_ci


covprob_final <- ggplot(coverage_results, aes(x = method, y = coverage)) +
    geom_hline(yintercept = 0.95, linetype = 2) +
    geom_boxplot() +
    facet_grid(model ~  N, scales = "free_y") +
    labs(x = "Method", y = "Empirical coverage probability",  title = "Empirical coverage probability") +
    #scale_color_brewer(palette = "Set2") +
    scale_x_discrete(labels = c("Stacked SDM", "sjSDM", "CBFM (25/5)", "CBFM (25/10)")) +
    scale_y_continuous(limits = c(0.5,1)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank())

covprob_final



#--------------------------
# Interval width for species-specific coefficients
#--------------------------
ciwidth_results <- results_nobs_ciwidth %>% 
    mutate(species = sapply(species, function(x) strsplit(x, "spp")[[1]][2] %>% as.numeric) %>% factor() %>% fct_inorder()) %>% 
    rename(covariate = name) %>% 
    mutate(covariate = factor(covariate) %>% fct_inorder()) %>% 
    group_by(N, model, species, covariate)
ciwidth_results <- ciwidth_results %>% 
    summarise(stackedsdm = mean(stackedsdm),
              sjsdm = mean(sjsdm),
              cbfm255 = mean(cbfm255),
              cbfm2510 = mean(cbfm2510)
              ) %>% 
    arrange(N, model, species) %>% 
    pivot_longer(stackedsdm:last_col(), values_to = "coverage", names_to = "method") %>% 
    mutate(method = factor(method) %>% fct_inorder())
levels(ciwidth_results$method) <- method_names


ciwidth_final <- ggplot(ciwidth_results, aes(x = method, y = coverage)) +
    geom_boxplot() +
    facet_grid(model ~  N, scales = "free_y") +
    labs(x = "Method", y = "Mean interval width", title = "Mean interval width") +
    scale_x_discrete(labels = c("Stacked SDM", "sjSDM", "CBFM (25/5)", "CBFM (25/10)")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ciwidth_final


p <- covprob_final | ciwidth_final
p


#-----------------------
# Tjur R-squared
#-----------------------
method_names <- c("Stacked SDM", "LVM (NNGP)", "LVM (GPP)", "sjSDM", "CBFM (25/5)", "CBFM (25/10)")
colnames(results_nobs_tjuR2)[-c(1:4)] <- method_names

tjurr2_results <- results_nobs_tjuR2 %>%  
    pivot_longer("Stacked SDM":last_col(), names_to = "method") %>% 
    mutate(method = factor(method) %>% fct_inorder()) %>% 
    mutate(species = factor(species) %>% fct_inorder()) 
levels(tjurr2_results$method) <- method_names
tjurr2_results <- tjurr2_results %>% 
    group_by(N, model, species, method) %>% 
    summarise(value = mean(value))


tjurr2_final <- ggplot(tjurr2_results, aes(x = method, y = value)) +
    geom_boxplot() +
    facet_grid(model ~ N, scales = "free_y") +
    labs(y = "Tjur R-squared", x = "Method", title = "Tjur R-squared") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank())

tjurr2_final


#-----------------------
# Log score
#-----------------------
colnames(results_nobs_log_score)[-c(1:4)] <- method_names

logscore_results <- results_nobs_log_score %>% 
    pivot_longer("Stacked SDM":last_col(), names_to = "method") %>% 
    mutate(method = factor(method) %>% fct_inorder()) %>% 
    mutate(species = factor(species) %>% fct_inorder())   
levels(logscore_results$method) <- method_names
logscore_results <- logscore_results %>% 
    group_by(N, model, species, method) %>% 
    summarise(value = mean(value))

    
logscore_final <- ggplot(logscore_results, aes(x = method, y = value)) +
    geom_boxplot() +
    facet_grid(model ~ N, scales = "free_y") +
    labs(y = "Log score", x = "Method", fill = "Method", title = "Log score") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank())

logscore_final


#-----------------------
# Brier score
#-----------------------
colnames(results_nobs_brier_score)[-c(1:4)] <- method_names

brierscore_results <- results_nobs_brier_score %>% 
    pivot_longer("Stacked SDM":last_col(), names_to = "method") %>% 
    mutate(method = factor(method) %>% fct_inorder()) %>% 
    mutate(species = factor(species) %>% fct_inorder())   
levels(brierscore_results$method) <- method_names
brierscore_results <- brierscore_results %>% 
    group_by(N, model, species, method) %>% 
    summarise(value = mean(value))


brierscore_final <- ggplot(brierscore_results, aes(x = method, y = value)) +
    geom_boxplot() +
    facet_grid(model ~ N, scales = "free_y") +
    labs(y = "Brier score", x = "Method", fill = "Method", title = "Brier score") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank())

brierscore_final


#-----------------------
# Computing time
#-----------------------
colnames(results_nobs_time)[-c(1:4)] <- method_names

comptimes_results <- results_nobs_time %>%  
    pivot_longer("Stacked SDM":last_col(), names_to = "method") %>% 
    mutate(method = factor(method) %>% fct_inorder()) 
levels(comptimes_results$method) <- method_names

computingtime_final <- ggplot(comptimes_results, aes(x = method, y = value/60)) +
    geom_boxplot() +
    facet_grid(model ~ N, scales = "free_y") +
    scale_y_log10() +
    labs(y = "Computing time (in minutes)", x = "Method", title = "Computing time") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank())

computingtime_final


p <- (brierscore_final | tjurr2_final) / (logscore_final | computingtime_final)
p


#-----------------------
# Residual spatial correlation between species (using multiple definitions)
#-----------------------
method_names <- c("LVM (NNGP/Omega)", "LVM (NNGP/Linpred)", "LVM (GPP/Omega)", "LVM (GPP/Linpred)", "sjSDM", 
                  "CBFM (25/5/G)", "CBFM (25/5/Linpred)", "CBFM (25/10/G)", "CBFM (25/10/Linpred)")
colnames(results_nobs_residualcor)[-c(1:3)] <- method_names

residualcor_results <- results_nobs_residualcor %>%  
   pivot_longer("LVM (NNGP/Omega)":last_col(), names_to = "method") %>% 
   mutate(method = factor(method) %>% fct_inorder()) 
levels(residualcor_results$method) <- method_names

residualcor_results_final <- ggplot(residualcor_results, aes(x = method, y = value)) +
   geom_boxplot() +
   facet_grid(model ~ N, scales = "free_y") +
   scale_y_log10() +
   labs(y = "Frobenius norm", x = "Method", title = "Residual correlation between species") +
   theme_bw() +
   theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank())

residualcor_results_final


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
