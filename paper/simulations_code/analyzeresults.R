##--------------------
## Template code for processing results from simulation study -- presence-absence responses. 
## It is assumed that the user has already run the script templatesimulationstudy_PAresponses.R several times using different numbers of spatial locations (N_total = 625, 1250, 2500), and consequently obtained results stored in simstudy_PAresponses_n625.rds, simstudy_PAresponses_n1250.rds, and simstudy_PAresponses_n2500.rds. 
## All results are assuming the two missing covariates are spatially stationary.
##--------------------
here::i_am("desired_directory/analyzeresults.R") # install.packages("here")

rm(list = ls())
library(tidyverse)
library(patchwork)

num_spp <- 30
nsims <- 200

simresults_N625 <- readRDS("simstudy_PAresponses_n625.rds")
simresults_N1250 <- readRDS("simstudy_PAresponses_n1250.rds")
simresults_N2500 <- readRDS("simstudy_PAresponses_n2500.rds")


#--------------------------
# Point performance on species-specific coefficients
#--------------------------
#' Results for HMSC are not calculated for presence-absence responses since they the link is misspecified by construction i.e., the spatial multivariate presence-absence data are generated assumed a logit link, and all methods use the logit link except for HMSC which can only use the probit link.

# Bias
method_names <- c("Stacked SDM", "sjSDM", "CBFM_25_5", "CBFM_25_10")

results_nobs <- NULL
for(k0 in 1:nsims) {
    cw_results500 <- simresults_N625[[k0]]$slopes_estimates %>%
        mutate(diff_stackedsdm = stackedsdm -  value,
               diff_sjsdm = sjsdm - value,
               diff_cbfm255 = cbfm255 - value,
               diff_cbfm2510 = cbfm2510 - value
               ) %>%
        dplyr::select(species, name, value, starts_with("diff"))
    
    cw_results1000 <- simresults_N1250[[k0]]$slopes_estimates %>%
        mutate(diff_stackedsdm = stackedsdm - value,
               diff_sjsdm = sjsdm - value,
               diff_cbfm255 = cbfm255 - value,
               diff_cbfm2510 = cbfm2510 - value
               ) %>%
        dplyr::select(species, name, value, starts_with("diff"))

    cw_results2000 <- simresults_N2500[[k0]]$slopes_estimates %>%
        mutate(diff_stackedsdm = stackedsdm - value,
               diff_sjsdm = sjsdm - value,
               diff_cbfm255 = cbfm255 - value,
               diff_cbfm2510 = cbfm2510 - value
               ) %>%
        dplyr::select(species, name, value, starts_with("diff"))

    results_nobs <- rbind(results_nobs, 
                          data.frame(dataset = k0, model = "Stationary covariates", N = "N = 500", cw_results500),
                          data.frame(dataset = k0, model = "Stationary covariates", N = "N = 1000", cw_results1000),
                          data.frame(dataset = k0, model = "Stationary covariates", N = "N = 2000", cw_results2000))
    
    }
    

biasresults <- results_nobs %>% 
    mutate(species = sapply(species, function(x) strsplit(x, "spp")[[1]][2] %>% as.numeric)) %>% 
    mutate(species = factor(species) %>% fct_inorder()) %>% 
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
    scale_x_discrete(labels = c("Stacked SDM", "sjSDM", "CBFM (d = 25; q = 5)", "CBFM (d = 25; q = 10)")) +
    labs(x = "Method", y = "Bias", fill = "Method", title = "Empirical bias") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")

bias_final


# RMSE
rmseresults <- results_nobs %>% 
    mutate(species = sapply(species, function(x) strsplit(x, "spp")[[1]][2] %>% as.numeric)) %>% 
    mutate(species = factor(species) %>% fct_inorder()) %>% 
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
    scale_x_discrete(labels = c("Stacked SDM", "sjSDM", "CBFM (d = 25; q = 5)", "CBFM (d = 25; q = 10)")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")

rmse_final

p <- bias_final | rmse_final
p


#--------------------------
# Coverage probability of species-specific coefficients 
#--------------------------
method_names <- c("Stacked_SDM", "sjSDM", "CBFM_25_5", "CBFM_25_10")

results_nobs <- results2_nobs <- NULL
for(k0 in 1:nsims) {
    cw_results500 <- cbind(simresults_N625[[k0]]$slopes_CIcoverage)
    cw_results1000 <- cbind(simresults_N1250[[k0]]$slopes_CIcoverage)
    cw_results2000 <- cbind(simresults_N2500[[k0]]$slopes_CIcoverage)
    results_nobs <- rbind(results_nobs, 
                          data.frame(dataset = k0, N = "N = 500", model = "Stationary covariates", cw_results500),
                          data.frame(dataset = k0, N = "N = 1000", model = "Stationary covariates", cw_results1000),
                          data.frame(dataset = k0, N = "N = 2000", model = "Stationary covariates", cw_results2000))
    
    
    cw_results500 <- cbind(simresults_N625[[k0]]$slopes_CIwidth)
    cw_results1000 <- cbind(simresults_N1250[[k0]]$slopes_CIwidth)
    cw_results2000 <- cbind(simresults_N2500[[k0]]$slopes_CIwidth)
    results2_nobs <- rbind(results2_nobs, 
                           data.frame(dataset = k0, N = "N = 500", model = "Stationary covariates", cw_results500),
                           data.frame(dataset = k0, N = "N = 1000", model = "Stationary covariates", cw_results1000),
                           data.frame(dataset = k0, N = "N = 2000", model = "Stationary covariates", cw_results2000))

    }


coverage_results <- results_nobs %>% 
    mutate(species = sapply(species, function(x) strsplit(x, "spp")[[1]][2] %>% as.numeric)) %>% 
    mutate(species = factor(species) %>% fct_inorder()) %>% 
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
levels(coverage_results$method) <- method_names


covprob_final <- ggplot(coverage_results, aes(x = method, y = coverage)) +
    geom_hline(yintercept = 0.95, linetype = 2) +
    geom_boxplot() +
    facet_grid(model ~  N, scales = "free_y") +
    labs(x = "Method", y = "Empirical coverage probability",  title = "Empirical coverage probability") +
    scale_x_discrete(labels = c("Stacked SDM", "sjSDM", "CBFM (d = 25; q = 5)", "CBFM (d = 25; q = 10)")) +
    scale_y_continuous(limits = c(0.5,1)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")

covprob_final



#--------------------------
# Interval width for species-specific coefficients
#--------------------------
ciwidth_results <- results2_nobs %>% 
    mutate(species = sapply(species, function(x) strsplit(x, "spp")[[1]][2] %>% as.numeric)) %>% 
    mutate(species = factor(species) %>% fct_inorder()) %>% 
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
    labs(x = "Method", y = "Empirical coverage probability", title = "Mean interval width") +
    scale_x_discrete(labels = c("Stacked SDM", "sjSDM", "CBFM (d = 25; q = 5)", "CBFM (d = 25; q = 10)")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")

ciwidth_final


p <- covprob_final | ciwidth_final
p



#-----------------------
# Tjur R-squared
#-----------------------
method_names <- c("Stacked SDM", "LVM (NNGP)", "LVM (GPP)", "sjSDM", "CBFM (d = 25; q = 5)", "CBFM (d = 25; q = 10)")

results_nobs <- NULL
for(k0 in 1:nsims) {
    cw_results500 <- simresults_N625[[k0]]$tjur_r2 %>% 
        rownames_to_column(var = "species")
    colnames(cw_results500)[-1] <- method_names
    cw_results1000 <- simresults_N1250[[k0]]$tjur_r2 %>% 
        rownames_to_column(var = "species")
    colnames(cw_results1000)[-1] <- method_names
    cw_results2000 <- simresults_N2500[[k0]]$tjur_r2 %>% 
        rownames_to_column(var = "species")
    colnames(cw_results2000)[-1] <- method_names
    results_nobs <- rbind(results_nobs, 
                          data.frame(dataset = k0, N = "N = 500", model = "Stationary covariates", cw_results500),
                          data.frame(dataset = k0, N = "N = 1000", model = "Stationary covariates", cw_results1000),
                          data.frame(dataset = k0, N = "N = 2000", model = "Stationary covariates", cw_results2000))                          
    }


tjurr2_results <- results_nobs %>% 
    pivot_longer(Stacked.SDM:last_col(), names_to = "method") %>% 
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
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")

tjurr2_final


#-----------------------
# AUC
#-----------------------
results_nobs <- NULL
for(k0 in 1:nsims) {
    cw_results500 <- simresults_N625[[k0]]$aucs %>% 
        rownames_to_column(var = "species")
    colnames(cw_results500)[-1] <- method_names
    cw_results1000 <- simresults_N1250[[k0]]$aucs %>% 
        rownames_to_column(var = "species")
    colnames(cw_results1000)[-1] <- method_names
    cw_results2000 <- simresults_N2500[[k0]]$aucs %>% 
        rownames_to_column(var = "species")
    colnames(cw_results2000)[-1] <- method_names
    results_nobs <- rbind(results_nobs, 
                          data.frame(dataset = k0, N = "N = 500", model = "Stationary covariates", cw_results500),
                          data.frame(dataset = k0, N = "N = 1000", model = "Stationary covariates", cw_results1000),
                          data.frame(dataset = k0, N = "N = 2000", model = "Stationary covariates", cw_results2000))                          
    }


auc_results <- results_nobs %>%  
    pivot_longer(Stacked.SDM:last_col(), names_to = "method") %>% 
    mutate(method = factor(method) %>% fct_inorder()) %>% 
    mutate(species = factor(species) %>% fct_inorder())   
levels(auc_results$method) <- method_names
auc_results <- auc_results %>% 
    group_by(N, model, species, method) %>% 
    summarise(value = mean(value))


auc_final <- ggplot(auc_results, aes(x = method, y = value)) +
    geom_boxplot() +
    facet_grid(model ~ N, scales = "free_y") +
    labs(y = "AUC", x = "Method", fill = "Method", title = "AUC") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")

auc_final


#-----------------------
# Predictive deviance
#-----------------------
results_nobs <- NULL
for(k0 in 1:nsims) {
    cw_results500 <- (simresults_N625[[k0]]$predlogLik/500) %>% 
        rownames_to_column(var = "species")
    colnames(cw_results500)[-1] <- method_names
    cw_results1000 <- (simresults_N1250[[k0]]$predlogLik/1000) %>% 
        rownames_to_column(var = "species")
    colnames(cw_results1000)[-1] <- method_names
    cw_results2000 <- (simresults_N2500[[k0]]$predlogLik/2000) %>% 
        rownames_to_column(var = "species")
    colnames(cw_results2000)[-1] <- method_names
    results_nobs <- rbind(results_nobs, 
                          data.frame(dataset = k0, N = "N = 500", model = "Stationary covariates", cw_results500),
                          data.frame(dataset = k0, N = "N = 1000", model = "Stationary covariates", cw_results1000),
                          data.frame(dataset = k0, N = "N = 2000", model = "Stationary covariates", cw_results2000))                          
    }


predlogLik_results <- results_nobs %>% 
    pivot_longer(Stacked.SDM:last_col(), names_to = "method") %>% 
    mutate(method = factor(method) %>% fct_inorder()) %>% 
    mutate(species = factor(species) %>% fct_inorder())   
levels(predlogLik_results$method) <- method_names
predlogLik_results <- predlogLik_results %>% 
    group_by(N, model, species, method) %>% 
    summarise(value = mean(value))


predlogLik_final <- ggplot(predlogLik_results, aes(x = method, y = -2*value)) +
    geom_boxplot() +
    facet_grid(model ~ N, scales = "free_y") +
    labs(y = "Predictive deviance", x = "Method", fill = "Method", title = "Predictive deviance per observational unit") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")

predlogLik_final


#-----------------------
# Computing time
#-----------------------
results_nobs <- NULL
for(k0 in 1:nsims) {
    cw_results500 <- simresults_N625[[k0]]$time
    colnames(cw_results500) <- method_names
    cw_results1000 <- simresults_N1250[[k0]]$time
    colnames(cw_results1000) <- method_names
    cw_results2000 <- simresults_N2500[[k0]]$time
    colnames(cw_results2000) <- method_names
    results_nobs <- rbind(results_nobs, 
                          data.frame(dataset = k0, N = "N = 500", model = "Stationary covariates", cw_results500),
                          data.frame(dataset = k0, N = "N = 1000", model = "Stationary covariates", cw_results1000),
                          data.frame(dataset = k0, N = "N = 2000", model = "Stationary covariates", cw_results2000))
    }


comptimes_results <- results_nobs %>%  
    pivot_longer(Stacked.SDM:last_col(), names_to = "method") %>% 
    mutate(method = factor(method) %>% fct_inorder()) 
levels(comptimes_results$method) <- method_names

computingtime_final <- ggplot(comptimes_results, aes(x = method, y = value/60)) +
    geom_boxplot() +
    facet_grid(model ~ N, scales = "free_y") +
    scale_y_log10() +
    labs(y = "Computing time", x = "Method", title = "Computing time (in minutes)") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")

computingtime_final


p <- (auc_final | tjurr2_final) / (predlogLik_final | computingtime_final)
p



