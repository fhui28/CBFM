## See github issues for outstanding things to do!!!

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

function() {
     y = simy_train
     data = dat_train
     formula = useformula 
     ziformula <- NULL
     weights = NULL
     B_space = train_basisfunctions
     B_time = NULL
     B_spacetime = NULL
     ncores = detectCores() - 4
     family = ztnb2() 
     start_params = NULL
     control = list(trace = 1)
     G_control = list(NULL)
     Sigma_control = list(NULL)
     knots = NULL
     ziknots = NULL
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
     TMB_directories = list(cpp = system.file("executables", package = "CBFM"), compile = system.file("executables", package = "CBFM"))
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
