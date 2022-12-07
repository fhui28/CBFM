#' ---
#' title: Testing out some additional prediction functions purely for Chris Haak's application
#' author: FKCH
#' date: Code started Dec 2022
#' ---
#' 
#' 
##------------------------------
# Setup
##------------------------------
rm(list = ls())
library(tidyverse)
library(autoFRK)
library(FRK)
library(MASS)
library(mvabund)
library(mvtnorm)
library(ROCR)
library(sp)
library(RandomFields)
library(gamlss)
source("additionalpredictfns.R")
library(CBFM)

sessionInfo()
# R version 3.6.3 (2020-02-29)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Linux Mint 20.3
# 
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
# LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0
# 
# locale:
#      [1] LC_CTYPE=en_AU.UTF-8       LC_NUMERIC=C               LC_TIME=en_AU.UTF-8        LC_COLLATE=en_AU.UTF-8     LC_MONETARY=en_AU.UTF-8   
# [6] LC_MESSAGES=en_AU.UTF-8    LC_PAPER=en_AU.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_AU.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#      [1] parallel  splines   stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#      [1] CBFM_0.1                TMB_1.9.1               gamlss.tr_5.1-7         gamlss_5.4-3            nlme_3.1-144            gamlss.dist_6.0-3      
# [7] gamlss.data_6.0-2       RandomFields_3.3.14     RandomFieldsUtils_1.2.5 sp_1.5-0                ROCR_1.0-11             mvtnorm_1.1-3          
# [13] mvabund_4.2.1           MASS_7.3-51.5           FRK_2.0.5               autoFRK_1.4.3           spam_2.8-0              forcats_0.5.1          
# [19] stringr_1.4.0           dplyr_1.0.9             purrr_0.3.4             readr_2.1.2             tidyr_1.2.0             tibble_3.1.7           
# [25] ggplot2_3.3.6           tidyverse_1.3.1        
# 
# loaded via a namespace (and not attached):
#      [1] readxl_1.4.0         backports_1.4.1      Hmisc_4.7-0          filehashSQLite_0.2-4 tdigest_0.3.0        plyr_1.8.7           usethis_2.1.6       
# [8] digest_0.6.29        gratia_0.7.3         foreach_1.5.2        htmltools_0.5.2      LatticeKrig_8.4      viridis_0.6.2        fansi_1.0.3         
# [15] magrittr_2.0.3       checkmate_2.1.0      memoise_2.0.1        doParallel_1.0.17    cluster_2.1.0        remotes_2.4.2        tzdb_0.3.0          
# [22] modelr_0.1.8         xts_0.12.1           prettyunits_1.1.1    jpeg_0.1-9           colorspace_2.0-3     blob_1.2.3           rvest_1.0.2         
# [29] haven_2.5.0          xfun_0.31            callr_3.7.0          crayon_1.5.1         jsonlite_1.8.3       iterators_1.0.14     survival_3.1-8      
# [36] zoo_1.8-10           glue_1.6.2           gtable_0.3.0         filematrix_1.3       car_3.1-0            pkgbuild_1.3.1       maps_3.4.0          
# [43] abind_1.4-5          scales_1.2.0         DBI_1.1.2            rstatix_0.7.0        miniUI_0.1.1.1       Rcpp_1.0.8.3         viridisLite_0.4.0   
# [50] xtable_1.8-4         htmlTable_2.4.0      foreign_0.8-75       bit_4.0.4            dotCall64_1.0-1      Formula_1.2-4        tweedie_2.3.3       
# [57] intervals_0.15.2     profvis_0.3.7        htmlwidgets_1.5.4    httr_1.4.4           FNN_1.1.3.1          RColorBrewer_1.1-3   ellipsis_0.3.2      
# [64] urlchecker_1.0.1     pkgconfig_2.0.3      nnet_7.3-13          dbplyr_2.2.1         utf8_1.2.2           tidyselect_1.1.2     rlang_1.0.4         
# [71] reshape2_1.4.4       later_1.3.0          munsell_0.5.0        cellranger_1.1.0     tools_3.6.3          cachem_1.0.6         cli_3.3.0           
# [78] generics_0.1.2       RSQLite_2.2.14       devtools_2.4.4       broom_1.0.0          fastmap_1.1.0        processx_3.5.3       knitr_1.39          
# [85] bit64_4.0.5          fs_1.5.2             filehash_2.4-3       mime_0.12            mvnfast_0.2.7        xml2_1.3.3           compiler_3.6.3      
# [92] rstudioapi_0.13      png_0.1-7            ggsignif_0.6.3       reprex_2.0.1         spacetime_1.2-7      statmod_1.4.36       stringi_1.7.6       
# [99] ps_1.7.0             desc_1.4.1           fields_14.0          lattice_0.20-40      Matrix_1.4-1         vctrs_0.4.1          pillar_1.7.0        
# [106] lifecycle_1.0.1      data.table_1.14.2    patchwork_1.1.1      httpuv_1.6.5         R6_2.5.1             latticeExtra_0.6-29  promises_1.2.0.1    
# [113] gridExtra_2.3        codetools_0.2-16     sessioninfo_1.2.2    pkgload_1.3.0        assertthat_0.2.1     rprojroot_2.0.3      withr_2.5.0         
# [120] mgcv_1.8-31          hms_1.1.1            grid_3.6.3           rpart_4.1-15         carData_3.0-5        sparseinv_0.1.3      ggpubr_0.4.0        
# [127] numDeriv_2016.8-1.1  shiny_1.7.1          lubridate_1.8.0      base64enc_0.1-3     

     

##------------------------------
## Pulled from help file for makeahurdle function 
## **Example 2: Fitting a CBFM to spatial multivariate hurdle negative binomial data**
## simulated from a spatial hurdle latent variable model
## Please note the data generation process (thus) differs from CBFM.
##------------------------------
set.seed(2021)
num_sites <- 1000 # 500 (units) sites for training set + 500 sites for testing.
num_spp <- 50 # Number of species
num_X <- 4 # Number of regression slopes

# Generate and combine latent variables model for presence-absence and zero truncated Poisson X
spp_slopes_pa <- matrix(runif(num_spp * num_X, -1, 1), nrow = num_spp)
spp_slopes_ztnb <- matrix(runif(num_spp * num_X, -1, 1), nrow = num_spp)
spp_intercepts_pa <- runif(num_spp, -2, 0)
spp_intercepts_ztnb <- runif(num_spp, -2, 0)
spp_dispersion_ztnb <- runif(num_spp, 0, 5)

# Simulate spatial coordinates and environmental covariate components
xy <- data.frame(x = runif(num_sites, 0, 5), y = runif(num_sites, 0, 5))
X <- rmvnorm(num_sites, mean = rep(0,4))
colnames(X) <- c("temp", "depth", "chla", "O2")
dat <- data.frame(xy, X)
mm <- model.matrix(~ temp + depth + chla + O2 - 1, data = dat) %>%
scale %>%
as.matrix

# Simulate latent variable components
true_lvs_pa <- RFsimulate(model = RMexp(var=1, scale=2),
x = xy$x, y = xy$y, n = 2)@data %>%
as.matrix
spp_loadings_pa <- matrix(runif(num_spp * 2, -1, 1), nrow = num_spp)
true_lvs_ztnb <- RFsimulate(model = RMexp(var=1, scale=2.5),
x = xy$x, y = xy$y, n = 2)@data %>%
as.matrix
spp_loadings_ztnb <- matrix(runif(num_spp * 2, -1, 1), nrow = num_spp)
set.seed(NULL)

# Simulate spatial multivariate abundance data
eta <- tcrossprod(cbind(1,mm), cbind(spp_intercepts_pa,spp_slopes_pa)) +
tcrossprod(true_lvs_pa, spp_loadings_pa)
simy_pa <- matrix(rbinom(num_sites * num_spp, size = 1,
prob = plogis(eta)), nrow = num_sites)

# Now simulate spatial count data from a truncated NB distribution
eta <- tcrossprod(cbind(1,mm), cbind(spp_intercepts_ztnb,spp_slopes_ztnb)) +
tcrossprod(true_lvs_ztnb, spp_loadings_ztnb)
ztNBR <- trun.r(par = 0, family = "NBI", type = "left")
simy_ztnb <- matrix(ztNBR(num_sites * num_spp, mu = exp(eta)), nrow = num_sites)

# Spatial multivariate count data from a hurdle NB model is then the product of the two
simy <- simy_pa *  simy_ztnb

# Form training and test sets
dat_train <- dat[1:500,]
dat_test <- dat[501:1000,]
simy_train <- simy[1:500,]
simy_test <- simy[501:1000,]

# Delete the "component" responses and present you only observe the final response
rm(eta, simy_pa, simy_ztnb, X, mm, spp_loadings_pa, spp_loadings_ztnb, true_lvs_pa, true_lvs_ztnb,
xy, simy, dat, ztNBR)



# Set up spatial basis functions for CBFM -- Most users will start here!
# We will use this basis functions in both CBFM fits
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

# Fit presence-absence CBFM
simy_pa_train <- (simy_train > 0)*1
useformula <- ~ temp + depth + chla + O2
fitcbfm_pa <- CBFM(y = simy_pa_train, formula = useformula, data = dat_train,
B_space = train_basisfunctions, family = binomial(), control = list(trace = 1))
rm(simy_pa_train)

# Now fit zero-truncated negative binomial CBFM
# This could take a while...
fitcbfm_ztnb <- CBFM(y = simy_train, formula = useformula, data = dat_train,
                     B_space = train_basisfunctions, family = ztnb2(), control = list(trace = 1))


# Finally, form the hurdle CBFM, and remove the component CBFMs to save space!
fitcbfm_hurdle <- makeahurdle(pa_object = fitcbfm_pa, count_object = fitcbfm_ztnb)
rm(fitcbfm_pa, fitcbfm_ztnb)

fitcbfm_hurdle



##------------------------------
# Comparing different forms of predictions
##------------------------------
oospred_cbfm_timespacetime_hurdle_mean <- predict(fitcbfm_hurdle, newdata_pa = dat_test, newdata_count = dat_test,
                                   new_B_space_pa = test_basisfunctions, new_B_space_count = test_basisfunctions,
                                   type = "response")


ptpred_ztnb_test1 <- predict(fitcbfm_hurdle$count_fit, newdata = dat_test,
                                    new_B_space = test_basisfunctions,
                                    type = "link")
ptpred_ztnb_test2 <- predict_ztcount_CBFM(fitcbfm_hurdle$count_fit, newdata = dat_test,
                                    new_B_space = test_basisfunctions,
                                    type = "y", parameter_uncertainty = FALSE)



ptpred_pa <- predict_PA_CBFM(fitcbfm_hurdle$pa_fit, newdata = dat_test,
                     new_B_space = test_basisfunctions, 
                     type = "response", parameter_uncertainty = TRUE)
ptpred_ztnb <- predict_ztcount_CBFM(fitcbfm_hurdle$count_fit, newdata = dat_test,
                       new_B_space = test_basisfunctions, 
                       type = "link", parameter_uncertainty = TRUE)



oospred_cbfm_timespacetime_hurdle_median <- array(NA, dim = dim(ptpred_pa))
for(k0 in 1:dim(ptpred_pa)[3]) {
     oospred_cbfm_timespacetime_hurdle_median[,,k0] <- gamlss.dist::qZANBI(0.5, 
                                                                     mu = exp(ptpred_ztnb[,,k0])+1e-12,  
                                                                     sigma = matrix(fitcbfm_hurdle$count_fit$dispparam, nrow = nrow(ptpred_pa), ncol = ncol(ptpred_pa), byrow = TRUE),
                                                                     nu = 1-ptpred_pa[,,k0]) %>% 
          matrix(., nrow = nrow(ptpred_pa), ncol = ncol(ptpred_pa))
     }




quantile_threshold <- 0.9
oospred_cbfm_timespacetime_hurdle_trim <- array(NA, dim = dim(ptpred_pa))
for(k0 in 1:dim(ptpred_pa)[3]) {
     oospred_cbfm_timespacetime_hurdle_trim[,,k0] <- pmin(
          oospred_cbfm_timespacetime_hurdle_mean,
          gamlss.dist::qZANBI(quantile_threshold, 
                         mu = exp(ptpred_ztnb[,,k0])+1e-12,  
                         sigma = matrix(fitcbfm_hurdle$count_fit$dispparam, nrow = nrow(ptpred_pa), ncol = ncol(ptpred_pa), byrow = TRUE),
                         nu = 1-ptpred_pa[,,k0]) %>% 
          matrix(., nrow = nrow(ptpred_pa), ncol = ncol(ptpred_pa))
          )
     }




linexpred_hurdlenb <- function(a = 0.01, p, phi, mu) { 
     out <- (1-p) + p * (1/(1+phi*mu))^(1/phi) * (1 - (1/(1+phi*mu))^(1/phi))^(-1) * ((1 - mu*exp(-a)/(mu+1/phi))^(-1/phi) - 1) # E(exp(-aY)) where Y ~ HurdleNB(p, mu, phi), and p is probability of presence; from https://etd.ohiolink.edu/apexprod/rws_etd/send_file/send?accession=osu1543573678017356&disposition=inline
     -1/a * log(out)
     }

oospred_cbfm_timespacetime_hurdle_linex <- array(NA, dim = dim(ptpred_pa))
for(k0 in 1:dim(ptpred_pa)[3]) {
     oospred_cbfm_timespacetime_hurdle_linex[,,k0] <- linexpred_hurdlenb(a = 0.001, 
                                                                         p = ptpred_pa[,,k0], mu = exp(ptpred_ztnb[,,k0])+1e-12,
                                                                         phi = matrix(fitcbfm_hurdle$count_fit$dispparam, nrow = nrow(ptpred_pa), ncol = ncol(ptpred_pa), byrow = TRUE))
     }




