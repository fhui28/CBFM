##-----------------------------
## Application to NOAA fall bottom trawl survey 
## Forms part of the manuscript Hui et al., Spatio-Temporal Joint Species Distribution Modelling: A Basis Function Approach
##-----------------------------
rm(list = ls())
library(tidyverse)
library(foreach)
library(doParallel)
library(mgcv) 
library(autoFRK) 
library(FRK) 
library(CBFM) #devtools::install_github("fhui28/CBFM")
library(Hmsc) #devtools::install_version("Hmsc", version = "3.0-9", repos = "http://cran.us.r-project.org") # NB: Need to use previous version since new version has bugs with setting up XData and constructing predictions
library(corrplot)
library(gganimate)
library(patchwork)

##-----------------------------
## Load in data and set up training and test splits
##-----------------------------
load(file = "Hui_et_al_2022_NMFS_BTS.RData") # Data from 2000 to 2019 of the NOAA NEFSC fall bottom trawl survey. The data are publicly available at <https://www.fisheries.noaa.gov/inport/item/22560>. Observations with missing values for any of the 5 covariates were removed.  For more information on the cleaning of the data to produce the RData above, please contact <chaak@umass.edu> 

# Basic dimensions of data
dim(Y)
dim(XData)

summary(XData)
summary(Y)

# Set years to use for training and test
train_years <- 2000:2015
test_years <- 2016:2020

sel_training_units <- which(XData$YEAR %in% train_years)
numtrain <- length(sel_training_units)
numtrain

sel_test_units <- which(XData$YEAR %in% test_years)
numtest <- length(sel_test_units)
numtest

# Wrangling the vessel effect
# There are some "PC"s in there, where the Pisces was used in place of the Henry Bigelow (HB) in the fall of 2017. The two vessels have equivalent efficiency, so we will treat PC as HB...
# Rename "PC" to "HB"
table(XData$SVVESSEL) 
XData$SVVESSEL[which(XData$SVVESSEL=="PC")]<-"HB"

# Standardize covariates
XData$DEPTH <- XData$DEPTH %>% scale #%>% log 
XData$SURFTEMP <- XData$SURFTEMP %>% scale
XData$BOTTEMP <- XData$BOTTEMP %>% scale
XData$SURFSALIN <- XData$SURFSALIN %>% scale
XData$BOTSALIN <- XData$BOTSALIN %>% scale

summary(XData)
summary(Y)



##-----------------------------
## Fit CBFM 
##-----------------------------
# Set up spatial basis functions (50 adaptive thin-plate splines) 
num_spbasisfunctions <- 50
train_sp_basisfunctions <- mrts(XData[sel_training_units,] %>% dplyr::select(UTM_X,UTM_Y), k = num_spbasisfunctions) %>%
  as.matrix %>%
  {.[,-(1)]}
  
test_sp_basisfunctions <- mrts(XData[sel_training_units,] %>% dplyr::select(UTM_X,UTM_Y), k = num_spbasisfunctions) %>% 
  predict(newx = XData[sel_test_units,] %>% dplyr::select(UTM_X,UTM_Y)) %>%
  as.matrix %>%
  {.[,-c(1)]} 


## Temporal basis functions (7 Gaussian basis functions)
time_dat <- data.frame(year = XData$YEAR, y = 0)
table(time_dat$year)
sp::coordinates(time_dat) <- ~ year + y
time_basisfunctions <- local_basis(manifold = real_line(), loc = matrix(seq(2000,2019,by=3), ncol = 1), 
                                   type = "Gaussian", scale = rep(5,seq(2000,2019,by=3) %>% length)) 
# show_basis(time_basisfunctions, ggplot()) + 
#   geom_point(data = data.frame(time = time_dat@coords[,1], y = 0), aes(x = time, y = 0)) +
#   geom_vline(xintercept = 2015)
train_time_basisfunctions <- eval_basis(time_basisfunctions, time_dat@coords) 
test_time_basisfunctions <- train_time_basisfunctions[time_dat@coords[,1] > 2015,]
train_time_basisfunctions <- train_time_basisfunctions[time_dat@coords[,1] <= 2015,]

# Check number of the columns are entirely a vector of zeros
colSums(train_time_basisfunctions)
colSums(test_time_basisfunctions)


# Fit CBFM
myformula <- ~ SVVESSEL + s(SURFTEMP) + s(BOTTEMP) + s(SURFSALIN) + s(BOTSALIN) + s(DEPTH)
fit_cbfmspacetime <- CBFM(y = Y[sel_training_units,], data = XData[sel_training_units,], 
                    formula_X = myformula,
                    B_space =  train_sp_basisfunctions, B_time = train_time_basisfunctions, 
                    family = binomial(), 
                    control = list(trace = 1), 
                    G_control = list(rank = c(10,10)), 
                    Sigma_control = list(rank = c(5,2))) 


# Basic residual diagnostics                                  
plot(fit_cbfmspacetime)

# Basic summary to examine statistical significance of fits
summary(fit_cbfmspacetime)

# Time taken in minutes
fit_cbfmspacetime$time_taken/60



##-----------------------------
## Explore the fitted CBFM in more detail -- Covariate importance and variance partitioning
##-----------------------------

s <- summary(fit_cbfmspacetime)
s

sapply(s$summary_tables, function(x) x$smooth_terms[,4]) 
(sapply(s$summary_tables, function(x) x$smooth_terms[,4]) < 0.05) %>% apply(1,table)
# Smooths appear to be statistically significant for a reasonable number of species, except for salinity and surface temperature

sapply(s$summary_tables, function(x) x$anova_terms[,3]) 
(sapply(s$summary_tables, function(x) x$anova_terms[,3]) %>% 
  unlist %>% 
  {. < 0.05}) %>% 
  table
# Vessel type is statistically significant for many species


# Construct heatmap for species-covariates, coloring them according to whether they were statistically significant or not
predictors_significance <- rbind(
  sapply(s$summary_tables, function(x) x$anova_terms[,3] < 0.05),
  sapply(s$summary_tables, function(x) x$smooth_terms[,4] < 0.05) 
  )
rownames(predictors_significance) <- c("Survey vessel", "Surface temp.", "Bottom temp.", "Surface salinity.", "Bottom salinity.", "Depth")
colnames(predictors_significance) <-  str_replace(colnames(predictors_significance), "\\.", " ")
predictors_significance <- predictors_significance %>% 
   t %>% 
   as.data.frame %>% 
   rownames_to_column(var = "Species") %>% 
   pivot_longer(-Species) %>% 
   mutate(name = fct_inorder(name), value = factor(value))
predictors_significance$value <- fct_infreq(predictors_significance$value)
statsignificance_plot <- ggplot(predictors_significance, aes(x = name, y = Species, fill = value)) +
   geom_tile() +
   scale_fill_viridis_d() +
   coord_flip() +
   theme_bw() +
   theme(axis.text.x = element_text(angle = 70, hjust = 1), legend.position = "bottom") +
   labs(x = "Covariate", y = "Species", fill = "Statistically\nsignificant", title = "(a)")


# Variance partitioning 
v <- varpart(fit_cbfmspacetime)
colnames(fit_cbfmspacetime$betas)
groupX <- c(1,2,rep(3:7, each = 9)) # Group columns in model matrix up based on smooths
rm(v)
length(groupX)
v <- varpart(fit_cbfmspacetime, groupX = groupX)

round(v$varpart_X, 4)
round(v$varpart_B_space, 4)
round(v$varpart_B_time, 4)

# Construct a nice stacked barplot representing the variance partitioning
v_pretty <- rbind(v$varpart_X[-1,], v$varpart_B_space, v$varpart_B_time) %>% 
  as.data.frame
rownames(v_pretty) <- c("Survey vessel", "Surface temp.", "Bottom temp.", "Surface salinity", "Bottom salinity", "Depth", "Spatial basis functions", "Temporal basis functions")
colnames(v_pretty) <-  str_replace(colnames(v_pretty), "\\.", " ")
v_pretty <- v_pretty %>% 
  t %>% 
  as.data.frame %>% 
  rownames_to_column(var = "Species") %>% 
  pivot_longer(-Species)
v_pretty$name <- fct_inorder(v_pretty$name)
varpartition_plot <-ggplot(v_pretty, aes(x = Species, y = value, fill = name)) +
  geom_bar(position = "stack", stat = "identity") +
  labs(x = "Species", y = "Proportion", fill = "Model component", title = "(b)") + #Percentage of variance explained
  theme_bw() +
  scale_fill_viridis_d() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 70, hjust = 1))
# For many species, the four smoothed environment explain a heck of a lot, while the spatial basis functions explain more while the temporal basis functions explain relative little. I am not sure how much I trust these results though given the nature of the problems with potential starting values


# NB: In more recent versions of the CBFM package, the output contains some more information that will assist with visualizing smooth and parametric effects [courtesy of the gratia package!]. Please see the output all_parametric_effects and all_smooth_estimates from a CBFM fit.
# One can also obtain information regarding deviance explained, both in total and on a per-species basis, from a CBFM fit.  Please see the output deviance_explained and deviance_explained_perspecies from a CBFM fit.
# One can also apply the AIC(), AICc() function to a CBFM fit to obtain information criterion for a CBFM fit.


##-----------------------------
## Explore the fitted CBFM in more detail -- Covariance and correlation matrices
##-----------------------------
# Baseline between species correlations (take with a grain of salt!) 
corrplot(fit_cbfmspacetime$G_space, title = "Baseline correlations attributable to spatial basis functions", method = "square", type = "lower", diag = FALSE, order = "hclust", mar = c(0,0,3,5))

corrplot(fit_cbfmspacetime$G_time, title = "Baseline correlations attributable to temporal basis functions", method = "square", type = "lower", diag = FALSE, order = "hclust", mar = c(0,0,3,1))
# Spatial correlation look OK, but temporal correlation looks pretty terrible. The latter is partly due to the nature and small amount of temporal basis functions.

##-----------------------------
## Explore the fitted CBFM in more detail -- Spatio-temporal ordination
##-----------------------------
time_dat <- data.frame(time = X[sel_training_units,"YEAR"]) %>% 
   group_by(time) %>% 
   summarise(numeric_time = median(time)) # For each month-year combination, find a representative value (median in this case) in time_numeric 
sel_subsample <- sample(1:numtrain, 1000) %>% sort # Take a sub-sample of the spatial BFs to speed up computation
ord_time_BFs <- eval_basis(time_basisfunctions, as.matrix(time_dat$time)) 
ord_time_BFs <- ord_time_BFs[rep(1:nrow(ord_time_BFs), each = 1000), ]
ord_sp_BFs <- train_sp_basisfunctions[sel_subsample,]
ord_sp_BFs <- ord_sp_BFs[rep(1:1000, nrow(time_dat)),]

LVs_dat <- data.frame(time = rep(time_dat$time, each = 1000),
                      X[sel_training_units,c("LON","LAT")][rep(sel_subsample, nrow(time_dat)),],
                      lvs = ordinate(fit_cbfmspacetime, num_comp = 4, new_B_space = ord_sp_BFs, new_B_time = ord_time_BFs)$scores
                      )

LVs_dat <- LVs_dat %>% 
   pivot_longer(lvs.Axis1:lvs.Axis4)
LVs_dat$name <- fct_inorder(LVs_dat$name) 
levels(LVs_dat$name) <- paste0("Axis.", 1:4)
lvs_plot <- ggplot(LVs_dat, aes(x = LON, y = LAT, color = value)) +
  borders("world", xlim = c(-82,-60), ylim = c(30,48), fill = "lightgreen") + # From <https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.7150>
  coord_map(xlim = c(-82,-60), ylim = c(30,48)) +
  geom_point(size = 0.8, alpha = 0.5) +
  facet_wrap(. ~ name, nrow = 2) +
  scale_color_viridis_c() +
  labs(x = "Longitude", y = "Latitude", title = "Ordination scores", color = "Value") +
  theme_bw()


graph1.animation <- lvs_plot +
  transition_manual(frames = time) +
  labs(subtitle = "Date: {current_frame}") 

gganimate::animate(graph1.animation, height = 8, width = 8, units = "in", res = 100, nframes = length(unique(LVs_dat$time)), fps = 5)
# Consistent with results above, we do a clear spatial trends in the ordination scores, but the scores are hardly hardly changing in time


##-----------------------------
## Alternative model 1: Stacked GAMs 
##-----------------------------
registerDoParallel(cores = detectCores()-2)

stackedgams_fn <- function(j, y, formula_X, data) {
   fit0 <- gam(tmp_formula, data = data.frame(response = y[,j], data), method = "ML", family = binomial()) # Tried method = "REML" but made relatively little difference, although note CBFM does use REML as a default
   return(fit0)
   }
    

myformula <- ~ SVVESSEL + s(SURFTEMP) + s(BOTTEMP) + s(SURFSALIN) + s(BOTSALIN) + s(DEPTH)
stackedgams <- foreach(j = 1:ncol(Y)) %dopar% stackedgams_fn(j = j, y = Y[sel_training_units,],
                                                             formula_X = myformula, 
                                                             data = X[sel_training_units,])


# NB: We also tried stacked GLMs where covariates included linear and quadratic instead of smooth terms for the covariates. Not surprisingly this performed a lot worse and so we omitted it here. 



##-----------------------------
## Alternative model 2: Additive spatio-temporal LVM 
##-----------------------------
# Note a separate script is created here as a lot of wrangling and preparation is needed to fit the model. Please see the details of that script for more details
source(file = "spatiotemporalVMfit.R")
# Result of the above is, aside from a bunch of plots, the creation of an oospred_spacetimeLV.rds file containing the performance metrics


##-----------------------------
## Alternative model 3: CBFM but with linear and quadratic terms for covariates
##-----------------------------

# Uses the information set up above for CBFM
# Note there is an argument that one could/should include all pairwise interaction terms as well, but we have omitted that here!
myformula_quad <- ~ SVVESSEL + poly(SURFTEMP, 2, raw = TRUE) + poly(BOTTEMP, 2, raw = TRUE) + poly(SURFSALIN, 2, raw = TRUE) + poly(BOTSALIN, 2, raw=TRUE) + poly(DEPTH, 2, raw=TRUE)
fit_cbfmspacetime_quad <- CBFM(y = Y[sel_training_units,], data = XData[sel_training_units,], 
                    formula_X = myformula_quad,
                    B_space =  train_sp_basisfunctions, B_time = train_time_basisfunctions, 
                    family = binomial(), 
                    control = list(trace = 1), 
                    G_control = list(rank = c(10,10)), 
                    Sigma_control = list(rank = c(5,2))) 


##------------------------
## Compare the two CBFM fits with the stacked GAM in terms of AIC 
##------------------------
sapply(stackedgams, function(x) AIC(x)) %>% 
   sum

AIC(fit_cbfmspacetime_quad, use_edf = TRUE)

AIC(fit_cbfmspacetime, use_edf = TRUE)

##------------------------
## Compare out-of-sample performance (predicted probabilities of presence) on test data of the four models fitted above
## Recall the training set was from 2000-2015, while performance is assessed on data from 2016-2019. There are 39 demersal fish species
##------------------------
length(sel_training_units)
length(sel_test_units)


oospred_cbfm_spacetime <- predict(fit_cbfmspacetime, newdata = X[sel_test_units,], 
                                  new_B_space = test_sp_basisfunctions, new_B_time = test_time_basisfunctions, type = "response")
oospred_stackedgam <- sapply(stackedgams, function(x) predict(x, newdata = X[sel_test_units,], type = "response"))
oospred_cbfm_spacetime_quad <- predict(fit_cbfmspacetime_quad, newdata = X[sel_test_units,], 
                                       new_B_space = test_sp_basisfunctions, new_B_time = test_time_basisfunctions, type = "response")


# Out of sample predictions spatio-temporal LVMs, constructed from Hmsc fits in spatiotemporalVMfit.R script
oospred_spacetimeLV <- readRDS(file = "oospred_spacetimeLV.rds")
summary(oospred_spacetimeLV)


# Three performance metrics
compare_tjurR2 <- data.frame(
  stackedGAM = sapply(1:ncol(Y), function(j) mean(oospred_stackedgam[,j][which(Y[sel_test_units,j]==1)]) - mean(oospred_stackedgam[,j][which(Y[sel_test_units,j]==0)])),
  HMSC_multiscale = oospred_spacetimeLV$tjur_r2,
  cbfm_spacetime_quad = sapply(1:ncol(Y), function(j) mean(oospred_cbfm_spacetime_quad[,j][which(Y[sel_test_units,j]==1)]) - mean(oospred_cbfm_spacetime_quad[,j][which(Y[sel_test_units,j]==0)])),
  cbfm_spacetime = sapply(1:ncol(Y), function(j) mean(oospred_cbfm_spacetime[,j][which(Y[sel_test_units,j]==1)]) - mean(oospred_cbfm_spacetime[,j][which(Y[sel_test_units,j]==0)]))
  )

compare_auc <- data.frame(
  stackedGAM = sapply(1:ncol(Y), function(j) {pred <- ROCR::prediction(oospred_stackedgam[,j], labels = Y[sel_test_units,j]) %>% ROCR::performance(measure = "auc"); pred@y.values[[1]] } ),
  HMSC_multiscale = oospred_spacetimeLV$aucs,
  cbfm_spacetime_quad = sapply(1:ncol(Y), function(j) {pred <- ROCR::prediction(oospred_cbfm_spacetime_quad[,j], labels = Y[sel_test_units,j]) %>% ROCR::performance(measure = "auc"); pred@y.values[[1]] } ),
  cbfm_spacetime = sapply(1:ncol(Y), function(j) {pred <- ROCR::prediction(oospred_cbfm_spacetime[,j], labels = Y[sel_test_units,j]) %>% ROCR::performance(measure = "auc"); pred@y.values[[1]] } )
  )

compare_deviance <- data.frame(
  stackedGAM = sapply(1:ncol(Y), function(j) -2*mean(dbinom(Y[sel_test_units,j], 1, prob = oospred_stackedgam[,j], log = TRUE))),
  HMSC_multiscale = oospred_spacetimeLV$pred_deviance,
  cbfm_spacetime_quad = sapply(1:ncol(Y), function(j) -2*mean(dbinom(Y[sel_test_units,j], 1, prob = oospred_cbfm_spacetime_quad[,j], log = TRUE))),
  cbfm_spacetime = sapply(1:ncol(Y), function(j) -2*mean(dbinom(Y[sel_test_units,j], 1, prob = oospred_cbfm_spacetime[,j], log = TRUE)))
  )

rownames(compare_tjurR2) <- rownames(compare_auc) <- rownames(compare_deviance) <- colnames(Y)
compare_tjurR2 <- compare_tjurR2 %>% 
  rownames_to_column(var = "Species") %>% 
  pivot_longer(-Species, names_to = "method") %>% 
  mutate(method = fct_inorder(method))
compare_auc <- compare_auc %>% 
  rownames_to_column(var = "Species") %>% 
  pivot_longer(-Species, names_to = "method") %>% 
  mutate(method = fct_inorder(method))
compare_deviance <- compare_deviance %>% 
  rownames_to_column(var = "Species") %>% 
  pivot_longer(-Species, names_to = "method") %>% 
  mutate(method = fct_inorder(method))


# Comparative boxplots of the metrics
p_tjur2 <- ggplot(compare_tjurR2, aes(x = method, y = value)) +
   geom_boxplot() +
   scale_x_discrete(labels = c("Stacked GAM", "LVM", "CBFM (parametric)", "CBFM (smooth)")) +
   labs(y = "Tjur R2", x = "Method", title = "Tjur R-squared") +
   theme_bw() 

p_auc <- ggplot(compare_auc, aes(x = method, y = value)) +
   geom_boxplot() +
   labs(y = "AUC", x = "Method", title = "AUC") +
   scale_x_discrete(labels = c("Stacked GAM", "LVM", "CBFM (parametric)", "CBFM (smooth)")) +
   theme_bw()

p_deviance <- ggplot(compare_deviance, aes(x = method, y = value)) +
   geom_boxplot() +
   labs(y = "Predictive deviance", x = "Method", title = "Predictive deviance per observational unit") +
   scale_x_discrete(labels = c("Stacked GAM", "LVM", "CBFM (parametric)", "CBFM (smooth)")) +
   theme_bw() 

(p_auc / p_tjur2 / p_deviance )

