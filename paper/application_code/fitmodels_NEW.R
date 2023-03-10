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
library(Hmsc)
library(corrplot)
library(ggforce)
library(gganimate)
library(gifski)
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
## EDA to explore how covariates and species responses in all data vary over space and time -- SKIP IF EDA ALREADY DONE PREVIOUSLY! 
##-----------------------------
bigdat <- data.frame(
  year = XData$YEAR[sel_all_units],
  time = XData$UTC_TOWDATE[sel_all_units] %>% round_date(unit = "month"),
  XData[sel_all_units,] %>% dplyr::select(UTM_X, UTM_Y, LON, LAT, DEPTH, SURFTEMP, BOTTEMP, SURFSALIN, BOTSALIN),
  Y[sel_all_units,]) %>%
  pivot_longer(-c(year,time,LON,LAT,UTM_X,UTM_Y), values_to = "value")

DEPTH_plot <- ggplot(bigdat %>% dplyr::filter(name == "DEPTH"), aes(x = LON, y = LAT, color = value)) +
  borders("world", xlim = c(-82,-60), ylim = c(30,48), fill = "lightgreen") + # From <https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.7150>
  coord_map(xlim = c(-82,-60), ylim = c(30,48)) +
  geom_point(size = 0.8, alpha = 0.5) +
  #facet_wrap(. ~ year, nrow = 4) +
  scale_color_viridis_c() +
  labs(x = "Longitude", y = "Latitude", title = "Depth", color = "Value") +
  theme_bw() +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

graph1.animation <- DEPTH_plot +
  transition_manual(frames = time) +
  labs(subtitle = "Date: {current_frame}") 

gganimate::animate(graph1.animation, height = 4, width = 4, units = "in", res = 100, nframes = length(unique(bigdat$time)), fps = 5)
anim_save("plots/DEPTH_sampling.gif")


SURFTEMP_plot <- ggplot(bigdat %>% dplyr::filter(name == "SURFTEMP"), aes(x = LON, y = LAT, color = value)) +
  borders("world", xlim = c(-82,-60), ylim = c(30,48), fill = "lightgreen") + # From <https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.7150>
  coord_map(xlim = c(-82,-60), ylim = c(30,48)) +
  geom_point(size = 0.8, alpha = 0.5) +
  #facet_wrap(. ~ year, nrow = 4) +
  scale_color_viridis_c() +
  labs(x = "Longitude", y = "Latitude", title = "Surface temp.", color = "Value") +
  theme_bw() +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

graph1.animation <- SURFTEMP_plot +
  transition_manual(frames = time) +
  labs(subtitle = "Date: {current_frame}") 

gganimate::animate(graph1.animation, height = 4, width = 4, units = "in", res = 100, nframes = length(unique(bigdat$time)), fps = 5)
anim_save("plots/SURFTEMP_sampling.gif")


BOTTEMP_plot <- ggplot(bigdat %>% dplyr::filter(name == "BOTTEMP"), aes(x = LON, y = LAT, color = value)) +
  borders("world", xlim = c(-82,-60), ylim = c(30,48), fill = "lightgreen") + # From <https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.7150>
  coord_map(xlim = c(-82,-60), ylim = c(30,48)) +
  geom_point(size = 0.8, alpha = 0.5) +
  #facet_wrap(. ~ year, nrow = 4) +
  scale_color_viridis_c() +
  labs(x = "Longitude", y = "Latitude", title = "Seabed temp.", color = "Value") +
  theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

graph1.animation <- BOTTEMP_plot +
  transition_manual(frames = time) +
  labs(subtitle = "Date: {current_frame}") 

gganimate::animate(graph1.animation, height = 4, width = 4, units = "in", res = 100, nframes = length(unique(bigdat$time)), fps = 5)
anim_save("plots/BOTTEMP_sampling.gif")


SURFSALIN_plot <- ggplot(bigdat %>% dplyr::filter(name == "SURFSALIN"), aes(x = LON, y = LAT, color = value)) +
  borders("world", xlim = c(-82,-60), ylim = c(30,48), fill = "lightgreen") + # From <https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.7150>
  coord_map(xlim = c(-82,-60), ylim = c(30,48)) +
  geom_point(size = 0.8, alpha = 0.5) +
  #facet_wrap(. ~ year, nrow = 4) +
  scale_color_viridis_c() +
  labs(x = "Longitude", y = "Latitude", title = "Surface Salinity", color = "Value") +
  theme_bw() + 
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

graph1.animation <- SURFSALIN_plot +
  transition_manual(frames = time) +
  labs(subtitle = "Date: {current_frame}") 

gganimate::animate(graph1.animation, height = 4, width = 4, units = "in", res = 100, nframes = length(unique(bigdat$time)), fps = 5)
anim_save("plots/SURFSALIN_sampling.gif")


BOTSALIN_plot <- ggplot(bigdat %>% dplyr::filter(name == "BOTSALIN"), aes(x = LON, y = LAT, color = value)) +
  borders("world", xlim = c(-82,-60), ylim = c(30,48), fill = "lightgreen") + # From <https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.7150>
  coord_map(xlim = c(-82,-60), ylim = c(30,48)) +
  geom_point(size = 0.8, alpha = 0.5) +
  #facet_wrap(. ~ year, nrow = 4) +
  scale_color_viridis_c() +
  labs(x = "Longitude", y = "Latitude", title = "Bottom Salinity", color = "Value") +
  theme_bw() +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

graph1.animation <- BOTSALIN_plot +
  transition_manual(frames = time) +
  labs(subtitle = "Date: {current_frame}") 

gganimate::animate(graph1.animation, height = 4, width = 4, units = "in", res = 100, nframes = length(unique(bigdat$time)), fps = 5)
anim_save("plots/BOTSALIN_sampling.gif")


## Sampling over time plot
time_dat <- bigdat %>% 
  dplyr::filter(name == "DEPTH") %>% 
  group_by(time) %>%
  summarise(n = length(value))
time_dat$time <- ymd(time_dat$time)

ggplot(time_dat, aes(x = time, y = n)) +
  scale_x_date(date_labels = "%Y") +
  geom_bar(stat = "identity") +
  theme_bw() +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(x = "Time", y = "Number of observations")
ggsave(file = "plots/samplingovertime.pdf", width = 8, height = 4)


rm(bigdat, BOTSALIN_plot, BOTTEMP_plot, DEPTH_plot, SURFSALIN_plot, SURFTEMP_plot, graph1.animation)


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
                    formula = myformula, 
                    B_space =  train_sp_basisfunctions, B_time = as.matrix(train_time_basisfunctions), 
                    ncores = 12, # (using 12 here to be equivalent with the HMSC setup), but this can adjusted.
                    family = binomial(link = "logit"), 
                    control = list(trace = 1), 
                    G_control = list(rank = c(10,10)), 
                    Sigma_control = list(rank = c(5,2))) 


##-----------------------------
## Explore the fitted CBFM in more detail -- Basic residual diagnostics
##-----------------------------
pdf(file = "plots/CBFM_residualplots.pdf", width = 10, height = 10)
plot(fit_cbfmspacetime, pch = ".", ask = FALSE)
dev.off()
 
 
# Basic summary to examine statistical significance of fits. By default, summary.CBFM runs in parallel using (detectCores() - 1) cores to construct the summaries more efficiently. However, the command below has been set up to use only two cores for simplicity.
summary(fit_cbfmspacetime, ncores = 2)


# Time taken in minutes
fit_cbfmspacetime$time_taken/60


##-----------------------------
## Explore the fitted CBFM in more detail -- Covariate importance and variance partitioning
##-----------------------------
s <- summary(fit_cbfmspacetime, cores = 2)
s

sapply(s$summary_tables, function(x) x$smooth_terms[,4]) 
(sapply(s$summary_tables, function(x) x$smooth_terms[,4]) < 0.05) %>% apply(1,table)
# Smooths appear to be statistically significant for a reasonable number of species, except for surface salinity and surface temperature. This is not too surprising given we dealing with pelagics.

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
# For many species, the four smoothed environment explain a lot of the variation, the spatial basis functions also explain a lot, while the temporal basis functions explain relative little. 


p <- (statsignificance_plot + coord_flip()) / varpartition_plot
ggsave(p, file = "plots/effectsplots.pdf", width = 10, height = 12)


# One can also obtain information regarding deviance explained, both in total and on a per-species basis, from a CBFM fit.  Please see the output deviance_explained and deviance_explained_perspecies from a CBFM fit.
# One can also apply the AIC(), AICc() function to a CBFM fit to obtain information criterion for a CBFM fit; see below for more details.


rm(p, predictors_significance, s, statsignificance_plot, v, v_pretty, varpartition_plot)


##-----------------------------
## Explore the fitted CBFM in more detail -- Covariance and correlation matrices
##-----------------------------
# Baseline between species correlations (take with a grain of salt!) 
pdf(file = "plots/G_space.pdf", width = 8, height = 8)
corrplot(fit_cbfmspacetime$G_space, title = "Baseline correlations attributable to spatial basis functions", method = "square", type = "lower", diag = FALSE, order = "hclust", mar = c(0,0,3,5))
dev.off()

pdf(file = "plots/G_time.pdf", width = 8, height = 8)
corrplot(fit_cbfmspacetime$G_time, title = "Baseline correlations attributable to temporal basis functions", method = "square", type = "lower", diag = FALSE, order = "hclust", mar = c(0,0,3,1))
dev.off()
# Spatial correlation look OK, but temporal correlation looks pretty terrible. The latter is partly due to the nature and small amount of temporal basis functions.


##-----------------------------
## Explore the fitted CBFM in more detail -- Spatio-temporal ordination
##-----------------------------
time_dat <- data.frame(time = XData[sel_training_units,"YEAR"]) %>% 
   group_by(time) %>% 
   summarise(numeric_time = median(time)) # For each month-year combination, find a representative value (median in this case) in time_numeric 
sel_subsample <- sample(1:numtrain, 1000) %>% sort # Take a sub-sample of the spatial BFs to speed up computation
ord_time_BFs <- eval_basis(time_basisfunctions, as.matrix(time_dat$time)) 
ord_time_BFs <- ord_time_BFs[rep(1:nrow(ord_time_BFs), each = 1000), ]
ord_sp_BFs <- train_sp_basisfunctions[sel_subsample,]
ord_sp_BFs <- ord_sp_BFs[rep(1:1000, nrow(time_dat)),]

LVs_dat <- data.frame(time = rep(time_dat$time, each = 1000),
                      XData[sel_training_units,c("LON","LAT")][rep(sel_subsample, nrow(time_dat)),],
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
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


graph1.animation <- lvs_plot +
  transition_manual(frames = time) +
  labs(subtitle = "Date: {current_frame}") 

gganimate::animate(graph1.animation, height = 8, width = 8, units = "in", res = 100, nframes = length(unique(LVs_dat$time)), fps = 5)
anim_save("plots/CBFM_ordination.gif")
# Consistent with results above, we do a clear spatial trends in the ordination scores, but the scores are hardly hardly changing in time

rm(LVs_dat, lvs_plot, ord_sp_BFs, ord_time_BFs, sel_subsample, graph1.animation, time_dat)


##-----------------------------
## Alternative model 1: Stacked GAMs 
##-----------------------------
registerDoParallel(cores = detectCores()-2)

stackedgams_fn <- function(j, y, formula_X, data) {
  tmp_formula <- as.formula(paste("response", paste(as.character(formula_X),collapse="") )) 
  
  fit0 <- gam(tmp_formula, data = data.frame(response = y[,j], data), method = "ML", family = binomial()) # Tried method = "REML" but made relatively little difference, although note CBFM does use REML as a default
  return(fit0)
  }


myformula <- ~ SVVESSEL + s(SURFTEMP) + s(BOTTEMP) + s(SURFSALIN) + s(BOTSALIN) + s(DEPTH)
stackedgams <- foreach(j = 1:ncol(Y)) %dopar% stackedgams_fn(j = j, y = Y[sel_training_units,],
                                                             formula_X = myformula, 
                                                             data = XData[sel_training_units,])
saveRDS(stackedgams, file = "fit_stackedgams.rds")

# Note we also tried stacked GLMs where covariates included linear and quadratic instead of smooth terms for the covariates. Not surprisingly this performed a lot worse and so we omitted it here. 


##-----------------------------
## Alternative model 2: Additive spatio-temporal LVM 
##-----------------------------
# Note a separate script is created here as a lot of wrangling and preparation is needed to fit the model. Please see the details of that script for more details
source(file = "spatiotemporalVMfit.R")
# Result of the above is, aside from a bunch of plots and saved RDS files, is the specific creation of an oospred_spacetimeLV.rds file containing the performance metrics.


##-----------------------------
## Alternative model 3: CBFM but with linear and quadratic terms for covariates
##-----------------------------
# Uses the information set up above for CBFM
# Note there is an argument that one could/should include all pairwise interaction terms as well, but we have omitted that here!
myformula_quad <- ~ SVVESSEL + poly(SURFTEMP, 2, raw = TRUE) + poly(BOTTEMP, 2, raw = TRUE) + poly(SURFSALIN, 2, raw = TRUE) + poly(BOTSALIN, 2, raw=TRUE) + poly(DEPTH, 2, raw=TRUE)
fit_cbfmspacetime_quad <- CBFM(y = Y[sel_training_units,], data = XData[sel_training_units,], 
                    formula = myformula_quad,
                    B_space =  train_sp_basisfunctions, B_time = as.matrix(train_time_basisfunctions), 
                    family = binomial(), ncores = 6,
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
## Compare out-of-sample performance (predicted probabilities of presence)  on test data of the four models fitted above
## Recall the training set was from 2000-2015, while performance is assessed on data from 2016-2019. There are 39 demersal fish species
##------------------------
length(sel_training_units)
length(sel_test_units)


oospred_stackedgam <- sapply(stackedgams, function(x) predict(x, newdata = XData[sel_test_units,], type = "response"))
oospred_cbfm_spacetime <- predict(fit_cbfmspacetime, newdata = XData[sel_test_units,],
                                  new_B_space = test_sp_basisfunctions, new_B_time = test_time_basisfunctions, type = "response")
oospred_cbfm_spacetime_quad <- predict(fit_cbfmspacetime_quad, newdata = XData[sel_test_units,],
                                       new_B_space = test_sp_basisfunctions, new_B_time = test_time_basisfunctions, type = "response")
oospred_spacetimeLV <- readRDS(file = "oospred_spacetimeLV.rds") # Out of sample predictions spatio-temporal LVMs, constructed from Hmsc fits in spatiotemporalVMfit.R script
summary(oospred_spacetimeLV)


# Performance metrics
compare_tjurR2 <- data.frame(
    stackedGAM = sapply(1:ncol(Y), function(j) mean(oospred_stackedgam[,j][which(Y[sel_test_units,j]==1)]) - mean(oospred_stackedgam[,j][which(Y[sel_test_units,j]==0)])),
    HMSC_multiscale = oospred_spacetimeLV$tjur_r2,
    cbfm_spacetime_quad = sapply(1:ncol(Y), function(j) mean(oospred_cbfm_spacetime_quad[,j][which(Y[sel_test_units,j]==1)]) - mean(oospred_cbfm_spacetime_quad[,j][which(Y[sel_test_units,j]==0)])),
    cbfm_spacetime = sapply(1:ncol(Y), function(j) mean(oospred_cbfm_spacetime[,j][which(Y[sel_test_units,j]==1)]) - mean(oospred_cbfm_spacetime[,j][which(Y[sel_test_units,j]==0)]))
    )

compare_logscore <- data.frame(
    stackedGAM = sapply(1:ncol(Y), function(j) mean(dbinom(Y[sel_test_units,j], 1, prob = oospred_stackedgam[,j], log = TRUE))),
    HMSC_multiscale = oospred_spacetimeLV$logscore,
    cbfm_spacetime_quad = sapply(1:ncol(Y), function(j) mean(dbinom(Y[sel_test_units,j], 1, prob = oospred_cbfm_spacetime_quad[,j], log = TRUE))),
    cbfm_spacetime = sapply(1:ncol(Y), function(j) mean(dbinom(Y[sel_test_units,j], 1, prob = oospred_cbfm_spacetime[,j], log = TRUE)))
    )


compare_brierscore <- data.frame(
    stackedGAM = sapply(1:ncol(Y), function(j) mean((Y[sel_test_units,j] - oospred_stackedgam[,j])^2)),
    HMSC_multiscale = oospred_spacetimeLV$model_error,
    cbfm_spacetime_quad = sapply(1:ncol(Y), function(j) mean((Y[sel_test_units,j] - oospred_cbfm_spacetime_quad[,j])^2)),
    cbfm_spacetime = sapply(1:ncol(Y), function(j) mean((Y[sel_test_units,j] - oospred_cbfm_spacetime[,j])^2))
    )


rownames(compare_tjurR2) <- rownames(compare_logscore) <- rownames(compare_brierscore) <-colnames(Y)
compare_tjurR2 <- compare_tjurR2 %>% 
  rownames_to_column(var = "Species") %>% 
  pivot_longer(-Species, names_to = "method") %>% 
  mutate(method = fct_inorder(method))
compare_logscore <- compare_logscore %>% 
  rownames_to_column(var = "Species") %>% 
  pivot_longer(-Species, names_to = "method") %>% 
  mutate(method = fct_inorder(method))
compare_brierscore <- compare_brierscore %>% 
  rownames_to_column(var = "Species") %>% 
  pivot_longer(-Species, names_to = "method") %>% 
  mutate(method = fct_inorder(method))


# Comparative boxplots of the metrics
p_tjur2 <- ggplot(compare_tjurR2 %>% 
                     dplyr::filter(method %in% c("stackedGAM", "HMSC_multiscale", "cbfm_spacetime_quad", "cbfm_spacetime")) %>% 
                     pivot_wider(names_from = method, values_from = value) %>% 
                     mutate(stackedGAM = stackedGAM - cbfm_spacetime, 
                            HMSC_multiscale = HMSC_multiscale - cbfm_spacetime, 
                            cbfm_spacetime_quad = cbfm_spacetime_quad - cbfm_spacetime,
                            cbfm_spacetime = cbfm_spacetime - cbfm_spacetime) %>% 
                     pivot_longer(-Species, names_to = "method") %>% 
                     dplyr::filter(method %in% c("stackedGAM", "HMSC_multiscale", "cbfm_spacetime_quad")) %>% 
                     mutate(method = fct_inorder(method)), 
                  aes(x = method, y = value)) +
   geom_boxplot() +
   geom_hline(yintercept = 0, linetype = 2) +
   scale_x_discrete(labels = c("Stacked GAM", "LVM", "CBFM (parametric)")) +
   labs(y = "Differences in Tjur R2 (relative to CBFM with smoothing terms)", x = "Method", title = "Tjur R-squared") +
   theme_bw() +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

p_logscore <- ggplot(compare_logscore %>%
                         dplyr::filter(method %in% c("stackedGAM", "HMSC_multiscale", "cbfm_spacetime_quad", "cbfm_spacetime")) %>% 
                         pivot_wider(names_from = method, values_from = value) %>% 
                         mutate(stackedGAM = stackedGAM - cbfm_spacetime, 
                                HMSC_multiscale = HMSC_multiscale - cbfm_spacetime, 
                                cbfm_spacetime_quad = cbfm_spacetime_quad - cbfm_spacetime,
                                cbfm_spacetime = cbfm_spacetime - cbfm_spacetime) %>% 
                         pivot_longer(-Species, names_to = "method") %>% 
                         dplyr::filter(method %in% c("stackedGAM", "HMSC_multiscale", "cbfm_spacetime_quad")) %>% 
                         mutate(method = fct_inorder(method)), 
                  aes(x = method, y = value)) +
   geom_boxplot() +
   geom_hline(yintercept = 0, linetype = 2) +
   labs(y = "Differences in Log score (relative to CBFM with smoothing terms)", x = "Method", title = "Log score per observational unit") +
   scale_x_discrete(labels = c("Stacked GAM", "LVM", "CBFM (parametric)")) +
   theme_bw() +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

p_brierscore <- ggplot(compare_brierscore %>% 
                           dplyr::filter(method %in% c("stackedGAM", "HMSC_multiscale", "cbfm_spacetime_quad", "cbfm_spacetime")) %>% 
                           pivot_wider(names_from = method, values_from = value) %>% 
                           mutate(stackedGAM = stackedGAM - cbfm_spacetime, 
                                  HMSC_multiscale = HMSC_multiscale - cbfm_spacetime, 
                                  cbfm_spacetime_quad = cbfm_spacetime_quad - cbfm_spacetime,
                                  cbfm_spacetime = cbfm_spacetime - cbfm_spacetime) %>% 
                           pivot_longer(-Species, names_to = "method") %>% 
                           dplyr::filter(method %in% c("stackedGAM", "HMSC_multiscale", "cbfm_spacetime_quad")) %>% 
                           mutate(method = fct_inorder(method)), 
                  aes(x = method, y = value)) +
  geom_boxplot() +
   geom_hline(yintercept = 0, linetype = 2) +
   scale_y_continuous(limits = c(-0.025, 0.04)) +
   labs(y = "Differences in Brier score (relative to CBFM with smoothing terms)", x = "Method", title = "Brier score") +
   scale_x_discrete(labels = c("Stacked GAM", "LVM", "CBFM (parametric)")) +
  #coord_flip() +
  theme_bw() +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

(p_brierscore | p_tjur2 | p_logscore)
ggsave(filename = "plots/oospredictiveperformance.pdf", width = 12, height = 8)


p_tjurR2_sppprevalence <- ggplot(compare_tjurR2 %>% 
                          dplyr::filter(method %in% c("stackedGAM", "HMSC_multiscale", "cbfm_spacetime_quad", "cbfm_spacetime")) %>% 
                          pivot_wider(names_from = method, values_from = value) %>% 
                          mutate(stackedGAM = stackedGAM - cbfm_spacetime, 
                                 HMSC_multiscale = HMSC_multiscale - cbfm_spacetime, 
                                 cbfm_spacetime_quad = cbfm_spacetime_quad - cbfm_spacetime,
                                 cbfm_spacetime = cbfm_spacetime - cbfm_spacetime) %>%
                          mutate(Prevalence = colMeans(Y)) %>% 
                          pivot_longer(-c(Species, Prevalence), names_to = "method") %>% 
                          dplyr::filter(method %in% c("stackedGAM", "HMSC_multiscale", "cbfm_spacetime_quad")) %>% 
                          mutate(method = fct_inorder(method)),
                       aes(x = Prevalence, y = value, color = method, shape = method, group = method)) +
   geom_point(size = 2) +
   geom_hline(yintercept = 0, linetype = 2) +
   scale_color_viridis_d(labels = c("Stacked GAM", "LVM", "CBFM (parametric)")) +
   labs(y = "Differences in Tjur R2", x = "Species prevalence", title = "Tjur R-squared", color = "Method") +
   theme_bw() +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
   guides(shape = FALSE)


p_logscore_sppprevalence <- ggplot(compare_logscore %>% 
                          dplyr::filter(method %in% c("stackedGAM", "HMSC_multiscale", "cbfm_spacetime_quad", "cbfm_spacetime")) %>% 
                          pivot_wider(names_from = method, values_from = value) %>% 
                          mutate(stackedGAM = stackedGAM - cbfm_spacetime, 
                                 HMSC_multiscale = HMSC_multiscale - cbfm_spacetime, 
                                 cbfm_spacetime_quad = cbfm_spacetime_quad - cbfm_spacetime,
                                 cbfm_spacetime = cbfm_spacetime - cbfm_spacetime) %>%
                          mutate(Prevalence = colMeans(Y)) %>% 
                          pivot_longer(-c(Species, Prevalence), names_to = "method") %>% 
                          dplyr::filter(method %in% c("stackedGAM", "HMSC_multiscale", "cbfm_spacetime_quad")) %>% 
                          mutate(method = fct_inorder(method)),
                       aes(x = Prevalence, y = value, color = method, shape = method, group = method)) +
   geom_point(size = 2) +
   geom_hline(yintercept = 0, linetype = 2) +
   scale_color_viridis_d(labels = c("Stacked GAM", "LVM", "CBFM (parametric)")) +
   labs(y = "Differences in Log score", x = "Species prevalence", title = "Log score per observational unit", color = "Method") +
   theme_bw() +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
   guides(shape = FALSE)


p_brierscore_sppprevalence <- ggplot(compare_brierscore %>% 
                          dplyr::filter(method %in% c("stackedGAM", "HMSC_multiscale", "cbfm_spacetime_quad", "cbfm_spacetime")) %>% 
                          pivot_wider(names_from = method, values_from = value) %>% 
                          mutate(stackedGAM = stackedGAM - cbfm_spacetime, 
                                 HMSC_multiscale = HMSC_multiscale - cbfm_spacetime, 
                                 cbfm_spacetime_quad = cbfm_spacetime_quad - cbfm_spacetime,
                                 cbfm_spacetime = cbfm_spacetime - cbfm_spacetime) %>%
                          mutate(Prevalence = colMeans(Y)) %>% 
                          pivot_longer(-c(Species, Prevalence), names_to = "method") %>% 
                          dplyr::filter(method %in% c("stackedGAM", "HMSC_multiscale", "cbfm_spacetime_quad")) %>% 
                          mutate(method = fct_inorder(method)),
                       aes(x = Prevalence, y = value, color = method, shape = method, group = method)) +
   geom_point(size = 2) +
   geom_hline(yintercept = 0, linetype = 2) +
   scale_color_viridis_d(labels = c("Stacked GAM", "LVM", "CBFM (parametric)")) +
   labs(y = "Differences in Brier score", x = "Species prevalence", title = "Brier score", color = "Method") +
   theme_bw() +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
   guides(shape = FALSE)


(p_brierscore_sppprevalence / p_tjurR2_sppprevalence / p_logscore_sppprevalence) + plot_layout(guides = "collect") 
ggsave(file = "plots/oospredictiveperformancebyprevalence.pdf", width = 8, height = 10)



##---------------------------------
## Plots of the covariate effects for each species from CBFM
##---------------------------------
fit_cbfmspacetime$all_smooth_estimates$species <- fit_cbfmspacetime$all_smooth_estimates$species %>%
   fct_inorder
fit_cbfmspacetime$all_smooth_estimates <- fit_cbfmspacetime$all_smooth_estimates %>% 
   gratia::add_confint()


p <- ggplot(data = fit_cbfmspacetime$all_smooth_estimates %>% subset(smooth == "s(DEPTH)")) +
   geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = DEPTH), alpha = 0.2) +
   geom_line(aes(x = DEPTH, y = est), show.legend = FALSE) +
   geom_rug(aes(x = DEPTH), data = XData[sel_training_units,], sides = "b", color = "black") +
   facet_wrap_paginate(. ~ species, nrow = 2, ncol = 2, scales = "free") +
   labs(x = "Sea depth", y = "Effect") +
   theme_bw()

pdf(file = "plots/DEPTH_covariateplots.pdf")
for(i in 1:n_pages(p)) {
   print(p + facet_wrap_paginate(. ~ species, nrow = 2, ncol = 2, page = i, scales = "free"))
   }
dev.off()


p <- ggplot(data = fit_cbfmspacetime$all_smooth_estimates %>% subset(smooth == "s(SURFTEMP)")) +
   geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = SURFTEMP), alpha = 0.2) +
   geom_line(aes(x = SURFTEMP, y = est), show.legend = FALSE) +
   geom_rug(aes(x = SURFTEMP), data = XData[sel_training_units,], sides = "b", color = "black") +
   facet_wrap_paginate(. ~ species, nrow = 2, ncol = 2, scales = "free") +
   labs(x = "Sea surface temperature", y = "Effect") +
   theme_bw()

pdf(file = "plots/SURFTEMP_covariateplots.pdf")
for(i in 1:n_pages(p)) {
   print(p + facet_wrap_paginate(. ~ species, nrow = 2, ncol = 2, page = i, scales = "free"))
   }
dev.off()


p <- ggplot(data = fit_cbfmspacetime$all_smooth_estimates %>% subset(smooth == "s(BOTTEMP)")) +
   geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = BOTTEMP), alpha = 0.2) +
   geom_line(aes(x = BOTTEMP, y = est), show.legend = FALSE) +
   geom_rug(aes(x = BOTTEMP), data = XData[sel_training_units,], sides = "b", color = "black") +
   facet_wrap_paginate(. ~ species, nrow = 2, ncol = 2, scales = "free") +
   labs(x = "Seabed temperature", y = "Effect") +
   theme_bw()

pdf(file = "plots/BOTTEMP_covariateplots.pdf")
for(i in 1:n_pages(p)) {
   print(p + facet_wrap_paginate(. ~ species, nrow = 2, ncol = 2, page = i, scales = "free"))
   }
dev.off()


p <- ggplot(data = fit_cbfmspacetime$all_smooth_estimates %>% subset(smooth == "s(SURFSALIN)")) +
   geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = SURFSALIN), alpha = 0.2) +
   geom_line(aes(x = SURFSALIN, y = est), show.legend = FALSE) +
   geom_rug(aes(x = SURFSALIN), data = XData[sel_training_units,], sides = "b", color = "black") +
   facet_wrap_paginate(. ~ species, nrow = 2, ncol = 2, scales = "free") +
   labs(x = "Sea surface salinity", y = "Effect") +
   theme_bw()

pdf(file = "plots/SURFSALIN_covariateplots.pdf")
for(i in 1:n_pages(p)) {
   print(p + facet_wrap_paginate(. ~ species, nrow = 2, ncol = 2, page = i, scales = "free"))
   }
dev.off()


p <- ggplot(data = fit_cbfmspacetime$all_smooth_estimates %>% subset(smooth == "s(BOTSALIN)")) +
   geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = BOTSALIN), alpha = 0.2) +
   geom_line(aes(x = BOTSALIN, y = est), show.legend = FALSE) +
   geom_rug(aes(x = BOTSALIN), data = XData[sel_training_units,], sides = "b", color = "black") +
   facet_wrap_paginate(. ~ species, nrow = 2, ncol = 2, scales = "free") +
   labs(x = "Seabed salinity", y = "Effect") +
   theme_bw()

pdf(file = "plots/BOTSALIN_covariateplots.pdf")
for(i in 1:n_pages(p)) {
   print(p + facet_wrap_paginate(. ~ species, nrow = 2, ncol = 2, page = i, scales = "free"))
   }
dev.off()





sessionInfo()
# R version 4.2.2 Patched (2022-11-10 r83330)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 20.04.5 LTS
# 
# Matrix products: default
# BLAS/LAPACK: /opt/intel/compilers_and_libraries_2020.4.304/linux/mkl/lib/intel64_lin/libmkl_rt.so
# 
# locale:
#     [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8   
# [6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#     [1] parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#     [1] patchwork_1.1.2   gganimate_1.0.8   corrplot_0.92     Hmsc_3.0-13       coda_0.19-4       CBFM_0.1          TMB_1.9.2        
# [8] FRK_2.1.3         autoFRK_1.4.3     spam_2.9-1        mgcv_1.8-41       nlme_3.1-161      doParallel_1.0.17 iterators_1.0.14 
# [15] foreach_1.5.2     lubridate_1.9.1   forcats_1.0.0     stringr_1.5.0     dplyr_1.1.0       purrr_1.0.1       readr_2.1.3      
# [22] tidyr_1.3.0       tibble_3.1.8      ggplot2_3.4.0     tidyverse_1.3.2  
# 
# loaded via a namespace (and not attached):
#     [1] utf8_1.2.2           tidyselect_1.2.0     RSQLite_2.2.20       htmlwidgets_1.6.1    tweedie_2.3.5        grid_4.2.2          
# [7] pROC_1.18.0          devtools_2.4.5       munsell_0.5.0        codetools_0.2-18     interp_1.1-3         statmod_1.5.0       
# [13] gifski_1.6.6-1       miniUI_0.1.1.1       withr_2.5.0          colorspace_2.1-0     knitr_1.42           rstudioapi_0.14     
# [19] ggsignif_0.6.4       gratia_0.8.1.11      MCMCpack_1.6-3       bit64_4.0.5          farver_2.1.1         rprojroot_2.0.3     
# [25] vctrs_0.5.2          generics_0.1.3       xfun_0.36            timechange_0.2.0     R6_2.5.1             fields_14.1         
# [31] filehashSQLite_0.2-6 filematrix_1.3       cachem_1.0.6         assertthat_0.2.1     promises_1.2.0.1     scales_1.2.1        
# [37] nnet_7.3-18          googlesheets4_1.0.1  gtable_0.3.1         processx_3.8.0       ggokabeito_0.1.0     mcmc_0.9-7          
# [43] rlang_1.0.6          MatrixModels_0.5-1   splines_4.2.2        rstatix_0.7.1        gargle_1.3.0         tdigest_0.4.1       
# [49] broom_1.0.3          checkmate_2.1.0      reshape2_1.4.4       abind_1.4-5          modelr_0.1.10        backports_1.4.1     
# [55] httpuv_1.6.8         Hmisc_4.7-2          tools_4.2.2          usethis_2.1.6        ellipsis_0.3.2       gamlss.data_6.0-2   
# [61] RColorBrewer_1.1-3   sessioninfo_1.2.2    gamlss_5.4-12        Rcpp_1.0.10          plyr_1.8.8           base64enc_0.1-3     
# [67] progress_1.2.2       ps_1.7.2             prettyunits_1.1.1    ggpubr_0.5.0         rpart_4.1.19         deldir_1.0-6        
# [73] viridis_0.6.2        urlchecker_1.0.1     zoo_1.8-11           haven_2.5.1          cluster_2.1.4        fs_1.6.0            
# [79] magrittr_2.0.3       data.table_1.14.6    spacetime_1.2-8      SparseM_1.81         reprex_2.0.2         mvnfast_0.2.7       
# [85] truncnorm_1.0-8      googledrive_2.0.0    matrixStats_0.63.0   pkgload_1.3.2        hms_1.1.2            mime_0.12           
# [91] xtable_1.8-4         jpeg_0.1-10          readxl_1.4.1         gridExtra_2.3        LatticeKrig_8.4      compiler_4.2.2      
# [97] maps_3.4.1           crayon_1.5.2         htmltools_0.5.4      BayesLogit_2.1       later_1.3.0          tzdb_0.3.0          
# [103] Formula_1.2-4        filehash_2.4-3       DBI_1.1.3            tweenr_2.0.2         dbplyr_2.3.0         MASS_7.3-58.2       
# [109] Matrix_1.5-1         car_3.1-1            cli_3.6.0            gamlss.tr_5.1-7      dotCall64_1.0-2      pkgconfig_2.0.3     
# [115] numDeriv_2016.8-1.1  foreign_0.8-82       sp_1.6-0             xml2_1.3.3           gamlss.dist_6.0-5    rvest_1.0.3         
# [121] callr_3.7.3          digest_0.6.31        pracma_2.4.2         cellranger_1.1.0     intervals_0.15.2     htmlTable_2.4.1     
# [127] sparseinv_0.1.3      curl_5.0.0           shiny_1.7.4          quantreg_5.94        lifecycle_1.0.3      jsonlite_1.8.4      
# [133] carData_3.0-5        desc_1.4.2           viridisLite_0.4.1    fansi_1.0.4          pillar_1.8.1         lattice_0.20-45     
# [139] fastmap_1.1.0        httr_1.4.4           pkgbuild_1.4.0       survival_3.4-0       glue_1.6.2           xts_0.12.2          
# [145] remotes_2.4.2        FNN_1.1.3.1          png_0.1-8            bit_4.0.5            stringi_1.7.12       profvis_0.3.7       
# [151] blob_1.2.3           latticeExtra_0.6-30  memoise_2.0.1        ape_5.6-2       
