##-----------------------------
## Application to NOAA fall bottom trawl survey 
## This separate script is just for fitting and constructing predictions from a spatio-tempora LVM, fitted using HMSC. It assumes that this is only run as part of fitmodels.R
## Forms part of the manuscript Hui et al., Spatio-Temporal Joint Species Distribution Modelling: A Basis Function Approach
##-----------------------------
library(Hmsc)
library(abind)
library(geoR)
library(leaderCluster)
library(sp)
library(sf)
library(corrplot)
library(RhpcBLASctl)


##------------------------
## Prepration for multi-scale LVs 
##--------------------------
# Divide by 1000 (to get it in km) and round down to 100m to reduce precision
XData$UTM_X <- round((XData$UTM_X/1000),1)
XData$UTM_Y <- round((XData$UTM_Y/1000),1)
plot(XData$UTM_X,XData$UTM_Y)

# Are there any duplicate values in coords? (HMSC does not like coincident points)
dup.coords(XData[,c("UTM_X","UTM_Y")])

# There were a few duplicates: add a small (.1 km) jitter to only the duplicates, holding one of them fixed
XData[,c("UTM_X","UTM_Y")]<-jitterDupCoords(XData[,c("UTM_X","UTM_Y")], max=.1, min=.1, fix.one=TRUE)

# Are there any duplicates now?
dup.coords(XData[,c("UTM_X","UTM_Y")])
# No...good!

# What are the min and max distances b/w points?
summary(dist(XData[,c("UTM_X","UTM_Y")])) 
# Note these span five orders of magnitude, which can pose estimation issues for HMSC (refer to https://github.com/hmsc-r/HMSC/issues/33)

# Following the advice of package maintainers, use a multi-scale approach: a broad-scale spatial LV (random field), assigning observations to 100km grid squares, and a fine-scale LV, clustering nearby points (within 1km) at their local centroids to maintain a min distance between obs



# Make clustered coords for fine scale LV: Group observations that are within a 1km radius of each other and then assign them to a shared location at their centroid.
# This will reduce the number of unique locations and give us a minimum spacing b/w observations
clusters <- leaderCluster(XData[,c("UTM_X","UTM_Y")], radius = 1)

# How many different locations are there now?
nrow(clusters$cluster_centroids)
str(clusters)
summary(clusters$cluster_centroids)

# Check for duplicates
dup.coords(clusters$cluster_centroids)

# Check min dist b/w clusters
summary(dist(clusters$cluster_centroids))
# 1km

# Get grid point coords for each haul
cluster_id_coords<-clusters$cluster_centroids[clusters$cluster_id,]
dim(cluster_id_coords)
plot(cluster_id_coords)

# Mutate cluster_id and cluster coords into XData
XData$cluster_id <- clusters$cluster_id
XData$UTMX_cluster <- cluster_id_coords[,1]
XData$UTMY_cluster <- cluster_id_coords[,2]

# Plot original locations, and then overlaying new (clustered) locations
plot(XData$UTM_X, XData$UTM_Y, col = XData$cluster_id) 
points(clusters$cluster_centroids[,1], clusters$cluster_centroids[,2], col="red", pch=3, cex=0.5)


# Now make coords for broad scale gridded LV
xycoords_sp <- SpatialPoints(XData[,c("UTM_X","UTM_Y")])
xycoords_sf <- st_as_sf(xycoords_sp)

# Create grid
grid <- makegrid(xycoords_sp, cellsize = 100) # 100 km grid
grid_sp <-SpatialPoints(grid)
grid_sf <- st_as_sf(grid_sp)

# Distances b/w grid points
summary(dist(grid))

# Get nearest point in grid to xycoords
grid_id <- st_nearest_feature(xycoords_sf,grid_sf)
str(grid_id)
length(grid_id) #5892 obs
length(unique(grid_id)) #48 cels

# Plot grid cell centers, and original coordinates colored according to their grid cel membership
plot(XData$UTM_X, XData$UTM_Y, col = XData$grid_id)
points(grid, col = "red", pch = 3, cex = 0.5)

# Get grid point coords for each haul
grid_id_coords <- grid[grid_id,]
dim(grid_id_coords)
plot(grid_id_coords)

# Add a grid id column in Xdata
XData$grid_id<-grid_id

# Add grid coords in Xdata
XData$UTMX_grid <- grid_id_coords[,1]
XData$UTMY_grid <- grid_id_coords[,2]

head(XData)


# Split data in training and fit [for HMSC script only]
X.train <- XData[sel_training_units,]
X.test <- XData[sel_test_units,]
Y.train <- Y[sel_training_units,]
Y.test <- Y[sel_test_units,]

summary(XData)
summary(X.train)
dim(X.train)
dim(X.test)


##-----------------------------
## Define study design
##-----------------------------
# Training data
studyDesign <- data.frame (Haul = factor(X.train$HAULID),
                           Cluster = factor(X.train$cluster_id),
                           Grid = factor(X.train$grid_id),
                           Year = factor(X.train$YEAR)
                           )
head(studyDesign)                


##---------------------------------------
## Set up spatial & temporal coordinates
##---------------------------------------
# For coarse scale spatial grid
Grid_coords <- cbind.data.frame (studyDesign$Grid, X.train[,c("UTMX_grid","UTMY_grid")]) %>%
  unique()
row.names(Grid_coords) <- Grid_coords[,1]
Grid_coords[1] <- NULL

# For fine scale space with clustering
Cluster_coords <- cbind.data.frame (studyDesign$Cluster, X.train[,c("UTMX_cluster","UTMY_cluster")]) %>%
  unique()
row.names(Cluster_coords) <- Cluster_coords[,1]
Cluster_coords[1] <- NULL

# For time
Year_coords <- cbind.data.frame(studyDesign$Year,X.train[,"YEAR"]) %>% 
  unique()
row.names(Year_coords) <- Year_coords[,1]
Year_coords[1] <- NULL


##-----------------------------
## Define spatial latent variables
##-----------------------------
# Broad scale spatial LVs using 100km grid coords -- exponential spatial structure with 3 LVs
rL.grid <- HmscRandomLevel(sData = Grid_coords)
rL.grid <- setPriors(rL.grid, nfMin = 3, nfMax = 3) 
rL.grid

# Constrain prior to broad scales (no smaller than the 100km grid resolution)
rL.grid$alphapw
#...make it uniform  across the remaining distances
rL.grid$alphapw <- rL.grid$alphapw[-1,] # drop first row of the prior
mindist <- 100 # set to min grid distance
maxdist <- max(rL.grid$alphapw[,1]) # use max distance in default prior (which is the max dist b/w any 2 points)
nrow(rL.grid$alphapw)

# Now replace first column with a vector of distances from mindist to maxdist in equal intervals
rL.grid$alphapw[,1] <- c(seq(from=mindist, to=maxdist, length.out=100))
# Replace probabilities, making it uniform
rL.grid$alphapw[,2] <- c(rep(1/100, 100))
sum(rL.grid$alphapw[,2])
rL.grid$alphapw

# Fine scale spatial random effect using clustered coordinates -- NNGP spatial structure 20 neighbors and 2 LVs
rL.cluster.nngp <- HmscRandomLevel(sData = Cluster_coords, sMethod = 'NNGP', nNeighbours = 20)
rL.cluster.nngp <- setPriors(rL.cluster.nngp, nfMin = 2, nfMax = 2) 
rL.cluster.nngp

# Constrain scale to be less than that of the broad-scale effect (< 100km)
rL.cluster.nngp$alphapw
nrow(rL.cluster.nngp$alphapw)
maxdist <- 100 # set maximum scale of 100km
# Replace first column with distances from 0-100, with equal intervals
rL.cluster.nngp$alphapw[,1] <- seq(from=0, to=maxdist, length.out=101)
sum(rL.cluster.nngp$alphapw[,2])

rL.cluster.nngp$alphapw
rL.cluster.nngp

 
##-----------------------------
## Define temporal latent variables
##-----------------------------
# Structured temporal effect at the level of years
rL.temporal.year <- HmscRandomLevel(sData = Year_coords)
rL.temporal.year <- setPriors(rL.temporal.year, nfMin = 2, nfMax = 2) # limit to 2 LVs


##-----------------------------
## Define and fit model 
##-----------------------------
#formula
myformula_quad <- ~ SVVESSEL + poly(SURFTEMP, 2, raw = TRUE) + poly(BOTTEMP, 2, raw = TRUE) + poly(SURFSALIN, 2, raw = TRUE) +  poly(BOTSALIN, 2, raw=TRUE) + poly(DEPTH, 2, raw=TRUE)


# NB: There is some weird error about non numeric variables in X.train, even though "str" says they are numeric. For now, duplicate needed columns into new object and force numeric
X.train2 <- X.train %>% 
  dplyr::select(SVVESSEL, DEPTH:BOTSALIN) %>%
  as.data.frame

head(X.train2)

X.train2[,c(2:6)] <- sapply(X.train2[,c(2:6)], as.numeric)


# With Cluster for finescale RE (min dist b/w obs = 1km)
lvmmod_spatiotemp_multiscale_clust <- Hmsc(Y = Y.train, XData = X.train2, 
                                     XFormula = myformula_quad, distr = "probit",
                                     studyDesign = studyDesign, 
                                     ranLevels = list("Grid" = rL.grid, "Cluster" = rL.cluster.nngp, "Year" = rL.temporal.year))


# Setting a few things to maximize speed: we used 4 openmp cores per chain - may need to be changed depending on your machine...
blas_set_num_threads(1)
omp_set_num_threads(4)


# Sample: 10,000 burnin and 1000 samples with thinning at 100, so each chain is 110,000 long
tic <- proc.time()
fit_spatiotemp_multiscale_clust <- sampleMcmc(lvmmod_spatiotemp_multiscale_clust, thin = 100, samples = 1000, transient = 10000,
                                        nChains = 3,  nParallel = 3, useSocket = FALSE, verbose = 10, updater = list(GammaEta = FALSE))
toc <- proc.time()
toc - tic


saveRDS(fit_spatiotemp_multiscale, file = "LVM_fit_multiscale.rds")


##-----------------------------
## Explore fitted spatio-temporal LVM
##-----------------------------
# Extract the posterior distribution from the model object and convert it into a coda object.
mpost <- convertToCodaObject(fit_spatiotemp_multiscale_clust)

# Gelman-Rubin statistics/rhats
gelman.diag(mpost$Beta, multivariate=FALSE)$psrf #betas 
gelman.diag(mpost$Alpha[[1]], multivariate=FALSE)$psrf # spatial scale parameter for broad
gelman.diag(mpost$Alpha[[2]], multivariate=FALSE)$psrf # spatial scale parameter for fine
gelman.diag(mpost$Alpha[[3]], multivariate=FALSE)$psrf # temporal scale parameter
# how about lambdas (loadings??)
gelman.diag(mpost$Lambda[[1]])$psrf # loadings for space-broad
gelman.diag(mpost$Lambda[[2]])$psrf # loadings for space-fine
gelman.diag(mpost$Lambda[[3]])$psrf # loadings for time
# Notes from Chris Haak's exploration 21/02/23: For the "cluster" fit, both fine-scale spatial LVs initialize, but the issue of bad mixing still persists...  that is, for some reason, in 1 (of the 3) chains, the estimates of alpha (scale parameter) for the 2 fine-scale spatial LVs look like they are flip-flopped (so that two of the chains agree, but the other chain has LV1 taking roughly the value of LV2 in the other 2 chains, and vice versa).  I still have no clue why this is happening, or how to fix it.  I could try removing the heavy zero-weighting on the priors for alpha (but that is their default).


# Some MCMC plots
par(mfrow = c(3,2))

psrf.beta <- gelman.diag(mpost$Beta, multivariate = FALSE)$psrf
ess.beta <- effectiveSize(mpost$Beta)
hist(ess.beta, breaks = 50)
hist(psrf.beta, breaks = 50)

ess.alpha.space <- effectiveSize(mpost$Alpha[[1]])
psrf.alpha.space <- gelman.diag(mpost$Alpha[[1]], multivariate = FALSE)$psrf
hist(ess.alpha.space, breaks = 50)
hist(psrf.alpha.space, breaks = 50)

ess.alpha.time <- effectiveSize(mpost$Alpha[[2]])
psrf.alpha.time <- gelman.diag(mpost$Alpha[[2]], multivariate = FALSE)$psrf
hist(ess.alpha.time, breaks = 50)
hist(psrf.alpha.time, breaks = 50)


# Trace-plots
plot(mpost$Alpha[[1]]) # alpha for broad-scale spatial LV
summary(mpost$Alpha[[1]]) 

plot(mpost$Alpha[[2]]) # alpha for fine scale LV 
summary(mpost$Alpha[[2]]) 

plot(mpost$Alpha[[3]]) #alpha for temporal LV 
summary(mpost$Alpha[[3]]) 

dev.off()


##-----------------------------
## Predict to test data
##-----------------------------
# Set up space and time coords for test set REs
Grid.test <- X.test[,c("UTMX_grid","UTMY_grid")] %>% as.matrix ()
Cluster.test <- X.test[,c("UTMX_cluster","UTMY_cluster")] %>% as.matrix ()
Haul.test <- X.test[,c("UTM_X","UTM_Y")] %>% as.matrix ()
Year.test <- X.test["YEAR"] %>% as.matrix ()


#Set up gradient for prediction
Gradient_clust <- prepareGradient(hM = fit_spatiotemp_multiscale_clust, XDataNew = X.test, 
                                  sDataNew = list(Grid = Grid.test, Cluster = Cluster.test, Year = Year.test)
                                  )


#Generate predictions using gradient, pooling all 1000 samples of each chain 
#note - we could use just the last n samples as well, by changing "start" value
testpreds_clust <- predict(object = fit_spatiotemp_multiscale_clust, Gradient = Gradient_clust, predictEtaMean = TRUE,expected = TRUE, useSocket = FALSE,
                           post = poolMcmcChains(postList = fit_spatiotemp_multiscale_clust$postList, start = 1))


saveRDS(testpreds_clust, file = "HMSC_oospreds_clust.rds")


##----------------------------
## Assess OOS Prediction Performance
##----------------------------
# Compute three performance metrics, the same as in the fitmodels.R script
getpreds <- testpreds_clust %>%
  abind(along = 3) %>%
  apply(., c(1,2), mean, na.rm = TRUE)

tjur_r2<- sapply(1:ncol(Y.test), function(j) { 
   m1 <- getpreds[which(Y.test[,j] == 1),j] %>% mean(na.rm = TRUE)
   m0 <- getpreds[which(Y.test[,j] == 0),j] %>% mean(na.rm = TRUE)
   m1 - m0     
   })

aucs <- sapply(1:ncol(Y.test), function(j) { 
   pred <- ROCR::prediction(getpreds[,j], labels = Y.test[,j]) %>% ROCR::performance(measure = "auc")
   pred@y.values[[1]]
   })

logscore <- sapply(1:ncol(Y.test), function(j) {
   mean(dbinom(Y.test[,j], size = 1 , prob = getpreds[,j], log = TRUE))
   })

model_error <- sapply(1:ncol(Y.test), function(j) {
   mean((Y.test[,j] - getpreds[,j])^2)
   })


oospred_spacetimeLV <- data.frame(tjur_r2, aucs, logscore, model_error)
summary(oospred_spacetimeLV)

saveRDS(oospred_spacetimeLV, file = "oospred_spacetimeLV.rds")


##----------------------
## Variance partitioning
## Using a conditional approach to be consistent with what CBFM does.
##----------------------
allBeta <- lapply(poolMcmcChains(fit_spatiotemp_multiscale$postList), function(a) a[["Beta"]]) %>% abind(along = 3)
allXBeta <- sapply(1:length(fit_spatiotemp_multiscale$postList), function(k) (fit_spatiotemp_multiscale$XScaled %*% allBeta[,,k]) %>% apply(., 2, var))
rm(allBeta)


fit_spatiotemp_multiscale$studyDesign$Haul %>% unique %>% length
fit_spatiotemp_multiscale$studyDesign$Cluster %>% unique %>% length
fit_spatiotemp_multiscale$studyDesign$Grid %>% unique %>% length
fit_spatiotemp_multiscale$studyDesign$Year %>% unique %>% length


allEta <- lapply(poolMcmcChains(fit_spatiotemp_multiscale$postList), function(a) a[["Eta"]][[1]]) %>% abind(along = 3)
allLambda <- lapply(poolMcmcChains(fit_spatiotemp_multiscale$postList), function(a) a[["Lambda"]][[1]]) %>% abind(along = 3)
allEtaLambda1 <- sapply(1:length(fit_spatiotemp_multiscale$postList), function(k) (allEta[fit_spatiotemp_multiscale$studyDesign$Grid,,k] %*% allLambda[,,k]) %>% apply(., 2, var))

allEta <- lapply(poolMcmcChains(fit_spatiotemp_multiscale$postList), function(a) a[["Eta"]][[2]]) %>% abind(along = 3)
allLambda <- lapply(poolMcmcChains(fit_spatiotemp_multiscale$postList), function(a) a[["Lambda"]][[2]]) %>% abind(along = 3)
allEtaLambda2 <- sapply(1:length(fit_spatiotemp_multiscale$postList), function(k) (allEta[fit_spatiotemp_multiscale$studyDesign$Haul,,k] %*% allLambda[,,k]) %>% apply(., 2, var))

allEta <- lapply(poolMcmcChains(fit_spatiotemp_multiscale$postList), function(a) a[["Eta"]][[3]]) %>% abind(along = 3)
allLambda <- lapply(poolMcmcChains(fit_spatiotemp_multiscale$postList), function(a) a[["Lambda"]][[3]]) %>% abind(along = 3)
allEtaLambda3 <- sapply(1:length(fit_spatiotemp_multiscale$postList), function(k) (allEta[fit_spatiotemp_multiscale$studyDesign$Year,,k] %*% allLambda[,,k]) %>% apply(., 2, var))
rm(allEta, allLambda)

var_Xbeta <- allXBeta/(allXBeta + allEtaLambda1 + allEtaLambda2 + allEtaLambda3)
var_LV <- (allEtaLambda1 + allEtaLambda2 + allEtaLambda3)/(allXBeta + allEtaLambda1 + allEtaLambda2 + allEtaLambda3)
var_Xbeta %>% rowMeans
var_LV %>% rowMeans

var_Xbeta %>% rowMeans %>% summary
var_LV %>% rowMeans %>% summary



##-----------------------
## Residual correlations
##----------------------
ResCorrs <- computeAssociations(fit_spatiotemp_multiscale_clust)
str(ResCorrs)

# Mean residual corrs for coarse-scale spatial RE
corrplot(ResCorrs[[1]]$mean, method = "color",
         col = colorRampPalette(c("blue","white","red"))(200),
         tl.cex = 0.6, 
         tl.col = "black",
         title = paste("Random effect level:", fit_spatiotemp_multiscale_clust$rLNames[1]), mar = c(0,0,1,0))

# Mean residual corrs for fine-scale spatial RE
corrplot(ResCorrs[[2]]$mean, method = "color",
         col = colorRampPalette(c("blue","white","red"))(200),
         tl.cex = 0.6, 
         tl.col = "black",
         title = paste("Random effect level:", fit_spatiotemp_multiscale_clust$rLNames[2]), mar=c(0,0,1,0))

# Mean residual corrs for temporal (year) RE
corrplot(ResCorrs[[3]]$mean, method = "color",
         col = colorRampPalette(c("blue","white","red"))(200),
         tl.cex = 0.6, 
         tl.col = "black",
         title = paste("Random effect level:", fit_spatiotemp_multiscale_clust$rLNames[3]), mar=c(0,0,1,0))





