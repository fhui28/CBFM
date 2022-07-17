##-----------------------------
## Application to NOAA fall bottom trawl survey 
## This separate script is just for fitting and constructing predictions from a spatio-tempora LVM, fitted using HMSC 
## Forms part of the manuscript Hui et al., Spatio-Temporal Joint Species Distribution Modelling: A Basis Function Approach
##-----------------------------
library(Hmsc)
library(abind)
library(geoR)
library(leaderCluster)
library(sp)
library(sf)

##------------------------
## Prepration for multi-scale LVs 
##--------------------------
# Divide by 1000 (to get it in km)
XData$UTM_X <- XData$UTM_X/1000
XData$UTM_Y <- XData$UTM_Y/1000
plot(XData$UTM_X,XData$UTM_Y)

# Are there any duplicate values in coords? (HMSC does not like coincident points)
dup.coords(XData[,c("UTM_X","UTM_Y")])

# What are the min and max distances b/w points?
summary(dist(XData[,c("UTM_X","UTM_Y")])) 
#Note these span roughly six orders of magnitude, which can pose estimation issues for HMSC (refer to https://github.com/hmsc-r/HMSC/issues/33)

## Following the advice of package maintainers, use a multi-scale approach: a broad-scale spatial LV (random field), assigning observations to 100km grid squares, and a fine-scale LV, clustering nearby points (within 1km) at their local centroids to maintain a min distance between obs




# Make clustered coords for fine scale LV: Group observations that are within a 1km radius of each other and then assign them to a shared location at their centroid 
# this will reduce the number of unique locations and give us a minimum spacing b/w observations
clusters <- leaderCluster(XData[,c("UTM_X","UTM_Y")], radius = 1)

# How many different locations are there now?
nrow(clusters$cluster_centroids)

summary(clusters$cluster_centroids)

# Check for duplicates
dup.coords(clusters$cluster_centroids)

# Check min dist b/w clusters
summary(dist(clusters$cluster_centroids))

# Add a cluster id column in Xdata
XData$cluster_id <- clusters$cluster_id
length(XData$cluster_id)
length(unique(XData$cluster_id))

# Plot original locations, and then overlaying new (clustered) locations
plot(XData$UTM_X, XData$UTM_Y, col = XData$cluster_id) 
points(clusters$cluster_centroids[,1], clusters$cluster_centroids[,2], col="red", pch=3, cex=0.5)


# Make coords for broad scale gridded LV
xycoords_sp <- SpatialPoints(XData[,c("UTM_X","UTM_Y")])
xycoords_sf <- st_as_sf(xycoords_sp)

# Create grid
grid <- makegrid(xycoords_sp, cellsize = 100) # 100 km grid
grid_sp <-SpatialPoints(grid)
grid_sf <- st_as_sf(grid_sp)

# Get nearest point in grid to xycoords
grid_id <- st_nearest_feature(xycoords_sf,grid_sf)
str(grid_id)
length(unique(grid_id))

# Distances b/w grid points
summary(dist(grid))

# Add a grid id column in Xdata
XData$grid_id<-grid_id

# Plot grid cell centers, and original coordinates colored according to their grid cel membership
plot(XData$UTM_X, XData$UTM_Y, col = XData$grid_id)
points(grid, col="red", pch=3, cex=0.5)


# Split data in training and fit [for HMSC script only]
X.train <- XData[sel_training_units,]
X.test <- XData[sel_test_units,]
Y.train <- Y[sel_training_units,]
Y.test <- Y[sel_test_units,]


##-----------------------------
## Define study design
##-----------------------------
# Training data
studyDesign <- data.frame(sample = as.factor(rownames(X.train)))
studyDesign$Year <- as.factor(X.train$YEAR)
studyDesign$gridfac <- as.factor(X.train$grid_id)
studyDesign$clusterfac <- as.factor(X.train$cluster_id)
str(studyDesign)


# Test data
studyDesign.test <- data.frame(sample = as.factor(rownames(X.test)))
studyDesign.test$Year <- as.factor(X.test$YEAR)
studyDesign.test$gridfac <- as.factor(X.test$grid_id)
studyDesign.test$clusterfac <- as.factor(X.test$cluster_id)
str(studyDesign.test)

##-----------------------------
## Make spatial coordinates
##-----------------------------
# Coarse scale spatial grid coordinates

# Get ids of clusters in training data
traingrids <- unique(studyDesign$grid)

# Pull corresponding centroids from clusters object and assign rownames
xycoords_grid <- grid[traingrids,]
rownames(xycoords_grid)<-traingrids
head(xycoords_grid)
nrow(xycoords_grid)

# Check distances and repeats to confirm
dup.coords(xycoords_grid)
summary(dist(xycoords_grid))
summary(xycoords_grid)

# Do the same for test data
testgrids <- unique(studyDesign.test$grid)
xycoords_grid.test <- grid[testgrids,]
rownames(xycoords_grid.test) <- testgrids



# Fine scale spatial cluster coordinates

# Get ids of clusters in training data
trainclusters <- unique(studyDesign$cluster)

# Pull corresponding centroids from clusters object and assign rownames
xycoords_clusters <- clusters$cluster_centroids[trainclusters,]
rownames(xycoords_clusters) <- trainclusters
head(xycoords_clusters)

# Check distances and repeats
dup.coords(xycoords_clusters)
summary(dist(xycoords_clusters))
summary(xycoords_clusters)

# Do the same for the test data
testclusters <- unique(studyDesign.test$cluster)
xycoords_clusters.test <- clusters$cluster_centroids[testclusters,]
rownames(xycoords_clusters.test) <- testclusters


##-----------------------------
## Make temporal coordinates
##-----------------------------
time_dat <- X.train$YEAR
table(time_dat)

# Make a matrix with integer values from 1 to length of years
yearcoords <- matrix(unique(time_dat)) 
yearcoords
rownames(yearcoords) <- levels(studyDesign$Year)
head(yearcoords)


# Do the same for test data
time_dat.test <- X.test$YEAR
table(time_dat.test)
yearcoords.test <- matrix(unique(time_dat.test))
rownames(yearcoords.test) <- levels(studyDesign.test$Year)


##-----------------------------
## Define spatial latent variables
##-----------------------------
# Broad scale spatial LVs using 100km grid coords -- exponential spatial structure with 3 LVs
rL.grid = HmscRandomLevel(sData = xycoords_grid)
rL.grid = setPriors(rL.grid, nfMin = 3, nfMax = 3) 
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

# Do the same For test data (but don't need to define priors for test set)
rL.grid.test <- HmscRandomLevel(sData = rbind(xycoords_grid, xycoords_grid.test))


# Fine scale spatial random effect using clustered coordinates -- NNGP spatial structure 20 neighbors and 2 LVs
rL.cluster.nngp <- HmscRandomLevel(sData = xycoords_clusters, sMethod = 'NNGP', nNeighbours = 20)
rL.cluster.nngp <- setPriors(rL.cluster.nngp, nfMin = 2, nfMax = 2) 
rL.cluster.nngp


# Constrain scale to be less than that of the broad-scale effect (< 100km)
rL.cluster.nngp$alphapw
nrow(rL.cluster.nngp$alphapw)
maxdist<-100 # set maximum scale of 100km
# Replace first column with distances from 0-100, with equal intervals
rL.cluster.nngp$alphapw[,1] <- seq(from=0, to=maxdist, length.out=101)
sum(rL.cluster.nngp$alphapw[,2])

rL.cluster.nngp$alphapw
rL.cluster.nngp

# Do the same For test data (but don't need to define priors for test set)
rL.cluster.nngp.test <- HmscRandomLevel(sData = rbind(xycoords_clusters, xycoords_clusters.test), sMethod = 'NNGP', nNeighbours = 20)


##-----------------------------
# Define temporal latent variables
##-----------------------------
# Structured temporal effect at the level of years
rL.temporal.year <- HmscRandomLevel(sData = yearcoords)
rL.temporal.year <- setPriors(rL.temporal.year, nfMin = 2, nfMax = 2) # limit to 2 LVs

# Do the same For test data (but don't need to define priors for test set)
rL.temporal.year.test <- HmscRandomLevel(sData = rbind(yearcoords, yearcoords.test))


##-----------------------------
# Define and fit model 
##-----------------------------
myformula_quad <- ~ SVVESSEL + poly(SURFTEMP, 2, raw = TRUE) + poly(BOTTEMP, 2, raw = TRUE) + poly(SURFSALIN, 2, raw = TRUE) + poly(BOTSALIN, 2, raw=TRUE) + poly(DEPTH, 2, raw=TRUE)

lvmmod_spatiotemp_multiscale <- Hmsc(Y = Y.train, XData = X.train %>% dplyr::select(SVVESSEL, DEPTH:BOTSALIN),
                                     XFormula = myformula_quad, distr = "probit",
                                     studyDesign = studyDesign, 
                                     ranLevels = list("gridfac" = rL.grid, "clusterfac" = rL.cluster.nngp, "Year" = rL.temporal.year))


# Sample
tic <- proc.time()
fit_spatiotemp_multiscale <- sampleMcmc(lvmmod_spatiotemp_multiscale, thin = 100, samples = 1000, transient = 10000,
                                      nChains = 3,  nParallel = 3, verbose = verbose,updater=list(GammaEta=FALSE))

toc <- proc.time()
toc - tic



##-----------------------------
# Explore fitted spatio-temporal LVM
##-----------------------------

#extract the posterior distribution from the model object and convert it into a coda object.
mpost <- convertToCodaObject(fit_spatiotemp_multiscale)

# Gelman-Rubin statistics/rhats
gelman.diag(mpost$Beta, multivariate=FALSE)$psrf #betas 
gelman.diag(mpost$Alpha[[1]], multivariate=FALSE)$psrf # spatial scale parameter for broad
gelman.diag(mpost$Alpha[[2]], multivariate=FALSE)$psrf # spatial scale parameter for fine
gelman.diag(mpost$Alpha[[3]], multivariate=FALSE)$psrf # temporal scale parameter
# how about lambdas (loadings??)
gelman.diag(mpost$Lambda[[1]])$psrf # loadings for space-broad
gelman.diag(mpost$Lambda[[2]])$psrf # loadings for space-fine
gelman.diag(mpost$Lambda[[3]])$psrf # loadings for time


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


# Residual correlations
ns <- ncol(myfit$Y) # Set NS to number of species
ns
sppairs <- matrix(sample(x = 1:ns^2, size = 100))
tmp <- mpost$Omega[[1]]
for (chain in 1:length(tmp)) {
   tmp[[chain]] <- tmp[[chain]][,sppairs]
   }
ess.omega <- effectiveSize(tmp)
psrf.omega <- gelman.diag(tmp, multivariate=FALSE)$psrf
hist(ess.omega, breaks=50)
hist(psrf.omega, breaks=50)


# Trace-plots
plot(mpost$Alpha[[1]]) # alpha for broad-scale spatial LV
summary(mpost$Alpha[[1]]) 

plot(mpost$Alpha[[2]]) # alpha for fine scale LV 
summary(mpost$Alpha[[2]]) 

plot(mpost$Alpha[[3]]) #alpha for temporal LV 
summary(mpost$Alpha[[3]]) 

dev.off()


##-----------------------------
# Predict to test data
##-----------------------------
testpreds <- predict(fit_spatiotemp_multiscale, XData = X.test, studyDesign = studyDesign.test, 
                     ranLevels = list("gridfac" = rL.grid.test, "clusterfac" = rL.cluster.nngp.test, "Year" = rL.temporal.year.test), 
                     expected = TRUE)


# Compute three performance metrics, the same as in the fitmodels.R script
getpreds <- testpreds %>%
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

pred_deviance <- sapply(1:ncol(Y.test), function(j) {
   -2*mean(dbinom(Y.test[,j], size = 1 , prob = getpreds[,j], log = TRUE))
   })


oospred_spacetimeLV <- data.frame(tjur_r2, aucs, pred_deviance)

saveRDS(oospred_spacetimeLV, file = "oospred_spacetimeLV.rds")


rm(list = ls())


