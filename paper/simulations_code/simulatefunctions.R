##------------------------
## Functions to simulate covariates and spatial multivariate abundance data
##--------------------------
sim_covariates <- function(n = 1250, missing_type = "stationary", longitude_range = c(145, 150), latitude_range = c(-40, -35), num_blocks = 25, seed = NULL) {
     
     set.seed(seed)
     xy <- data.frame(x = runif(n, longitude_range[1], longitude_range[2]), y = runif(n, latitude_range[1], latitude_range[2]))
     
     # Form spatial blocks based on spatial range, based on checkerboard style
     if(num_blocks > 0) {
          actual_num_blocks <- ceiling(sqrt(num_blocks))^2
          split_longitude <- seq(longitude_range[1], longitude_range[2], length = sqrt(actual_num_blocks)+1)
          split_latitude <- seq(latitude_range[1], latitude_range[2], length = sqrt(actual_num_blocks)+1)
          index_longitude <- cut(xy$x, breaks = split_longitude, labels = 1:(length(split_longitude)-1)) %>%
               as.numeric
          index_latitude <- cut(xy$y, breaks = split_latitude, labels = 1:(length(split_latitude)-1)) %>%
               as.numeric
          xy$block_number <- (index_latitude-1)*sqrt(actual_num_blocks) + index_longitude
          }
 
     xy$temp <- rnorm(n) %>%
          scales::rescale(to = range(16:20)) %>%
          scale

     xy$depth <- rchisq(n, df = 3) %>%
          scales::rescale(to = range(10:100)) %>%
          scale
               
     if(missing_type == "stationary") {
          xy$chla <- RFsimulate(model = RMexp(var=1, scale=0.1*n), x = xy$x, y = xy$y, n = 1)@data[,1] %>%
               scale
          
          xy$O2 <- RFsimulate(model = RMexp(var=1, scale=0.5*n), x = xy$x, y = xy$y, n = 1)@data[,1] %>%
               scale
          }
     if(missing_type == "nonstationary") {
          makecovs <- mrts(xy[,1:2], k = 15)[,-1]
     
          xy$chla <- rnorm(n, mean = makecovs[,3]-makecovs[,5]+makecovs[,10], sd = 0.1) %>%
               scale
     
          xy$O2 <- rnorm(n, mean = makecovs[,2]-makecovs[,4]+makecovs[,8], sd = 0.1) %>%
               scale          
          }
     
     
     
     set.seed(NULL)
     return(list(data = xy, split_longitude = split_longitude, split_latitude = split_latitude))
     }



create_life <- function(env, formula_X, family, spp_slopes, spp_intercepts = NULL, seed = NULL) {     
     
     #--------------------
     # Checks and balances
     #--------------------
     if(class(env)[1] != "SpatialPointsDataFrame")
          stop("env should be an object for class \"SpatialPointsDataFrame\".")
     if(!(family$family %in% c("binomial","poisson")))
          stop("Family is currently not supported...sorry!")
     
     num_units <- nrow(env@coords)
     num_spp <- nrow(spp_slopes)
     sim_dat <- as.matrix(model.matrix(formula_X, data = env@data))
     
     
     #--------------------
     # Species-sepcific intercept values (overall prevelance of specie), based on values observed for species in Kerguelen PLateau RCP models)
     #--------------------
     if(is.null(spp_intercepts)) {
          #set.seed(2021)
          #betamean <- 0.2
          #betabeta <- 15
          #betaalpha <- betamean/(1-betamean) * betabeta
          #prevalences <- rbeta(num_spp, betaalpha, betabeta) #prevalences with mean of betaalpha/(betaalpha+betabeta)
          #rm(betamean, betabeta, betaalpha)
          #spp_intercepts <- family$linkfun(prevalences) 
          #if(family$family == "poisson")
          spp_intercepts <- rnorm(num_spp, 0, 1)      
          #if(family$family == "binomial" & family$link == "logit")
          #      spp_intercepts <- spp_intercepts / 1.6 # To make data generation more comparable to probit
          #set.seed(NULL)
          }

     
     #--------------------
     # Form component due to environmental covariates
     #--------------------
     spp_slopes <- as.matrix(spp_slopes)
     spp_coefficients <- cbind(spp_intercepts, spp_slopes)
     rownames(spp_coefficients) <- paste0("spp", 1:num_spp)
     all_eta <- all_resp <- tcrossprod(sim_dat, spp_coefficients)
     
     
     #--------------------
     # Simulate responses
     #--------------------
     set.seed(seed)
     if(family$family == "binomial") 
          all_resp <- matrix(rbinom(num_units*num_spp, size = 1, prob = family$linkinv(all_eta)), nrow = num_units)
     if(family$family == "poisson") 
          all_resp <- matrix(rpois(num_units*num_spp, lambda = family$linkinv(all_eta)), nrow = num_units)        
          
          
     # Bling up!
     rownames(all_resp) <- paste0("unit", 1:num_units)
     colnames(all_resp) <- rownames(spp_coefficients)
     out <- list(resp = all_resp, eta = all_eta, env = env, spp_coefficients = spp_coefficients)
     
     set.seed(NULL)
     return(out)
     }

     
