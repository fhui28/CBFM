#' @title Basic plot residual diagnostics from a CBFM fit
#' 
#' @description 
#' `r lifecycle::badge("experimental")`
#' 
#' Five potential plots are currently available for some basic residual diagnostics for a fitted \code{CBFM} project: 1) a plot of residuals against
#' estimated linear predictors; 2) a normal probability or quantile-quantile plot of residuals with simulated point-wise 95\% confidence interval envelope; 3) plot of residuals against observational unit index; 4) a plot of residuals again column index; 5) scale-location plot.
#' 
#' @param x An object of class "CBFM".
#' @param type The type of residuals to be used in constructing the plots. Currently the options available are: "response" (default), "pearson", "PIT", "dunnsmyth", and "partial". Can be abbreviated.
#' @param which_plot If a subset of the plots is desired, then a vector containing subset of the integers 1,2,3,4,5.
#' @param titles Titles to appear above each plot.
#' @param species_colors Either a scalar is supplied for a a vector with length of number of species in the spatio-temporal multivariate abundance data. If the former than all species use this color in the plots. If the latter then the vector specified the colors to use for each species. Defaults to \code{NULL}, which results in each species having a unique color based on the [grDevices::rainbow()] palette.
#' @param smooth Should a smoother be added to each plot?
#' @param envelope_col A vector of length 2, specifying the colors to use for the lines and shade respectively for the simulation envelopes.
#' @param which_species A vector indexing the species to plot, if the residual plots should be constructed for only a subset of species. Defaults to \code{NULL}, in which case all species are plotted. This may be useful if the number of species is quite large.
#' @param seed This can be used set the seed when constructing the PIT and Dunn-Smyth residuals, which for discrete responses involve some degree of jittering.  
#' @param ... Additional graphical arguments.

#' 
#' @details 
#' This function is heavily adapted from [gllvm::plot.gllvm()] and [boral::plot.boral()]. A lot of credit goes to the authors of the \code{gllvm} package, especially Jenni Niku, for the code!
#' 
#' As basic residual diagnostics, these plots should behave a follows: 1) the plot of residuals versus linear predictors should not exhibit any noticeable pattern such as (inverse) fan-shape or a trend; 2) the normal probability plot should have the residuals lying approximately on a straight line and almost all residuals lying within the simulation envelopes; 3) a plot of residuals against observational unit index should not exhibit any noticeable pattern such as (inverse) fan-shape or a trend. The plot can also be used to look for potential outlying observational units; 4) a plot of residuals against species index should not exhibit any noticeable pattern such as (inverse) fan-shape or a trend. The plot can also be used to look for potential outlying species; 5) the scale-location plot should not exhibit any trend. 
#' 
#' If the above does not occur then it may suggest one or more modeling assumptions such as the assumed response distribution, or the model used for the measured covariates, may ot sufficiently satisfied. 
#' 
#' # Warning:
#' This function only provides basic diagnostic plots. For spatio-temporal data, there are also more specialized methods and plots for assessing the assumptions made in relation to the spatial and/or temporal correlation, and the practitioner is encouraged to examine these. Basic plots may include plotting the residuals as a function of the spatio-temporal coordinates, and constructing autocorrelation plots and covariograms of the residuals; see for example Li (2003), Hyndman and Athanasopoulos (2018), Plant (2018), and Wilke et al. (2019) for some *general* references which contain sections on residual analysis for spatial and time-series data. There are also more (currently) bespoke methods e.g., Bose et al., (2018), although they may be hard to implement and may also not be too useful for the CBFM approach to analyzing spatio-temporal data in general.  
#' 
#' @author Francis K.C. Hui <fhui28@gmail.com>, Chris Haak
#' 
#' @references 
#' Bose, M., Hodges, J. S., and Banerjee, S. (2018). Toward a diagnostic toolkit for linear models with Gaussian‐process distributed random effects. Biometrics, 74, 863-873.
#' 
#' Li, W. K. (2003). Diagnostic checks in time series. CRC Press.
#' 
#' Hyndman, R. J., and Athanasopoulos, G. (2018). Forecasting: principles and practice. OTexts.
#' 
#' Plant, R. E. (2018). Spatial data analysis in ecology and agriculture using R. CRC Press.
#' 
#' Wikle, C. K., Zammit-Mangion, A., and Cressie, N. (2019). Spatio-temporal Statistics with R. Chapman and Hall/CRC 
#' 
#' @seealso [CBFM()] for fitting CBFMs, [fitted.values.CBFM()] for calculating fitted values from a CBFM fit, and [residuals.CBFM()] for calculating various types of residuals.
#' 
#' @examples
#' \donttest{
#' library(autoFRK)
#' library(FRK)
#' library(MASS)
#' library(mvabund)
#' library(mvtnorm)
#' library(ROCR)
#' library(sp)
#' library(RandomFields)
#' library(tidyverse)
#' 
#' ##------------------------------
#' ## Example 1: Fitting a CBFM to spatial multivariate presence-absence data 
#' ## simulated from a spatial latent variable model
#' ## Please note the data generation process (thus) differs from CBFM.
#' ##------------------------------
#' set.seed(2021)
#' num_sites <- 500 # 500 (units) sites 
#' num_spp <- 50 # Number of species
#' num_X <- 4 # Number of regression slopes
#' 
#' spp_slopes <- matrix(runif(num_spp * num_X, -1, 1), nrow = num_spp)
#' spp_intercepts <- runif(num_spp, -2, 0)
#' 
#' # Simulate spatial coordinates and environmental covariate components
#' xy <- data.frame(x = runif(num_sites, 0, 5), y = runif(num_sites, 0, 5))
#' X <- rmvnorm(num_sites, mean = rep(0,4)) 
#' colnames(X) <- c("temp", "depth", "chla", "O2")
#' dat <- data.frame(xy, X)
#' mm <- model.matrix(~ temp + depth + chla + O2 - 1, data = dat) %>% 
#' scale %>% 
#' as.matrix
#' 
#' # Simulate latent variable component
#' true_lvs <- RFsimulate(model = RMexp(var=1, scale=2), 
#' x = xy$x, y = xy$y, n = 2)@data %>% 
#' as.matrix
#' spp_loadings <- matrix(runif(num_spp * 2, -1, 1), nrow = num_spp) 
#' set.seed(NULL)
#' 
#' # Simulate spatial multivariate abundance data (presence-absence)
#' eta <- tcrossprod(cbind(1,mm), cbind(spp_intercepts,spp_slopes)) + 
#' tcrossprod(true_lvs, spp_loadings)
#' simy <- matrix(rbinom(num_sites * num_spp, size = 1, 
#' prob = binomial()$linkinv(eta)), nrow = num_sites)
#' rm(X, mm, spp_loadings, true_lvs, xy, eta)
#' 
#' 
#' # Set up spatial basis functions for CBFM -- Most practitioners will start here! 
#' num_basisfunctions <- 25 # Number of spatial basis functions to use
#' basisfunctions <- mrts(dat[,c("x","y")], num_basisfunctions) %>% 
#' as.matrix %>%
#' {.[,-(1)]} # Remove the first intercept column
#' 
#' # Fit CBFM 
#' useformula <- ~ temp + depth + chla + O2
#' fitcbfm <- CBFM(y = simy, formula_X = useformula, data = dat, 
#' B_space = basisfunctions, family = binomial(), control = list(trace = 1))
#' 
#' plot(fitcbfm, ask = TRUE)
#' 
#' plot(fitcbfm, ask = TRUE, pch = ".") # To speed up plotting...
#' }
#' 
#' @export
#' @importFrom graphics abline boxplot lines panel.smooth par plot points polygon
#' @importFrom grDevices rainbow
#' @importFrom mgcv gam predict.gam
#' @importFrom stats qnorm qqnorm qqline quantile rnorm
#' @importFrom tdigest tdigest tquantile

plot.CBFM <- function(x, which_plot = 1:5, type = "dunnsmyth", titles = c("Residuals vs. linear predictors", "Normal probability plot", "Residuals vs. unit index", "Residuals vs. species index","Scale-Location plot"), species_colors = NULL, smooth = TRUE, envelope_col = c("blue","lightblue"), which_species = NULL, seed = NULL, ...) {
     num_units <- nrow(x$y)
     num_spp <- ncol(x$y)
        
     sppind <- 1:num_spp
     if(!is.null(which_species))
          sppind <- sort(which_species)
          
     if(length(sppind) > num_spp)
          stop("Length of which_species exceeded the number of species on x$y.")
     if(any(which_species > num_spp))
          stop("which_species should be a vector of integers ranging from 1 to ncol(x$y) i.e., the number of species.")
     
     if(any(which_plot > 5))
          stop("which_plot should be a vector of integers ranging from 1 to 5. There are only five possible plots that this function currently offers.")
     
     
     # Form plot titles
     mains <- rep("", 5)
     mains[which_plot] <- titles[which_plot]

     res <- residuals(object = x, type = type, seed = seed)
     dsres <- res[, sppind]
     etamat <- x$linear_predictor[,sppind]
     xxx <- boxplot(c(etamat), outline = FALSE, plot = FALSE)$stats     
     

     # Form colors for species - done by prevalence
     csum <- order(colSums(as.matrix(x$y))[sppind])
     if(!is.null(species_colors)) {
          col <- rep(1, num_spp)
          col[1:num_spp] <- species_colors
          } 
     if(is.null(species_colors)) {
          if(num_spp < 8)
               col <- (1:num_spp)[csum]
          else
               col <- rainbow(num_spp + 1)[2:(num_spp + 1)][csum]
          }
          

     gr.pars <- list(...)
     par(...)

     
     # Residuals versus linear predictors
     if(1 %in% which_plot) {
          if(is.null(gr.pars$xlim)) {
               plot(etamat, dsres, xlab = "Linear predictors", ylab = "Residuals", type = "n", col = rep(col, each = num_units), 
                    main = mains[1], xlim = c(min(xxx), max(xxx)))
               abline(0, 0, col = "grey", lty = 3)
               } 
          else {
               plot(etamat, dsres, xlab = "Linear predictors", ylab = "Residuals", type = "n", col = rep(col, each = num_units), 
                    main = mains[1], ...)
               abline(0, 0, col = "grey", lty = 3)
               }
 
          if(smooth) 
               .gamEnvelope(etamat, dsres, col = rep(col, each = num_units), envelopes = TRUE, envelope.col = envelope_col, ...)
          
          }
          

     # Normal probability or quantile-quantile plot of residuals with simulated point-wise 95\% confidence interval envelope          
     if(2 %in% which_plot) {
          qq.x <- qqnorm(c(dsres), main = mains[2], ylab = "Dunn-Smyth residuals", col = rep(col, each = num_units), cex = 0.5, xlab = "theoretical quantiles")
          #qqline(c(dsres), col = envelope_col[1])
          
         K <- 200
         yy <- quantile(dsres, c(0.25, 0.75), names = FALSE, type = 7, na.rm = TRUE)
         xx <- qnorm(c(0.25, 0.75))
         slope <- diff(yy) / diff(xx)
         int <- yy[1] - slope * xx[1]
         all_ris <- matrix(rnorm(num_units * num_spp * K, mean = int, sd = slope), ncol = K)
         Ym <- apply(all_ris, 2, sort)
         Xm <- sort(qq.x$x)
         #cis <- apply(Ym, 1, quantile, probs = c(0.025, 0.975))
         cis <- apply(Ym, 1, function(x) { tquantile(tdigest(x, 1000), probs = c(0.025, 0.975)) })
 
         num_obs <- num_units * num_spp
         polygon(Xm[c(1:num_obs,num_obs:1)], c(cis[1,],cis[2, num_obs:1]), col = envelope_col[2], border = NA)
         points(qq.x, col = rep(col, each = num_units), cex = 0.5)
         qqline(c(dsres), col = envelope_col[1])
         }
          
          
     # Residuals against observational unit index          
     if(3 %in% which_plot) {
          plot(rep(1:num_units, num_spp), dsres, xlab = "Unit index", ylab = "Residuals", col = rep(col, each = num_units), main = mains[3], ...);
          abline(0, 0, col = "grey", lty = 3)
          if(smooth) 
               panel.smooth(rep(1:num_units, num_spp), dsres, col = rep(col, each = num_units), col.smooth = envelope_col[1], ...)
          }

          
     # Residuals against species index          
     if(4 %in% which_plot) {
          plot(rep(1:num_spp, each = num_units), dsres, xlab = "Species index", ylab = " Residuals", col = rep(col[csum], each = num_units), main = mains[4], ...) 
          abline(0, 0, col = "grey", lty = 3)
          if(smooth) 
               panel.smooth(rep(1:num_spp, each = num_units), dsres, col = rep(col[csum], each = num_units), col.smooth = envelope_col[1], ...)
          }

          
     # Scale-location plot
     if(5 %in% which_plot) {
          sqres <- sqrt(abs(dsres))
          yl <- as.expression(substitute(sqrt(abs(YL)), list(YL = as.name("Residuals"))))
          
          if(is.null(gr.pars$xlim)) {
               plot(etamat, sqres, xlab = "Linear predictors", ylab = yl, col = rep(col, each = num_units), main = mains[5], xlim = c(min(xxx), max(xxx)), ...)
               } 
          else {
               plot(etamat, sqres, xlab = "Linear predictors", ylab = yl, col = rep(col, each = num_units), main = mains[5], ...)
               }

          if(smooth) 
               panel.smooth(etamat, sqres, col = rep(col, each = num_units), col.smooth = envelope_col[1], ...)
          }
     
     }

     
## Pulled straight from gllvm package. Thanks to Jenni for this function!
.gamEnvelope <- function(x, y, line.col = "red", envelope.col = c("blue","lightblue"), col = 1, envelopes = TRUE, ...) {
     xSort <- sort(x, index.return = TRUE)
     gam.yx <- gam(y[xSort$ix] ~ xSort$x)
     pr.y <- predict.gam(gam.yx, se.fit = TRUE)
     
     n.obs <- length(xSort$ix)
     prHi <- pr.y$fit + 1.96*pr.y$se.fit
     prLow <- pr.y$fit - 1.96*pr.y$se.fit
     
     if(envelopes) 
          polygon(xSort$x[c(1:n.obs,n.obs:1)], c(prHi,prLow[n.obs:1]), col = envelope.col[2], border = NA)
     
     lines(xSort$x, pr.y$fit, col = envelope.col[1])
     abline(h = 0, col = 1)
     points(x, y, col = col, ...)
     }
     
     
