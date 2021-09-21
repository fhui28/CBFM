# ## Truncated negative binomial with log-link
# ## Adapted from the ztpoisson family function on countreg and negbin family function in mgcv 
## SCRAPPED AS IT EITHER DOES NOT WORK OR IS TOO SLOW!
# .ztnb2_lambda_to_mean <- function(lambda, phi) {
#     ifelse(lambda <= 0, 1, lambda / (1 - (1 + phi*lambda)^(-1/phi)))
#     }
#     
# .ztnb2_mean_to_lambda <- function(mean, phi, interval = c(-1e-5,5e5)) {
#     f <- function(x, phi, mu) 
#         mu - x/(1 - (1 + phi * x)^(-1/phi)) 
#     
#     out <- NULL
#     for(i in 1:length(mean)) {
#         rootfind <- stats::uniroot(f, interval = interval, phi = phi, mu = mean[i])
#         out <- c(out, ifelse(rootfind$root > 0, rootfind$root, 1e-6))
#         }
#     return(out)
#     }     
# 
# 
# 
# ztnb2 <- function(theta = stop("'theta' must be specified")) {
#     # log(lambda) = eta => lambda = exp(eta)
#     # mean = lambda / P(Y=0) = lambda / 1 - (1+phi*lambda)^(-1/phi)
# 
#     ##--------------------------
#     ## Some sub functions first
#     ##--------------------------
#     # See <https://rdrr.io/cran/mgcv/man/fix.family.link.html>. Basically, d2link is the second derivative of the link function w.r.t. the mean of the family supplied; similar for d3link and d4ink (4th derivatives needed if Fisher scoring is not used).
#     # Using a combination of inverse function theorem <https://en.wikipedia.org/wiki/Inverse_functions_and_differentiation> and Wolfram alpha
#     # + derivative of log(lambda(mu))/mu = 1/lambda * dlambda/dmu
#     # + second derivative is = -1/lambda^2 (dlambda/dmu)^2 + 1/lambda * d^2lambda/dmu^2
#     # + third derivative is = 2/lambda^3 (dlambda/dmu)^3 - 1/lambda^2 * 2 * dlambda/dmu * d^2lambda/dmu^2 - 1/lambda^2 * dlambda/dmu * d^2lambda/dmu^2 + 1/lambda * d^3lambda/dmu^3
#     # Second and third derivatives below confirmed numercially
#     ztnb2_dlink <- function(mu) {
#         phi <- 1/get(".Theta")[1]
#         lambda <- .ztnb2_mean_to_lambda(mean = mu, phi = phi)
#         
#         dmu_dlambda <- -( (1+phi*lambda)^(-1/phi) + lambda*(1+phi*lambda)^(-(phi+1)/phi) - 1 )/((1+phi*lambda)^(-1/phi) - 1)^2
#         return(1/lambda * 1/dmu_dlambda)
#         }
#     
#     ztnb2_d2link <- function(mu) {
#         phi <- 1/get(".Theta")[1]
#         lambda <- .ztnb2_mean_to_lambda(mean = mu, phi = phi)
#     
#         dmu_dlambda <- -( (1+phi*lambda)^(-1/phi) + lambda*(1+phi*lambda)^(-(phi+1)/phi) - 1 )/((1+phi*lambda)^(-1/phi) - 1)^2
#         d2mu_dlambda2 <- lambda * ( (2*(phi*lambda+1)^(-2/phi-2))/(1-(phi*lambda+1)^(-1/phi))^3 - ((-1/phi-1)*phi*(phi*lambda+1)^(-1/phi-2))/(1-(phi*lambda+1)^(-1/phi))^2) - (2*(phi*lambda+1)^(-1/phi-1))/(1-(phi*lambda+1)^(-1/phi))^2
#     
#         d2lambda_dmu2 <- -d2mu_dlambda2 * (1/dmu_dlambda)^3
#         
#         return(-1/lambda^2 * (1/dmu_dlambda)^2 + 1/lambda * d2lambda_dmu2)
#         }
#      
#     ztnb2_d3link <- function(mu) {
#         phi <- 1/get(".Theta")[1]
#         lambda <- .ztnb2_mean_to_lambda(mean = mu, phi = phi)
#     
#         dmu_dlambda <- -( (1+phi*lambda)^(-1/phi) + lambda*(1+phi*lambda)^(-(phi+1)/phi) - 1 )/((1+phi*lambda)^(-1/phi) - 1)^2
#         d2mu_dlambda2 <- lambda * ( (2*(phi*lambda+1)^(-2/phi-2))/(1-(phi*lambda+1)^(-1/phi))^3 - ((-1/phi-1)*phi*(phi*lambda+1)^(-1/phi-2))/(1-(phi*lambda+1)^(-1/phi))^2) - (2*(phi*lambda+1)^(-1/phi-1))/(1-(phi*lambda+1)^(-1/phi))^2
#         d3mu_dlambda3 <- ((phi*lambda+1)^(1/phi-3) * (-2*((phi^2+2)*lambda+3*phi)*(phi*lambda+1)^(1/phi) + (phi+1)*((phi-1)*lambda+3)*(phi*lambda+1)^(2/phi) + (phi-1)*(phi*lambda+lambda+3))) / ((phi*lambda+1)^(1/phi)-1)^4
#         
#         d2lambda_dmu2 <- -d2mu_dlambda2 * (1/dmu_dlambda)^3
#         d3lambda_dmu3 <- -d3mu_dlambda3 * (1/dmu_dlambda)^4 + 3 * d2mu_dlambda2^2 * (1/dmu_dlambda)^5
#         
#         return(2/lambda^3 * (1/dmu_dlambda)^3 - 3/lambda^2 * 1/dmu_dlambda * d2lambda_dmu2 + 1/lambda * d3lambda_dmu3)
#         }
#     
#     ztnb2_d4link <- function(mu) {
#         phi <- 1/get(".Theta")[1]
#         
#         innerztnb2_d3link <- function(mu) {
#             lambda <- .ztnb2_mean_to_lambda(mean = mu, phi = phi)
#         
#             dmu_dlambda <- -( (1+phi*lambda)^(-1/phi) + lambda*(1+phi*lambda)^(-(phi+1)/phi) - 1 )/((1+phi*lambda)^(-1/phi) - 1)^2
#             d2mu_dlambda2 <- lambda * ( (2*(phi*lambda+1)^(-2/phi-2))/(1-(phi*lambda+1)^(-1/phi))^3 - ((-1/phi-1)*phi*(phi*lambda+1)^(-1/phi-2))/(1-(phi*lambda+1)^(-1/phi))^2) - (2*(phi*lambda+1)^(-1/phi-1))/(1-(phi*lambda+1)^(-1/phi))^2
#             d3mu_dlambda3 <- ((phi*lambda+1)^(1/phi-3) * (-2*((phi^2+2)*lambda+3*phi)*(phi*lambda+1)^(1/phi) + (phi+1)*((phi-1)*lambda+3)*(phi*lambda+1)^(2/phi) + (phi-1)*(phi*lambda+lambda+3))) / ((phi*lambda+1)^(1/phi)-1)^4
#             
#             d2lambda_dmu2 <- -d2mu_dlambda2 * (1/dmu_dlambda)^3
#             d3lambda_dmu3 <- -d3mu_dlambda3 * (1/dmu_dlambda)^4 + 3 * d2mu_dlambda2^2 * (1/dmu_dlambda)^5
#             
#             return(2/lambda^3 * (1/dmu_dlambda)^3 - 3/lambda^2 * 1/dmu_dlambda * d2lambda_dmu2 + 1/lambda * d3lambda_dmu3)
#             }
# 
#         out <- numDeriv::grad(func = innerztnb2_d3link, x = mu)
#         return(out)
#         }
#     
#     ## dvar is the derivative of the variance function w.r.t. the mean of the family supplied
#     ## d2var the 2nd derivative of the variance function w.r.t. the mean of the family supplied
#     ## variance = mu * (1+phi*lambda - p0 mu^2) = mu * (1 - mu + lambda*(phi+1)) once you substitute in p0 = 1-lambda/mu
#     dvar <- function(mu) {
#         phi <- 1/get(".Theta")[1]
#         lambda <- .ztnb2_mean_to_lambda(mean = mu, phi = phi)
#         
#         dmu_dlambda <- -( (1+phi*lambda)^(-1/phi) + lambda*(1+phi*lambda)^(-(phi+1)/phi) - 1 )/((1+phi*lambda)^(-1/phi) - 1)^2
#         return(1 - 2*mu + lambda*(phi+1) + mu*(phi+1)*(1/dmu_dlambda))
#         }
#     
#     d2var <- function(mu) {
#         phi <- 1/get(".Theta")[1]
#         lambda <- .ztnb2_mean_to_lambda(mean = mu, phi = phi)
#         
#         dmu_dlambda <- -( (1+phi*lambda)^(-1/phi) + lambda*(1+phi*lambda)^(-(phi+1)/phi) - 1 )/((1+phi*lambda)^(-1/phi) - 1)^2
#         d2mu_dlambda2 <- lambda * ( (2*(phi*lambda+1)^(-2/phi-2))/(1-(phi*lambda+1)^(-1/phi))^3 - ((-1/phi-1)*phi*(phi*lambda+1)^(-1/phi-2))/(1-(phi*lambda+1)^(-1/phi))^2) - (2*(phi*lambda+1)^(-1/phi-1))/(1-(phi*lambda+1)^(-1/phi))^2
#         
#         d2lambda_dmu2 <- -d2mu_dlambda2 * (1/dmu_dlambda)^3
#         
#         return(-2 + 2*(phi+1)*(1/dmu_dlambda) + mu*(phi+1)*d2lambda_dmu2)
#         }
#     
#     d3var <- function(mu) {
#         phi <- 1/get(".Theta")[1]
#         lambda <- .ztnb2_mean_to_lambda(mean = mu, phi = phi)
#         
#         dmu_dlambda <- -( (1+phi*lambda)^(-1/phi) + lambda*(1+phi*lambda)^(-(phi+1)/phi) - 1 )/((1+phi*lambda)^(-1/phi) - 1)^2
#         d2mu_dlambda2 <- lambda * ( (2*(phi*lambda+1)^(-2/phi-2))/(1-(phi*lambda+1)^(-1/phi))^3 - ((-1/phi-1)*phi*(phi*lambda+1)^(-1/phi-2))/(1-(phi*lambda+1)^(-1/phi))^2) - (2*(phi*lambda+1)^(-1/phi-1))/(1-(phi*lambda+1)^(-1/phi))^2
#         d3mu_dlambda3 <- ((phi*lambda+1)^(1/phi-3) * (-2*((phi^2+2)*lambda+3*phi)*(phi*lambda+1)^(1/phi) + (phi+1)*((phi-1)*lambda+3)*(phi*lambda+1)^(2/phi) + (phi-1)*(phi*lambda+lambda+3))) / ((phi*lambda+1)^(1/phi)-1)^4
#             
#         d2lambda_dmu2 <- -d2mu_dlambda2 * (1/dmu_dlambda)^3
#         d3lambda_dmu3 <- -d3mu_dlambda3 * (1/dmu_dlambda)^4 + 3 * d2mu_dlambda2^2 * (1/dmu_dlambda)^5
#         
#         return(3*(phi+1)*d2lambda_dmu2 + mu*(phi+1)*d3lambda_dmu3)
#         }
# 
#     # See <https://data.princeton.edu/wws509/notes/countmoments>
#     variance <- function(mu) {
#         phi <- 1/get(".Theta")[1]
#         lambda <- .ztnb2_mean_to_lambda(mean = mu, phi = phi)
#         
#         p0 <- (1 + phi*lambda)^(-1/phi)
#         return(mu * (1 + phi * lambda - mu * p0))
#         }
#     
#     validmu <- function(mu) 
#         all(mu > 1)
# 
#     dev.resids <- function(y, mu, wt) {
#         phi <- 1/get(".Theta")[1]
#         -2 * wt * (dztnb2(y, mean = mu, phi = phi, log = TRUE) - dztnb2(y, mean = y, phi = phi, log = TRUE))
#         }
#      
#     aic <- function(y, n, mu, wt, dev) {
#         phi <- 1/get(".Theta")[1]
#         -2 * sum(dztnb2(y, mean = mu, phi = phi, log = TRUE) * wt)
#         }
#     
#     getTheta <- function() 
#         get(".Theta")    
#           
#     ls <- function(y, w, n, scale) { 
#         phi <- 1/get(".Theta")[1]
#         c(sum(dztnb2(y, lambda = y + 1e-8, phi = phi, log = TRUE) * w), 0, 0) 
#         }
#     
#     ##--------------------------
#     ## Feed the family!
#     ##--------------------------
#     link <- "log"
#     env <- new.env(parent = .GlobalEnv)
#     assign(".Theta", theta, envir = env)
#     
#     environment(ztnb2_dlink) <- environment(ztnb2_d2link) <- environment(ztnb2_d3link) <- environment(ztnb2_d4link) <- environment(dvar) <- environment(d2var) <- environment(d3var) <- environment(variance) <- environment(dev.resids) <- environment(aic) <- environment(getTheta) <- env
#     
#     stats <- structure(list(
#         linkfun = function(mu) log(.ztnb2_mean_to_lambda(mean = mu, phi = 1/theta)),
#         linkinv = function(eta) .ztnb2_lambda_to_mean(lambda = exp(eta), phi = 1/theta),
#         mu.eta = function(eta) {
#             lambda <- exp(eta)
#             phi <- 1/theta
#             mu <- .ztnb2_lambda_to_mean(lambda = lambda, phi = phi)
#             mu * (1 - mu * (1+phi*lambda)^(1/phi-1))
#             },
#         valideta = function(eta) TRUE,
#         name = "log"), 
#         class = "link-glm"
#         )
# 
#     initialize <- expression({
#         if(any(y < 1)) 
#             stop("zero or negative values not allowed for the zero-truncated negative binomial family.")
#         n <- rep.int(1, nobs)
#         mustart <- y + 0.1
#         })     
#      
#      
#     structure(list(
#         family = "ztnegative.binomial", 
#         link = link,
#         linkfun = stats$linkfun,
#         linkinv = stats$linkinv,
#         variance = variance,
#         dvar = dvar,
#         d2var = d2var,
#         d3var = d3var,
#         dev.resids = dev.resids,
#         aic = aic,
#         mu.eta = stats$mu.eta,
#         initialize = initialize,
#         ls = ls, 
#         validmu = validmu,
#         valideta = stats$valideta,
#         getTheta = getTheta,
#         d2link = ztnb2_d2link,
#         d3link = ztnb2_d3link,
#         d4link = ztnb2_d4link
#         ),
#         class = "family")
#      }
# 
# 
# dztnb2 <- function(x, lambda, mean, phi, log = FALSE) {
#     if(!missing(mean)) 
#         lambda <- .ztnb2_mean_to_lambda(mean, phi = phi)
#         
#     rval <- dnbinom(x, mu = lambda, size = 1/phi, log = TRUE) - pnbinom(0, mu = lambda, size = 1/phi, lower.tail = FALSE, log.p = TRUE)
#     rval[x < 1] <- -Inf
#     rval[lambda <= 0] <- 0
#     if(log) 
#         rval 
#     else exp(rval)
#     }



#-----------------------------
#-----------------------------
library(mvtnorm)
library(countreg)
library(mgcv)

data("CrabSatellites", package = "countreg")
cs <- CrabSatellites[, c("satellites", "width", "color")]
cs$color <- as.numeric(cs$color)
cs <- subset(cs, subset = satellites > 0)

# zt_p <- zerotrunc(satellites ~ ., data = cs)
# summary(zt_p)
# fit_ztp <- gam(satellites ~ s(width) + color, data = cs, family = ztpoisson())
# summary(fit_ztp)

zt_nb <- zerotrunc(satellites ~ ., data = cs, dist = "negbin")
summary(zt_nb)

fit_ztnb <- gam(satellites ~ width + color, data = cs, family = ztnb2(theta = 4.6))
                     



##-------------------------------
# ## Beta distribution with probit-link
# betaprobitfam <- function() {v
#      link <- "probit"
#      linkfun <- function(mu) 
#           return(qnorm(mu))
#      linkinv <- function(eta) 
#           return(pnorm(eta))
#      mu.eta <- function(eta) 
#           return(pmax(dnorm(eta), .Machine$double.eps))
#      variance <- function(mu, phi) 
#           return(mu*(1-mu)/(1+phi))
# 
#      structure(list(family = "Beta", linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta, variance = variance, link = link), class = "family")
#      }

