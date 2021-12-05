#' @title Family object for zero-truncated negative binomial distribution 
#' 
#' @description 
#' `r lifecycle::badge("stable")`
#' 
#' A specified family object for the zero-truncated negative binomial distribution, using the log link function for the negative binomial distribution component.
#'
#' @details 
#' This family object was created specifically for fitting CBFMs. 
#' 
#' This family object was created specifically for fitting CBFMs. However it relies heavily on [gamlss.tr::trun()], so big shout outs to the Mikis and the maintainers of that package.
#' 
#' @return An object of class "gamlss.family" and "family"
#' 
#' @author Francis K.C. Hui <fhui28@gmail.com>, Chris Haak
#' 
#' @export 
# #' @importFrom gamlss.dist NBI
# #' @importFrom gamlss.tr trun
#' @md

ztnb2 <- function() {
    #maketrunfam <- trun(par = 0, family = "NBI", type = "left")()
    maketrunfam <- list(
      family = c("NBItr", "left truncated Negative Binomial type I"), 
      parameters = list(mu = TRUE, sigma = TRUE), 
      nopar = 2, 
      type = "Discrete", 
      mu.link = "log", 
      #sigma.link = "log", 
      mu.linkfun = function (mu) log(mu), 
      #sigma.linkfun = function (mu) log(mu), 
      mu.linkinv = function (eta) pmax(exp(eta), .Machine$double.eps), 
      #sigma.linkinv = function (eta) pmax(exp(eta), .Machine$double.eps), 
      mu.dr = function (eta) pmax(exp(eta), .Machine$double.eps), 
      #sigma.dr = function (eta) pmax(exp(eta), .Machine$double.eps), 
      #dldm = function (y, mu, sigma) 
      #  attr(gamlss::numeric.deriv(dNBItr(y, mu, sigma, log = TRUE), "mu", delta = NULL), "gradient"), 
      #d2ldm2 = function (mu, sigma) { -1/(mu * (1 + mu * sigma)) }, 
      #dldd = function (y, mu, sigma) 
      #  attr(gamlss::numeric.deriv(dNBItr(y, mu, sigma, log = TRUE), "sigma", delta = NULL), "gradient"), 
      #d2ldd2 = function (y, mu, sigma) {
      #  dldd <- -((1/sigma)^2) * (digamma(y + (1/sigma)) - digamma(1/sigma) - log(1 + mu * sigma) - (y - mu) * sigma/(1 + mu * sigma))
      #  d2ldd2 <- -dldd^2
      #  d2ldd2 <- ifelse(d2ldd2 < -1e-15, d2ldd2, -1e-15)
      #  d2ldd2
      #  }, 
      #d2ldmdd = function (y) rep(0, length(y)), 
      #G.dev.incr = function (y, mu, sigma, ...)  -2 * dNBItr(y, mu = mu, sigma = sigma, log = TRUE), 
      #rqres = structure(expression(rqres(pfun = "pNBItr", type = "Discrete", ymin = 1, y = y, mu = mu, sigma = sigma)), 
      #                  srcfile = <environment>, wholeSrcref = structure(c(1L, 0L, 3L, 0L, 0L, 0L, 1L, 3L), 
      #                                                                   srcfile = <environment>, class = "srcref")), 
      mu.initial = expression(mu <- (y + mean(y))/2), 
      #sigma.initial = expression(sigma <- rep(max(((var(y) - mean(y))/(mean(y)^2)), 0.1), length(y))), 
      mu.valid = function (mu) all(mu > 0), 
      #sigma.valid = function (sigma) all(sigma > 0), 
      y.valid = function (y) all(y >= 0),
      mean = function (mu, sigma) mu, 
      variance = function (mu, sigma) mu + sigma * mu^2
      )
    maketrunfam$family <- "ztnegative.binomial"
    maketrunfam$link <- "log"
    maketrunfam$actual_variance <- function(lambda, phi) {
      mu <- lambda / pnbinom(0, mu = lambda, size = 1/phi, lower.tail = FALSE)
      mu * (1 + phi*lambda + lambda - mu)
      }
    maketrunfam
     }



.dztnbinom <- function(x, mu, size, log = FALSE) {
  rval <- dnbinom(x, mu = mu, size = size, log = TRUE) - pnbinom(0, mu = mu, size = size, lower.tail = FALSE, log.p = TRUE)
  rval[x < 1] <- -Inf
  rval[mu <= 0] <- 0
  if(log)
    rval
  else
    exp(rval)
    }


## Negative score of log zero truncated NB distribution wrt to eta, the linear predictor of the NB distribution component, where mu = exp(eta) 
## Adapted heavily from and acknowledgements go to the authors of the the countreg package
.sc_ztnbinom <- function(eta, y, size) {
  mu <- exp(eta)
        
  # Score of NB distribution wrt mu
  snbinom <- function(y, mu, size) {
    s <- y/mu - (y + size)/(mu + size) 
    s[(y < 0) | (abs(y - round(y)) > sqrt(.Machine$double.eps))] <- 0
    s
    }
        
  s <- snbinom(y, mu = mu, size = size)
  logratio <- pnbinom(0, mu = mu, size = size, log.p = TRUE) - pnbinom(0, mu = mu, size = size, lower.tail = FALSE, log.p = TRUE)
  s <- s - exp(logratio + log(size) - log(mu + size))
  s[(y < 1) | (abs(y - round(y)) > sqrt(.Machine$double.eps))] <- 0
  
  return(-s*mu) # Score of ztnb wrt to mu * eta(eta)
  }


## Negative second derivative of log zero truncated NB distribution wrt to eta, the linear predictor of the NB distribution component, where mu = exp(eta) 
## Adapted heavily from and acknowledgements go to the authors of the the countreg package
.hess_ztnbinom <- function(eta, y, size) {
  mu <- exp(eta)
  mu[mu < .Machine$double.eps] <- .Machine$double.eps

  # Score of NB distribution wrt mu
  snbinom <- function(y, mu, size) {
    s <- y/mu - (y + size)/(mu + size) 
    s[(y < 0) | (abs(y - round(y)) > sqrt(.Machine$double.eps))] <- 0
    s
    }
        
  # Hessian of NB wrt mu
  hnbinom <- function(y, mu, size) {
    h <- -y/mu^2 + (y + size)/(mu + size)^2
    h[(y < 0) | (abs(y - round(y)) > sqrt(.Machine$double.eps))] <- 0
    h
    }        
        
  s <- snbinom(y, mu = mu, size = size)
  logratio <- pnbinom(0, mu = mu, size = size, log.p = TRUE) - pnbinom(0, mu = mu, size = size, lower.tail = FALSE, log.p = TRUE)
  s <- s - exp(logratio + log(size) - log(mu + size))
  s[(y < 1) | (abs(y - round(y)) > sqrt(.Machine$double.eps))] <- 0
        
  h <- hnbinom(y = y, mu = mu, size = size)
  h <- h + exp(logratio - 2*log(mu + size) + log(size) + log(1+size)) + exp(2*logratio + 2*log(size) - 2*log(mu+size))
  h[(y < 1) | (abs(y - round(y)) > sqrt(.Machine$double.eps))] <- 0
        
  return(-mu^2*h - mu*s) # exp(eta) * score ztnb wrt mu + exp(2*eta) + hessian ztnb wrt mu
  }
