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
#' This family object was created specifically for fitting CBFMs. However it follows heavily from the [gamlss.tr::trun()] so big shout outs to the Mikis and the maintainers of that package.
#' 
#' @return An object of class "gamlss.family" and "family"
#' 
#' @author Francis K.C. Hui <fhui28@gmail.com>, Chris Haak
#' 
#' @export 
#' 
#' @importFrom gamlss.dist NBI
#' @importFrom gamlss.tr trun
#' @md

ztnb2 <- function() {
    maketrunfam <- gamlss.tr::trun(par = 0, family = "NBI", type = "left")()
    maketrunfam$family <- "ztnegative.binomial"
    maketrunfam$link <- "log"
    maketrunfam
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
