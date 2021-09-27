#' @title Family object for zero-truncated Poisson distribution 
#' 
#' @description 
#' `r lifecycle::badge("stable")`
#' 
#' A specified family object for the zero-truncated Poisson distribution, using the log link function for the mean of the Poisson distribution component.
#'
#' @details 
#' This family object was created specifically for fitting CBFMs. However it follows heavily from the [gamlss.tr::trun()] so big shout outs to the Mikis and the maintainers of that package.
#' 
#' @return An object of class "gamlss.family" and "family"
#' 
#' @author Francis K.C. Hui <fhui28@gmail.com>, Chris Haak
#' 
#' @export 
#' 
#' @importFrom gamlss.dist PO
#' @importFrom gamlss.tr trun
#' @md

ztpoisson <- function() {
    maketrunfam <- trun(par = 0, family = "PO", type = "left")()
    maketrunfam$family <- "ztpoisson"
    maketrunfam$link <- "log"
    maketrunfam$actual_variance <- function(mu, lambda) {
        mu * (1 + lambda - mu) 
        }
    maketrunfam
    }


.dztpois <- function(x, lambda, log = FALSE) {
    rval <- dpois(x, lambda, log = TRUE) - ppois(0, lambda, lower.tail = FALSE, log.p = TRUE)
    rval[x < 1] <- -Inf
    rval[lambda <= 0] <- 0
    if(log)
        rval
        else
                exp(rval)
        }


.hess_ztpoisson <- function(eta, y) {
    lambda <- exp(eta)
    lambda[lambda < .Machine$double.eps] <- .Machine$double.eps
    
    # First derivative of ztpoisson wrt to lambda
    s <- y/lambda - 1 - exp(-lambda)/(1 - exp(-lambda))
    s[(y < 1) | (abs(y - round(y)) > sqrt(.Machine$double.eps))] <- 0
    
    # Second derivative of ztpoisson wrt to lambda
    h <- - y/lambda^2 + exp(-lambda)/(1 - exp(-lambda))^2
    h[(y < 1) | (abs(y - round(y)) > sqrt(.Machine$double.eps))] <- 0
    
    return(-lambda^2*h - lambda*s) # exp(eta) * score ztnb wrt mu + exp(2*eta) + hessian ztnb wrt mu
    }
