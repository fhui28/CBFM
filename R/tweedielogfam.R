#' @title Family object for Tweedie distribution
#' 
#' @description 
#' `r lifecycle::badge("stable")`
#' 
#' A specified family object for the Tweedie distribution using the log link function and a power mean-variance relationship
#'
#' @details 
#' This family object was cretated specifically for fitting CBFMs to non-negative continuous data with a "spike" at zero. 
#' 
#' @return An object of class "family" (which has a concise print method). This is a list with the following elements: 
#' \item{family }{The family name i.e., "tweedie".}

#' \item{link }{The link function name i.e., "log".}

#' \item{linkfun }{The link function i.e., the log link.}

#' \item{linkinv }{The inverse link i.e., the exponential function.}

#' \item{mu.eta }{The derivative of the inverse-link function with respect to the linear predictor. For the log link, this is equivalent to the variance function itself i.e., \eqn{d\mu/d\eta = \mu}.}

#' \item{variance }{The variance as a function of the mean. This itself is a function which two arguments, the mean \code{mu}, the dispersion parameter \code{phi}, and the power parameter \code{power}. The mean-variance relationship is given by \eqn{V = \phi\mu^\rho} where \eqn{\mu} denotes the mean,  \eqn{\phi} denotes the dispersion parameter, and \eqn{\rho} is the power parameter.}
#' 
#' @author Francis K.C. Hui <fhui28@gmail.com>, Chris Haak
#' 
#' @export
#' @md

tweedielogfam <- function() {
     link <- "log"
     linkfun <- function(mu) 
          return(log(mu))
     linkinv <- function(eta) 
          return(pmax(exp(eta), .Machine$double.eps))
     mu.eta <- function(eta) 
          return(pmax(exp(eta), .Machine$double.eps))
     variance <- function(mu, power, phi) 
          return(pmax(phi * mu^power, .Machine$double.eps))

     structure(list(family = "tweedie", link = link, linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta, variance = variance), class = "family")
     }

