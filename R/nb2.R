#' @title Family object for negative binomial distribution
#' 
#' @description 
#' `r lifecycle::badge("stable")`
#' 
#' A specified family object for the negative binomial distribution using the log link function and a quadratic mean-variance relationship.
#'
#' @details 
#' This family object was cretated specifically for fitting CBFMs to overdispersed count data 
#' 
#' @return An object of class "family" (which has a concise print method). This is a list with the following elements: 
#' \item{family }{The family name i.e., "negative.binomial".}
#' \item{link }{The link function name i.e., "log".}
#' \item{linkfun }{The link function i.e., the log link.}
#' \item{linkinv }{The inverse link i.e., the exponential function.}
#' \item{mu.eta }{The derivative of the inverse-link function with respect to the linear predictor. For the log link, this is equivalent to the variance function itself i.e., \eqn{d\mu/d\eta = \mu}.}
#' \item{variance }{The variance as a function of the mean. This itself is a function which two arguments, the mean \code{mu} and the dispersion parameter \code{phi}. The mean-variance relationship is given by \eqn{V = \mu + \phi\mu^2} where \eqn{\mu} denotes the mean and \eqn{\phi} is the dispersion parameter.}
#' 
#' @author Francis K.C. Hui <fhui28@gmail.com>, Chris Haak
#' 
#' @export
#' @md

nb2 <- function() {
     ##--------------------------
     ## Feed the family!
     ##--------------------------
	link <- "log"
	linkfun <- function(mu) 
          return(log(mu))
	linkinv <- function(eta) 
          return(pmax(exp(eta), .Machine$double.eps))
	mu.eta <- function(eta) 
          return(pmax(exp(eta), .Machine$double.eps))
	variance <- function(mu, phi) 
          return(pmax(mu+phi*mu^2, .Machine$double.eps))
  
	structure(list(family = "negative.binomial", link = link, linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta, variance = variance), class = "family")
	}
	
	
