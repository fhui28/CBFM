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
#' @return An object of class "family".
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
	
	
