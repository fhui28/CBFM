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
#' @return An object of class "family". 
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

