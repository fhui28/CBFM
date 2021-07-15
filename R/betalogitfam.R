#' @title Family object for beta distribution 
#' 
#' @description 
#' `r lifecycle::badge("stable")`
#' 
#' A specified family object for the beta distribution using the logit link function. 
#'
#' @details 
#' This family object was cretated specifically for fitting CBFMs to continuous data between but excluding zero and one. 
#' 
#' @return A matrix of fitted values, with the same dimensions as the response matrix used in the fitting process.
#' 
#' @author Francis K.C. Hui <fhui28@gmail.com>, Chris Haak
#' 
#' @export
#' @md

betalogitfam <- function() {
     link <- "logit"
     linkfun <- function(mu) 
           return(log(mu/(1-mu)))
     linkinv <- function(eta) 
          return(exp(eta)/(1+exp(eta)))
     mu.eta <- function(eta) {
             ifelse(abs(eta)>30,.Machine$double.eps, exp(eta)/(1+exp(eta))^2) 
             }
     variance <- function(mu, phi) 
          return(mu*(1-mu)/(1+phi))

     structure(list(family = "Beta", link = link, linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta, variance = variance), class = "family")
     }

     
