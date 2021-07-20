#' @title Family object for beta distribution 
#' 
#' @description 
#' `r lifecycle::badge("stable")`
#' 
#' A specified family object for the beta distribution using the logit link function. 
#'
#' @details 
#' This family object was created specifically for fitting CBFMs. 
#' 
#' @return An object of class "family".
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

     
