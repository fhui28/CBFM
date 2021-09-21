#' @title Family object for zero-inflated negative binomial distribution 
#' 
#' @description 
#' `r lifecycle::badge("stable")`
#' 
#' A specified family object for the zero-inflated negative binomial distribution, using the log link function for the negative binomial distribution component.
#'
#' @details 
#' This family object was created specifically for fitting CBFMs. 
#' 
#' @return An object of class "family".
#' 
#' @author Francis K.C. Hui <fhui28@gmail.com>, Chris Haak
#' 
#' @export zinb2
#' @md

zinb2 <- function() {
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
        actual_variance <- function(mu, zeroinfl_prob, phi) 
                  return((1-zeroinfl_prob) * mu * (1 + mu * (zeroinfl_prob + phi)))

        structure(list(family = "zinegative.binomial", link = link, linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta, variance = variance,
                       actual_variance = actual_variance), class = "family")
        }

     

.dzinegativebinomial_log <- function(eta, y, zeroinfl_prob, phi) {
        logp <- log(1-zeroinfl_prob) + dnbinom(y, mu = exp(eta), size = 1/phi, log=TRUE)
        logp[y == 0] <- log(exp(logp[y == 0]) + zeroinfl_prob)
        
        return(sum(logp))
        }

## Adapted from the VGAM package, and acknowledgement goes to its authors
## Assumes q, lambda, zeroinfl_prob, phi are all of the same dimension 
.pzinegativebinomial <- function(q, lambda, zeroinfl_prob, phi) {
        ans <- zeroinfl_prob + (1 - zeroinfl_prob) * pnbinom(q, mu = lambda, size = 1/phi)
        ans[q < 0] <- 0

        prob <- (1/phi) / (1/phi + lambda)
        prob0 <- prob^(1/phi)
        deflat.limit <- -prob0 / (1 - prob0)
        ans[zeroinfl_prob < deflat.limit] <- NaN
        ans[zeroinfl_prob > 1] <- NaN        
        return(ans)
        }