#' @title Family object for zero-inflated Poisson distribution 
#' 
#' @description 
#' `r lifecycle::badge("stable")`
#' 
#' A specified family object for the zero-inflated Poisson distribution, using the log link function for the Poisson distribution component.
#'
#' @details 
#' This family object was created specifically for fitting CBFMs. 
#' 
#' @return An object of class "family".
#' 
#' @author Francis K.C. Hui <fhui28@gmail.com>, Chris Haak
#' 
#' @export zipoisson
#' @md

zipoisson <- function() {
        link <- "log"
        linkfun <- function(mu) 
                return(log(mu))
        linkinv <- function(eta) 
                return(exp(eta))
        mu.eta <- function(eta) {
                pmax(exp(eta), .Machine$double.eps)
                }
        variance <- function(mu) 
                  return(mu)
        actual_variance <- function(mu, zeroinfl_prob) 
                  return(mu + zeroinfl_prob/(1-zeroinfl_prob) * mu^2)

        structure(list(family = "zipoisson", link = link, linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta, variance = variance,
                       actual_variance = actual_variance), class = "family")
        }

     
# .dzipoisson <- function(y, lambda, zeroinfl_prob) {
#         p <- (1-zeroinfl_prob) * dpois(y, lambda)
#         p[y == 0] <- exp(logp[y == 0]) + (zeroinfl_prob[y == 0])
#         
#         return(sum(logp))
#         }

## Adapted from the countreg/pscl package, and acknowledgement goes to them
.dzipoisson_log <- function(eta, y, zeroinfl_prob) {
        logp <- log(1-zeroinfl_prob) + dpois(y, lambda = exp(eta), log=TRUE)
        logp[y == 0] <- log(exp(logp[y == 0]) + zeroinfl_prob)
        
        return(sum(logp))
        }

## Adapted from the VGAM package, and acknowledgement goes to them
.pzipois <- function(q, lambda, zeroinfl_prob = 0) {
        LLL <- max(length(zeroinfl_prob), length(lambda), length(q))
        if (length(zeroinfl_prob)  != LLL) 
                zeroinfl_prob  <- rep_len(zeroinfl_prob,  LLL)
        if (length(lambda) != LLL) 
                lambda <- rep_len(lambda, LLL)
        if (length(q)      != LLL) 
                q      <- rep_len(q,      LLL)
        
        ans <- ppois(q, lambda)
        ans <- ifelse(q < 0, 0, zeroinfl_prob + (1 - zeroinfl_prob) * ans)

        deflat.limit <- -1 / expm1(lambda)
        ans[zeroinfl_prob < deflat.limit] <- NaN
        ans[zeroinfl_prob > 1] <- NaN
        
        ans
        }