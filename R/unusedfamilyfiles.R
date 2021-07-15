# ## Truncated negative binomial with log-link
# ztnegative.binomial <- function() {
#      link <- "log"
#      linkinv <- function(eta, phi) {
#           return(mean_ztnbinom(mu = pmax(exp(eta), .Machine$double.eps), size = 1/phi))
#           }
#      variance <- function(mu, phi) 
#           return(var_ztnbinom(mu = mu, size = 1/phi))
# 
#      structure(list(family = "ztnegative.binomial", linkinv = linkinv, variance = variance, link = link), class = "family")
#      }


# ## Beta distribution with probit-link
# betaprobitfam <- function() {
#      link <- "probit"
#      linkfun <- function(mu) 
#           return(qnorm(mu))
#      linkinv <- function(eta) 
#           return(pnorm(eta))
#      mu.eta <- function(eta) 
#           return(pmax(dnorm(eta), .Machine$double.eps))
#      variance <- function(mu, phi) 
#           return(mu*(1-mu)/(1+phi))
# 
#      structure(list(family = "Beta", linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta, variance = variance, link = link), class = "family")
#      }

