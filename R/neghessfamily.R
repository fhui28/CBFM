## Negative second derivative for a bunch of distributions
## Not exported

.neghessfamily <- function(family, eta, y, phi = NULL, trial_size, return_matrix = FALSE) {
     if(family$family[1] %in% c("Beta")) {
          .sbetaprobit <- function(eta, y, phi) {
               ystar <- qlogis(y)
               mu <- betar(link = "logit")$linkinv(eta)   
               mustar <- digamma(mu * phi) - digamma((1 - mu) * phi)
               
               rval <- cbind(phi * (ystar - mustar) * betar(link = "logit")$mu.eta(eta))
               return(rval)
               }
          
          out <- -grad(.sbetaprobit, x = eta, y = y, phi = phi)    
          }
     if(family$family[1] %in% c("binomial")) {
#           if(family$link == "probit") {
#                mu <- binomial(link = "probit")$linkinv(eta)
# 
#                d2mu_deta2 <- grad(dnorm, x = eta)
#                if(is.matrix(eta))
#                     d2mu_deta2 <- matrix(d2mu_deta2, nrow = nrow(y), ncol = ncol(y))
#                
#                out <- (-y/mu^2-(trial_size-y)/(1-mu)^2)*(binomial(link = "probit")$mu.eta(eta)^2) + (y/mu-(trial_size-y)/(1-mu))*d2mu_deta2
#                rm(mu, d2mu_deta2)
#                }
          if(family$link == "logit") {
               mu <- binomial(link = "logit")$linkinv(eta)
               out <- -trial_size*binomial(link = "logit")$var(mu)
               }
          }
     if(family$family[1] %in% c("Gamma")) {
          mu <- Gamma(link = "log")$linkinv(eta)
          out <- -(1/phi) * y / mu
          }
     if(family$family[1] %in% c("gaussian")) {
          out <- -1/phi
          }
     if(family$family[1] == "negative.binomial") {
          mu <- exp(eta)
          out <- -(phi*y + 1) * mu / (1 + phi*mu)^2               
          }
     if(family$family[1] %in% c("poisson")) {
          out <- -exp(eta) 
          }
#      if(family$family[1] == "ztpoisson") {
#           out <- ztpoisson()$mu.eta(eta)
#           }
#      if(family$family[1] == "ztnegative.binomial") {
#           mu <- exp(eta)
#           out <- sztnbinom(y, mu = mu, size = 1/phi, parameter = "mu", drop=FALSE)*mu + hztnbinom(y, mu = mu, size = 1/phi, parameter = "mu", drop=FALSE)*mu^2
#           rm(mu)
#           }

          
     if(return_matrix) {
          if(is.matrix(y))
               return(matrix(-out, nrow = nrow(y), ncol = ncol(y)))
          if(!is.matrix(y))
               return(as.vector(-out))
          }
     if(!return_matrix)
          return(as.vector(-out))
     }
     
     
