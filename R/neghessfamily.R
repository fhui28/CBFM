## Negative second derivative for a bunch of distributions
## Hidden and not exported

.neghessfamily <- function(family, eta, y, phi = NULL, powerparam = NULL, zeroinfl_prob_intercept = NULL, trial_size, 
                           return_matrix = FALSE, domore = FALSE, tol = 1e-5) {
     
        if(family$family[1] %in% c("Beta")) {
                .sbetalogit <- function(eta, y, phi) {
                ystar <- qlogis(y)
                mu <- betar(link = "logit")$linkinv(eta)   
                mustar <- digamma(mu * phi) - digamma((1 - mu) * phi)
               
                rval <- cbind(phi * (ystar - mustar) * betar(link = "logit")$mu.eta(eta))
                return(rval)
                }
          
                out <- grad(.sbetalogit, x = eta, y = y, phi = phi)    
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
                        mu <- plogis(eta)
                        out <- trial_size*binomial()$var(mu)
                        }
                }
        if(family$family[1] %in% c("Gamma")) {
                mu <- Gamma(link = "log")$linkinv(eta)
                out <- (1/phi) * y / mu
                }
        if(family$family[1] %in% c("gaussian")) {
                out <- 1/phi
                }
        if(family$family[1] == "negative.binomial") {
                mu <- exp(eta)
                out <- (phi*y + 1) * mu / (1 + phi*mu)^2               
                }
        if(family$family[1] %in% c("poisson")) {
                out <- exp(eta) 
                }
#      if(family$family[1] == "ztpoisson") {
#           out <- ztpoisson()$mu.eta(eta)
#           }
#      if(family$family[1] == "ztnegative.binomial") {
#           mu <- exp(eta)
#           out <- sztnbinom(y, mu = mu, size = 1/phi, parameter = "mu", drop=FALSE)*mu + hztnbinom(y, mu = mu, size = 1/phi, parameter = "mu", drop=FALSE)*mu^2
#           rm(mu)
#           }
     if(family$family[1] %in% c("tweedie")) {
             two_minus_powerparam <- 2-powerparam
             exp_two_minus_powerparam_linpred <- exp(two_minus_powerparam * eta)
             exp_one_minus_powerparam_linpred <- exp((two_minus_powerparam-1) * eta)
             out <-  1/phi * (two_minus_powerparam*exp_two_minus_powerparam_linpred - y*(two_minus_powerparam-1)*exp_one_minus_powerparam_linpred)    
             }
     if(family$family[1] %in% c("zipoisson")) {
                lambda <- exp(eta)
                rhat <- numeric(length(y))
                rhat[y == 0] <- as.vector(plogis(zeroinfl_prob_intercept + lambda))[y == 0]
                out <- lambda * (1-rhat) * (1-lambda*rhat)
         
                if(domore) {
                 # out is already the collection of weights for betasbetas and basiseffectsbasiseffects. So we need the other terms involving the zero-inflation component...
                        dhat <- exp(zeroinfl_prob_intercept) / (exp(zeroinfl_prob_intercept) + exp(-lambda))
                        phat <- plogis(zeroinfl_prob_intercept)
                        out_zeroinflzeroinfl <- phat * (1-phat) - ((y == 0) * 1) * dhat * (1-dhat)
                
                        out_zeroinflbetas <- -(exp(zeroinfl_prob_intercept) * lambda * exp(-lambda)) / (exp(zeroinfl_prob_intercept) + exp(-lambda))^2
                        out_zeroinflbetas[y > 0] <- 0
                        }
                }

        
        out[out < tol] <- tol ## At the moment, needed primarily for ZIP models where weights have be negative (by design?!)

        if(!domore)
                return(as.vector(out))
        if(domore) {
                if(family$family[1] != "zipoisson")
                        return(list(out = as.vector(out)))
                if(family$family[1] == "zipoisson")
                        return(list(out = as.vector(out), 
                                    out_zeroinflzeroinfl = out_zeroinflzeroinfl, # matrix of the same dimension as y
                                    out_zeroinflbetas = out_zeroinflbetas # matrix of the same dimension as y)
                                    )) 
                }
        }
     
     

## E-step functions for zero-inflated distributions -- calculate the posterior probability of being in the zero-inflation component
## Hidden and not exported
.estep_fn <- function(family, cwfit, y, X, B) {
        num_units <- nrow(y)
        num_spp <- ncol(y)
        out <- Matrix(0, nrow = num_units, ncol = num_spp, sparse = TRUE)
        if(family$family == "zipoisson") {
                fitvals <- exp(tcrossprod(X, cwfit$betas) + tcrossprod(B, cwfit$basis_effects_mat))
                zeroinfl_prob <- plogis(cwfit$zeroinfl_prob_intercept)
                
                for(j in 1:num_spp) {
                        sel_zerospp <- which(y[,j] == 0)
                        out[sel_zerospp,j] <- zeroinfl_prob[j] / (zeroinfl_prob[j] + (1-zeroinfl_prob[j])*dpois(0, fitvals[sel_zerospp,j]))
                        }
                }
        
        return(out)
        }


