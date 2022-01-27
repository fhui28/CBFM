## Negative second derivative for a bunch of distributions
## Hidden and not exported
## Some help from Wolfram alpha differentiation online!

.neghessfamily <- function(family, eta, y, phi = NULL, powerparam = NULL, zeroinfl_prob_intercept = NULL, trial_size, return_matrix = FALSE, 
                           domore = FALSE, tol = 1e-6) {
        if(family$family[1] %in% c("Beta")) {
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
        if(family$family[1] == "ztpoisson") {
                out <- .hess_ztpoisson(eta = eta, y = y)
                # lambda <- exp(eta)
                # out <- exp(eta) + exp(eta)/(exp(lambda)-1) + exp(2*eta) * exp(lambda) / (exp(lambda)-1)^2
                # rm(lambda)
                }
        if(family$family[1] == "ztnegative.binomial") {
                out <- .hess_ztnbinom(eta = eta, y = y, size = 1/phi)
                }
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
        if(family$family[1] %in% c("zinegative.binomial")) { 
             lambda <- exp(eta)
             phat <- plogis(zeroinfl_prob_intercept)
             
             score_beta <- function(x, phi, phat) {
                exp(x) * (1 + phi * exp(x))^(-1) / (1 + phat * (1 + phi * exp(x))^(1/phi))     
                }
             out <- grad(score_beta, x = eta, phi = phi, phat = phat)  * as.numeric(y == 0) # Being lazy here!
             out <- out + (lambda * (1 + phi * y) / (1 + phi * lambda)^2) * as.numeric(y > 0)

             if(domore) {
                 # out is already the collection of weights for betasbetas and basiseffectsbasiseffects. So we need the other terms involving the zero-inflation component...
                     dhat <- exp(zeroinfl_prob_intercept) * (1+phi*lambda)^(-1/phi) / (exp(zeroinfl_prob_intercept) + (1+phi*lambda)^(-1/phi))^2
                     phat <- plogis(zeroinfl_prob_intercept)
                     out_zeroinflzeroinfl <- phat * (1-phat) - ((y == 0) * 1) * dhat
                
                     out_zeroinflbetas <- -(exp(zeroinfl_prob_intercept)*lambda*(1+phi*lambda)^(-1/phi-1)) / (exp(zeroinfl_prob_intercept) + (1+phi*lambda)^(-1/phi))^2
                     out_zeroinflbetas[y > 0] <- 0
                     }
             }

        out[out < 0] <- 0 ## At the moment, needed primarily for zero-inflated and zero-truncated models where weights can be negative (by design?!)
        out[!is.finite(out)] <- 0 ## Needed when for zero-truncated models where weights on the very rare occasional be stupid?

        if(!domore)
                return(as.vector(out))
        if(domore) {
                if(!(family$family[1] %in% c("zipoisson","zinegative.binomial")))
                        return(list(out = as.vector(out)))
                if(family$family[1] %in% c("zipoisson","zinegative.binomial"))
                        return(list(out = as.vector(out), 
                                    out_zeroinflzeroinfl = out_zeroinflzeroinfl, # matrix of the same dimension as y
                                    out_zeroinflbetas = out_zeroinflbetas # matrix of the same dimension as y)
                                    )) 
                }
        }
     
     

