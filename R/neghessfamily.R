## Negative second derivative for a bunch of distributions
## Hidden and not exported

.neghessfamily <- function(family, eta, y, phi = NULL, powerparam = NULL,
                           zieta = NULL, trial_size,
                           return_matrix = FALSE, domore = FALSE, tol = 1e-6) {
     if(family$family[1] %in% c("binomial")) {
          if(family$link == "logit") {
               mu <- plogis(eta)
               out <- trial_size*binomial()$var(mu)
               }
          }
     if(family$family[1] %in% c("poisson")) {
                out <- exp(eta)
                }
     if(family$family[1] == "negative.binomial") {
                mu <- exp(eta)
                out <- (phi*y + 1) * mu / (1 + phi*mu)^2
                }
     if(family$family[1] %in% c("gaussian")) {
                out <- 1/phi
                }
     if(family$family[1] %in% c("Gamma")) {
                out <- y * exp(-eta) / phi
                #mu <- Gamma(link = "log")$linkinv(eta)
                #out <- (1/phi) * y / mu
                }
     if(family$family[1] %in% c("Beta")) {
          out <- grad(.sbetalogit, x = eta, y = y, phi = phi)
          }
     if(family$family[1] == "ztpoisson") {
                out <- .hess_ztpoisson(eta = eta, y = y)
                }
     if(family$family[1] == "ztnegative.binomial") {
                out <- .hess_ztnbinom(eta = eta, y = y, size = 1/phi)
                }
     if(family$family[1] %in% c("tweedie")) {
             two_minus_p <- 2 - powerparam
             one_minus_p <- 1 - powerparam
             exp_2mp_eta <- exp(two_minus_p * eta)
             exp_1mp_eta <- exp(one_minus_p * eta)
             out <- (two_minus_p * exp_2mp_eta - y * one_minus_p * exp_1mp_eta) / phi
             }
     if(family$family[1] %in% c("zipoisson")) {
          lambda <- exp(eta)
          y_zero <- y == 0

          # Pre-allocate and compute rhat only for y==0 cases
          rhat <- numeric(n)
          if(any(y_zero)) {
               rhat[y_zero] <- plogis(zieta[y_zero] + lambda[y_zero])
               }

          out <- lambda * (1 - rhat) * (1 - lambda * rhat)

          if(domore) {
               # out is already the collection of weights for betasbetas and basiseffectsbasiseffects. So we need the other terms involving the zero-inflation component...
               expetalambda <- exp(zieta + lambda)
               expetalambda[expetalambda > .Machine$double.xmax] <- .Machine$double.xmax
               dhat <- expetalambda / (expetalambda + 1)

               phat <- exp(zieta)
               phat <- exp(zieta) / (1 + exp(zieta))
               out_zeroinflzeroinfl <- phat * (1 - phat) - y_zero * dhat * (1 - dhat)
               out_zeroinflzeroinfl[out_zeroinflzeroinfl < 0] <- 0 ## At the moment, needed primarily for zero-inflated and zero-truncated models where weights can be negative (by design?!)
                  
               out_zeroinflbetas <- -(expetalambda * lambda) / (expetalambda + 1)^2
               out_zeroinflbetas[!y_zero] <- 0
               }
             }
        if(family$family[1] %in% c("zinegative.binomial")) { 
             lambda <- exp(eta)
             phat <- plogis(zieta)
             
             # Probability of being in the count component
             f0 <- (1 + phi * lambda)^(-1/phi)
             prob_zero <- phat + (1 - phat) * f0
             w <- ifelse(y > 0, 1, ((1 - phat) * f0) / prob_zero)

             # Derivative of the weight w.r.t eta (only needed for y=0)
             df0 <- -lambda * (1 + phi * lambda)^(-1/phi - 1)
             dw <- ifelse(y > 0, 0, ((1 - phat) * df0 / prob_zero) - ((1 - phat) * f0 * (1 - phat) * df0) / (prob_zero^2))

             term1 <- w * ((-lambda * (1 + phi * y)) / (1 + phi * lambda)^2)
             term2 <- dw * ((y - lambda) / (1 + phi * lambda))
             out <- -(term1 + term2)

             if(domore) {
                 # out is already the collection of weights for betasbetas and basiseffectsbasiseffects. So we need the other terms involving the zero-inflation component...

                 f0 <- dnbinom(0, size = 1/phi, mu = lambda)
                 # w = P(Structural | Y=0)
                 w <- ifelse(y == 0, phat / (phat + (1 - phat) * f0), 0)

                 # For y > 0, only the -p(1-p) term remains
                 out_zeroinflzeroinfl <- -((w * (1 - w)) - (phat * (1 - phat)))
                 #dhat <- exp(zieta) * (1+phi*lambda)^(-1/phi) / (exp(zieta) + (1+phi*lambda)^(-1/phi))^2
                 #phat <- plogis(zieta)
                 #out_zeroinflzeroinfl <- phat * (1-phat) - ((y == 0) * 1) * dhat
                 out_zeroinflzeroinfl[out_zeroinflzeroinfl < 0] <- 0 ## At the moment, needed primarily for zero-inflated and zero-truncated models where weights can be negative (by design?!)


                 # Posterior probability of coming from the count process
                 w_count <- ifelse(y == 0, ((1 - phat) * f0) / (phat + (1 - phat) * f0), 0)
                 grad_eta_part <- -lambda / (1 + phi * lambda )

                 # The cross derivative: d/dzieta [w_count * grad_eta_part]
                 # d(w_count)/d(zieta) = -w_count * (1 - w_count)
                 out_zeroinflbetas <- w_count * (1 - w_count) * grad_eta_part
                 #out_zeroinflbetas <- -(exp(zieta) * lambda * (1+phi*lambda)^(-1/phi-1)) / (exp(zieta) + (1+phi*lambda)^(-1/phi))^2
                 #out_zeroinflbetas[y > 0] <- 0
             }
        }

     out[out < 0] <- 0 ## At the moment, needed primarily for zero-inflated and zero-truncated models where weights can be negative (by design?!)
     out[!is.finite(out)] <- 0 ## Needed when for zero-truncated models where weights on the very rare occasional be stupid?
     if(!domore)
          return(as.vector(out))
     if(domore) {
          if(!(family$family[1] %in% c("zipoisson","zinegative.binomial")))
               return(list(out = as.vector(out)))
          if(family$family[1] %in% c("zipoisson","zinegative.binomial")) {
               out_zeroinflzeroinfl[!is.finite(out_zeroinflzeroinfl)] <- 0 
               out_zeroinflbetas[!is.finite(out_zeroinflbetas)] <- 0 
               return(list(out = as.vector(out), 
                           out_zeroinflzeroinfl = out_zeroinflzeroinfl, # matrix of the same dimension as y
                           out_zeroinflbetas = out_zeroinflbetas # matrix of the same dimension as y)
                           )) 
               }
          }
     }



