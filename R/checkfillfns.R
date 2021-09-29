## Check and fill functions
## Hidden and not exported

.check_B_forms <- function(B_space, B_time, B_spacetime, Sigma = NULL, G = NULL, extra_check = FALSE) {
     if(is.null(B_space) & is.null(B_time) & is.null(B_spacetime))
          stop("At least one of B_space/B_time/B_spacetime must be supplied.")
     
     if(extra_check) {
          if(!is.list(Sigma))
               stop("Sigma should be a list containing at least one of following three elements: 1) space, which is the covariance matrix for the random slopes corresponding to B_space; 2) time, which is the analogous covariance matrix corresponding to B_time; 3) spacetime, which is the analogous covariance matrix corresponding to B_spacetime.")
          if(!is.list(G))
               stop("G should be a list containing at least one of following three elements: 1) space, which is the baseline between-response correlation matrix for the random slopes corresponding to B_space; 2) time, which is the analogous covariance matrix correlation to B_time; 3) spacetime, which is the analogous correlation matrix corresponding to B_spacetime.")
          
          if(!is.null(B_space) & (is.null(G$space) | is.null(Sigma$space)))
               stop("Because B_space is supplied, then Sigma$space and G$space must also be supplied.")
          if(!is.null(B_time) & (is.null(G$time) | is.null(Sigma$time)))
               stop("Because B_time is supplied, then Sigma$time and G$time must also be supplied.")
          if(!is.null(B_spacetime) & (is.null(G$spacetime) | is.null(Sigma$spacetime)))
               stop("Because B_spacetime is supplied, then Sigma$spacetime and G$spacetime must also be supplied.")
          }
     
     }
     

.check_BX <- function(B, X, tol = 0.9) {     
     if(ncol(X) > 1) {
          sel_interceptcols <- which(apply(X, 2, function(x) length(unique(x)) == 1) == TRUE)
          if(length(sel_interceptcols) > 0)
               X2 <- X[,-sel_interceptcols, drop = FALSE]
          if(length(sel_interceptcols) == 0)
               X2 <- X
               
          sel_interceptcols <- which(apply(B, 2, function(x) length(unique(x)) == 1) == TRUE)
          if(length(sel_interceptcols) > 0)
               B2 <- B[,-sel_interceptcols, drop = FALSE]
          if(length(sel_interceptcols) == 0)
               B2 <- B

          getcrosscor <- stats::cor(cbind(as.matrix(B2),as.matrix(X2)))[1:ncol(B2), -(1:ncol(B2))]
          if(any(abs(getcrosscor) > tol)) 
               warning("There are columns of X and B that are highly correlated. This may complicate interpretaton of the regression coefficients from X (i.e., betas), due to potential concerns regarding spatial confounding.")          
          
#           getBcor <- cor(as.matrix(B2))
#           getBcor <- getBcor[lower.tri(getBcor)] 
#           if(any(abs(getBcor) > tol)) 
#                warning("There are columns of formed B that are highly correlated. This could cause potential issues and instability in the fitted model.")          
          }
     }
     

.check_family <- function(family, y, trial_size) {
    if(!(family$family[1] %in% c("gaussian", "Gamma", "negative.binomial", "poisson", "binomial", "tweedie", "beta", 
                                 "zipoisson", "zinegative.binomial", "ztpoisson"))) #"ztnegative.binomial"     
        stop("Family currently not supported. Sorry!")
    if((family$family %in% c("ztpoisson", "ztnegative.binomial")) & any(y < 1))
        message("For zero truncated distributions, observed zero counts will be (obviously) ignored in the model.")
     
    if(family$family[1] == "gaussian" & family$link != "identity")
        stop("Currently Gaussian family only permits the identity link.")
    if(family$family[1] %in% c("Gamma","negative.binomial","poisson", "tweedie", "zipoisson", "zinegative.binomial", "ztpoisson", "ztnegative.binomial") & family$link != "log")
        stop("Supplied family currently only permits the log link function to be used.")
    if(family$family[1] %in% c("beta","binomial") & family$link != "logit")
        stop("Supplied family currently only permits the logit link function to be used.")
    if(family$family[1] %in% c("binomial")) {
        if(!(length(trial_size) %in% c(1,length(y))))
            stop("If the binomial family, is used, then trial_size must equal either a scalar or a matrix of the same dimensions as y.")
        }
     }

     
.check_offset <- function(offset = NULL, y) {
     if(!is.null(offset)) { 
          if(!is.matrix(offset)) 
               stop("offset should be a matrix with the same dimensions as y.")
          if(nrow(offset) != nrow(y)) 
               stop("offset should be a matrix with the same dimensions as y.")
          if(ncol(offset) != ncol(y)) 
               stop("offset should be a matrix with the same dimensions as y.")
          } 
    }

        
.check_ranks <- function(num_spp, rank_G, num_basisfns, rank_Sigma) {
     if(num_spp <= 2)
          message("rank_G ignored for models containing two or less responses.")
     if(num_spp > 2) {
          dof <- 0.5 * ((num_spp - rank_G)^2 - num_spp - rank_G)
          if(dof < 0 | (num_spp <= rank_G)) 
          stop("rank_G is too many factors for the number of responses.")
          }
     
     
     if(num_basisfns <= 2) 
        stop("There should be at least two basis functions included in the model.")
     if(num_basisfns == 2) {
          message("rank_Sigma ignored for models containing two basis functions.")
          }
     if(num_basisfns > 2) {
          dof <- 0.5 * ((num_basisfns - rank_Sigma)^2 - num_basisfns - rank_Sigma)
          if(dof < 0 | (num_basisfns <= rank_Sigma)) 
               stop("rank_Sigma is too many factors for the number of basis functions.")
          }
     }
     

.check_ranks2 <- function(num_spp, which_B_used, G_control, vec_num_basisfns, Sigma_control) {
     for(k0 in 1:length(which_B_used)) { 
          if(which_B_used[k0] == 1) {
               .check_ranks(num_spp = num_spp, rank_G = G_control$rank[sum(which_B_used[1:k0])], 
                            num_basisfns = vec_num_basisfns[k0], rank_Sigma = Sigma_control$rank[sum(which_B_used[1:k0])])
               }
          }
     }

     
.check_start_params <- function(start_params, num_spp, num_basisfns, num_X) {
     if(!is.null(start_params$betas)) {
          if(!is.matrix(start_params$betas))
               stop("start_params$betas should be a matrix of starting values for the regression coefficients corresponding to formula_X.")
          if(nrow(start_params$betas) != num_spp | ncol(start_params$betas) != num_X)
               stop("The dimensions of start_params$betas are not appropriate. Please check and amend.")
          }
     if(!is.null(start_params$basis_effects_mat)) {
          if(!is.matrix(start_params$basis_effects_mat))
               stop("start_params$basis_effects_mat should be a matrix of starting values for the random slopes corresponding to basis functions.")
          if(nrow(start_params$basis_effects_mat) != num_spp | ncol(start_params$basis_effects_mat) != num_basisfns)
               stop("The dimensions of start_params$basis_effects_mat are not appropriate. Please check and amend.")
          }
     if(!is.null(start_params$dispparam)) {
          if(length(start_params$dispparam) != num_spp)
               stop("The length of start_params$dispparam is not appropriate. Please check and amend.")
          }
     if(!is.null(start_params$powerparam)) {
          if(length(start_params$powerparam) != num_spp)
               stop("The length of start_params$powerparam is not appropriate. Please check and amend.")
          }          
     if(!is.null(start_params$zeroinfl_prob)) {
          if(length(start_params$zeroinfl_prob) != num_spp)
               stop("The length of start_params$zeroinfl_prob is not appropriate. Please check and amend.")
          }          
     }

     
.check_X_formula <- function(formula_X, data) {
     formulaX <- as.formula(formula_X)
          
     termsinformula <- as.character(formula_X)
     if(length(termsinformula) == 3)
          termsinformula <- termsinformula[-2]
     formula_X <- as.formula(termsinformula)
          
     return(formula_X)
     }     
     
     
.family_to_counter <- function(family) {
     if(family$family[1] == "beta")
          family_counter <- 1
     if(family$family[1] == "binomial")
          family_counter <- 2
     if(family$family[1] == "Gamma")
          family_counter <- 3
     if(family$family[1] == "negative.binomial")
          family_counter <- 4
     if(family$family[1] == "gaussian")
          family_counter <- 5
     if(family$family[1] == "poisson")
          family_counter <- 6
     if(family$family[1] == "tweedie")
          family_counter <- 7
     if(family$family[1] == "ztnegative.binomial")
          family_counter <- 8
     if(family$family[1] == "ztpoisson")
          family_counter <- 9
     if(family$family[1] == "zipoisson")
          family_counter <- 10
     if(family$family[1] == "zinegative.binomial")
          family_counter <- 11
     
     return(family_counter)
     }
          
.fill_control <- function(control) {
     if(is.null(control$maxit))
          control$maxit <- 100
     if(is.null(control$optim_lower))
          control$optim_lower <- -10
     if(is.null(control$optim_upper))
          control$optim_upper <- 10
     if(is.null(control$tol))
          control$tol <- 1e-4
     if(is.null(control$initial_betas_dampen))
          control$initial_betas_dampen <- 1
     if(is.null(control$subsequent_betas_dampen))
          control$subsequent_betas_dampen <- 0.25
     if(is.null(control$convergence_type))
          control$convergence_type <- "parameters"
     if(is.null(control$ridge))
          control$ridge <- 0
     if(is.null(control$trace))
          control$trace <- 0

     control$convergence_type <- match.arg(control$convergence_type, choices = c("parameters","linear_predictor", "logLik"))
     return(control)
     }
     
     
.fill_G_control <- function(control, which_B_used) {
     if(is.null(control$rank))
          control$rank <- rep(5, sum(which_B_used))
     if(sum(which_B_used) != length(control$rank))
          stop("G_control$rank should be a vector with length depending on whether B_space/B_time/B_spacetime are supplied. Each element corresponds to the rank of G to use for B_space/B_time/B_spacetime. For example, if B_space and B_spacetime are both supplied, then G_control$rank should be a vector with length 2.")
     if(is.null(control$nugget_profile))          
          control$nugget_profile = seq(0.05, 0.95, by = 0.05)
     if(is.null(control$maxit))
          control$maxit <- 100
     if(is.null(control$tol))
          control$tol <- 1e-4
     if(is.null(control$method))
          control$method <- "LA"
     if(is.null(control$inv_method))
          control$inv_method <- "chol2inv"
     if(is.null(control$trace))
          control$trace <- 0
          
     control$method<- match.arg(control$method, choices = c("LA","simple"))
     #control$inv_method <- match.arg(control$inv_method, choices = c("chol2inv","schulz"))
     
     return(control)
     }
     

.fill_Sigma_control <- function(control, which_B_used) {
     if(is.null(control$rank))
          control$rank <- rep(5, sum(which_B_used))
     if(sum(which_B_used) != length(control$rank)) {
          stop("Sigma_control$rank should be a vector with length depending on whether B_space/B_time/B_spacetime are supplied. Each element corresponds to the rank of Sigma to use for B_space/B_time/B_spacetime. For example, if B_space and B_spacetime are both supplied, then Sigma_control$rank should be a vector with length 2.")
          }
     if(is.null(control$maxit))
          control$maxit <- 100
     if(is.null(control$tol))
          control$tol <- 1e-4
     if(is.null(control$method))
          control$method <- "LA"
     if(is.null(control$inv_method))
          control$inv_method <- "chol2inv"
     if(is.null(control$trace))
          control$trace <- 0
          
     control$method<- match.arg(control$method, choices = c("LA","simple"))
     #control$inv_method <- match.arg(control$inv_method, choices = c("chol2inv","schulz"))

     return(control)
     }
     
                

.ifelse_size <- function(trial_size, trial_size_length, j, num_units, family = binomial()) {
     if(family$family[1] != "binomial")
          return(rep(trial_size,2)) # This can be set to anything. But it is set to a vector for compatability with TMB
     if(trial_size_length)
          return(rep(trial_size, num_units))
     if(!trial_size_length)
          return(trial_size[,j])
     }


#Rcpp::cppFunction(depends = "RcppArmadillo", code = '
#     Rcpp::List cppmatinv(const arma::mat &Am) {
#     arma::mat B = inv(Am);
#     return Rcpp::List::create( Rcpp::Named("Imp") = B);
#     }')
