## Check and fill functions
## Hidden and not exported

# Only used to create_CBFM_life
.check_B_forms <- function(B_space, B_time, B_spacetime, Sigma = NULL, G = NULL, extra_check = FALSE) {
     if(is.null(B_space) & is.null(B_time) & is.null(B_spacetime))
          stop("At least one of B_space/B_time/B_spacetime must be supplied.")
     
     if(!is.null(B_space)) {
          if(!is.matrix(B_space))
               stop("B_space must be a matrix.")
          }
     if(!is.null(B_time)) {
          if(!is.matrix(B_time))
               stop("B_time must be a matrix.")
     }
     if(!is.null(B_spacetime)) {
          if(!is.matrix(B_spacetime))
               stop("B_spacetime must be a matrix.")
     }
     
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


     
.check_customSigma_Gstructure <- function(Sigma_control, G_control, which_B_used) {
     if(is.null(Sigma_control[["custom_space"]]) & which_B_used[1]) {
          if(G_control$structure[1] != "unstructured")
               stop("If Sigma_control$custom_space is not supplied i.e., it is estimated, then the corresponding element in G_control$structure can only be set to \"unstructured\" i.e., an unrestructured, possibly rank-reduced correlation matrix.")
          }
     if(is.null(Sigma_control[["custom_time"]]) & which_B_used[2]) {
          if(G_control$structure[sum(which_B_used[1:2])] != "unstructured")
               stop("If Sigma_control$custom_time is not supplied i.e., it is estimated, then the corresponding element in G_control$structure can only be set to \"unstructured\" i.e., an unrestructured, possibly rank-reduced correlation matrix.")
          }
     if(is.null(Sigma_control[["custom_spacetime"]]) & which_B_used[3]) {
          if(G_control$structure[sum(which_B_used[1:3])] != "unstructured")
               stop("If Sigma_control$custom_spacetime is not supplied i.e., it is estimated, then the corresponding element in G_control$structure can only be set to \"unstructured\" i.e., an unrestructured, possibly rank-reduced correlation matrix.")
          }
     
     
     if(any(G_control$structure %in% c("identity", "homogeneous"))) {
          warning("The default choice of G_control$structure = \"unstructured\" is not being used. Using these other options should not be done unless you know what are doing, especially as very little checks are made to ensure parameter identifiability of the model in these settings!")
          }
     
     
     if(which_B_used[1]) {
          if(G_control$structure[1] %in% c("identity", "homogeneous")) { if(is.null(Sigma_control[["custom_space"]])) {
               stop("For any G_control$structure set to \"identity\" or \"homogeneous\", the corresponding element of Sigma_control$custom_xxx must be supplied to ensure parameter identifiability of the model in these settings")
               } }
          }
     
     if(which_B_used[2]) {
          if(G_control$structure[sum(which_B_used[1:2])] %in% c("identity", "homogeneous")) { if(is.null(Sigma_control[["custom_time"]])) {
                    stop("For any G_control$structure set to \"identity\" or \"homogeneous\", the corresponding element of Sigma_control$custom_xxx must be supplied to ensure parameter identifiability of the model in these settings")
               } }
          }
          
     if(which_B_used[3]) {
          if(G_control$structure[sum(which_B_used[1:3])] %in% c("identity", "homogeneous")) { if(is.null(Sigma_control[["custom_spacetime"]])) {
               stop("For any G_control$structure set to \"identity\" or \"homogeneous\", the corresponding element of Sigma_control$custom_xxx must be supplied to ensure parameter identifiability of the model in these settings")
               } }
          }
     
     }



.check_family <- function(family, y, trial_size) {
    if(!(family$family[1] %in% c("gaussian", "Gamma", "negative.binomial", "poisson", "binomial", "tweedie", "beta", 
                                 "zipoisson", "zinegative.binomial", "ztpoisson", "ztnegative.binomial")))     
        stop("Family currently not supported. Sorry!")
    if((family$family %in% c("ztpoisson", "ztnegative.binomial")) & any(y < 1))
        message("For zero truncated distributions, observed zero counts will be ignored in the model.")
     
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


.check_G_correlation <- function(custom_Sigma, G_structure) {
     if(is.null(custom_Sigma)) #' If no custom Sigma is supplied, then G must be estimated as a correlation matrix
          out <- 1
     if(!is.null(custom_Sigma)) { #' If custom Sigma is supplied, then G can be estimated as correlation or covariance depending on structure
          if(G_structure == "unstructured")
               out <- 0
          if(G_structure == "identity")
               out <- 1
          if(G_structure == "homogeneous")
               out <- 1
          }
     
     return(out)
     } 
     

.check_nonzeromeans <- function(nonzeromean_B_space, nonzeromean_B_time, nonzeromean_B_spacetime) {
     if(nonzeromean_B_space)
          message("A non-zero mean vector is being used for the distribution of the spatial basis function coefficients. Please check this is what you want!")
     if(nonzeromean_B_time)
          message("A non-zero mean vector is being used for the distribution of the temporal basis function coefficients. Please check this is what you want!")
     if(nonzeromean_B_spacetime)
          message("A non-zero mean vector is being used for the distribution of the spatio-temporal basis function coefficients. Please check this is what you want!")
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
     #if(num_spp <= 2)
     #     message("rank_G ignored for models containing two or less responses.")
     if(num_spp > 2) {
         # Full rank does nothing
         
         if(rank_G != "full") {
             rank_G <- as.numeric(rank_G)
         
             dof <- 0.5 * ((num_spp - rank_G)^2 - num_spp - rank_G)
             if(dof < 0 | (num_spp <= rank_G)) 
                stop("rank_G is too many factors for the number of responses.")
            }
         }
     
     
     #if(num_basisfns <= 2) 
     #   stop("There should be at least two basis functions included in the model.")
     if(num_basisfns == 2) {
          message("rank_Sigma ignored for models containing two basis functions.")
          }
     if(num_basisfns > 2) {
         # Full rank does nothing
         
         if(rank_Sigma != "full") {
             rank_Sigma <- as.numeric(rank_Sigma)
          
             dof <- 0.5 * ((num_basisfns - rank_Sigma)^2 - num_basisfns - rank_Sigma)
             if(dof < 0 | (num_basisfns <= rank_Sigma)) 
                 stop("rank_Sigma is too many factors for the number of basis functions.")
             }
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
               stop("start_params$betas should be a matrix of starting values for the regression coefficients corresponding to formula.")
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

     
.check_X_formula <- function(formula, data) {
     formulaX <- as.formula(formula)
          
     termsinformula <- as.character(formula)
     if(length(termsinformula) == 3)
          termsinformula <- termsinformula[-2]
     formula <- formula(paste(termsinformula, collapse = " "))
          
     return(formula)
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
          
.fill_control <- function(control, num_spp, which_B_used, nonzeromean_B_space, nonzeromean_B_time, nonzeromean_B_spacetime) {
     if(is.null(control$maxit))
          control$maxit <- 100
     if(is.null(control$inner_maxit))
          control$inner_maxit <- 1
     if(is.null(control$optim_lower))
          control$optim_lower <- -50
     if(is.null(control$optim_upper))
          control$optim_upper <- 50
     if(is.null(control$tol))
          control$tol <- 1e-4
     if(is.null(control$final_maxit))
          control$final_maxit <- 100
     if(is.null(control$initial_betas_dampen))
          control$initial_betas_dampen <- 1
     if(!is.null(control$initial_betas_dampen)) {
        if(!(length(control$initial_betas_dampen) %in% c(1)))
            stop("control$initial_betas_dampen should either be a scalar.") # or a vector equal to the number of species i.e., ncol(y).
        }
     if(is.null(control$subsequent_betas_dampen))
          control$subsequent_betas_dampen <- 0.25
     if(!is.null(control$subsequent_betas_dampen)) {
        if(!(length(control$subsequent_betas_dampen) %in% c(1)))
            stop("control$subsequent_betas_dampen should either be a scalar.") #or a vector equal to the number of species i.e., ncol(y).
        }
     
     if(nonzeromean_B_space == TRUE & which_B_used[1] == 0)
         stop("If B_space is not supplied, then nonzeromean_B_space can not be set to TRUE.")
     if(nonzeromean_B_time == TRUE & which_B_used[2] == 0)
         stop("If B_time is not supplied, then nonzeromean_B_time can not be set to TRUE.")
     if(nonzeromean_B_spacetime == TRUE & which_B_used[3] == 0)
         stop("If B_spacetime is not supplied, then nonzeromean_B_spacetime can not be set to TRUE.")
     
     if(is.null(control$convergence_type))
          control$convergence_type <- "parameters_MSE"
     if(is.null(control$gam_method))
          control$gam_method <- "REML"
     if(is.null(control$initial_ridge))
          control$initial_ridge <- 0
     if(is.null(control$ridge))
          control$ridge <- 0
     if(is.null(control$initial_ziridge))
          control$initial_ziridge <- 0
     if(is.null(control$ziridge))
          control$ziridge <- 0
     if(is.null(control$trace))
          control$trace <- 0

     control$convergence_type <- match.arg(control$convergence_type, choices = c("parameters_MSE", "parameters_norm", "parameters_relative", "logLik_relative"))
     control$gam_method <- match.arg(control$gam_method, choices = c("ML","REML", "GCV.Cp", "GACV.Cp"))
     return(control)
     }
     
     
.fill_G_control <- function(control, which_B_used, num_spp, Sigma_control) {
     if(is.null(control$rank))
          control$rank <- rep(5, sum(which_B_used))
     if(length(control$rank) == 1)
          control$rank <- rep(control$rank, sum(which_B_used))
     
     if(sum(which_B_used) != length(control$rank))
          stop("G_control$rank should be a vector with length depending on whether B_space/B_time/B_spacetime are supplied. Each element corresponds to the rank of G to use for B_space/B_time/B_spacetime. For example, if B_space and B_spacetime are both supplied, then G_control$rank should be a vector with length 2.
               Please note ranks still needs to be supplied even when custom Gs are used (although the corresponding element is ignored in such case).")
     
     if(is.null(control$structure))
          control$structure <- rep("unstructured", sum(which_B_used))
     if(length(control$structure) == 1)
          control$structure <- rep(control$structure, sum(which_B_used))
     if(sum(which_B_used) != length(control$structure))
          stop("G_control$structure should be a vector with length depending on whether B_space/B_time/B_spacetime are supplied. Each element corresponds to the structure of G to use for B_space/B_time/B_spacetime. For example, if B_space and B_spacetime are both supplied, then G_control$structure should be a vector with length 2.
               Please note structure still needs to be supplied even when custom Gs are used (although the corresponding element is ignored in such case).")
     # if(is.null(control$min_sp))
     #      control$min_sp <- rep(0, sum(which_B_used))
     # if(length(control$min_sp) == 1)
     #      control$min_sp <- rep(control$min_sp, sum(which_B_used))
     # if(sum(which_B_used) != length(control$min_sp))
     #      stop("G_control$min_sp should be a vector with length depending on whether B_space/B_time/B_spacetime are supplied. Each element corresponds to the minimum smoothing parameter to apply to an (identity structured) G for B_space/B_time/B_spacetime. For example, if B_space and B_spacetime are both supplied, then G_control$min_sp should be a vector with length 2, irrespective of how many of G_control$structure was set to an identity structure.
     #           Put simply, please note min_sp still needs to be supplied even when an unstructured G is used, and even when custom Gs are used (although the corresponding element is ignored in such cases).")
     
     if(is.null(control$nugget_profile))          
          control$nugget_profile = seq(0.05, 0.95, by = 0.05)
     if(is.null(control$maxit))
          control$maxit <- 100
     if(is.null(control$tol))
          control$tol <- 1e-4
     if(is.null(control$method))
          control$method <- "REML"
     if(is.null(control$inv_method))
          control$inv_method <- "chol2inv"
     if(is.null(control$trace))
          control$trace <- 0
          
     control$which_custom_G_used <- c(0,0,0)
     if(!is.null(control[["custom_space"]])) {
        if(which_B_used[1] == 0)
            stop("Please do not supply G_control$custom_space if B_space is also not supplied.")
        if(nrow(control[["custom_space"]]) != num_spp | ncol(control[["custom_space"]]) != num_spp)
            stop("G_control$custom_space should be a square matrix with the same dimensions as ncol(y).") 
         control$which_custom_G_used[1] <- 1
         }
     if(!is.null(control[["custom_time"]])) {
         if(which_B_used[2] == 0)
            stop("Please do not supply G_control$custom_time if B_time is also not supplied.")
        if(nrow(control[["custom_time"]]) != num_spp | ncol(control[["custom_time"]]) != num_spp)
            stop("G_control$custom_time should be a square matrix with the same dimensions as ncol(y).") 
        control$which_custom_G_used[2] <- 1
        }
     if(!is.null(control[["custom_spacetime"]])) {
         if(which_B_used[3] == 0) 
            stop("Please do not supply G_control$custom_spacetime if B_spacetime is also not supplied.")
        if(nrow(control[["custom_spacetime"]]) != num_spp | ncol(control[["custom_spacetime"]]) != num_spp)
            stop("G_control$custom_spacetime should be a square matrix with the same dimensions as ncol(y).") 
        control$which_custom_G_used[3] <- 1
        }
          
     control$method <- match.arg(control$method, choices = c("REML", "simple", "ML"))
    #control$inv_method <- match.arg(control$inv_method, choices = c("chol2inv","schulz"))
     
    return(control)
    }
     

.fill_Sigma_control <- function(control, 
                                which_B_used, 
                                num_spacebasisfns, num_timebasisfns, num_spacetimebasisfns, 
                                G_control) {
    if(is.null(control$rank))
          control$rank <- rep(5, sum(which_B_used))
     if(length(control$rank) == 1)
          control$rank <- rep(control$rank, sum(which_B_used))
     if(sum(which_B_used) != length(control$rank)) {
          stop("Sigma_control$rank should be a vector with length depending on whether B_space/B_time/B_spacetime are supplied. Each element corresponds to the rank of Sigma to use for B_space/B_time/B_spacetime. For example, if B_space and B_spacetime are both supplied, then Sigma_control$rank should be a vector with length 2. 
               Please note ranks still needs to be supplied even when custom Sigmas are used (although the corresponding rank is ignored in such case).")
          }
     
     if(is.null(control$maxit))
          control$maxit <- 100
     if(is.null(control$tol))
          control$tol <- 1e-4
     if(is.null(control$method))
          control$method <- "REML"
     if(is.null(control$inv_method))
          control$inv_method <- "chol2inv"
     if(is.null(control$trace))
          control$trace <- 0
    
     control$which_custom_Sigma_used <- c(0,0,0)
     if(!is.null(control[["custom_space"]])) {
        if(which_B_used[1] == 0)
            stop("Please do not supply Sigma_control$custom_space if B_space is also not supplied.")
        if(is.matrix(control[["custom_space"]])) {
            if(nrow(control[["custom_space"]]) != num_spacebasisfns | ncol(control[["custom_space"]]) != num_spacebasisfns)
                stop("Sigma_control$custom_space should be a square matrix with the same dimensions as ncol(B_space).")             
            }
        if(is.list(control[["custom_space"]])) {
               if(length(control[["custom_space"]]) == 1)
                    stop("If the list Sigma_control$custom_space is of length 1 i.e., only contains one matrix, please reformat Sigma_control$custom_space to just be a single matrix instead of a list.")
               for(j in 1:length(control[["custom_space"]])) {
                    if(nrow(control[["custom_space"]][[j]]) != num_spacebasisfns | ncol(control[["custom_space"]][[j]]) != num_spacebasisfns)
                         stop("Each element in the list Sigma_control$custom_space should be a square matrix with the same dimensions as ncol(B_space).")
                    } 
               }
          if(is.list(control[["custom_space"]])) {
               if(G_control$structure[sum(which_B_used[1])] != "homogeneous") #' I am almost certain you can set this identity as well and it would work. But anyway...
                    stop("If multiple (a list of) matrices are supplied to Sigma_control$custom_space, then the corresponding element in G_control$structure must be set to \"homogeneous\"...this is a current constraint of CBFM. Sorry!")
          }
         control$which_custom_Sigma_used[1] <- 1
         }
     if(!is.null(control[["custom_time"]])) {
         if(which_B_used[2] == 0)
            stop("Please do not supply Sigma_control$custom_time if B_time is also not supplied.")
          if(is.matrix(control[["custom_time"]])) {
               if(nrow(control[["custom_time"]]) != num_timebasisfns | ncol(control[["custom_time"]]) != num_timebasisfns)
                    stop("Sigma_control$custom_time should be a square matrix with the same dimensions as ncol(B_time).") 
               }
          if(is.list(control[["custom_time"]])) {
               if(length(control[["custom_time"]]) == 1)
                    stop("If the list Sigma_control$custom_time is of length 1 i.e., only contains one matrix, please reformat Sigma_control$custom_time to just be a single matrix instead of a list.")
               for(j in 1:length(control[["custom_time"]])) {
                    if(nrow(control[["custom_time"]][[j]]) != num_timebasisfns | ncol(control[["custom_time"]][[j]]) != num_timebasisfns)
                         stop("Each element in the list Sigma_control$custom_time should be a square matrix with the same dimensions as ncol(B_time).")
                    } 
               }
          if(is.list(control[["custom_time"]])) {
               if(G_control$structure[sum(which_B_used[1:2])] != "homogeneous") #' I am almost certain you can set this identity as well and it would work. But anyway...
                    stop("If multiple (a list of) matrices are supplied to Sigma_control$custom_time, then the corresponding element in G_control$structure must be set to \"homogeneous\"...this is a current constraint of CBFM. Sorry!")
          }
          control$which_custom_Sigma_used[2] <- 1
        }
     if(!is.null(control[["custom_spacetime"]])) {
         if(which_B_used[3] == 0) 
            stop("Please do not supply Sigma_control$custom_spacetime if B_spacetime is also not supplied.")
          if(is.matrix(control[["custom_spacetime"]])) {
               if(nrow(control[["custom_spacetime"]]) != num_spacetimebasisfns | ncol(control[["custom_spacetime"]]) != num_spacetimebasisfns)
                    stop("Sigma_control$custom_spacetime should be a square matrix with the same dimensions as ncol(B_spacetime).")
               }
          if(is.list(control[["custom_spacetime"]])) {
               if(length(control[["custom_spacetime"]]) == 1)
                    stop("If the list Sigma_control$custom_spacetime is of length 1 i.e., only contains one matrix, please reformat Sigma_control$custom_spacetime to just be a single matrix instead of a list.")
               for(j in 1:length(control[["custom_spacetime"]])) {
                    if(nrow(control[["custom_spacetime"]][[j]]) != num_spacetimebasisfns | ncol(control[["custom_spacetime"]][[j]]) != num_spacetimebasisfns)
                    stop("Each element in the list Sigma_control$custom_spacetime should be a square matrix with the same dimensions as ncol(B_spacetime).")
                    } 
               }
          if(is.list(control[["custom_spacetime"]])) {
               if(G_control$structure[sum(which_B_used[1:3])] != "homogeneous") #' I am almost certain you can set this identity as well and it would work. But anyway...
                    stop("If multiple (a list of) matrices are supplied to Sigma_control$custom_spacetime, then the corresponding element in G_control$structure must be set to \"homogeneous\"...this is a current constraint of CBFM. Sorry!")
               }
        control$which_custom_Sigma_used[3] <- 1
        }

     control$method <- match.arg(control$method, choices = c("REML","simple","ML"))
     #control$inv_method <- match.arg(control$inv_method, choices = c("chol2inv","schulz"))
     
     if(is.null(control$lower_space_lambdas))
          control$lower_space_lambdas <- 0
     if(is.null(control$upper_space_lambdas))
          control$upper_space_lambdas <- Inf
     if(is.null(control$lower_time_lambdas))
          control$lower_time_lambdas <- 0
     if(is.null(control$upper_time_lambdas))
          control$upper_time_lambdas <- Inf
     if(is.null(control$lower_spacetime_lambdas))
          control$lower_spacetime_lambdas <- 0
     if(is.null(control$upper_spacetime_lambdas))
          control$upper_spacetime_lambdas <- Inf
     
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
