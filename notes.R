# Could use remotes to do installaton for countreg, but then likely can not officially put it up on CRAN
Remotes: svn:://svn.r-forge.r-project.org/svnroot/countreg/

# Use document::document to oxygenise a single R script to help file e.g.,      
document::document(file_name = "R/createlife.R", check_package = FALSE, output_directory = "R/tmpman/")

#Removed from DESCRIPTION since we now turn off compilation of the TMB C++ files [learning from how James Thorson does it with VAST!]
LinkingTo: 
    TMB,
    RcppEigen

    
y = simy_train
formula_X = useformula
data = dat_train
B_space = train_basisfunctions
family = nb2()
B_time = NULL
B_spacetime = NULL
offset = NULL
ncores = NULL
trial_size = 1
dofit = TRUE
stderrors = TRUE
start_params = list(betas = NULL, basis_effects_mat = NULL, dispparam = NULL, powerparam = NULL)
TMB_directories = list(cpp = system.file("executables", package = "CBFM"), compile = system.file("executables", package = "CBFM"))
control = list(maxit = 1000, optim_lower = -5, optim_upper = 5, convergence_type = "parameters", tol = 1e-4, seed = NULL, trace = 1, ridge = 0)
Sigma_control = list(rank = 5, maxit = 1000, tol = 1e-4, method = "LA", trace = 0)
G_control = list(rank = 5, nugget_profile = seq(0.05, 0.95, by = 0.05), maxit = 1000, tol = 1e-4, method = "LA", trace = 0)

    
    
Ginv = new_LoadingnuggetG_space$covinv
basis_effects_mat = new_fit_CBFM_ptest$basis_effects_mat[,1:num_spacebasisfns,drop=FALSE]
Sigmainv = new_LoadingnuggetSigma_space$covinv
B = B_space
y_vec = as.vector(y)
linpred_vec = c(new_fit_CBFM_ptest$linear_predictor)
dispparam = new_fit_CBFM_ptest$dispparam
powerparam = new_fit_CBFM_ptest$powerparam
return_correlation = TRUE



#' compare_betas <- data.frame(
#' spp_slopes %>% 
#' as.data.frame %>% 
#' rownames_to_column(var = "spp") %>% 
#' pivot_longer(temp:O2, values_to = "true"))
#' compare_betas$stacked <- coef(fitstacked)[-1,] %>% 
#' as.vector
#' compare_betas$cbfm <- coef(fitcbfm)[,-1] %>% 
#' t %>% 
#' as.vector
#' compare_betas$name <- fct_inorder(compare_betas$name)
#' 
#' ggplot(compare_betas, aes(x = true, y = cbfm)) +
#' geom_point() +
#' facet_wrap(name ~ ., nrow = 2) +
#' geom_point(aes(x = true, y = stacked), data = compare_betas, color = "blue") +
#' geom_abline(intercept = 0, slope = 1, linetype = 2) +
#' theme_bw()
#' 
