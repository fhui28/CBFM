# Could use remotes to do installaton for countreg, but then likely can not officially put it up on CRAN
Remotes: svn:://svn.r-forge.r-project.org/svnroot/countreg/

#Removed from DESCRIPTION since we now turn off compilation of the TMB C++ files [learning from how James Thorson does it with VAST!]
LinkingTo: 
    TMB,
    RcppEigen


##--------------------------------------
y = simy_train
formula_X = useformula
data = dat_train
B_space = train_basisfunctions
family =  zipoisson()
B_time = NULL
B_spacetime = NULL
offset = NULL
ncores = NULL
trial_size = 1
dofit = TRUE
stderrors = TRUE
start_params = list(betas = NULL, basis_effects_mat = NULL, dispparam = NULL, powerparam = NULL, zeroinfl_prob = NULL)
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
zeroinfl_prob_intercept = new_fit_CBFM_ptest$zeroinfl_prob_intercept
return_correlation = TRUE



object = fitcbfm
newdata = NULL 
manualX = NULL
new_B_space = NULL
new_B_time = NULL
new_B_spacetime = NULL
type = "response"
se_fit = TRUE
coverage = 0.95
ncores = NULL

