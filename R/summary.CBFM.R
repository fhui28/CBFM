#' @title Summary for a CBFM bit
#' 
#' @description 
#' `r lifecycle::badge("experimental")`
#'
#' Takes a fitted \code{CBFM} object and produces some useful summaries from it.
#'
#' @param object An object of class "CBFM".
#' @param coverage The coverage probability of the confidence intervals for the regression coefficients. Defaults to 0.95, which corresponds to 95% confidence intervals.
#' @param digits The number of significant figures to use when printing.
#' @param ... Not used.
#'
#' @details 
#' Currently, the function is pretty basic in set up, returning the estimated matrix of species-specific regression coefficients corresponding to the model matrix created, the estimated matrix of species-specific regression coefficients corresponding to the combined matrix of basis functions, and the estimated vector of species-specific probabilities of zero-inflation (on the logit scale), for distributions which require one. If the \code{object$stderrors == TRUE}, then standard errors, Wald-tests, and Wald confidence intervals are also returned. As discussed in [CBFM()], the Wald tests are based on the Bayesian posterior covariance matrix of the coefficients.   
#' 
#' Please note all P-values are computed without considering uncertainty in the smoothing parameter estimates.
#'  
#' This summary function is currently **extremely** basic, and in future versions of the package we are hoping to include much more (hopefully useful!) information. Apologies for the delay! 
#'
#' @return An object of class "summary.CBFM" which includes the following components, not necessarily in the order below (and as appropriate):
#' \describe{
#' \item{call: }{The matched function call of \code{object}.}
#' \item{betas: }{The estimated matrix of species-specific regression coefficients corresponding to the model matrix created, rounded.}
#' \item{basis_effects_mat: }{The estimated matrix of species-specific regression coefficients corresponding to the combined matrix of basis functions. }
#' \item{zeroinfl_prob_intercept: }{The estimated vector of species-specific probabilities of zero-inflation, for distributions which require one, rounded. *Note this is presented on the logit scale*, that is the model returns \eqn{log(\pi_j/(1-\pi_j))} where \eqn{\pi_j} is the probability of zero-inflation. This is the same as the intercept term of a logistic regression model for the probabilities of zero-inflation, hence the name. } 
#' \item{betas_results: }{If the \code{object$stderrors == TRUE}, then a data frame containing the point estimates, standard errors, corresponding Wald statistics and P-values, and the lower and upper limit of Wald confidence intervals for the regression coefficients corresponding to the model matrix created.}
#' \item{zeroinf_prob__intercept_resultstab: }{If the \code{object$stderrors == TRUE}, then a data frame containing the point estimates, standard errors, corresponding Wald statistics and P-values, and the lower and upper limit of Wald confidence intervals for the probabilities of zero-inflation (if appropriate). Please note that the Wald-test is a test of whether the intercept is statistically different from i.e., whether the probability of zero-inflation is statistically different from 0.5. This may not be that useful in practice. }
#' }
#'
#' @details # Warning
#' Note that if the model matrix created includes smoothing terms as in a generalized additive model or GAM, the Wald P-values and corresponding Wald confidence intervals are approximate and neglect corresponding smoothing parameter uncertainty; please see [mgcv::summary.gam()] for more details. Besides, with the current primitive set up of this function, the currently Wald-type inference for individual smoothing coefficients is unlikely to be all that useful anyway!
#' 
#' @author Francis K.C. Hui <fhui28@gmail.com>, Chris Haak
#' 
#' @seealso [CBFM()] for fitting CBFMs and [coef.CBFM()] for extracting regression coefficients related to the covariates
#' 
#' @examples
#' \donttest{
#' # Please see examples in the help file for the main CBFM function 
#' }
#'
#' @export
#' @importFrom stats qnorm pnorm plogis
#' @md


summary.CBFM <- function(object, coverage = 0.95, digits = max(3L, getOption("digits") - 3L), ...) {
    num_spp <- nrow(object$betas)
    num_basisfns <- nrow(object$basis_effects_mat)
     
    summary_output <- list(call = object$call, 
                           betas = round(object$betas, digits), 
                           basis_effects_mat = round(object$basis_effects_mat, digits))
    if(object$family$family[1] %in% c("zipoisson","zinegative.binomial"))
        summary_output$zeroinfl_prob_intercept <- round(object$zeroinfl_prob_intercept, digits)
     
    if(object$stderrors) {
        get_std_errs <- sqrt(diag(object$covar_components$topleft))
        ci_alpha <- qnorm((1-coverage)/2, lower.tail = FALSE)
        
        if(object$family$family[1] %in% c("zipoisson","zinegative.binomial")) {
            sel_zeroinfl <- seq(1, length = nrow(summary_output$betas), by = ncol(summary_output$betas)+1)
            zeroinfl_prob_resultstab <- data.frame(
                estimate = object$zeroinfl_prob_intercept, 
                std_err = get_std_errs[sel_zeroinfl]
                )
            zeroinfl_prob_resultstab$z_value <- zeroinfl_prob_resultstab$estimate / zeroinfl_prob_resultstab$std_err
            zeroinfl_prob_resultstab$p_value <- 2*pnorm(abs(zeroinfl_prob_resultstab$z_value), lower.tail = FALSE)
            zeroinfl_prob_resultstab$lower <- zeroinfl_prob_resultstab$estimate - ci_alpha * zeroinfl_prob_resultstab$std_err
            zeroinfl_prob_resultstab$upper = zeroinfl_prob_resultstab$estimate + ci_alpha * zeroinfl_prob_resultstab$std_err
            zeroinfl_prob_resultstab <- round(zeroinfl_prob_resultstab, digits)
            colnames(zeroinfl_prob_resultstab) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)", "Lower CI", "Upper CI")
            
            summary_output$zeroinf_prob__intercept_resultstab <- zeroinfl_prob_resultstab
            get_std_errs <- get_std_errs[-sel_zeroinfl]
            }
        
        betas_resultstab <- data.frame(
            as.data.frame.table(t(object$betas)), 
            std_err = get_std_errs
            )
        betas_resultstab$z_value <- betas_resultstab[,3] / betas_resultstab$std_err
        betas_resultstab$p_value <- 2*pnorm(abs(betas_resultstab$z_value), lower.tail = FALSE)
        betas_resultstab$lower <- betas_resultstab[,3] - ci_alpha * get_std_errs
        betas_resultstab$upper <- betas_resultstab[,3] + ci_alpha * get_std_errs
        betas_resultstab[,-(1:2)] <- round(betas_resultstab[,-(1:2)], digits)
        
        colnames(betas_resultstab) <- c("Response", "Covariate", "Estimate", "Std. Error", "z value", "Pr(>|z|)", "Lower CI", "Upper CI")
          
        summary_output$betas_results <- betas_resultstab
        }
          
    class(summary_output) <- "summary.CBFM"
    return(summary_output) 
    }	


# # This version is heavily based on and inspired a lot more by the summary.gam function...acknowledgements go to Simon Wood for his mgcv package!
# summary2.CBFM <- function(object, coverage = 0.95, digits = max(3L, getOption("digits") - 3L), ...) {
#     num_spp <- nrow(object$betas)
#     num_basisfns <- nrow(object$basis_effects_mat)
#      
#     summary_output <- list(call = object$call, 
#                            betas = round(object$betas, digits), 
#                            basis_effects_mat = round(object$basis_effects_mat, digits))
#     if(object$family$family[1] %in% c("zipoisson","zinegative.binomial"))
#         summary_output$zeroinfl_prob_intercept <- round(object$zeroinfl_prob_intercept, digits)
#      
#     if(object$stderrors) {
#         ci_alpha <- qnorm((1-coverage)/2, lower.tail = FALSE)
#         get_std_errs <- sqrt(diag(object$covar_components$topleft))
#         
#         # Probabilities of zero-inflation 
#         if(object$family$family[1] %in% c("zipoisson","zinegative.binomial")) {
#             sel_zeroinfl <- seq(1, length = nrow(summary_output$betas), by = ncol(summary_output$betas)+1)
#             zeroinfl_prob_resultstab <- data.frame(
#                 estimate = object$zeroinfl_prob_intercept, 
#                 std_err = get_std_errs[sel_zeroinfl]
#                 )
#             zeroinfl_prob_resultstab$z_value <- zeroinfl_prob_resultstab$estimate / zeroinfl_prob_resultstab$std_err
#             zeroinfl_prob_resultstab$p_value <- 2*pnorm(abs(zeroinfl_prob_resultstab$z_value), lower.tail = FALSE)
#             zeroinfl_prob_resultstab$lower <- zeroinfl_prob_resultstab$estimate - ci_alpha * zeroinfl_prob_resultstab$std_err
#             zeroinfl_prob_resultstab$upper = zeroinfl_prob_resultstab$estimate + ci_alpha * zeroinfl_prob_resultstab$std_err
#             zeroinfl_prob_resultstab <- round(zeroinfl_prob_resultstab, digits)
#             colnames(zeroinfl_prob_resultstab) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)", "Lower CI", "Upper CI")
#             
#             summary_output$zeroinf_prob__intercept_resultstab <- zeroinfl_prob_resultstab
#             get_std_errs <- get_std_errs[-sel_zeroinfl]
#             }
#         
#         
#          ## A local pseudo-inverse function -- straight from summary.gam
#          pinv <- function(V, M, rank.tol = 1e-6) {
#              D <- eigen(V,symmetric=TRUE)
#              M1 <- length(D$values[D$values > rank.tol * D$values[1]])
#              if(M>M1) 
#                  M<-M1 # avoid problems with zero eigen-values
#              if(M+1 <= length(D$values)) 
#                  D$values[(M+1):length(D$values)]<-1
#              D$values<- 1/D$values
#              if(M+1 <= length(D$values)) D$values[(M+1):length(D$values)]<-0
#              res <- D$vectors %*% tcrossprod(D$values)  ##D$u%*%diag(D$d)%*%D$v
#              attr(res,"rank") <- M
#              res
#             }
#           
#         get_std_errs <- matrix(get_std_errs, nrow = num_spp, byrow = TRUE)
#         tmp_formula <- as.formula(paste("response", paste(as.character(object$formula_X),collapse="") ) )
#         spp_results_fn <- function(j) {
#             nullfit <- gam(tmp_formula, data = data.frame(response = object$y[,j], object$data), fit = TRUE, control = list(maxit = 1))
#             
#             # Individual parametric coefficient p-values...        
#             if(sum(nullfit$nsdf) > 0) {
#                 pstart <- 1
#                 ind <- 1:nullfit$nsdf
#                 p.coeff <- object$betas[j, ind]
#                 p.se <- get_std_errs[j, ind]
#                 p.t <- p.coeff / p.se
#                 p_table <- cbind(p.coeff, p.se, p.t, 2*pnorm(abs(p.t),lower.tail = FALSE), 
#                                  p.coeff - ci_alpha * p.se, p.coeff + ci_alpha * p.se)   
#                 dimnames(p_table) <- list(names(p.coeff), c("Estimate", "Std. Error", "z value", "Pr(>|z|)", "Lower CI", "Upper CI"))
#                 p_table <- as.data.frame(p_table)
#                 }    
#             if(sum(nullfit$nsdf) == 0) {
#                 p_table <- NULL
#                 }    
#             
#             
#             # Parametric terms, with factors treated whole... 
#             pterms <- if(is.list(nullfit$pterms)) nullfit$pterms else list(object$pterms)
#             if(!is.list(nullfit$assign)) object$assign <- list(object$assign)
#             npt <- length(unlist(lapply(pterms, attr, "term.labels")))
#             if(npt > 0)  
#                 pTerms.df <- pTerms.chi.sq <- pTerms.pv <- array(0, npt)
#             term.labels <- rep("",0)
#             k <- 0 ## total term counter        
#             
#             for(k1 in 1:length(pterms)) {
#                 tlj <- attr(pterms[[k1]], "term.labels") 
#                 nt <- length(tlj)
#                 if (k1 > 1 && nt > 0) 
#                     tlj <- paste(tlj, k1-1, sep = ".")
#                 term.labels <- c(term.labels, tlj)
#                 
#                 if(nt>0) { # individual parametric terms
#                     np <- length(nullfit$assign[[k1]])
#                     ind <- pstart[k1] - 1 + 1:np 
#                     Vb <- grep(paste0("response",j,"$"), rownames(object$covar_components$topleft))
#                     Vb <- object$covar_components$topleft[covmat_j, covmat_j, drop = FALSE]
#                     Vb <- covmat[ind, ind, drop = FALSE] 
#                     bp <- array(object$betas[j, ind], np)
#     
#                     for(i in 1:nt) { 
#                         k <- k + 1
#                         ind <- object$assign[[k1]] == i
#                         b <- bp[ind]
#                         V <- Vb[ind,ind]
#                         ## pseudo-inverse needed in case of truncation of parametric space 
#                         if(length(b) == 1) { 
#                             V <- 1/V 
#                             pTerms.df[k] <- nb <- 1      
#                             pTerms.chi.sq[k] <- V * b * b
#                             } 
#                         else {
#                             V <- pinv(V, length(b), rank.tol = 0.5*.Machine$double.eps)
#                             pTerms.df[k] <- nb <- attr(V,"rank")      
#                             pTerms.chi.sq[k] <- crossprod(b, V) %*% b
#                             }
#                         pTerms.pv[k] <- pchisq(pTerms.chi.sq[k], df = nb, lower.tail = FALSE)
#                         } ## for (i in 1:nt)
#                     }  ## if (nt>0)
#                 }
#             
#             if(npt) {
#                 attr(pTerms.pv,"names") <- term.labels
#                 pTerms_table <- cbind(pTerms.df, pTerms.chi.sq, pTerms.pv)   
#                 dimnames(pTerms_table) <- list(term.labels, c("df", "Chi.sq", "p-value"))
#                 } 
#             if(!npt) { 
#                 pTerms_table <- NULL
#                 }
#             
#             
#             # Smooth terms...
#             m <- length(nullfit$smooth) # number of smooth terms
#             df <- edf1 <- edf <- s.pv <- chi.sq <- array(0, m)
#             if (m > 0) { # form test statistics for each smooth
#                 X <- nullfit$R # NEED TO CHANGE THIS! 
#                 
#                 ii <- 0            
#                 for (i in 1:m) { ## loop through smooths
#                     start <- nullfit$smooth[[i]]$first.para
#                     stop <- nullfit$smooth[[i]]$last.para
# 
#                     V <- grep(paste0("response",j,"$"), rownames(object$covar_components$topleft))
#                     V <- object$covar_components$topleft[covmat_j, covmat_j, drop = FALSE]
#                     V <- V[start:stop, start:stop, drop = FALSE]
#                       
#                     p <- object$betas[j, start:stop]  # params for smooth
#                     edf1i <- edfi <- sum(nullfit$edf[start:stop]) # edf for this smooth
#                     # Extract alternative edf estimate for this smooth, if possible...  
#                     if(!is.null(object$edf1)) 
#                         edf1i <-  sum(nullfit$edf1[start:stop])
#                     Xt <- X[,start:stop, drop = FALSE]  
#                     fx <- if (inherits(object$smooth[[i]],"tensor.smooth")&& !is.null(object$smooth[[i]]$fx)) 
#                         all(object$smooth[[i]]$fx) else object$smooth[[i]]$fixed
#                     if(!fx&&object$smooth[[i]]$null.space.dim == 0 && !is.null(object$R)) { ## random effect or fully penalized term
#                         res <- if (re.test) reTest(object,i) else NULL
#                         } 
#                     else { ## Inverted Nychka interval statistics
#                         if (est.disp) rdf <- residual.df else rdf <- -1
#                         res <- testStat(p,Xt,V,min(ncol(Xt),edf1i),type=0,res.df = rdf)
#                     }
#                     
#                     if(!is.null(res)) {
#                         ii <- ii + 1
#                         df[ii] <- res$rank
#                         chi.sq[ii] <- res$stat
#                         s.pv[ii] <- res$pval 
#                         edf1[ii] <- edf1i 
#                         edf[ii] <- edfi 
#                         names(chi.sq)[ii]<- object$smooth[[i]]$label
#                         }
#                     }
#                 
#                 if(ii==0) 
#                     df <- edf1 <- edf <- s.pv <- chi.sq <- array(0, 0) 
#                 else {
#                     df <- df[1:ii]
#                     chi.sq <- chi.sq[1:ii]
#                     edf1 <- edf1[1:ii]
#                     edf <- edf[1:ii]
#                     s.pv <- s.pv[1:ii]
#                     }
#                 s_table <- cbind(edf, df, chi.sq, s.pv)      
#                 dimnames(s_table) <- list(names(chi.sq), c("edf", "Ref.df", "Chi.sq", "p-value"))
#                 }
#                     
#                     
#             out <- list(p_table = p_table, pTerms_table = pTerms_table, s_table = s_table)
#             return(out)
#             }
#         }
#     
#     class(summary_output) <- "summary.CBFM"
#     return(summary_output) 
#     }	


## Old stuff 
        # num_spp <- nrow(object$betas)
        # cov_names <- colnames(object$betas)[-1]
        # betas_resultstab$std_err_plusintercept <- betas_resultstab$std_err
        # betas_resultstab$lower_plusintercept <- betas_resultstab$lower
        # betas_resultstab$upper_plusintercept <- betas_resultstab$upper
        # for(k0 in 1:length(cov_names)) {          
        #     sel_cov <- c(1,grep(cov_names[k0],colnames(object$betas)))
        #     tmpX <- matrix(0, nrow = 1, ncol = ncol(object$betas))
        #     tmpX[,sel_cov] <- 1
        #     tmpB_space <- tmpB_time <- tmpB_spacetime <- NULL
        #     if(object$num_B_space > 0)
        #         tmpB_space <- Matrix(0, nrow = 1, ncol = object$num_B_space, sparse = TRUE)
        #     if(object$num_B_time > 0)
        #         tmpB_time <- Matrix(0, nrow = 1, ncol = object$num_B_time, sparse = TRUE)
        #     if(object$num_B_spacetime > 0)
        #         tmpB_spacetime <- Matrix(0, nrow = 1, ncol = object$num_B_spacetime, sparse = TRUE)
        #             
        #     tmp_pred <- suppressWarnings(predict.CBFM(object, newdata = object$data, manualX = tmpX, 
        #                                               new_B_space = tmpB_space, new_B_time = tmpB_time, new_B_spacetime = tmpB_spacetime,
        #                                               type = "link", se_fit = TRUE, coverage = coverage))
        #        
        #     find_indices <- grep(paste0(cov_names[k0],"$"), betas_resultstab$Var2)
        #     find_indices2 <- grep(paste0(cov_names[k0],"$"), colnames(object$betas))
        #     betas_resultstab$std_err_plusintercept[find_indices] <- round(tmp_pred$stderr, digits)
        #     betas_resultstab$lower_plusintercept[find_indices] <-  round(c(object$betas[,find_indices2]) - ci_alpha * tmp_pred$stderr, digits)
        #     betas_resultstab$upper_plusintercept[find_indices] <- round(c(object$betas[,find_indices2]) + ci_alpha * tmp_pred$stderr, digits)
        #     rm(tmpX, tmpB_space, sel_cov)
        #     }
        #   
