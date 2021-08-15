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
            as.data.frame.table(object$betas), 
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

