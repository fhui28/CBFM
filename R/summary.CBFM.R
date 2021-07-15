#' @title Summary for a CBFM bit
#' 
#' @description 
#' `r lifecycle::badge("experimental")`
#'
#' Takes a fitted \code{CBFM} object and produces some useful summaries from it.
#'
#' @param object An object of class "CBFM".
#' @param coverage The coverage probability of the confidence intervals for the regression coefficients. Defaults to 0.95, which corresponds to 95% confidence intervals.
#' @param digits The number of significant figures to use when printing
#' @param ... Not used.
#'
#' @details 
#' TBC
#'
#' @return TBC
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
#' @importFrom stats qnorm 
#' @md


summary.CBFM <- function(object, coverage = 0.95, digits = max(3L, getOption("digits") - 3L), ...) {
     num_spp <- nrow(object$betas)
     num_basisfns <- nrow(object$basis_effects_mat)
     
     summary_output <- list(call = object$call, betas = round(object$betas, digits), basis_effects = round(object$basis_effects_mat, digits))
    
     if(object$stderrors) {          
          get_std_errs <- sqrt(diag(object$covar_components$topleft))
          ci_alpha <- qnorm((1-coverage)/2, lower.tail = FALSE)
          betas_resultstab <- data.frame(as.data.frame.table(round(object$betas, digits)), 
               std_err = round(get_std_errs, digits),
               lower = round(c(object$betas) - ci_alpha * get_std_errs, digits), 
               upper = round(c(object$betas) + ci_alpha * get_std_errs, digits)
               )

          num_spp <- nrow(object$betas)
          cov_names <- colnames(object$betas)[-1]
          betas_resultstab$std_err_plusintercept <- betas_resultstab$std_err
          betas_resultstab$lower_plusintercept <- betas_resultstab$lower
          betas_resultstab$upper_plusintercept <- betas_resultstab$upper
          for(k0 in 1:length(cov_names)) {          
               sel_cov <- c(1,grep(cov_names[k0],colnames(object$betas)))
               tmpX <- matrix(0, nrow = 1, ncol = ncol(object$betas))
               tmpX[,sel_cov] <- 1
               tmpB_space <- tmpB_time <- tmpB_spacetime <- NULL
               if(object$num_B_space > 0)
                    tmpB_space <- Matrix(0, nrow = 1, ncol = object$num_B_space, sparse = TRUE)
               if(object$num_B_time > 0)
                    tmpB_time <- Matrix(0, nrow = 1, ncol = object$num_B_time, sparse = TRUE)
               if(object$num_B_spacetime > 0)
                    tmpB_spacetime <- Matrix(0, nrow = 1, ncol = object$num_B_spacetime, sparse = TRUE)
                    
               tmp_pred <- suppressWarnings(predict.CBFM(object, newdata = object$data, manualX = tmpX, 
                    new_B_space = tmpB_space, new_B_time = tmpB_time, new_B_spacetime = tmpB_spacetime,
                    type = "link", se_fit = TRUE, coverage = coverage))
               
               find_indices <- grep(paste0(cov_names[k0],"$"), betas_resultstab$Var2)
               find_indices2 <- grep(paste0(cov_names[k0],"$"), colnames(object$betas))
               betas_resultstab$std_err_plusintercept[find_indices] <- round(tmp_pred$stderr, digits)
               betas_resultstab$lower_plusintercept[find_indices] <-  round(c(object$betas[,find_indices2]) - ci_alpha * tmp_pred$stderr, digits)
               betas_resultstab$upper_plusintercept[find_indices] <- round(c(object$betas[,find_indices2]) + ci_alpha * tmp_pred$stderr, digits)
               rm(tmpX, tmpB_space, sel_cov)
              }
          
          colnames(betas_resultstab)[1:3] <- c("Response", "Predictor", "Estimate")          
          summary_output$betas_results <- betas_resultstab
          #summary_output$basis_effects_results <- basis_effects_resultstab
          }
          
    
    class(summary_output) <- "summary.CBFM"
    summary_output 
    }	

