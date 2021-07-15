#' @title Print a CBFM object
#' 
#' @description 
#' `r lifecycle::badge("stable")
#'
#' The default print method for a \code{CBFM} object.
#' 
#' @param x An object of class "CBFM".
#' @param ... Not used.
#'
#' @details 
#' Print out details such as the function call, assumed family/response distribution, number of observational units and species, response-environment relationship fitted as given by the formula, which sets of basis functions are used, and the ranks.`
#'
#' @author Francis K.C. Hui <fhui28@gmail.com>, Chris Haak
#'
#' @seealso [CBFM()] for fitting CBFMs and [fitted.values.CBFM()] for calculating fitted values from a CBFM fit.
#'
#' @export
#' @md

print.CBFM <- function(x, ...) {
     message("Call:")
     print(x$call)
     message()

     message("Family: ", x$family$family[1], "\nNo. of units: ", nrow(x$fitted), "\nNo. of responses: ", ncol(x$fitted)) 
     message("Responses-environment relationship fitted: ", x$formula_X) 
     message("Number of columns in model matrix induced by formula_X: ", ncol(x$betas))
     
     B_names <- c("B_space", "B_time", "B_spacetime")
     message("Basis functions included: ", B_names[x$which_B_used==1])
     message("Total number of basis functions: ", x$num_B)
     
     if(x$which_B_used[1])
          message("B_space:", "\n\tRank of baseline between-response correlation matrix, G: ", ncol(x$Loading_G_space), "\n\tRank of basis function covariance matrix, Sigma: ", ncol(x$Loading_Sigma_space)) 
     if(x$which_B_used[2])
          message("B_time:", "\n\tRank of baseline between-response correlation matrix, G: ", ncol(x$Loading_G_time), "\n\tRank of basis function covariance matrix, Sigma: ", ncol(x$Loading_Sigma_time)) 
     if(x$which_B_used[3])
          message("B_spacetime:", "\n\tRank of baseline between-response correlation matrix, G: ", ncol(x$Loading_G_spacetime), "\n\tRank of basis function covariance matrix, Sigma: ", ncol(x$Loading_Sigma_spacetime)) 
          
     invisible(x)
     }

