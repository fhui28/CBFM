#' @title Print a (hurdle) CBFM object
#' 
#' @description 
#' `r lifecycle::badge("stable")
#'
#' The default print method for a \code{CBFM} or \code{CBFM_hurdle} object.
#' 
#' @param x An object of class \code{CBFM} or \code{CBFM_hurdle}.
#' @param ... Not used.
#'
#' @details 
#' Print out details such as the function call, assumed family/response distribution, number of observational units and species, response-environment relationship fitted as given by the formula, and which sets of basis functions are used.`
#'
#' For a hurdle CBFM, details are provided for both of the component CBFMs separately. 
#' 
#' 
#' @author Francis K.C. Hui <fhui28@gmail.com>, Chris Haak
#'
#' @seealso [CBFM()] for fitting CBFMs and [makeahurdle()] for forming a hurdle CBFM based on two component CBFMs (a presence-absence component and a zero-truncated count component).
#'
#' @examples
#' \donttest{
#' # Please see the examples in the help file for the CBFM and makeahurdle function.s
#' }
#' 
#' @aliases print.CBFM print.CBFM_hurdle
#' @method print CBFM 
#' @export
#' @md

print.CBFM <- function(x, ...) {
        message("Call:")
        print(x$call)
        message()

        if(!(x$family$family[1] %in% c("ztpoisson")))
                message("Family: ", x$family$family[1], "\nNo. of units: ", nrow(x$fitted), "\nNo. of responses: ", ncol(x$fitted)) 
        if(x$family$family[1] %in% c("ztpoisson","ztnegative.binomial"))
                message("Family: ", x$family$family[1], "\nNo. of units: ", nrow(x$fitted), " (note zero counts in the response matrix are ignored in the model) \nNo. of responses: ", ncol(x$fitted)) 
        message("Responses-environment relationship fitted: ", x$formula) 
        message("Number of columns in model matrix induced by formula: ", ncol(x$betas))
        if(x$family$family[1] %in% c("zipoisson","zinegative.binomial")) {
             message("Responses-environment relationship fitted for modeling the probability of zero-inflation: ", x$ziformula) 
             message("Number of columns in model matrix induced by ziformula: ", ncol(x$zibetas))
             }
        message("Were standard errors calculated? (TRUE/FALSE): ", x$stderrors) 
     
        B_names <- c("B_space ", "B_time ", "B_spacetime ")
        message("\nBasis functions included: ", B_names[x$which_B_used==1])
        message("Total number of basis functions: ", x$num_B)
     
        if(sum(x$which_custom_Sigma_used) > 0)
                message("One or more custom basis function covariance matrices Sigma supplied.")
        
        if(sum(x$which_custom_G_used) > 0)
                message("One or more custom between species correlation matrices G supplied.")

        # if(x$which_B_used[1]) 
        #         message("Spatial component:", "\n\tRank of baseline between-response correlation matrix, G: ", ncol(x$Loading_G_space), "\n\tRank of basis function covariance matrix, Sigma: ", ifelse(x$which_custom_Sigma_used[1], NA, ncol(x$Loading_Sigma_space))) 
        # if(x$which_B_used[2])
        #         message("Temporal component:", "\n\tRank of baseline between-response correlation matrix, G: ", ncol(x$Loading_G_time), "\n\tRank of basis function covariance matrix, Sigma: ", ifelse(x$which_custom_Sigma_used[2], NA, ncol(x$Loading_Sigma_time))) 
        # if(x$which_B_used[3])
        #         message("Spatio-temporal component:", "\n\tRank of baseline between-response correlation matrix, G: ", ncol(x$Loading_G_spacetime), "\n\tRank of basis function covariance matrix, Sigma: ", ifelse(x$which_custom_Sigma_used[3], NA, ncol(x$Loading_Sigma_spacetime))) 
        
        if(any(x$which_nonzeromean_B > 0))
                message("\nOne or more of the normal distributions for the basis effect coefficients had non-zero mean vectors.")
          
        invisible(x)
        } 


#' @rdname print.CBFM
#' @method print CBFM_hurdle
#' @export
print.CBFM_hurdle <- function(x, ...) {
     message("A hurdle count CBFM with the following two components:")

     message("---------------------------")
     message("Presence-absence component:")
     message("---------------------------")
     print(x$pa_fit)
     
     message("\n\n---------------------------")
     message("Zero-truncated count component:")
     message("---------------------------")
     print(x$count_fit)

     invisible(x)
     }
