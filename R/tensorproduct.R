#' @title Tensor product
#' 
#' @description 
#' `r lifecycle::badge("stable")`
#' 
#' A simple function used to calculate the tensor product of two matrices. 
#' 
#' @param B1 A matrix.
#' @param B2 A matrix.
#' 
#' @details 
#' For the purposes of this package, this function may be used for constructing a matrix basis functions that is formed from the tensor product of two individual basis functions. A prime example of this is the construction of \code{B_spacetime} for use in the main function [CBFM()], where the spatio-temporal basis functions may be formed as a tensor product of spatial and temporal basis functions. 
#' 
#' The tensor product is formed by multiplying the first column of \code{B2} with \code{B1}, then the second column of \code{B2} with \code{B1}, and so on, working through all the columns of \code{B2}. 
#' 
#' @return A matrix resulting from the tensor product, with the number of columns equal to \code{ncol(B1)*ncol(B2)}. 
#' 
#' @details # Warnings
#' No attempt is made by the function to made the resulting matrix sparse, unless the suppled \code{B1} and \code{B2} are themselves sparse. Therefore, please be careful of how much memory the resulting object uses!
#' 
#' @author Francis K.C. Hui <fhui28@gmail.com>, Chris Haak
#' 
#' @seealso [CBFM()] for fitting CBFMs, with one example that using a tensor product formed set of basis functions.
#' 
#' @examples
#' \dontrun{
#' # Please see Example 3 in the help file for the main CBFM function 
#' }
#' 
#' @export
#' @md

tensorproduct <- function(B1, B2) {
     if(is.null(colnames(B1)))
          colnames(B1) <- paste0("B1_",1:ncol(B1))
     if(is.null(colnames(B2)))
          colnames(B2) <- paste0("B2_",1:ncol(B2))

     B <- lapply(1:ncol(B2), function(i)  (B2[,i] * B1)) 
     B <- do.call(cbind, B) 
     colnames(B) <-  sapply(1:ncol(B2), function(i) paste0(colnames(B1), "-", colnames(B2)[i])) 
     
     return(B)
     }
     
