#' @title Family object for zero-truncated Poisson distribution 
#' 
#' @description 
#' `r lifecycle::badge("stable")`
#' 
#' A specified family object for the zero-truncated Poisson distribution, using the log link function for the mean of the Poisson distribution component.
#'
#' @details 
#' This family object was created specifically for fitting CBFMs. However it follows heavily from the [gamlss.tr::trun()] so big shout outs to the Mikis and the maintainers of that package.
#' 
#' @return An object of class "gamlss.family" and "family"
#' 
#' @author Francis K.C. Hui <fhui28@gmail.com>, Chris Haak
#' 
#' @export 
#' 
#' @importFrom gamlss.dist PO
#' @importFrom gamlss.tr trun
#' @md

ztpoisson <- function() {
    maketrunfam <- trun(par = 0, family = "PO", type = "left")()
    maketrunfam$family <- "ztpoisson"
    maketrunfam$link <- "log"
    maketrunfam
    }

