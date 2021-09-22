#' @title Family object for zero-truncated negative binomial distribution 
#' 
#' @description 
#' `r lifecycle::badge("stable")`
#' 
#' A specified family object for the zero-truncated negative binomial distribution, using the log link function for the negative binomial distribution component.
#'
#' @details 
#' This family object was created specifically for fitting CBFMs. 
#' 
#' This family object was created specifically for fitting CBFMs. However it follows heavily from the [gamlss.tr::trun()] so big shout outs to the Mikis and the maintainers of that package.
#' 
#' @return An object of class "gamlss.family" and "family"
#' 
#' @author Francis K.C. Hui <fhui28@gmail.com>, Chris Haak
#' 
#' @export 
#' 
#' @importFrom gamlss.tr trun
#' @md

ztnb2 <- function() {
    maketrunfam <- gamlss.tr::trun(par = 0, family = "NBI", type = "left")()
    maketrunfam$family <- "ztnegative.binomial"
    maketrunfam$link <- "log"
    maketrunfam
     }

