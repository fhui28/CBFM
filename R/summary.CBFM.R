#' @title Summary for a CBFM bit
#' 
#' @description 
#' `r lifecycle::badge("experimental")`
#'
#' Takes a fitted \code{CBFM} object and produces some useful summaries from it.
#'
#' @param object An object of class "CBFM".
#' @param coverage The coverage probability for any confidence intervals of the regression coefficients which are calculated. Defaults to 0.95, which corresponds to 95% confidence intervals.
#' @param digits The number of significant figures to use when printing.
#' @param ncores To speed up calculation of the standard error estimates and summary tables, parallelization can be performed, in which case this argument can be used to supply the number of cores to use in the parallelization. Defaults to \code{detectCores()-1}.
#' @param ... Not used.
#'
#' @details 
#' Currently, the function returns estimated species-specific regression coefficients, and estimated vector of species-specific probabilities of zero-inflation (on the logit scale), for distributions which require one. If the \code{object$stderrors == TRUE}, then summary tables are also produced containing standard errors, hypothesis tests for parametric coefficients and smoothing terms, plus Wald confidence intervals (when possible). The set up of these summary tables, in particular the way the hypothesis tests for smoothing terms are performed, is adapted from that of [mgcv::summary.gam()]; we refer the reader to their help file as well as Wood (2017) for more details on these tests.
#' 
#' As discussed in [CBFM()], all tests are based on the Bayesian posterior covariance matrix of the coefficients. Please note all P-values are computed without considering uncertainty in the smoothing parameter estimates.
#' 
#'  
#' @return An object of class "summary.CBFM" which includes the following components, not necessarily in the order below (and as appropriate):
#' \describe{
#' \item{call: }{The matched function call of \code{object}.}

#' \item{betas: }{The estimated matrix of species-specific regression coefficients corresponding to the model matrix created, rounded.}

#' \item{basis_effects_mat: }{The estimated matrix of species-specific regression coefficients corresponding to the combined matrix of basis functions. }

#' \item{zeroinfl_prob_intercept: }{The estimated vector of species-specific probabilities of zero-inflation, for distributions which require one, rounded. *Note this is presented on the logit scale*, that is the model returns \eqn{log(\pi_j/(1-\pi_j))} where \eqn{\pi_j} is the probability of zero-inflation. This is the same as the intercept term of a logistic regression model for the probabilities of zero-inflation, hence the name. } 

#' \item{zeroinf_prob_intercept_resultstab: }{If the \code{object$stderrors == TRUE}, then a data frame containing the point estimates, standard errors, corresponding Wald statistics and P-values, and the lower and upper limit of Wald confidence intervals for the probabilities of zero-inflation (if appropriate). Please note that the Wald-test is a test of whether the intercept is statistically different from i.e., whether the probability of zero-inflation is statistically different from 0.5. This may not be that useful in practice. }

#' \item{summary_tables: }{If the \code{object$stderrors == TRUE}, then a list with length equal to the number of species i.e., \code{ncol(object$y)}. Each element in the list may contain up to three summary tables: 
#' 1) \code{parametric_coefs}, which is a summary table corresponding to (any) strictly parametric coefficients included in the model; 
#' 2) \code{anova_terms}, which is a summary table for (any) parametric terms included in the model. The difference between \code{anova_terms} and \code{parametric_coefs} matters when it comes to categorical predictors (most commonly), with \code{parametric_coefs} returning summaries for individual coefficients involved in parameterizing the factor, and \code{anova_terms} returning a single summary for the entire factor. 
#' 3) \code{smooth_terms}, which is s summary table for (any) smoothing terms included in the model. 
#'  
#' Each summary table itself may consists of (and as appropriate) the estimated coefficient, the estimated standard error, the estimated degrees of freedom and rank, the value of the test statistic, the resulting P-values, and the lower and upper limit of Wald confidence intervals.}
#' }
#'
#' @details # Warning
#' Note that if the model matrix created includes smoothing terms as in a generalized additive model or GAM, the Wald P-values and corresponding Wald confidence intervals are approximate and neglect corresponding smoothing parameter uncertainty; please see [mgcv::summary.gam()] for more details. 
#' 
#' The current summary function is pretty basic (apologies!), and in the future we hope to add some useful information to the summary output. 
#' 
#' Testing of the species-specific coefficients associated with spatial and/or temporal coefficients is *not* performed.
#' 
#' @author Francis K.C. Hui <fhui28@gmail.com>, Chris Haak
#' 
#' @references
#' Wood, S. N. (2017). Generalized additive models: An introduction with R. CRC press.
#' 
#' @seealso [CBFM()] for fitting CBFMs and [coef.CBFM()] for extracting regression coefficients related to the covariates
#' 
#' @examples
#' \donttest{
#' # Please see examples in the help file for the main CBFM function 
#' }
#'
#' @export
#' @import foreach  
#' @importFrom stats pchisq pf qnorm pnorm plogis
#' @importFrom mgcv model.matrix.gam
#' @md

# This version is heavily based on and inspired a lot more by the summary.gam function...acknowledgements go to Simon Wood for his mgcv package.
# Interesting note that in simulations, most of the time, the results from this are not too far from (but very slightly more conservative than) just ad-hoc applying summary.gam to the last iteration of the PQL estimation algorithm in CBFM! 
summary.CBFM <- function(object, coverage = 0.95, digits = max(3L, getOption("digits") - 3L), ncores = NULL, ...) {
    if(is.null(ncores))
        registerDoParallel(cores = detectCores()-1)
    if(!is.null(ncores))
        registerDoParallel(cores = ncores)
    
    num_spp <- nrow(object$betas)
    num_basisfns <- nrow(object$basis_effects_mat)

    summary_output <- list(call = object$call,
                           betas = round(object$betas, digits),
                           basis_effects_mat = round(object$basis_effects_mat, digits))
    if(object$family$family[1] %in% c("zipoisson","zinegative.binomial"))
        summary_output$zeroinfl_prob_intercept <- round(object$zeroinfl_prob_intercept, digits)

    if(object$stderrors) {
        ci_alpha <- qnorm((1-coverage)/2, lower.tail = FALSE)
        get_std_errs <- sqrt(diag(object$covar_components$topleft))

        # Probabilities of zero-inflation
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
            colnames(zeroinfl_prob_resultstab) <- c("Estimate", "Std. Error", "z value", "P-value", "Lower CI", "Upper CI")

            summary_output$zeroinf_prob_intercept_resultstab <- zeroinfl_prob_resultstab
            get_std_errs <- get_std_errs[-sel_zeroinfl]
            }


         ## A local pseudo-inverse function -- straight from summary.gam
         pinv <- function(V, M, rank.tol = 1e-6) {
             D <- eigen(V, symmetric = TRUE)
             M1 <- length(D$values[D$values > rank.tol * D$values[1]])
             if(M > M1)
                 M<-M1 # avoid problems with zero eigen-values
             if(M+1 <= length(D$values))
                 D$values[(M+1):length(D$values)]<-1
             D$values<- 1/D$values
             if(M+1 <= length(D$values)) D$values[(M+1):length(D$values)]<-0
             res <- D$vectors %*% tcrossprod(D$values)  ##D$u%*%diag(D$d)%*%D$v
             attr(res,"rank") <- M
             res
            }
         
         get_std_errs <- matrix(get_std_errs, nrow = num_spp, byrow = TRUE)
         tmp_formula <- as.formula(paste("response", paste(as.character(object$formula_X),collapse="") ) )
         spp_results_fn <- function(j) {
            nullfit <- gam(tmp_formula, data = data.frame(response = object$y[,j], object$data), fit = TRUE, control = list(maxit = 1))

            # Individual parametric coefficient p-values...
            if(sum(nullfit$nsdf) > 0) {
                pstart <- 1
                ind <- 1:nullfit$nsdf
                p.coeff <- object$betas[j, ind]
                p.se <- get_std_errs[j, ind]
                p.t <- p.coeff / p.se
                p_table <- cbind(p.coeff, p.se, p.t, 2*pnorm(abs(p.t),lower.tail = FALSE),
                                 p.coeff - ci_alpha * p.se, p.coeff + ci_alpha * p.se)
                dimnames(p_table) <- list(names(p.coeff), c("Estimate", "Std. Error", "z value", "P-value", "Lower CI", "Upper CI"))
                p_table <- as.data.frame(p_table)
                p_table <- round(p_table, digits)
                }
            if(sum(nullfit$nsdf) == 0) {
                p_table <- NULL
                }


            # Parametric terms, with factors treated whole...
            pterms <- if(is.list(nullfit$pterms)) nullfit$pterms else list(nullfit$pterms)
            if(!is.list(nullfit$assign)) nullfit$assign <- list(nullfit$assign)
            npt <- length(unlist(lapply(pterms, attr, "term.labels")))
            if(npt > 0)
                pTerms.df <- pTerms.chi.sq <- pTerms.pv <- array(0, npt)
            term.labels <- rep("",0)
            k <- 0 ## total term counter

            for(k1 in 1:length(pterms)) {
                tlj <- attr(pterms[[k1]], "term.labels")
                nt <- length(tlj)
                if (k1 > 1 && nt > 0)
                    tlj <- paste(tlj, k1-1, sep = ".")
                term.labels <- c(term.labels, tlj)

                if(nt>0) { # individual parametric terms
                    np <- length(nullfit$assign[[k1]])
                    ind <- pstart[k1] - 1 + 1:np
                    Vb <- grep(paste0(colnames(object$y)[j],"$"), rownames(object$covar_components$topleft))
                    Vb <- object$covar_components$topleft[Vb, Vb, drop = FALSE]
                    Vb <- Vb[ind, ind, drop = FALSE]
                    bp <- array(object$betas[j, ind], np)

                    for(i in 1:nt) {
                        k <- k + 1
                        ind <- nullfit$assign[[k1]] == i
                        b <- bp[ind]
                        V <- Vb[ind,ind]
                        ## pseudo-inverse needed in case of truncation of parametric space
                        if(length(b) == 1) {
                            V <- 1/V
                            pTerms.df[k] <- nb <- 1
                            pTerms.chi.sq[k] <- V * b * b
                            }
                        else {
                            V <- pinv(V, length(b), rank.tol = 0.5*.Machine$double.eps)
                            pTerms.df[k] <- nb <- attr(V,"rank")
                            pTerms.chi.sq[k] <- crossprod(b, V) %*% b
                            }
                        pTerms.pv[k] <- pchisq(pTerms.chi.sq[k], df = nb, lower.tail = FALSE)
                        } ## for (i in 1:nt)
                    }  ## if (nt>0)
                }

            if(npt) {
                attr(pTerms.pv,"names") <- term.labels
                pTerms_table <- cbind(pTerms.df, pTerms.chi.sq, pTerms.pv)
                dimnames(pTerms_table) <- list(term.labels, c("DF", "Chi-squared", "P-value"))
                pTerms_table <- as.data.frame(pTerms_table)
                pTerms_table <- round(pTerms_table, digits)
                }
            if(!npt) {
                pTerms_table <- NULL
                }


            # Smooth terms...
            m <- length(nullfit$smooth) # number of smooth terms
            df <- edf1 <- edf <- s.pv <- chi.sq <- array(0, m)
            s_table <- NULL
            if (m > 0) { # form test statistics for each smooth
                X <- model.matrix.gam(nullfit)
                ii <- 0
                for (i in 1:m) { ## loop through smooths
                    start <- nullfit$smooth[[i]]$first.para
                    stop <- nullfit$smooth[[i]]$last.para

                    V <- grep(paste0(colnames(object$y)[j],"$"), rownames(object$covar_components$topleft))
                    V <- object$covar_components$topleft[V, V, drop = FALSE]
                    V <- V[start:stop, start:stop, drop = FALSE]

                    p <- object$betas[j, start:stop]  # params for smooth
                    edf1i <- edfi <- sum(object$edf[start:stop,j]) # edf for this smooth
                    if(!is.null(object$edf1)) # Extract alternative edf estimate for this smooth, if possible
                        edf1i <-  sum(object$edf1[start:stop,j])
                    
                    Xt <- X[,start:stop, drop = FALSE]
                    #fx <- if(inherits(nullfit$smooth[[i]],"tensor.smooth")&& !is.null(nullfit$smooth[[i]]$fx))
                    #    all(nullfit$smooth[[i]]$fx) else nullfit$smooth[[i]]$fixed
                    #if(!fx && nullfit$smooth[[i]]$null.space.dim == 0 && !is.null(object$R)) { ## random effect or fully penalized term
                    #    res <- if (re.test) reTest(nullfit,i) else NULL
                    #    }
                    #else { ## Inverted Nychka interval statistics
                    res <- .testStat(p = p, X = Xt, V = V, rank = min(ncol(Xt), edf1i), type = 0, res.df = -1)
                
                    if(!is.null(res)) {
                        ii <- ii + 1
                        df[ii] <- res$rank
                        chi.sq[ii] <- res$stat
                        s.pv[ii] <- res$pval
                        edf1[ii] <- edf1i
                        edf[ii] <- edfi
                        names(chi.sq)[ii] <- nullfit$smooth[[i]]$label
                        }
                    }

                if(ii == 0)
                    s_table <- NULL
                if(ii > 0) {
                    df <- df[1:ii]
                    chi.sq <- chi.sq[1:ii]
                    edf1 <- edf1[1:ii]
                    edf <- edf[1:ii]
                    s.pv <- s.pv[1:ii]
                    
                    s_table <- cbind(edf, df, chi.sq, s.pv)
                    dimnames(s_table) <- list(names(chi.sq), c("EDF", "Rank", "Chi-squared", "P-value"))
                    s_table <- as.data.frame(s_table)
                    s_table <- round(s_table, digits)
                    }
                }

            out <- list(parametric_coefs = p_table, anova_terms = pTerms_table, smooth_terms = s_table)
            return(out)
            }
        
         spp_summaries <- foreach(j = 1:num_spp) %dopar% spp_results_fn(j = j)
         names(spp_summaries) <- colnames(object$y)
         
         summary_output$summary_tables <- spp_summaries
         }

    
    class(summary_output) <- "summary.CBFM"
    return(summary_output)
    }


# This function is take directly mgcv.R in the mgcv package. Full credit goes to Simon Wood for this.
.testStat <- function(p, X, V, rank = NULL, type = 0, res.df= -1) {
    qrx <- qr(X, tol = 0)
    R <- qr.R(qrx)
    V <- R %*% V[qrx$pivot, qrx$pivot, drop=FALSE] %*% t(R)
    V <- (V + t(V))/2
    ed <- eigen(V,symmetric=TRUE)
  
    ## remove possible ambiguity from statistic...
    siv <- sign(ed$vectors[1,])
    siv[siv==0] <- 1
    ed$vectors <- sweep(ed$vectors,2,siv,"*")

    k <- max(0,floor(rank)) 
    nu <- abs(rank - k) ## fractional part of supplied edf
    if (type==1) { ## round up is more than .05 above lower
        if (rank > k + .05||k==0) 
            k <- k + 1
        nu <- 0;rank <- k
        }

    if (nu>0) k1 <- k+1 else k1 <- k

    ## check that actual rank is not below supplied rank+1
    r.est <- sum(ed$values > max(ed$values)*.Machine$double.eps^.9)
    if (r.est < k1) { 
        k1 <- k <- r.est
        nu <- 0
        rank <- r.est 
        }

    ## Get the eigenvectors...
    # vec <- qr.qy(qrx,rbind(ed$vectors,matrix(0,nrow(X)-ncol(X),ncol(X))))
    vec <- ed$vectors
    if(k1<ncol(vec)) 
        vec <- vec[, 1:k1, drop = FALSE]

    ## deal with the fractional part of the pinv...
    if(nu>0&&k>0) {
        if(k>1) 
            vec[,1:(k-1)] <- t(t(vec[,1:(k-1)])/sqrt(ed$val[1:(k-1)]))
        b12 <- .5*nu*(1-nu)
        if(b12<0) 
            b12 <- 0
        b12 <- sqrt(b12)
        B <- matrix(c(1,b12,b12,nu),2,2)
        ev <- diag(ed$values[k:k1]^-.5,nrow=k1-k+1)
        B <- ev %*% B %*% ev
        eb <- eigen(B, symmetric = TRUE)
        rB <- eb$vectors%*%diag(sqrt(eb$values))%*%t(eb$vectors)
        vec1 <- vec
        vec1[,k:k1] <- t(rB%*%diag(c(-1,1))%*%t(vec[,k:k1]))
        vec[,k:k1] <- t(rB%*%t(vec[,k:k1]))
        } 
    else {
        vec1 <- vec <- if (k==0) t(t(vec)*sqrt(1/ed$val[1])) else t(t(vec)/sqrt(ed$val[1:k]))
        if (k==1) 
            rank <- 1
        }


    ## there is an ambiguity in the choise of test statistic, leading to slight differences in the p-value computation depending on which of 2 alternatives. is arbitrarily selected. Following allows both to be computed and p-values averaged (can't average test stat as dist then unknown) 
    d <- t(vec) %*% (R%*%p)
    d <- sum(d^2) 
    d1 <- t(vec1) %*% (R%*%p)
    d1 <- sum(d1^2)
    ##d <- d1 ## uncomment to avoid averaging

    rank1 <- rank ## rank for lower tail pval computation below

    ## note that for <1 edf then d is not weighted by EDF, and instead is 
    ## simply refered to a chi-squared 1

    if (nu>0) { ## mixture of chi^2 ref dist
        if (k1==1) 
            rank1 <- val <- 1 
        else { 
            val <- rep(1,k1) ##ed$val[1:k1]
            rp <- nu+1
            val[k] <- (rp + sqrt(rp*(2-rp)))/2
            val[k1] <- (rp - val[k])
            }
   
        if (res.df <= 0) 
            pval <- 0.5*(.psum.chisq(q = d, lb = val) + .psum.chisq(q = d1, lb = val))
        else {  ## (liu2(d,val) + liu2(d1,val))/2 else
            k0 <- max(1,round(res.df))
            pval <- 0.5*(.psum.chisq(q = 0, lb = c(val,-d/k0), df = c(rep(1,length(val)),k0)) + .psum.chisq(q = 0, lb = c(val,-d1/k0), df = c(rep(1,length(val)),k0)) ) 
            }
        }  
    else { pval <- 2 }
    
    ## integer case still needs computing, 
    ## OLD: also liu/pearson approx only good in upper tail. In lower tail, 2 moment approximation is better (Can check this by simply plotting the whole interesting range as a contour plot!)
    if (pval > 1) {
        if (res.df <= 0) 
            pval <- (pchisq(d,df=rank1,lower.tail=FALSE)+pchisq(d1,df=rank1,lower.tail=FALSE))/2 
        else
            pval <- (pf(d/rank1,rank1,res.df,lower.tail=FALSE)+pf(d1/rank1,rank1,res.df,lower.tail=FALSE))/2
        }
    
    list(stat = d, pval = min(1,pval), rank = rank)
    } 



# This function is adapted from mgcv.R in the mgcv package. Full credit goes to Simon Wood for this.
.psum.chisq <- function(q, lb, df = rep(1,length(lb)), nc = rep(0,length(lb)), sigz = 0, lower.tail = FALSE, 
                       tol = 2e-5, nlim = 100000, trace = FALSE) {
    ## compute Pr(q>\sum_j lb[j] X_j + sigz Z) where X_j ~ chisq(df[j],nc[j]), Z~N(0,1) and nc is
    ## a vector of non-centrality parameters. lb can be either sign. df should be integer. 
    p <- q
    r <- length(lb)
    if (length(df) == 1) 
      df <- rep(df,r)
    if (length(nc) == 1) 
      nc <- rep(nc,r)
    if (length(df)!=r||length(nc)!=r) 
      stop("lengths of lb, df and nc must match")
    df <- round(df)
    if (any(df<1)) 
      stop("df must be positive integers")
    if (all(lb == 0)) 
      stop("at least one element of lb must be non-zero")
    if (sigz<0) sigz <- 0
  
    for(i in 1:length(q)) {  
        p[i] <- if(all(nc==0)) .liu2(x = q[i], lambda = lb, h = df) else NA
        } 
    p
    } ##psum.chisq



# This function is take directly mgcv.R in the mgcv package. Full credit goes to Simon Wood for this.
.liu2 <- function(x, lambda, h = rep(1,length(lambda)), lower.tail = FALSE) {
    # Evaluate Pr[sum_i \lambda_i \chi^2_h_i < x] approximately.
    # Code adapted from CompQuadForm package of Pierre Lafaye de Micheaux and directly from....
    # H. Liu, Y. Tang, H.H. Zhang, A new chi-square approximation to the distribution of non-negative definite quadratic forms in non-central normal variables, Computational Statistics and Data Analysis, Volume 53, (2009), 853-856. Actually, this is just Pearson (1959) given that the chi^2 variables are central. 
    # Note that this can be rubbish in lower tail (e.g. lambda=c(1.2,.3), x = .15). But I won't think we ever consider this so we stick with it...
      
    if (length(h) != length(lambda)) 
        stop("lambda and h should have the same length!")
    
    ## Try Davies exact method in place of Liu et al/ Pearson approx.
    # res <- rep(NA, length(x))
    # for (i in 1:length(x)) {
    #         do_davies <-  try(davies(q = x[i], lambda = lambda, h = h, acc = 2e-5, lim = 100000), silent = TRUE)
    #         if(!inherits(do_davies, "try-error")) {
    #             if(do_davies$ifault == 2)
    #                 warning("Danger of round-off error.")
    #             res[i] <- min(do_davies$Qq, 1)
    #             }
    #         }            
    # if(all(!is.na(res)))
    #     return(res)
    
    ## Else use Liu approximation
    lh <- lambda*h
    muQ <- sum(lh)
  
    lh <- lh*lambda
    c2 <- sum(lh)
  
    lh <- lh*lambda
    c3 <- sum(lh)
  
    xpos <- x > 0
    res <- 1 + 0 * x
    if(sum(xpos)==0 || c2 <= 0) 
        return(res)

    s1 <- c3/c2^1.5
    s2 <- sum(lh*lambda)/c2^2

    sigQ <- sqrt(2*c2)

    t <- (x[xpos]-muQ)/sigQ

    if (s1^2>s2) {
        a <- 1/(s1-sqrt(s1^2-s2))
        delta <- s1*a^3-a^2
        l <- a^2-2*delta
        } 
    else {
        a <- 1/s1
        delta <- 0
        if (c3==0) return(res)
        l <- c2^3/c3^2
        }

    muX <- l+delta
    sigX <- sqrt(2)*a
    res[xpos] <- pchisq(t*sigX+muX, df = l, ncp = delta, lower.tail = lower.tail)
    res
    } ## liu2



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

# OLD_summary.CBFM <- function(object, coverage = 0.95, digits = max(3L, getOption("digits") - 3L), ...) {
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
#         get_std_errs <- sqrt(diag(object$covar_components$topleft))
#         ci_alpha <- qnorm((1-coverage)/2, lower.tail = FALSE)
#         
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
#             summary_output$zeroinf_prob__intercept_results <- zeroinfl_prob_resultstab
#             get_std_errs <- get_std_errs[-sel_zeroinfl]
#             }
#         
#         betas_resultstab <- data.frame(
#             as.data.frame.table(t(object$betas)), 
#             std_err = get_std_errs
#             )
#         betas_resultstab$z_value <- betas_resultstab[,3] / betas_resultstab$std_err
#         betas_resultstab$p_value <- 2*pnorm(abs(betas_resultstab$z_value), lower.tail = FALSE)
#         betas_resultstab$lower <- betas_resultstab[,3] - ci_alpha * get_std_errs
#         betas_resultstab$upper <- betas_resultstab[,3] + ci_alpha * get_std_errs
#         betas_resultstab[,-(1:2)] <- round(betas_resultstab[,-(1:2)], digits)
#         
#         colnames(betas_resultstab) <- c("Covariate", "Response", "Estimate", "Std. Error", "z value", "Pr(>|z|)", "Lower CI", "Upper CI")
#           
#         summary_output$betas_results <- betas_resultstab
#         }
#           
#     class(summary_output) <- "summary.CBFM"
#     return(summary_output) 
#     }	
# 
