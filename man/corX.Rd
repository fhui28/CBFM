% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/corX.R
\name{corX}
\alias{corX}
\title{Construct (cross-)correlations due to measured covariates for a CBFM fit}
\usage{
corX(
  object,
  newdata = NULL,
  newdata2 = NULL,
  se_cor = FALSE,
  coverage = 0.95,
  ncores = NULL,
  num_sims = 500
)
}
\arguments{
\item{object}{An object of class \code{CBFM}.}

\item{newdata}{A data frame containing the values of the covariates at which correlations are to be calculated. If this is not provided, then correlations corresponding to the original data are returned. If \code{newdata} is provided then it should contain all the variables needed for constructing correlations, that is, it can construct a model matrix from this as \code{object$formula}.}

\item{newdata2}{A second data frame containing the values of the covariates at which cross-correlations are to be calculated. If this is supplied, then \code{newdata} must also be supplied, as the function assumes then the user desires calculation of cross-correlations.}

\item{se_cor}{Should uncertainty intervals also be produced for the correlations? Defaults to \code{FALSE}.}

\item{coverage}{The coverage probability of the uncertainty intervals for the correlations. Defaults to 0.95, which corresponds to 95\% uncertainty intervals.}

\item{ncores}{To speed up calculation of the uncertainty estimates, parallelization can be performed, in which case this argument can be used to supply the number of cores to use in the parallelization. Defaults to \code{detectCores()-1}.}

\item{num_sims}{The number of Monte-Carlo examples to simulate.}
}
\value{
If \code{se_cor = FALSE}, then a (cross-)correlation matrix is returned. If \code{se_cor = TRUE}, then a list with the following components is returned:
\item{correlation: }{A matrix of (cross-)correlation values.}
\item{lower: }{A matrix of the lower bound of the uncertainty intervals for the correlations}
\item{upper: }{A matrix of the upper bound of the uncertainty intervals for the correlations}
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#experimental}{\figure{lifecycle-experimental.svg}{options: alt='[Experimental]'}}}{\strong{[Experimental]}}

Takes a fitted \code{CBFM} object calculates the between-species correlation matrix due to the measured covariates, along with corresponding uncertainty intervals if desired (via a simulation-based approach; see Hui et al., 2023). Similar to \code{\link[=predict.CBFM]{predict.CBFM()}}, this correlation matrix can be calculated based on a different sets of covariates to those used to actually fit the model. Additionally, the user can supplied two sets of covariates (data frames), in which case the function calculates cross-correlations (between and within species) between these two sets of covariates.
}
\details{
This function is adapted from and behaves somewhat similarly to \code{\link[boral:get.enviro.cor]{boral::get.enviro.cor()}}, in calculating a between-species correlation matrix due to the measured covariates i.e., shared environmental response, along with corresponding uncertainty intervals if desired. Recall the general form of the mean regression model for the CBFM is given by

\deqn{g(\mu_{ij}) = \eta_{ij} = x_i^\top\beta_j + b_i^\top a_j,}

where \eqn{x_i} denotes a vector of predictors for unit \eqn{i} i.e., the \eqn{i}-th row from the created model matrix, \eqn{\beta_j} denotes the corresponding regression coefficients for species \eqn{j}.

The covariance and hence correlation between two species \eqn{j} and \eqn{j'} that can be attributed to the measured covariates is then based on examining the components \eqn{x_i^\top\beta_j} and \eqn{x_i^\top\beta_{j'}} across the observational units; see equation 4 in Pollock et al., (2014) for the precise formula, as well as Warton et al., (2015), Hui (2016), and Hui et al., (2023) among others. A point estimate for this residual correlation can be constructed, along with a corresponding uncertainty interval by simulation if desired. For the latter, species-specific regression coefficients \eqn{\beta_j} are sampled from their approximate large sample normal distribution (i.e., basically a Gaussian approximation to the posterior distribution of the parameters; see \code{\link[=CBFM]{CBFM()}} and the section on estimation and inference), which are then used to calculate the correlations. This sampling and calculation is then performed a large number of times (as governed by \code{num_sims}) after which uncertainty intervals can be constructed by taking sample quantiles.

Note because this function calculates correlations as based on component of the linear predictor \eqn{x_i^\top\beta_j}, then it can not be applied to  \code{CBFM_hurdle} object (which by construction contains two linear predictors, so the user has to decide which component of the hurdle model they are interested in). Analogously, for zero-inflated CBFMs this function currently only calculates the between-species correlation matrix due to the measured covariates in \code{object$formula} i.e., the mean of the count component. It is \emph{not} able to calculate correlations due to measured covariates in \code{object$ziformula} i.e., in modeling the probability of zero-inflation.

With the above definition of the between-species correlation, note the predictors on which the correlations are constructed need not be the same those used in fitting the original CBFM i.e., the \eqn{x_i}'s can be different to those of \code{object$data}. This is handled via the \code{newdata} argument. Additionally, it is possible to calculate within and between species \emph{cross-correlations} across two different sets of measured predictors. That is, correlations are calculated between \eqn{x_i^\top\beta_j} and \eqn{x_{i2}^\top\beta_{j'}}, where \eqn{x_i} and \eqn{x_{i2}} can be different sets of measured predictors. This is handled by supplying both \code{newdata} and \code{newdata2} arguments simultaneously. Cross-correlations may be useful, say, if the two sets of measurement predictors reflect two different sets of sampling units, and we are interested in how similar (or lack of) the species communities are in terms of their environmental response across these two sets (Ovakainen et al., 2017). Another example if is the same set of observational units are visited at two different time points, and we are interested in how similarity (or lack of) the environmental responses within and between species are between these two time points.

NOTE: A cross-correlation matrix is not going to be a standard correlation matrix in the sense of having the ones alone the diagonal. This is because even for the same species \eqn{j = j'}, the correlation is not guaranteed to be equal to one as the covariates being considered can be different.
}
\examples{
\dontrun{
library(autoFRK)
library(FRK)
library(MASS)
library(mvabund)
library(mvtnorm)
library(ROCR)
library(sp)
library(geoR)
library(tidyverse)
library(corrplot)

##------------------------------
## **Example 1: Fitting a CBFM to spatial multivariate presence-absence data** 
## simulated from a spatial latent variable model
## Please note the data generation process (thus) differs from CBFM.
##------------------------------
set.seed(2021)
num_sites <- 1000 # 500 (units) sites for training set + 500 sites for external calculation
num_spp <- 50 # Number of species
num_X <- 4 # Number of regression slopes

spp_slopes <- matrix(runif(num_spp * num_X, -1, 1), nrow = num_spp)
spp_intercepts <- runif(num_spp, -2, 0)

# Simulate spatial coordinates and environmental covariate components
# We will use this information in later examples as well
xy <- data.frame(x = runif(num_sites, 0, 5), y = runif(num_sites, 0, 5))
X <- mvtnorm::rmvnorm(num_sites, mean = rep(0,4))
colnames(X) <- c("temp", "depth", "chla", "O2")
dat <- data.frame(xy, X)
mm <- model.matrix(~ temp + depth + chla + O2 - 1, data = dat) \%>\% 
scale \%>\% 
as.matrix

# Simulate latent variable component
# We will use this information in later examples as well
true_lvs <- grf(grid = cbind(xy$x, xy$y), nsim = 2, cov.model = "exponential",
cov.pars = c(1, 2))$data \%>\%
     as.matrix
spp_loadings <- matrix(runif(num_spp * 2, -1, 1), nrow = num_spp)
set.seed(NULL)

# Simulate spatial multivariate abundance data (presence-absence)
# We will use this information in later examples as well
eta <- tcrossprod(cbind(1,mm), cbind(spp_intercepts,spp_slopes)) + 
tcrossprod(true_lvs, spp_loadings)
simy <- matrix(rbinom(num_sites * num_spp, size = 1, 
prob = plogis(eta)), nrow = num_sites)

# Form training and test sets
dat_train <- dat[1:500,]
dat_external <- dat[501:1000,]
simy_train <- simy[1:500,]
rm(X, mm, spp_loadings, true_lvs, xy, simy, dat)



# Set up spatial basis functions for CBFM -- Most users will start here! 
# We will also use this basis functions in some later examples
num_basisfunctions <- 25 # Number of spatial basis functions to use
# Training set basis functions
train_basisfunctions <- mrts(dat_train[,c("x","y")], num_basisfunctions) \%>\% 
as.matrix \%>\%
{.[,-(1)]} # Remove the first intercept column

# Fit CBFM 
tic <- proc.time()
useformula <- ~ temp + depth + chla + O2
fitcbfm <- CBFM(y = simy_train, formula = useformula, data = dat_train, 
B_space = train_basisfunctions, family = binomial(), control = list(trace = 1))
toc <- proc.time()
toc - tic


# Calculate between-species correlations based on measured covariates in training data
getcor <- corX(fitcbfm)
corrplot(getcor, method = "square", type = "lower", order = "hclust")


# Calculate between-species correlations based on measured covariates in external (new) data
getcor <- corX(fitcbfm, newdata = dat_external)
corrplot(getcor, method = "square", type = "lower", order = "hclust")


# Calculate species cross-correlations between measured covariates in training and 
# external (new) data
# This may be useful, for example, if the training and external data were the same sites visited 
# at two points in time, and the user is interested in the similarity (or lack of) 
# of the environmental responses within and between species are between these two time points. 
# Note the resulting cross-correlation matrix is not strictly a correlation matrix 
# in the sense of having ones on the diagonals; all elements, including diagonals, 
# will lie between -1 and 1.  
getcrosscor <- corX(fitcbfm, newdata = dat_train, newdata2 = dat_external)
corrplot(getcrosscor, method = "square", type = "lower", is.corr = FALSE)
}

}
\references{
Hui, F. K. C. (2016). boral-Bayesian ordination and regression analysis of multivariate abundance data in R. Methods in Ecology and Evolution, 7, 744-750.

Hui, F. K. C., Warton, D. I., Foster, S. D., & Haak, C. R. (2023). Spatiotemporal joint species distribution modelling: A basis function approach. Methods in Ecology and Evolution, 14(8), 2150-2164.

Ovaskainen, O., Tikhonov, G., Norberg, A., Guillaume Blanchet, F., Duan, L., Dunson, D., and Abrego, N. (2017). How to make more out of community data? A conceptual framework and its implementation as models and software. Ecology letters, 20, 561-576.

Pollock, L. J., Tingley, R., Morris, W. K., Golding, N., O'Hara, R. B., Parris, K. M., Vesk, P. A., and McCarthy, M. A. (2014). Understanding co‐occurrence by modelling species simultaneously with a Joint Species Distribution Model (JSDM). Methods in Ecology and Evolution, 5, 397-406.

Warton, D. I., Blanchet, F. G., O'Hara, R. B., Ovaskainen, O., Taskinen, S., Walker, S. C., and Hui, F. K. C. (2015). So many variables: joint modeling in community ecology. Trends in Ecology and Evolution, 30, 766-779.
}
\seealso{
\code{\link[=CBFM]{CBFM()}} for fitting CBFMs, and \code{\link[=corB]{corB()}} for calculating residual between-species (cross-)correlations due to the basis functions.
}
\author{
Francis K.C. Hui \href{mailto:fhui28@gmail.com}{fhui28@gmail.com}, Chris Haak
}
