\name{reglogit}
\alias{reglogit}
\alias{regmlogit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Gibbs sampling for regularized logistic regression
}
\description{
Regularized (multinomial) logistic regression
by Gibbs sampling implementing subtly different 
MCMC schemes with varying efficiency depending on the data type 
(binary v. binomial, say) and the desired estimator (regularized maximum
likelihood, or Bayesian maximum a posteriori/posterior mean, etc.) through a 
unified interface.
}
\usage{
reglogit(T, y, X, N = NULL, flatten = FALSE, sigma = 1, nu = 1,
      kappa = 1, icept = TRUE, normalize = TRUE, zzero = TRUE, 
      powerprior = TRUE, kmax = 442, bstart = NULL, lt = NULL, 
      nup = list(a = 2, b = 0.1), save.latents = FALSE, verb = 100)
regmlogit(T, y, X, flatten = FALSE, sigma = 1, nu = 1, kappa = 1, 
      icept=TRUE, normalize = TRUE, zzero = TRUE, powerprior = TRUE, 
      kmax = 442, bstart = NULL, lt = NULL, nup = list(a=2, b=0.1),
      save.latents = FALSE, verb=100)
}
\arguments{
  \item{T}{
  a positive integer scalar specifying the number of MCMC rounds
}
  \item{y}{
  \code{reglogit} requires \code{logical} classification labels for Bernoulli 
  data, or countsfor Binomial data; for the latter, \code{N} must also be specified.
  \code{regmlogit} requires positive integer class labeels in \code{1:C} where
  \code{C} is the number of classes.
}
  \item{X}{
  a design \code{matrix} of predictors; can be a typical (dense) \code{matrix}
  or a sparse \code{\link{Matrix}} object.  When the design matrix
  is sparse (and is stored sparsely), this can produce a ~3x-faster execution via
  a more efficient update for the beta parameter.  But when it is not sparse
  (but is stored sparsely) the execution could be much slower
}
  \item{N}{
  an optional integer vector of total numbers of replicate trials 
  for each \code{X}-\code{y}, i.e., for Binomial data instead of Bernoulli
}
  \item{flatten}{
  a scalar \code{logical} that is only specified for Binomial data.  It
  indicates if pre-processing code should flatten the Binomial
  likelihood into a Bernoulli likelihood
}
  \item{sigma}{
  weights on the regression coefficients in the lasso penalty.  The
  default of \code{1} is sensible when \code{normalize = TRUE} since
  then the estimator for \code{beta} is equivariant under rescaling
}
  \item{nu}{
  a non-negative scalar indicating the initial value of the penalty
  parameter
}
  \item{kappa}{
  a positive scalar specifying the multiplicity;  \code{kappa = 1}
  provides samples from the Bayesian posterior distribution.  Larger
  values of \code{kappa} facilitates a simulated annealing approach
  to obtaining a regularized point estimator
}
  \item{icept}{
  a scalar \code{logical} indicating if an (implicit) intercept should
  be included in the model
}
  \item{normalize}{
   a scalar logical which, if \code{TRUE}, causes each variable is standardized
    to have unit L2-norm, otherwise it is left alone 
}
  \item{zzero}{
  a scalar \code{logical} indicating if the latent \code{z} variables to be
  sampled.  Therefore this indicator specifies if the cdf
  representation (\code{zzero = FALSE}) or pdf representation
  (otherwise) should be used
}
  \item{powerprior}{
  a scalar \code{logical} indicating if the prior should be powered up
  with multiplicity parameter \code{kappa} as well as the likelihood
}
  \item{kmax}{
  a positive integer indicating the number replacing infinity in the
  sum for mixing density in the generative expression for
  \code{lambda}
}
  \item{bstart}{
  an optional vector of length \code{p = ncol(X)} specifying initial
  values for the regression coefficients \code{beta}.   Otherwise
  standard normal deviates are used
}
  \item{lt}{
  an optional vector of length \code{n = nrow(X)} of initial values
  for the \code{lambda} latent
  variables.  Otherwise a vector of ones is used.
}
  \item{nup}{
  prior parameters \code{=list(a, b)} for the inverse Gamma distribution
  prior for \code{nu}, or \code{NULL}, which causes \code{nu} to be fixed
}
\item{save.latents}{ a scalar \code{logical} indicating wether or not
  a trace of latent \code{z}, \code{lambda} and \code{omega} values should be saved
  for each iteration.  Specify \code{save.latents=TRUE} for very large \code{X}
  in order to reduce memory swapping on low-RAM machines }
  \item{verb}{
  A positive integer indicating the number of MCMC rounds after which
  a progress statement is printed.  Giving \code{verb = 0} causes no
  statements to be printed
}
}
\details{

  These are the main functions in the package.  They support an omnibus
  framework for simulation-based regularized logistic regression.  The
  default arguments invoke a Gibbs sampling algorithm to sample from the
  posterior distribution of a logistic regression model with
  lasso-type (double-exponential) priors.  See the paper by Gramacy \&
  Polson (2012) for details.  Both cdf and pdf implementations are
  provided, which use slightly different latent variable
  representations, resulting in slightly different Gibbs samplers.  These
  methods extend the un-regularized methods of Holmes \& Held (2006)

  The \code{kappa} parameter facilitates simulated annealing (SA)
  implementations in order to help find the MAP, and other point
  estimators.  The actual SA algorithm is not provided in the package.
  However, it is easy to string calls to this function, using the
  outputs from one call as inputs to another, in order to establish a SA
  schedule for increasing kappa values.

  The \code{regmlogit} function is a wrapper around the Gibbs sampler
  inside \code{reglogit}, invoking \code{C-1} linked chains for \code{C}
  classes, extending the polychotomous regression scheme outlined by 
  Holmes \& Held (2006).  For an example with \code{regmlogit}, see
  \code{\link{predict.regmlogit}}
}
\value{
The output is a \code{list} object of type \code{"reglogit"} or 
\code{"regmlogit"} containing a subset of the following fields;
for \code{"regmlogit"} everyhing is expanded by one dimension into
an \code{array} or \code{matrix} as appropriate.

 \item{X }{ the input design \code{matrix}, possible adjusted by
 normalization or intercept }
 \item{y }{ the input response variable }
 \item{beta }{ a \code{matrix} of \code{T} sampled regression
 coefficients on the original input scale }
 \item{z }{ if \code{zzero = FALSE} a \code{matrix} of latent
 variables for the hierarchical cdf representation of the likelihood }
 \item{lambda }{ a \code{matrix} of latent variables for the
 hierarchical (cdf or pdf) representation of the likelihood }
 \item{lpost }{ a vector of log posterior probabilities of the parameters }
 \item{map }{ the \code{list} containing the maximum a' posterior
 parameters; \code{out$map$beta} is on the original scale of the data  }
 \item{kappa }{ the input multiplicity parameter }
 \item{omega }{ a \code{matrix} of latent variables for the
 regularization prior}

}
\references{
R.B. Gramacy, N.G. Polson. \dQuote{Simulation-based regularized
  logistic regression}. (2012) Bayesian Analysis, 7(3), p567-590; 
  arXiv:1005.3430; \url{http://arxiv.org/abs/1005.3430}

C. Holmes, K. Held (2006). \dQuote{Bayesian Auxilliary Variable Models for
   Binary and Multinomial Regression}. Bayesian Analysis, 1(1), p145-168. 
}
\author{
Robert B. Gramacy \email{rbg@vt.edu}
}

\seealso{
\code{\link{predict.reglogit}}, \code{\link{predict.regmlogit}}, 
\code{\link[monomvn]{blasso}} and \code{\link[monomvn]{regress}}
}
\examples{
## load in the pima indian data
data(pima)
X <- as.matrix(pima[,-9])
y <- as.numeric(pima[,9])

## pre-normalize to match the comparison in the paper
one <- rep(1, nrow(X))
normx <- sqrt(drop(one \%*\% (X^2)))
X <- scale(X, FALSE, normx)

## compare to the GLM fit
fit.logit <- glm(y~X, family=binomial(link="logit"))
bstart <- fit.logit$coef

## do the Gibbs sampling
T <- 300 ## set low for CRAN checks; increase to >= 1000 for better results
out6 <- reglogit(T, y, X, nu=6, nup=NULL, bstart=bstart, normalize=FALSE)

## plot the posterior distribution of the coefficients
burnin <- (1:(T/10)) 
boxplot(out6$beta[-burnin,], main="nu=6, kappa=1", ylab="posterior",
        xlab="coefficients", bty="n", names=c("mu", paste("b", 1:8, sep="")))
abline(h=0, lty=2)

## add in GLM fit and MAP with legend
points(bstart, col=2, pch=17)
points(out6$map$beta, pch=19, col=3)
legend("topright", c("MLE", "MAP"), col=2:3, pch=c(17,19))

## simple prediction
p6 <- predict(out6, XX=X)
## hit rate
mean(p6$c == y)

##
## for a polychotomous example, with prediction, 
## see ? predict.regmlogit
##

\dontrun{
## now with kappa=10
out10 <- reglogit(T, y, X, kappa=10, nu=6, nup=NULL, bstart=bstart, 
                            normalize=FALSE)

## plot the posterior distribution of the coefficients
par(mfrow=c(1,2))
boxplot(out6$beta[-burnin,], main="nu=6, kappa=1",  ylab="posterior",
        xlab="coefficients", bty="n",  names=c("mu", paste("b", 1:8, sep="")))
abline(h=0, lty=2) 
points(bstart, col=2, pch=17)
points(out6$map$beta, pch=19, col=3)
legend("topright", c("MLE", "MAP"), col=2:3, pch=c(17,19))
boxplot(out10$beta[-burnin,], main="nu=6, kappa=10",  ylab="posterior",
        xlab="coefficients", bty="n",  names=c("mu", paste("b", 1:8, sep="")))
abline(h=0, lty=2)
## add in GLM fit and MAP with legend
points(bstart, col=2, pch=17)
points(out10$map$beta, pch=19, col=3)
legend("topright", c("MLE", "MAP"), col=2:3, pch=c(17,19))
}

##
## now some binomial data
##

\dontrun{
## synthetic data generation
library(boot)
N <- rep(20, 100)
beta <- c(2, -3, 2, -4, 0, 0, 0, 0, 0)
X <- matrix(runif(length(N)*length(beta)), ncol=length(beta))
eta <- drop(1 + X \%*\% beta)
p <- inv.logit(eta)
y <- rbinom(length(N), N, p)

## run the Gibbs sampler for the logit -- uses the fast Binomial
## version; for a comparison, try flatten=FALSE
out <- reglogit(T, y, X, N)

## plot the posterior distribution of the coefficients
boxplot(out$beta[-burnin,], main="binomial data",  ylab="posterior", 
       xlab="coefficients", bty="n",
       names=c("mu", paste("b", 1:ncol(X), sep="")))
abline(h=0, lty=2)

## add in GLM fit, the MAP fit, the truth, and a legend
fit.logit <- glm(y/N~X, family=binomial(link="logit"), weights=N)
points(fit.logit$coef, col=2, pch=17)
points(c(1, beta), col=4, pch=16)
points(out$map$beta, pch=19, col=3)
legend("topright", c("MLE", "MAP", "truth"), col=2:4, pch=c(17,19,16))

## also try specifying a larger kappa value to pin down the MAP
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
\keyword{ classif }% __ONLY ONE__ keyword per line
