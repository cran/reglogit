\name{predict.reglogit}
\alias{predict.reglogit}
\alias{predict.regmlogit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Prediction for regularized (polychotomous) logistic regression models
}
\description{
  Sampling from the posterior predictive distribution of a regularized (multinomial)
  logistic regression fit, including entropy information for variability assessment
}
\usage{
\method{predict}{reglogit}(object, XX, burnin = round(0.1 * nrow(object$beta)), ...)
\method{predict}{regmlogit}(object, XX, burnin = round(0.1 * dim(object$beta)[1]), ...) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    a \code{"reglogit"}-class object or a \code{"regmlogit"}-class object, depending
    on whether binary or polychotomous methods were used for fitting
}
  \item{XX}{
  a \code{matrix} of predictive locations where \code{ncol(XX) == object$ncol(XX)}.
}
  \item{burnin}{
  a scalar positive \code{integer} indicate the number of samples of \code{object$beta}
  to discard as burn-in;  the default is 10\% of the number of samples
}
\item{...}{ For compatibility with generic \code{\link{predict}} method; not used}
}
\details{
  Applies the logit transformation (\code{reglogit}) or multinomial logit (\code{regmlogit})
  to convert samples of the linear predictor at \code{XX} into a samples from a predictive
  posterior probability distribution.  The raw probabilties, averages (posterior means), 
  entropies, and posterior mean casses (arg-max of the average probabilities) are returned.
}
\value{
  The output is a \code{list} with components explained below. For
  \code{predict.regmlogit} everyhing (except entropy) is expanded by one
  dimension into an \code{array} or \code{matrix} as appropriate.

  \item{p }{ a \code{nrow(XX) x (T-burnin)} sized \code{matrix} of
  probabilities (of class 1) from the posterior predictive distribution. }

  \item{mp }{ a vector of average probablities calculated over the rows of \code{p}}

  \item{pc }{ class labels formed by rouding (or arg max for \code{predict.regmlogit}) 
  the values in \code{mp}}

  \item{ent }{ The posterior mean entropy given the probabilities in \code{mp} }
  }

\references{
R.B. Gramacy, N.G. Polson. \dQuote{Simulation-based regularized
  logistic regression}. (2012) Bayesian Analysis, 7(3), p567-590; 
  arXiv:1005.3430; \url{http://arxiv.org/abs/1005.3430}

C. Holmes, K. Held (2006). \dQuote{Bayesian Auxilliary Variable Models for
   Binary and Multinomial Regression}. Bayesian Analysis, 1(1), p145-168. 
}

\author{
Robert B. Gramacy \email{rbgramacy@chicagobooth.edu}
}

\seealso{
  \code{\link{reglogit}} and \code{\link{regmlogit}}
}
\examples{
## see reglogit for a full example of binary classifiction complete with
## sampling from the posterior predictive distribution.  

## the example here is for polychotomous classification and prediction

\dontrun{
library(plgp)
x <- seq(-2, 2, length=40)
X <- expand.grid(x, x)
C <- exp2d.C(X)
xx <- seq(-2, 2, length=100)
XX <- expand.grid(xx, xx)
CC <- exp2d.C(XX)

## build cubically-expanded design matrix (with interactions)
Xe <- cbind(X, X[,1]^2, X[,2]^2, X[,1]*X[,2],
            X[,1]^3, X[,2]^3, X[,1]^2*X[,2], X[,2]^2*X[,1],
            (X[,1]*X[,2])^2)

## perform MCMC
T <- 1000
out <- regmlogit(T, C, Xe, nu=6, normalize=TRUE)

## create predictive (cubically-expanded) design matrix
XX <- as.matrix(XX)
XXe <- cbind(XX, XX[,1]^2, XX[,2]^2, XX[,1]*XX[,2],
             XX[,1]^3, XX[,2]^3, XX[,1]^2*XX[,2], XX[,2]^2*XX[,1],
             (XX[,1]*XX[,2])^2)

## predict class labels
p <- predict(out, XXe)

## make an image of the predictive surface
cols <- c(gray(0.85), gray(0.625), gray(0.4))
par(mfrow=c(1,3))
image(xx, xx, matrix(CC, ncol=length(xx)), col=cols, main="truth")
image(xx, xx, matrix(p$c, ncol=length(xx)), col=cols, main="predicted")
image(xx, xx, matrix(p$ent, ncol=length(xx)), col=heat.colors(128),
      main="entropy")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ classif }% __ONLY ONE__ keyword per line
