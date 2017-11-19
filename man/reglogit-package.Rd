\name{reglogit-package}
\alias{reglogit-package}
\docType{package}
\title{
Simulation-based Regularized Logistic Regression
}
\description{
Regularized (polychotomous) logistic regression
by Gibbs sampling. The package implements subtly different 
MCMC schemes with varying efficiency depending on the data type 
(binary v. binomial, say) and the desired estimator (regularized maximum
likelihood, or Bayesian maximum a posteriori/posterior mean, etc.) through a 
unified interface.
}
\details{
\tabular{ll}{
Package: \tab reglogit\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2011-08-05\cr
License: \tab LGPL \cr
LazyLoad: \tab yes\cr
}

See the documentation for the \code{\link{reglogit}} function
}
\author{
Robert B. Gramacy \email{rbg@vt.edu}


Maintainer: Robert B. Gramacy \email{rbg@vt.edu}
}
\references{
R.B. Gramacy, N.G. Polson. \dQuote{Simulation-based regularized
logistic regression}. (2010); arXiv:1005.3430; \url{http://arxiv.org/abs/1005.3430}
}

\keyword{ package }
\seealso{
\code{\link{reglogit}}, \code{\link[monomvn]{blasso}} and \code{\link[monomvn]{regress}}
}
\examples{
## see the help file for the reglogit function
}
