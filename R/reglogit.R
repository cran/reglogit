## functions for Gibbs sampling the parameters to the
## regularized logit regression model

## draw.beta:
##
## draw the regression coefficients conditional on the
## latent zs, lambdas, and omegas, and the penalty
## parameter nu, and thermo parameter kappa

draw.beta <- function(X, z, lambda, omega, nu, sigma,
                      kappa, kp, method)
  {
    ## negative nu and no sigma indicates L2 prior if sigma not null
    if(nu < 0 && !is.null(sigma)) n <- 1
    
    ## covariance of the conditional
    ## tXLi <- t(X) %*% diag(1/lambda)
    tXLi <- t(X * (1/lambda))
    Bi <- tXLi %*% X
    if(nu > 0) {
      bdi <- (kp/nu)^2 * 1/(sigma^2 * omega)
      if(ncol(X) == length(sigma) + 1) bdi <- c(0, bdi)
      diag(Bi) <- diag(Bi) + bdi
    } ## else Bi <- 0
    B <- solve(Bi) # + tXLi %*% X)

    ## mean of the conditional
    if(is.null(z)) { ## z=0 PDF representation

      ## choices between different (a,b) parameterizations
      if(method=="vaduva") { a <- 1; b <- kappa - 1 }
      else { a  <- 0.5; b <- kappa - 0.5 }
      kappa <- (a - 0.5*(a-b))

      ## special handling in Binomial case
      if(length(kappa == 1)) uk <- apply(X*kappa, 2, sum)
      else uk <- t(kappa) %*% X
      b <- B %*% uk

    } else ## CDF representation
      b <- B %*% tXLi %*% (z - 0.5*(1-kappa)*lambda)

    ## draw from the MVN
    return(rmvnorm(1, b, B))
  }


## draw.nu:
##
## draw from the full conditional distribution of the
## penalty parameter nu

draw.nu <- function(beta, sigma, kp, nup=list(a=0, b=0))
  {
    ## intercept?
    p <- length(sigma)
    if(p == length(beta) - 1) beta <- beta[-1]

    ## temper the prior
    a <- kp*(nup$a+1)-1
    b <- kp*nup$b

    ## update using beta
    a <- a + kp*p
    ## a <- a + p
    b <- b + kp*sum(abs(beta)/sigma)

    ## return a sample from an inverse gamma
    nu <- 1.0/rgamma(1, shape=a, scale=1/b)

    ## return sampled nu
    return(nu)
  }


## draw.z:
##
## draw the latent z-values given beta, y, and the
## latent lambdas

draw.z <- function(X, beta, lambda, kappa)
{
  n <- nrow(X)
  xbeta <- X %*% beta
  
  return(.C("draw_z_R",
            n = as.integer(n),
            xbeta = as.double(xbeta),
            beta = as.double(beta),
            lambda = as.double(lambda),
            kappa = as.double(kappa),
            kl = as.integer(length(kappa)),
            z = double(n),
            PACKAGE = "reglogit")$z)
}


## draw.omega:
##
## draw the latent omega variables given the betas
## and sigmas

draw.omega <- function(beta, sigma, nu, kp)
  {
    if(nu < 0) stop("nu must be positive")

    ## deal wth intercept in prior?
    if(length(sigma) == length(beta) - 1) beta <- beta[-1]
    m <- length(beta)
    
    return(.C("draw_omega_R",
              m = as.integer(m),
              beta = as.double(beta),
              sigma = as.double(sigma),
              nu = as.double(nu/kp),
              omega = double(m),
              PACKAGE = "reglogit")$omega)
  }


## draw.logit.prior
##
## testing function for the draw.lambda function
## to make sure that the powered-up proposals from the
## prior are working

draw.logit.prior <- function(n, kappa, kmax=442)
  {
     k <- 0:kmax
     lami <- 1 /(0.5*(1+k)*(kappa + k))
     lambda <- draw.lambda.prior(n, lami)
     x <- rnorm(length(lambda), sd=sqrt(lambda))
     return(x)
  }


## draw.lambda:
##
## draw the latent lambda variables that help us
## implement the logit link, via Metropolis-Hastings

draw.lambda <- function(X, beta, lambda, kappa, kmax,
                        zzero=TRUE, thin=kappa, method)
  {
    ## method ignored by the C-side if zzero = FALSE
    if(method=="vaduva") method <- 1
    else if(method=="slice") method <- 2
    else method <- 3

    ## calculate X * beta
    xbeta <- X %*% beta
    n <- length(xbeta)

    ## call C code
    return(.C("draw_lambda_R",
              n = as.integer(n),
              xbeta = as.double(xbeta),
              kappa = as.double(kappa),
              kl = as.integer(length(kappa)),
              kmax = as.integer(kmax),
              zzero = as.integer(zzero),
              thin = as.integer(thin),
              tl = as.integer(length(thin)),
              method = as.integer(method),
              lambda = as.double(lambda),
              PACKAGE = "reglogit")$lambda)
  }


## draw.lambda.prior:
##

draw.lambda.prior <- function(n, psii)
  {
    return(.C("draw_lambda_prior_R",
              n = as.integer(n),
              psii = as.double(psii),
              kmax = as.integer(length(psii)-1),
              lambda = double(n),
              PACKAGE = "reglogit")$lambda)
  }



## my.rinvgauss:
##
## draw from an Inverse Gaussian distribution following
## Gentle p193

my.rinvgauss <- function(n, mu, lambda)
  {
    return(.C("rinvgauss_R",
              n = as.integer(n),
              mu = as.double(mu),
              lambda = as.double(lambda),
              x = double(n),
              PACKAGE = "reglogit")$x)
  }


## calc.lpost:
##
## calculate the log posterior probability of the
## parameters in argument

calc.lpost <- function(yX, beta, nu, kappa, kp, sigma,
                       nup=list(a=0, b=0))
  {
    ## likelihood part
    llik <- - sum(kappa*log(1 + exp(-drop(yX %*% beta))))

    ## prior part
    
    ## deal with intercept in prior
    p <- length(sigma)
    if(p == length(beta) - 1) beta <- beta[-1]

    ## for beta
    if(nu > 0) lprior <- - kp*p*log(nu) - (kp/nu) * sum(abs(beta/sigma))
    ## if(nu > 0) lprior <- - p*log(nu) - (kp/nu) * sum(abs(beta/sigma))
    ## else if(!is.null(sigma)) lprior <- - sum((kp*beta/sigma)^2)
    else lprior <- 0

    ## inverse gamma prior for nu
    if(!is.null(nup))
      ## lprior <- lprior - kp*(2*(nup$a+1)*log(nu) + nup$b/(nu^2))
      lprior <- lprior - kp*((nup$a+1)*log(nu) + nup$b/nu)

    ## return log posterior
    return(llik + lprior)
  }


## preprocess:
##
## pre-process the X and y and kappa data depending on
## how binomial responses might be processed

preprocess <- function(X, y, N, flatten, kappa)
  {
    ## process y's and set up design matrices
    if(is.null(N) || flatten == TRUE) { 
      
      if(flatten) { ## flattened binomial case
        if(is.null(N)) stop("flatten only applies for N != NULL")
        ye <- NULL
        Xe <- NULL
        for(i in 1:length(N)) {
          Xe <- rbind(Xe, matrix(rep(X[i,], N[i]), nrow=N[i], byrow=TRUE))
          ye <- c(ye, rep(1, y[i]), rep(0, N[i]-y[i]))
        }
        y <- ye; X <- Xe;
      }

      ## create y in {-1,+1}
      ypm1 <- y
      ypm1[y == 0] <- -1
      
      ## create y *. X
      ## yX <- diag(ypm1) %*% X
      yX <- X * ypm1 ## diag(ypm1) %*% X
      
    } else {  ## unflattened binomial case
      
      ## sanity checks
      if(length(N) != length(y))
        stop("should have length(N) = length(y)")
      if(any(N < y)) stop("should have all N >= y")
      
      ## construct yX and re-use kappa
      n <- length(y)
      kappa <- kappa*c(y, N-y)
      knz <- kappa != 0
      yX <- rbind(X, -X)[knz,]
      y <- c(rep(1,n),rep(0,n))[knz]
      kappa <- kappa[knz]
      n <- length(kappa)
    }

    ## return the data we're going to work with
    return(list(yX=yX, y=y, kappa=kappa))
  }



## reglogit:
##
## function for Gibbs sampling from the a logistic
## regression model, sigma=NULL and nu < 0 indicates
## no prior

reglogit <- function(T, y, X, N=NULL, flatten=FALSE, sigma=1, nu=1,
                  kappa=1, icept=TRUE, normalize=TRUE,
                  zzero=TRUE, powerprior=TRUE, kmax=442,
                  bstart=NULL, lt=NULL, nup=list(a=2, b=0.1),
                  method=c("MH", "slice", "vaduva"), verb=100)
{
  ## getting and checking data size
  m <- ncol(X)
  n <- length(y)
  if(n != nrow(X)) stop("dimension mismatch")

  ## check the method argument
  method <- match.arg(method)
  
  ## design matrix processing
  X <- as.matrix(X)
  if(normalize) {
    one <- rep(1, n)
    normx <- sqrt(drop(one %*% (X^2)))
    if(any(normx == 0)) stop("degenerate X cols")
    X <- as.matrix(as.data.frame(scale(X, FALSE, normx)))
  } else normx <- rep(1, ncol(X))

  ## init sigma
  if(length(sigma) == 1) sigma <- rep(sigma, ncol(X))

  ## add on intecept term?
  if(icept) {
    X <- cbind(1, X)
    normx <- c(1, normx)
  }

  ## put the starting beta value on the right scale
  if(!is.null(bstart) && normalize) bstart <- bstart*normx
  else if(is.null(bstart)) bstart <- rnorm(m+icept)

  ## allocate beta, and set starting position
  beta <- matrix(NA, nrow=T, ncol=m+icept)
  beta[1,] <- bstart
  map <- list(beta=bstart)
  
  ## check if we are inferring nu
  if(!is.null(nup)) {
    nus <- rep(NA, T)
    if(nu <= 0) stop("starting nu must be positive")
    nus[1] <- nu
  } else nus <- NULL
  
  ## check for agreement between nu and sigma
  if(is.null(sigma) && nu > 0)
    stop("must define sigma for regularization")
  else if(nu == 0) sigma <- NULL

  ## initial values of the regularization latent variables
  if(nu > 0) {  ## omega lasso/L1 latents
    omega <- matrix(NA, nrow=T, ncol=m)
    ot <- omega[1,] <- 1
    omega[,1] <- 1
  } else { omega <- NULL; ot <- rep(1, m) }

  ## save the original kappa
  if(powerprior) kp <- kappa
  else kp <- 1

  ## process y's and set up design matrices
  nd <- preprocess(X, y, N, flatten, kappa)
  yX <- nd$yX; y <- nd$y; kappa <- nd$kappa
  n <- nrow(yX)
 
  ## initialize lambda latent variables
  lambda <- matrix(NA, nrow=T, ncol=n)
  if(is.null(lt)) lt <- rep(1, n)
  map$lambda <- lambda[1,] <- lt

  ## initial values for the logit latents
  if(!zzero) {
    z <- matrix(NA, nrow=T, ncol=n)
    z[1,] <- zt <- y
  } else { z <- zt <- NULL }

  ## allocate and initial log posterior calculation
  lpost <- rep(NA, T)
  map$lpost <- lpost[1] <-
    calc.lpost(yX, map$beta, nu, kappa, kp, sigma, nup)
  if(!is.null(nus)) map$nu <- nu
  
  ## Gibbs sampling rounds
  for(t in 2:T) {

    ## progress meter
    if(t %% verb == 0) cat("round", t, "\n")
    
    ## if regularizing, draw the latent omega values
    if(nu > 0) ot <- omega[t,] <- draw.omega(beta[t-1,], sigma, nu, kp)

    ## if logistic, draw the latent lambda values
    lt <- lambda[t,] <-
      draw.lambda(yX, beta[t-1,], lt, kappa, kmax, zzero, method=method)

    ## draw the latent z values
    if(!zzero) zt <- z[t,] <- draw.z(yX, beta[t-1,], lt, kappa)

    ## draw the regression coefficients
    beta[t,] <- draw.beta(yX, zt, lt, ot, nu, sigma, kappa, kp, method=method)

    ## maybe draw samples from nu
    if(!is.null(nus)) nu <- nus[t] <- draw.nu(beta[t,], sigma, kp, nup)
    
    ## calculate the posterior probability to find the map
    lpost[t] <- calc.lpost(yX, beta[t,], nu, kappa, kp, sigma, nup)

    ## update the map
    if(lpost[t] > map$lpost) {
      map <- list(beta=beta[t,], lpost=lpost[t], lambda=lambda[t,])
      if(!is.null(nus)) map$nu <- nus[t]
    }
      
  }

  ## un-normalize
  if(normalize) {
    beta <- as.matrix(as.data.frame(scale(beta, FALSE, scale=normx)))
    map$beta <- map$beta/normx
  }     

  ## construct the return object
  r <- list(X=X, y=y, beta=beta, lambda=lambda, lpost=lpost,
            map=map, kappa=kappa)
  if(nu > 0) r$omega <- omega
  if(normalize) r$normx <- normx
  if(!zzero) r$z <- z
  if(!is.null(nus)) r$nu <- nus

  ## assign a class to the object
  class(r) <- "reglogit"
  
  return(r)
}
