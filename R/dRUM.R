z.dRUM <- function(X, y, beta)
  {
    lambda <- exp(drop(X %*% beta))
    u <- runif(length(y))
    z <- log(lambda*u + y) - log(1-u+lambda*(1-y))
    return(z)
  }

beta.dRUM <- function(beta, X, z, b0, B0i)
  {
    ## proposal
    if(!is.null(B0i)) {
      B <- solve(B0i + (3/(pi^2)) * t(X) %*% X)
      b <- B %*% (B0i %*% b0 + (3/(pi^2)) * t(X) %*% z)
    } else {
      B <- solve((3/(pi^2)) * t(X) %*% X)
      b <- B %*% ((3/(pi^2)) * t(X) %*% z)
    }
    bnew <- drop(rmvnorm(1, b, B))

    ## calculate accpetance probability
    ## likelihood
    lalpha <- sum(dlogis(z - X %*% bnew, log=TRUE) - dlogis(z - X %*% beta, log=TRUE))
    ## prior
    if(!is.null(B0i)) {
      B0 <- solve(B0i)
      lalpha <- lalpha + dmvnorm(bnew, b0, B0, log=TRUE)
      lalpha <- lalpha - dmvnorm(beta, b0, B0, log=TRUE)
    }
    ## proposal
    lalpha <- lalpha + dmvnorm(beta, b, B, log=TRUE)
    lalpha <- lalpha - dmvnorm(bnew, b, B, log=TRUE)

    ## accept or reject
    if(runif(1) < exp(lalpha)) return(bnew)
    else return(beta)
  }


gibbs.dRUM <- function(T, X, y, beta.start=NULL, z.start=NULL, b0=NULL, B0=NULL,
                       verb=100)
  {
    ## storage
    beta <- matrix(NA, nrow=T, ncol=ncol(X))
    z <- matrix(NA, nrow=T, ncol=nrow(X))
    
    ## initial values
    if(is.null(beta.start)) beta[1,] <- rep(0, ncol(X))
    else beta[1,] <- beta.start
    if(is.null(z.start)) { z[1,] <- y; z[1,y == 0] <- -1
    } else beta[1,] <- beta.start
    
    ## Gibbs sampling
    for(t in 2:T) {

      if(t %% verb == 0) cat("round", t, "\n")
      
      z[t,] <- z.dRUM(X, y, beta[t-1,])
      beta[t,] <- beta.dRUM(beta[t-1,], X, z[t,], b0, B0)
      
    }

    return(list(beta=beta, z=z))
  }
  
