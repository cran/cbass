#' Draw samples of independent normals (matrix) given previous sample, and maximal values
#' 
#' @param mu n by d matrix of latent means
#' @param y n-length vector of maximal indices
#' @param z n by d matrix of latent random variables
#' @returns A new sample of n by d matrix of latent random variables
#' @export
#' @examples
#' set.seed(1)
#' n <- 100;  d <- 3
#' mu <- matrix(rnorm(n*d), n, d)
#' bound <- qnorm(1/d^(1/(d-1)))
#' mu[,1] <- bound
#' z <- mu
#' z[,-1] <- rnorm(length(mu[,-1]), mu[,-1], 1)
#' y <- apply(z, 1, which.max)
#' z.new <- sample.z(mu, y, z)
#' all(apply(z.new, 1, which.max) == y)
sample.z <- function(mu,    # n \times d matrix of latent means
                       y,     # n-length vector of maximal indices
                       z)     # n \times d matrix of latent random variables
{
  y <- as.numeric(y)
  d <- ncol(mu)
  n <- nrow(mu)
  
  
  # sd.star <- res$sd
  bound <- mu[1,1]
  
  z.out <- z
  z.out[,1] <- bound  # enforce this constraint
  # Sample those where all z's less than zero
  i0 <- which(y==1)
  if(length(i0) > 0){
    beta <- bound-mu[i0,2:d]
    # alpha <- -Inf
    z.out[i0,2:d] <- qnorm(runif(length(mu[i0,2:d]))*(pnorm(beta))) + mu[i0,2:d]
  }
  
  # Sample the non-max z's (and max z's, but max z's will be wrong)
  i1 <- which(y!=1)
  if(length(i1) > 0){
    ii <- cbind(i1, y[i1])  # max entry locations
    beta <- z[ii]-mu[i1,2:d]   # upper bound (transformed) for non-max z's
    # alpha <- -Inf
    z.out[i1,2:d] <- qnorm(runif(length(mu[i1,2:d]))*(pnorm(beta))) + mu[i1,2:d]
    
    # Now correctly sample max z's
    # z.temp <- z
    z.temp <- z.out   # swap to z.out, i.e., conditional on previous sample!
    z.temp[ii] <- NA
    z.max.2 <- apply(z.temp[i1,], 1, max, na.rm=TRUE)  # lower bound for max z's
    
    # beta <- zmax[i1]-mu[i1,]
    alpha <- z.max.2 - mu[ii]
    pa <- pnorm(alpha)
    z.out[ii] <- qnorm(pa + runif(length(i1))*(1-pa)) + mu[ii]
  }
  
  return(z.out)
}




# Slower, but guaranteed sampler, no numerical issues
#   used for initialization only 
sample.z.while <- function(mu,     # n \times d matrix of latent means
                           y,      # n-length vector of maximal indices
                           z,      # n \times d matrix of latent random variables
                           maxit=1e3)   # maximum number of iterations for while loops
{
  d <- ncol(mu)
  n <- nrow(mu)
  
  z.out <- z
  z.out[,1] <- 0
  z.out[,-1] <- matrix(rnorm(n*(d-1), c(mu[,-1])), n, d-1)
  
  i.check <- which(apply(z.out, 1, which.max) != y)
  count <- 0
  while(count < maxit & length(i.check) > 0){
    count <- count + 1
    
    z.out[i.check,-1] <- matrix(rnorm(length(i.check)*(d-1), c(mu[i.check,-1])), 
                                length(i.check), d-1)
    i.check <- i.check[which(apply(z.out[i.check,,drop=FALSE], 1, which.max) != y[i.check])]
  }
  
  return(z.out)
}



