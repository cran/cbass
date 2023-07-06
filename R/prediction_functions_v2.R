#' Generate chain of latent normal random variables for a given X, for values saved in 'mod'
#' 
#' @param mod CBASS model list
#' @param X matrix of covariates of same size / makeup as that used to create mod. If matrix not scaled to the unit interval, then it will be
#' @param nburn Number of samples to burn from the chain in mod, default 0
#' @param nsub Number of samples to subset to, default to those stored in mod
#' @returns An array of latent variables, nsub by nrow(X) by d
#' @export
#' @examples
#' set.seed(1)
#' n <- 100;  d <- 3
#' X <- matrix(runif(n*2, 0, 1), ncol=2)
#' mu <- scale(X)
#' bound <- qnorm(1/d^(1/(d-1)))
#' mu <- cbind(bound, mu)
#' z <- mu
#' z[,-1] <- rnorm(length(mu[,-1]), mu[,-1], 1)
#' y <- apply(z, 1, which.max)
#' mod <- fit.cbass(X, y, max.int=1, max.basis=10, nmcmc=1e3, nburn=500, nthin=10)
#' pred.chain <- pred.cbass(mod, X)
#' mu.hat <- apply(pred.chain, 2:3, mean)
#' round(p.mu(mu.hat[1,]), 3)
pred.cbass <- function(mod,X,
                          nburn=0,nsub=NULL)
{ 
  max.int=mod$max.int
  nmcmc<-nrow(mod$nbasis)
  d <- ncol(mod$nbasis)
  
  X <- as.matrix(X)
  if(max(X) > 1 | min(X) < 0){
    X <- t( (t(X) - mod$mins) / mod$ranges)
    if(max(X) > 1 | min(X) < 0){
      warning("After scaling, new X is outside ")
    }
  }
  Xt<-t(X)
  
  i.pred <- (nburn+1):nmcmc
  if(!is.null(nsub)){
    if( nsub < length(i.pred)){
      i.pred <- round(seq(nburn+1,nmcmc, length.out = as.numeric(nsub)))
    }
  }
  pred<-array(NA, c(length(i.pred),nrow(X),d))
  n <- nrow(X)
  
  bound <- qnorm(1/d^(1/(d-1)))
  k.range <- 2:d
  pred[,,1] <- bound
  
  
  count <- 0   # storage index
  for(i in i.pred){
    count <- count + 1
      for(k in k.range){
        if(mod$nbasis[i,k] > 0){
          B<-matrix(nrow=nrow(X),ncol=mod$nbasis[i,k]+1)
          B[,1]<-1
          for(j in 1:mod$nbasis[i,k]){
            ltemp <- mod$levels[i,j,,,k]
            if(max.int==1){
              ltemp <- matrix(ltemp, nrow=1)
            }
            B[,j+1]<-makeBasis.2(mod$signs[i,j,1:mod$nint[i,j,k],k],mod$vars[i,j,1:mod$nint[i,j,k],k],mod$knots[i,j,1:mod$nint[i,j,k],k],Xt,
                                 ordinal=mod$ordinal, levels=ltemp[1:mod$nint[i,j,k],,drop=FALSE])
          }
        } else {
          B <- matrix(1, nrow(X), 1)
        }
        
        pred[count,,k]<-B%*%mod$beta[i,1:(mod$nbasis[i,k]+1),k]
      }
  }
  return(pred)
}





#' Predict vector of probabilities from vector of latent means
#' 
#' @param mu d-length vector of latent means
#' @param d input to avoid computing length of mu
#' @param bound input of mu[1] to avoid computation
#' @param npts number of integration points, default 100
#' @param rel.tol number of integration points, default 1e-4 
#' @returns A d-length numeric vector of probabilities given input latent means
#' @export
#' @examples
#' set.seed(1)
#' mu <- rnorm(5)
#' p.mu(mu)
p.mu <- function(mu,  # d-length vector of latent means
                 d=NULL,   # optional input to avoid computing length
                 bound=NULL,  # optional input of mu[1] to avoid computing
                   npts=100,   # number of integration points
                 rel.tol=1e-4)   # relative tolerance for integration
{
  if(is.null(d)){
    d <- length(mu)
    bound <- qnorm(1/d^(1/(d-1)))
  }
  eps <- rel.tol/d
  
  mu[1] <- bound
  p.out <- sapply(1:d, function(x) prod(pnorm(mu[x]-mu[-x])))
  p.out[-1] <- p.out[-1] / sum(p.out[-1]) * (1-p.out[1])
  
  j.fix <- setdiff(2:(d-1), which(p.out< eps | p.out > 1-eps))
  if(length(j.fix) > 0){
    for(j in j.fix){    
      p.out[j] <- my.integrate(mu, j, npts, rel.tol)
    }
    p.out[d] <- 1-sum(p.out[-d])
  }
  
  return(p.out)
}



# Integration wrapper of f.mu with automatic bounds
my.integrate <- function(mu,  # d-length vector of latent means    
                         j,   # latent index to compute for
                         npts=100,    # number of integration points
                         rel.tol=1e-4)    # relative tolerance for integration
{
  lb <- max(c(mu[1], mu[j]-10))
  ub <- max(c(mu[1]+10, mu[j]+10))
  
  integrate(f.mu, lb, ub, mu=mu, jj=j,
            subdivisions = npts,
            rel.tol = rel.tol)$val   #/ (1-pnorm(bound-mu[j]))^(d-1) 
  
}


# Function to integrate
f.mu <- function(z,   # integration coordinate
                 mu,  # d-length vector of latent means   
                 jj)  # latent index to compute for
{
  # d <- length(mu)
  # prod(sapply(mu.vec[-c(1,jj)], function(x) pnorm(z-x)))*dnorm(z, mu.vec[jj],1)
  out <- dnorm(z, mu[jj],1)
  for(k in setdiff(2:length(mu), jj)){
    out <- out*pnorm(z-mu[k])
  }
  out #*prod(pnorm(vec))
}
