#' Augment X for missing data approach for MNAR
#' 
#' @param X matrix of covariates, including some missing values (NAs)
#' @returns Matrix same size as X, with augmented columns and zeros in the missing spots
#' @export
#' @examples
#' set.seed(1)
#' n <- 100
#' X <- matrix(runif(n*2, 0, 1), ncol=2)
#' X[sample(1:length(X), round(.1*length(X)))] <- NA
#' X.new <- augment.X(X)
#' sum(is.na(X.new))
augment.X <- function(X)
{
  X.aug <- as.matrix(X)
  
  if(sum(is.na(X.aug)) > 0){
    for(j in 1:ncol(X)){
      if(sum(is.na(X[,j])) > 0){
        X.aug[which(is.na(X.aug[,j])),j] <- 0
        X.aug <- cbind(X.aug, 1*is.na(X[,j]))
      }
    }
  } else {
    warning("No NAs found in X")
  }
  
  X.aug
}



# Compute symmetric square root of A, assuming it is real, symmetric, positive definite
my_symm_square_root <- function(A, invert=TRUE, prec=round(log10(1/.Machine$double.eps)))
{
  if(length(c(A)) == 1){
    if(as.logical(invert)){   # if scalar
      r1 <- 1/sqrt(c(A))
    } else {
      r1 <- sqrt(c(A))
    }
  } else {
    e1 <- eigen(A)
    vals <- round(Re(e1$val), prec)   # get rid of barely negative zeros
    vals[vals<0] <- 0
    if(as.logical(invert)){
      r1 <- e1$vec %*% as.matrix(diag(1/sqrt(vals))) %*% t(e1$vec)   # symmetric square root 
    } else {
      r1 <- e1$vec %*% as.matrix(diag(sqrt(vals))) %*% t(e1$vec)   # symmetric square root 
    }
  }
  r1
}



# Master MakeBasis function that includes categorial and ordinal
# Ordinal is indicator of size of data indicating which columns of X are ordinal
#   all inputs are vectors, 
#      except datat which is t(X)
#      except levels which is a ragged matrix of nint by max # factor levels matrix, indicating levels used in categorical predictor
makeBasis.2 <- function(signs,vars,knots,datat,degree=1,ordinal=NULL,levels=NULL)
{
  if(is.null(ordinal)){   # if categorical not indicated
    return(makeBasis.continuous(signs,vars,knots,datat,degree))
  } else {
    ordinal.vars <- as.logical(ordinal[vars])
    if(sum(ordinal.vars) > 0){
      temp1 <- makeBasis.continuous(signs[ordinal.vars],vars[ordinal.vars],knots[ordinal.vars],datat,degree)
    } else {
      temp1 <- 1
    }
    if(sum(!ordinal.vars) > 0){
      temp2 <- makeBasis.categorical(vars[!ordinal.vars],datat,levels[!ordinal.vars,,drop=FALSE], signs=signs[!ordinal.vars])
    } else {
      temp2 <- 1
    }
    return(temp1*temp2)
  }
}




# Categorical version of makebasis function
makeBasis.categorical <- function(vars,datat,levels,signs=rep(1, length(vars)))
{ 
  # function to make a basis function given signs, variables, knots, and data
  temp1<-rep(1, ncol(datat))
  for(pp in 1:length(vars)){ # faster than apply
    temp0 <- 1*(datat[vars[pp],] %in% levels[pp,])
    # if(as.numeric(signs[pp]) < 0){
    #   temp0 <- 1 - temp0  # swap 1's and zero's
    # }
    temp1<-temp1*temp0
  }
  
  return(temp1)
}



# Continuous version of makebasis function
makeBasis.continuous <- function(signs,vars,knots,datat,degree=1)
{ # function to make a basis function given signs, variables, knots, and data
  temp1<-pos(signs*(datat[vars,,drop=F]-knots))^degree # this only works for t(data)...
  if(length(vars)==1){
    return(c(temp1))
  } else{
    temp2<-1
    for(pp in 1:length(vars)){ # faster than apply
      temp2<-temp2*temp1[pp,]
    }
    return(temp2)
  }
}



# function to zero out any negative values (there are many ways to do this)
pos<-function(vec)
{ 
  (abs(vec)+vec)/2
}


# differences from probs given mu, and p
prob.wrapper <- function(mu, p)
{
  max(abs(p.mu(c(0,mu)) - p))
}

# wrapper to find mu for given p
invert.p.mu <- function(p, maxit=1e4)
{
  mu.start <- -qnorm(1-p)   # init at very rough guess
  mu.start <- mu.start[-1] - mu.start[1]/2
  # alphas <- seq(.0001, .9999, length.out = maxit)
  # delta <- 0*alphas
  # for(i in 1:maxit){
  #   delta[i] <- prob.wrapper(mu.start[-1]-mu.start[1]*alphas[i], p)
  # }
  # mu.start <- mu.start[-1] - mu.start[1]*alphas[which.min(delta)]
  
  # mu.start <- mu.start[-1]
  # mu.start <- mu.start[-1] - mu.start[1]*.5
  
  res <- optim(par=mu.start, 
               fn=prob.wrapper, method="BFGS",
               p=p, control=list(maxit=maxit))
  
  res$par <- c(0, res$par)
  res$mu <- res$par
  res
}


# Propose a birth and compute its probability
propose.birth <- function(m,  # number of variables to select
                          n,  # possible number of variables
                          prob=NULL, # probability of inclusion of each variable
                          type="uniform")
{
  if(length(prob) != n){
    # cat("i = ", i, "\n")
    # cat("move.type = ", move.type, "\n")
    # cat("n=", n, "; m=", m,  "; prob=", prob, "\n")
    # cat("nbasis=", nbasis, "\n")
    stop("Length of probabilities must match number of variables")
  }
  
  if(type == "uniform"){
    out <- sample(1:n, m, replace=FALSE, prob=NULL)
    # p.out <- 1/choose(n,m)
    p.out <- NA
  } else {
    prob.in <- prob/sum(prob)   # save input probability for later, modify prob as going
    
    out <- NULL  # set of variables selected
    for(j in 1:m){
      temp <- sample(1:n, 1, prob=prob)
      out <- c(out, temp)
      prob[temp] <- 0
    }
    
    p.out <- compute.prob.set(out, prob.in)  # compute probability
  }
  
  return(list(vars=out, p=p.out))
}



# Probability of choosing indices "set" from vector of probabilities "prob"
compute.prob.set <- function(set, prob)
{
  prob <- prob/sum(prob)  # renormalize just in case
  
  n <- length(prob)
  m <- length(set)
  
  if(m == 1){
    pout <- prob[set]
  } else {
    num <- prod(prob[set])   # probability of numerator, i.e. product of probabilities of set
    perms <- matrix(set[permutations(m)], ncol=m)   # possible orderings of set
    
    denom <- rep(1, nrow(perms)) # normalizing constant
    c.prob <- 0
    for(i in 1:(m-1)){
      c.prob <- c.prob + prob[perms[,i]]   # probability of already-chosen balls
      denom <- denom/(1-c.prob)   
    }
    pout <- num*sum(denom)
  }
  
  return(pout)
}


# generate all permutations of 1:n 
permutations <- function(n)
{
  if(n==1){
    return(matrix(1))
  } else {
    sp <- permutations(n-1)
    p <- nrow(sp)
    A <- matrix(nrow=n*p,ncol=n)
    for(i in 1:n){
      A[(i-1)*p+1:p,] <- cbind(i,sp+(sp>=i))
    }
    return(A)
  }
}


# Lambda sampler for truncated poisson distribution and gamma hyperprior
#  uses MH, but proposes from unconstrained posterior
sample.lambda.mh <- function(lambda.old, nbasis, max.basis,
                             h1, h2)
{
  lambda.new <- rgamma(1, h1 + nbasis, h2 + 1)
  
  if(lambda.new > 0){   # there is zero density at lambda.new == 0, which can happen numerically
    prop.fwd <- dgamma(lambda.new, h1 + nbasis, h2 + 1, log=TRUE)
    prop.rev <- dgamma(lambda.old, h1 + nbasis, h2 + 1, log=TRUE)
    
    d.new <- dpois(nbasis, lambda.new, log=TRUE) + dgamma(lambda.new, h1, h2, log=TRUE) - ppois(max.basis, lambda.new, lower.tail=TRUE, log.p=TRUE)
    d.old <- dpois(nbasis, lambda.old, log=TRUE) + dgamma(lambda.old, h1, h2, log=TRUE)  - ppois(max.basis, lambda.old, lower.tail=TRUE, log.p=TRUE)
    
    alpha <- d.new - d.old + prop.rev - prop.fwd
    # if(is.na(alpha)){
    #   cat("lambda.old=", lambda.old, "\n")
    #   cat("lambda.new=", lambda.new, "\n")
    #   cat("nbasis=", nbasis, "\n")
    #   cat("max.basis=", max.basis, "\n")
    #   cat("h1,h2=", c(h1, h2), "\n")
    #   stop("NA alpha in sample lambda")
    # }
    if(log(runif(1,0,1)) < alpha){
      lambda.out <- lambda.new
      accept <- TRUE
    } else {
      lambda.out <- lambda.old
      accept=FALSE
    }
  } else {   # zero density at lambda.new = 0, so automatically reject
    lambda.out <- lambda.old
    accept <- FALSE
  }
  
  return(list(lambda=lambda.out, accept=accept))
}

# Scale vector to unit interval
my.scale <- function(x)
{
  c(x - min(x))/c(max(x) - min(x))
}


