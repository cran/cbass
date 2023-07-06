#' Fit CBASS model using reversible jump MCMC
#' 
#' @param X n by p matrix of inputs on unit invteral
#' @param y n-length factor factor of categories
#' @param max.int maximum number of interactions, default 3
#' @param max.basis maximum number of basis functions for each latent mean function, default ncol(X)*10
#' @param tau2 prior variance of basis regression coefficients, default 10
#' @param nmcmc number of MCMC samples, default 1e4
#' @param nburn number of MCMC samples to burn, default nmcmc/2
#' @param nthin number of samples by which to thin, default 10
#' @param h1 first parameter for Gamma hyperprior on tau2, default 4
#' @param h2 second parameter for Gamma hyperprior on tau2, default 20(d-1)/n
#' @param p.int.prior prior for number of interactions, default 1/(1:max.int)
#' @param verbose should progress be printed out? default false
#' @param print.interval how often should progress be printed out? default every 1\%
#' @param init1 should model be initialized with single interaction model? default FALSE
#' @param ordinal indicator of ordinal predictors (non-categorical), usually computed automatically
#' @param writeout should samples be written out to text file? default FALSE
#' @param writedir where should text files be written? default tempdir()
#' @param mod initial / previous model fit, default NULL
#' @param restart should initial input model be used for starting chain? default FALSE
#' @importFrom stats dgamma dnorm dpois integrate optim pnorm ppois qnorm rgamma rnorm runif
#' @importFrom utils tail write.table
#' @returns A list of CBASS model parameters. LIST THEM.
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
#' mean(abs(mu - mu.hat))
#' plot(c(mu), c(mu.hat))
fit.cbass <- function(X, # training input matrix
                      y, # training classification vector
                      max.int=3, # maximum degree of interaction (Jmax)
                      max.basis=10*ncol(X), # maximum nuber of basis functions
                      tau2=10, # prior variance for basis coefficients
                      nmcmc=1e4, # number of MCMC iterations
                      nburn=round(nmcmc/2),   # number of MCMC iterations to burn
                      nthin=10,   # number of samples to thin by
                      h1=4,h2=20*(length(unique(y))-1)/nrow(X), # shape and rate for Gamma prior on mean number of basis functions (number of basis functions nbasis has Poisson(lam) prior where lam has a Gamma hyperprior => marginal prior for nbasis is negative binomial)
                      p.int.prior=1/(1:max.int),   # prior on number of interactions
                      verbose=FALSE,     # print out progress?
                      print.interval=round(nmcmc/100),    # how often to print out progress
                      init1=FALSE,   # initialize 
                      ordinal=NULL,    # indicator of ordinal predictors (non-categorical)
                      writeout=FALSE,   # write out MCMC samples as progress made?
                      writedir=tempdir(),  # where to write out MCMC samples
                      mod=NULL,      # input initial model
                      restart=FALSE   # use input initial model?
)
{
  
  # Fixed legacy imput
  mu.init=NULL   # initial mean
  lambda.init=1  # initial lambda value
  birth.prop="uniform" # should birth proposals be uniform or weighted?
  w1=10; w2=10    # weighted birth proposal parameters
  knot.prop="discrete"   # knot proposals
  tau.knot=NULL    # knots
  save.burn=FALSE   # save the burned iterations?
  
  X <- as.matrix(X)
  if(max(X) > 1 | min(X) < 0){
    mins.X <- apply(X, 2, min)
    ranges.X <- apply(X, 2, function(x) diff(range(x)))
    X <- apply(X, 2, my.scale)
  } else {
    mins.X <- rep(0, ncol(X))
    ranges.X <- rep(1, ncol(X))
  }
  
  
  #### Input processing
  Xt<-t(X)
  n<-length(y)
  p<-ncol(X)
  tau2.0 <- tau2
  knot.prop="discrete"
  tau.knot=NULL
  save.burn=FALSE
  if(is.null(nburn)){
    nburn=round(nmcmc/2)
  }
  max.vars <- max.basis*max.int
  max.prop <- 10
  
  if(is.null(ordinal)){
    ordinal <- rep(1, ncol(X))
    ordinal[apply(X, 2, function(z) length(unique(z)) <= 10 )] <- 0
  }
  
  min.nonzero <- FALSE
  if(is.null(min.nonzero)){
    min.nonzero <- min(10, round(.05*n))
  } else {
    min.nonzero <- min(min.nonzero, round(.05*n))
  }
  eps <- .Machine$double.eps  # threshold for condition number
  
  
  
  # Classes
  d <- length(table(y))
  y <- as.factor(y)
  classes <- levels(y) #sort(unique(y))
  y.num <- match(y, classes)   # numeric position of each class
  bound <- qnorm(1/d^(1/(d-1)))
  
  
  # Make X numeric if not
  X <- as.matrix(X)
  for(j in 1:ncol(X)){
    X[, j] <- as.numeric(X[, j])
  }
  
  # Possible knot locations and number of knot locations
  vars.list <- lapply(1:ncol(X), function(z) sort(unique(X[,z])))
  vars.len <- sapply(vars.list, length)
  ####
  
  
  
  #### storage
  i.save <- (nburn+1):nmcmc
  i.save <- i.save[which(i.save %% nthin == 0)]
  n.save <- length(i.save)
  knots<-signs<-vars<-array(dim=c(n.save,max.basis,max.int,d)) # this is the largest possible, our filling will be ragged
  nint<-array(dim=c(n.save,max.basis,d)) # degree of interaction J, again filling will be ragged
  beta<-array(dim=c(n.save,max.basis+1,d)) # +1 for intercept, again filling will be ragged
  nbasis<-nvars<-lam<-tau2<-matrix(NA, n.save, d) # error variance, poisson hyperprior, and number of basis functions
  z <- matrix(0, n, d)   # latent unobservable variables
  # lam<-rep(NA, n.save)
  valid.prop <- n.prop <- rep(NA, n.save)
  
  knots.curr<-signs.curr<-vars.curr<-array(dim=c(max.basis,max.int,d))
  nint.curr<-array(dim=c(max.basis,d))
  
  
  
  if(is.null(p.int.prior)){
    p.int.prior <- rep(1/max.int, max.int)
    int.prior="uniform"
  } else {
    if(is.numeric(p.int.prior) & length(p.int.prior) == max.int){
      p.int.prior <- p.int.prior/sum(p.int.prior)   # renormalize
    } else {
      stop("p.int.prior must be numeric vector of length equal to max.int")
    }
    int.prior="weighted"
  }
  
  
  
  
  
  # If restarting chain, initialize from model, otherwise, follow initialization plan
  if(as.logical(restart)){
    ii <- nrow(mod$nbasis)
    
    knots.curr <- mod$knots[ii,,,]
    signs.curr <- mod$signs[ii,,,]
    vars.curr <- mod$vars[ii,,,]
    nint.curr <- mod$nint[ii,,]
    beta.curr<-mod$beta[ii,,]
    
    if(sum(ordinal) < length(ordinal)){
      m1 <- max(apply(X[,!ordinal,drop=FALSE], 2, function(z) length(unique(z))))   # max number of levels
    } else {
      m1 = 2   # arbitrary
    }
    levels <- array(dim=c(n.save,max.basis,max.int,m1,d))   # levels for smoothed categorical preictors
    levels.curr <- mod$levels[ii,,,,]
    
    lam.curr <- mod$lam[ii,]
    tau2.curr <- mod$tau2[ii,]
    nbasis.curr <- mod$nbasis[ii,]
    nvars.curr <- mod$nvars[ii,]
    
    
    z <- mod$z
    p.int <- mod$p.int
    p.vars <- mod$p.vars
    X.curr <- mod$X
    Vinv.curr <- mod$V
    d.curr <- mod$d
    log.det.curr <- mod$log.det
    
    
    rm(mod)  ;  gc()
    
    
    levels[1,,,,] <- levels.curr
    nint[1,,] <- nint.curr
    knots[1,,,] <- knots.curr
    signs[1,,,] <- signs.curr
    vars[1,,,] <- vars.curr
    nbasis[1,] <- nbasis.curr
    nvars[1,] <- nvars.curr
    lam[1,] <- lam.curr
    tau2[1,] <- tau2.curr
    beta[1,,] <- beta.curr
    
    
  } else {  # generate initialization here
    if(birth.prop!= "uniform"){
      p.int <- matrix(w1, max.int, d)
      p.vars <- matrix(w2, p, d)
    } else {
      p.int <- matrix(NA, max.int, d)
      p.vars <- matrix(NA, p, d)
    }
    
    if(sum(ordinal) < length(ordinal)){
      m1 <- max(apply(X[,!ordinal,drop=FALSE], 2, function(z) length(unique(z))))   # max number of levels
    } else {
      m1 = 2   # arbitrary
    }
    levels <- array(dim=c(n.save,max.basis,max.int,m1,d))   # levels for smoothed categorical preictors
    levels.curr <- array(dim=c(max.basis,max.int,m1,d)) 
    
    beta.curr<-matrix(NA, max.basis+1, d)
    
    
    if(init1){   # initialize with single-interaction model
      mod.init <- init.mod.1(X, y, 
                                   tau2=tau2.0, 
                                   max.basis = min(max.basis, ncol(X)*2),
                                   h1=h1, h2=h2,
                                   nmcmc=5e3, nthin=100, nburn=3e3,
                                   birth.prop=birth.prop, 
                                   ordinal=ordinal)
      z <- mod.init$z
      X.curr <- mod.init$X
      Vinv.curr <- mod.init$Vinv
      d.curr <- mod.init$d
      log.det.curr <- mod.init$log.det
      
      knots.curr[1:dim(mod.init$knots)[1],1,] <- mod.init$knots
      signs.curr[1:dim(mod.init$signs)[1],1,] <- mod.init$signs
      vars.curr[1:dim(mod.init$vars)[1],1,] <- mod.init$vars
      nint.curr[1:nrow(mod.init$nint),] <- mod.init$nint
      levels.curr[1:dim(mod.init$levels)[1], 1, , ] <- mod.init$levels
      nbasis.curr <- mod.init$nbasis
      nvars.curr <- mod.init$nvars
      lam.curr <- mod.init$lam
      beta.curr[1:nrow(mod.init$beta),] <- mod.init$beta
      # tau2.curr <- rep(tau2.0,d)
      tau2.curr <- mod.init$tau2
      
      # h1=mod.init$h1
      # h2=1
      
      rm(mod.init)
      gc()
      
      levels[1,,,,] <- levels.curr
      nint[1,,] <- nint.curr
      knots[1,,,] <- knots.curr
      signs[1,,,] <- signs.curr
      vars[1,,,] <- vars.curr
      nbasis[1,] <- nbasis.curr
      nvars[1,] <- nvars.curr
      lam[1,] <- lam.curr
      tau2[1,] <- tau2.curr
      beta[1,,] <- beta.curr
      
    } else {
      tau2.curr <- rep(tau2.0,d)
      
      #### Initialize mean parameters
      if(is.null(mu.init)){
        frac <- table(y.num) / n  #  1/d
        # mu <- mu.init <- qnorm(frac) - qnorm(1/d)
        # mu <- outer(rep(1,n), mu)
        mu <- invert.p.mu(frac)$mu
        mu <- outer(rep(1,n), mu)
      } else {
        mu <- mu.init
      }
      
      # FIXED INIT!
      z <- bound
      for(j in 2:d){
        z <- cbind(z, rnorm(n, mean=mu[1,j]))
      }
      # fix y
      i.fix <- which(apply(z, 1, which.max) != y.num)
      for(i in i.fix){
        # eps <- runif(1)
        lb <- max(z[i, -y.num[i]])
        alpha <- lb - mu[i,y.num[i]]
        pa <- pnorm(alpha)
        z[i,y.num[i]] <- qnorm(pa + runif(1)*(1-pa)) +mu[i,y.num[i]]
      }
      
      # Good starting point for z
      z <- sample.z.while(mu, y.num, z, maxit=1e4)
      z[,1] <- bound  # ensure this is true
      # z <- sample.z.debug(mu, y.num, z)
      
      # initialize beta
      nbasis.curr<-nvars.curr<-rep(0,d)
      lam.curr <- rep(min(lambda.init,max.basis-.Machine$double.eps), d)
      
      X.curr <- Vinv.curr <- d.curr <- log.det.curr <- list()
      bhat <- rep(0,d)
      for(j in 1:d){
        X.curr[[j]]<-matrix(rep(1,n)) # matrix of current basis functions, so that yhat = X.curr %*% beta
        Vinv.curr[[j]]<-crossprod(X.curr[[j]])+1/tau2.curr[j] # V^(-1) from DMS
        log.det.curr[[j]] <- 1/2*determinant(Vinv.curr[[j]])$mod   # log determinant to avoid recalculating
        Xy <- crossprod(X.curr[[j]],z[,j])
        bhat[j] <- solve(Vinv.curr[[j]], Xy)
        d.curr[[j]] <- 1/2*sum(Xy*bhat[j])
      }
      beta.curr[1,]<-bhat
      beta[1,,] <- beta.curr
      
      nbasis[1,] <- nbasis.curr
      nvars[1,] <- nvars.curr
      lam[1,] <- lam.curr
      tau2[1,] <- tau2.curr
    }
    
  }  # end if restart
  
  
  if(writeout & is.character(writedir)){
    write.table(matrix(c(nint.curr, nrow=1)),
                file=file.path(writedir, "nint.txt"), append=TRUE, col.names = FALSE)
    write.table(matrix(c(beta.curr, nrow=1)),
                file=file.path(writedir, "beta.txt"), append=TRUE, col.names = FALSE)
    write.table(matrix(c(knots.curr, nrow=1)),
                file=file.path(writedir, "knots.txt"), append=TRUE, col.names = FALSE)
    write.table(matrix(c(signs.curr, nrow=1)),
                file=file.path(writedir, "signs.txt"), append=TRUE, col.names = FALSE)
    write.table(matrix(c(levels.curr, nrow=1)),
                file=file.path(writedir, "levels.txt"), append=TRUE, col.names = FALSE)
    write.table(matrix(c(vars.curr, nrow=1)),
                file=file.path(writedir, "vars.txt"), append=TRUE, col.names = FALSE)
  }
  
  
  # mu <- sapply(2:d, function(z) X.curr[[z]]%*%beta.curr[1:ncol(X.curr[[z]]),z])   # BMARS mean estimate
  # mu <- cbind(0, mu)
  
  basis.samples<-array(0, c(n.save, d, 3))
  count<-matrix(0, d, 3) # count how many times we accept birth, death, change
  colnames(count) <- dimnames(basis.samples)[[3]] <- c("birth", "death", "change")
  accept.lambda <- 1
  valid.prop.curr <- TRUE
  n.prop.curr <- 1
  inf.flags <- NULL
  
  
  # start MCMC
  for(i in 2:nmcmc){
    
    # Update z's
    mu <- sapply(2:d, function(z) X.curr[[z]]%*%beta.curr[1:ncol(X.curr[[z]]),z])   # BMARS mean estimate
    mu <- cbind(bound, mu)
    
    inf.check <- TRUE
    count0 <- 0
    while(inf.check){
      temp <- sample.z(mu, y.num, z)
      inf.check <- sum(is.infinite(temp)) != 0
      count0 <- count0 + 1
    }
    z <- temp
    if(count0 > 1){
      inf.flags <- rbind(inf.flags,
                         c(i,j,count0))
    }
    
    
    # randomize sampling order
    j.range <- 2:d  #sample(2:d)
    for(j in j.range){
      accept <- FALSE
      
      # Vinv.curr[[j]]<-crossprod(X.curr[[j]])+1/tau2[i,j]*diag(ncol(X.curr[[j]])) # V^(-1) from DMS
      Xy <- crossprod(X.curr[[j]],z[,j])
      bhat.cand<-solve(Vinv.curr[[j]], Xy)
      d.curr[[j]] <- 1/2*sum(Xy*bhat.cand)
      
      ## Reversible jump step
      if(nbasis.curr[j]==0){
        move.type<-'birth'
      } else if(nbasis.curr[j]==max.basis-1){
        move.type<-sample(c('death','change'),1)
      } else {
        move.type<-sample(c('birth','death','change'),1)
      }
      
      
      if(move.type=='birth'){
        valid.proposal.temp <- FALSE
        count.temp <- 0
        
        while(!valid.proposal.temp & count.temp < max.prop){
          count.temp <- count.temp + 1
          
          if(birth.prop=="uniform"){   # uniform proposal probability
            nint.cand <- sample(max.int,1) # sample degree of interaction for new basis function (different from DMS)
            
            vars.cand <- sample(p,nint.cand,replace = FALSE) # variables to use in new basis function (different from DMS)
            ll.birth.vars <- ll.birth.int <- 0 # cancels when proposing from same distr as prior
            if(int.prior != "uniform"){
              ll.birth.int <- log(p.int.prior[nint.cand]) + log(max.int)  # prior - uniform proposal
            }
          } else {
            nint.cand <- sample(max.int, 1, prob=p.int)
            
            res <- propose.birth(nint.cand,  # number of variables to select
                                 p,  # possible number of variables
                                 prob=p.vars, # probability of each variable
                                 type="weighted")
            vars.cand <- res$vars
            ll.birth.vars <- log(1/choose(p,nint.cand)) - log(res$p)   # uniform prior  - weighted proposal
            ll.birth.int <- log(p.int.prior) - log(sum(p.int)/p.int[nint.cand])   #  prior  - weighted proposal
          }
          
          
          if(knot.prop == "discrete"){
            knots.cand<-sapply(vars.list[vars.cand], sample, size=1) # sample knots for new basis function from uniform (discrete) distribution
          } else {
            knots.cand<-runif(nint.cand,0,1) # sample knots for new basis function from uniform (continuous) distribution
          }
          
          signs.cand<-sample(c(-1,1), nint.cand, replace = TRUE) # signs for new basis function
          levels.cand <- matrix(NA, nint.cand, dim(levels.curr[,,,j])[length(dim(levels.curr[,,,j]))]) #0*levels
          
          if(sum(!ordinal[vars.cand]) > 0){
            for(k in which(!ordinal[vars.cand])){
              # if(!ordinal[vars.cand[k]]){  # if categorical, sample a subset
              cat.poss <- vars.list[[vars.cand[k]]]  #sort(unique(tX[vars.cand[k], ]))
              if(length(cat.poss) <= 3){
                prob.temp <- rep(1, length(cat.poss)-1)
              } else {
                # prob.temp <- 1/sqrt(choose(length(cat.poss), 1:(length(cat.poss)-1)))   # wtd prob
                prob.temp <- rep(1, length(cat.poss)-1)
              }
              u <- sample(1:(length(cat.poss)-1), 1, prob=prob.temp)  # at least 1 and not all 
              levels.temp <- sample(cat.poss, u, replace=FALSE)
              levels.cand[k,1:u] <- levels.temp
              knots.cand[k] <- NA
              signs.cand[k] <- 1
              
              prob.temp <- prob.temp / sum(prob.temp)
              prob.prior <- 1/sqrt(choose(length(cat.poss), 1:(length(cat.poss)-1)))
              prob.prior <- prob.prior / sum(prob.prior)
              # prob.prior = prob.temp
              ll.birth.vars <- ll.birth.vars - log(prob.temp[u]) + log(prob.prior[u])
            }
          }
          
          basis.cand<-makeBasis.2(signs.cand,vars.cand,knots.cand,Xt,ordinal=ordinal,levels=levels.cand) # make the new basis function
          X.cand<-cbind(X.curr[[j]],basis.cand) # add the new basis function to the basis functions we already have
          
          Vinv.cand <- Vinv.curr[[j]]
          x.new <- c(crossprod(basis.cand,X.curr[[j]]))
          Vinv.cand <- rbind(Vinv.cand, x.new)
          Vinv.cand <- cbind(Vinv.cand, c(x.new, sum(basis.cand^2) + 1/tau2.curr[j]))
          chol.cand <- chol.default(Vinv.cand)
          diag.cand <- diag(chol.cand)
          cond.cand <- (min(diag.cand)/max(diag.cand))^2  # approx condition number
          
          nonzero <- sum(abs(basis.cand) > eps) >= min.nonzero
          nonzero2 <- cond.cand > eps   # computationally nonsingular
          
          valid.proposal.temp <- nonzero & nonzero2
        }
        
        if(valid.proposal.temp){  # continue to check if there are sufficient nonzero entries
          log.det.cand <- sum(log(diag.cand))   # 1/2 log determinant of Vinv
          Xy <- crossprod(X.cand,z[,j])
          bhat.cand <- chol2inv(chol.cand)%*%Xy
          d.cand <- sum(Xy*bhat.cand)/2
          llik.alpha <- .5*log(1/tau2.curr[j]) - log.det.cand + log.det.curr[[j]] + d.cand - d.curr[[j]]
          # llam <- log(lam.curr)-log(nvars.curr[j]+length(vars.cand)) 
          llam <- log(lam.curr[j])-log(nbasis.curr[j]+1)
          
          # only lambda prior enters into acceptance ratio when knots and vars are uniform
          alpha <- llik.alpha + llam + ll.birth.vars + ll.birth.int
        } else {  # else for failed test of interesting function to explore (computationally nonsingular e.g.)
          alpha= -Inf
        }

        
        
        if(log(runif(1))<alpha){ # accept, update current values
          accept <- TRUE
          X.curr[[j]]<-X.cand
          Vinv.curr[[j]]<-Vinv.cand
          d.curr[[j]]<-d.cand
          log.det.curr[[j]] <- log.det.cand
          nbasis.curr[j] <- nbasis.curr[j] + 1
          nvars.curr[j] <- nvars.curr[j] + length(vars.cand)
          nint.curr[nbasis.curr[j],j]<-nint.cand
          
          if(max.int == 1){
            knots.curr[nbasis.curr[j],1,j] <- knots.cand
            signs.curr[nbasis.curr[j],1,j] <- signs.cand
            vars.curr[nbasis.curr[j],1,j] <- vars.cand
            levels.curr[nbasis.curr[j],1,,j] <- levels.cand
          } else {
            knots.curr[nbasis.curr[j],1:nint.cand,j] <- knots.cand
            signs.curr[nbasis.curr[j],1:nint.cand,j] <- signs.cand
            vars.curr[nbasis.curr[j],1:nint.cand,j] <- vars.cand
            levels.curr[nbasis.curr[j],1:nint.cand,,j] <- levels.cand
          }
          
          # update probabilities if weighted
          if(birth.prop != "uniform"){
            p.int[nint.cand] <- p.int[nint.cand] + 1
            p.vars[vars.cand] <- p.vars[vars.cand] + 1
          }
        }
        
      } else if(move.type=='death'){
        count.temp <- 1
        valid.proposal.temp <- TRUE
        
        tokill<-sample(nbasis.curr[j],1) # which basis function we will delete
        X.cand<-X.curr[[j]][,-(tokill+1),drop=FALSE] # +1 to skip the intercept
        if(max.int==1){
          vars.kill <- vars.curr[tokill,1,j]
        } else {
          vars.kill <- vars.curr[tokill, 1:nint.curr[tokill,j],j]
        }
        
        # assume that this is well-conditioned
        Vinv.cand <- Vinv.curr[[j]][-(tokill+1), -(tokill+1), drop=FALSE]   
        Xy <- crossprod(X.cand,z[,j])
        chol <- chol.default(Vinv.cand)
        bhat.cand <-chol2inv(chol)%*%Xy
        d.cand <- sum(Xy*bhat.cand)/2
        log.det.cand <- sum(log(diag(chol)))
        
        llik.alpha <- -.5*log(1/tau2.curr[j]) - log.det.cand + log.det.curr[[j]] + d.cand - d.curr[[j]]
        # llam <- -log(lam.curr)+log(nvars.curr[j]) # nbasis
        llam <- -log(lam.curr[j])+log(nbasis.curr[j]) # nbasis
        
        # Birth proposal distr
        if(birth.prop=="uniform"){   # uniform proposal probability
          ll.birth.vars <- ll.birth.int <- 0 # cancels when proposing from prior
          if(int.prior != "uniform"){
            ll.birth.int <- -log(p.int.prior[nint.curr[tokill,j]]) - log(max.int)  # -prior + uniform proposal
          }
        } else {
          p.vars.temp <- p.vars
          p.vars.temp[vars.kill] <- p.vars.temp[vars.kill] - 1   # remove variables 
          res <- compute.prob.set(vars.kill, p.vars.temp)   # probability of choosing these variables again
          
          p.int.temp <- p.int
          p.int.temp[nint.curr[tokill,j]] <- p.int.temp[nint.curr[tokill,j]] - 1   # remove one interaction
          
          ll.birth.vars <- -(log(1/choose(p,nint.curr[tokill,j])) - log(res))   # uniform prior  - weighted proposal
          ll.birth.int <- -(log(p.int.prior[nint.curr[tokill,j]]) - log(sum(p.int.temp)/p.int.temp[nint[tokill]]))    # uniform prior  - weighted proposal
        }
        
        if(any(!ordinal[vars.kill])){
          k.ordinal <- which(!ordinal[vars.kill])
          for(k in k.ordinal){
            cat.poss <- vars.list[[vars.kill[k]]]  #sort(unique(tX[vars.cand[k], ]))
            prob.temp <- rep(1, length(cat.poss)-1)
            prob.temp <- prob.temp / sum(prob.temp)
            
            u <- sum(!is.na(levels.curr[tokill,k,,j]))
            prob.prior <- 1/sqrt(choose(length(cat.poss), 1:(length(cat.poss)-1)))
            prob.prior <- prob.prior / sum(prob.prior)
            # prob.prior = prob.temp
            ll.birth.vars <- ll.birth.vars + log(prob.temp[u]) - log(prob.prior[u])   # prior and proposal for number of levels
          }
        }
        
        alpha <- llik.alpha + llam + ll.birth.vars + ll.birth.int  # only lambda when knots, int, and vars are uniform
      
        
        if(log(runif(1))<alpha){ # accept, update
          accept <- TRUE
          X.curr[[j]]<-X.cand
          # bhat<-bhat.cand
          Vinv.curr[[j]]<-Vinv.cand
          d.curr[[j]]<-d.cand
          log.det.curr[[j]] <- log.det.cand
          nbasis.curr[j]<-nbasis.curr[j]-1
          nvars.curr[j] <- nvars.curr[j] - length(vars.kill)
          if(nbasis.curr[j]==0){
            nint.curr[,j]<-rep(NA, length(nint.curr[,j]))
            if(max.int==1){
              knots.curr[,1,j]<-signs.curr[,1,j]<-vars.curr[,1,j]<-rep(NA, length(knots.curr[,1,j]))
            } else {
              knots.curr[,,j]<-signs.curr[,,j]<-vars.curr[,,j]<-matrix(NA, nrow(knots.curr[,,j]), ncol(knots.curr[,,j]))
            }
            # update probabilities if weighted
            if(birth.prop != "uniform"){
              p.int <- rep(w1, max.int)  # NO BASIS FUNCTIONS == UNIFORM PROBABILITY
              p.vars <- rep(w2, p)
            }
          } else{
            nint.curr[tokill,j] <- NA
            nint.curr[1:(nbasis.curr[j]+1),j] <- nint.curr[c(setdiff(1:(nbasis.curr[j]+1),tokill),tokill),j]
            if(max.int==1){
              knots.curr[tokill,1,j] <- NA
              knots.curr[1:(nbasis.curr[j]+1),1,j] <- knots.curr[c(setdiff(1:(nbasis.curr[j]+1),tokill),tokill),1,j]
              signs.curr[tokill,1,j] <- NA
              signs.curr[1:(nbasis.curr[j]+1),1,j] <- signs.curr[c(setdiff(1:(nbasis.curr[j]+1),tokill),tokill),1,j]
              vars.curr[tokill,1,j] <- NA
              vars.curr[1:(nbasis.curr[j]+1),1,j] <- vars.curr[c(setdiff(1:(nbasis.curr[j]+1),tokill),tokill),1,j]
              levels.curr[tokill,1,,j] <- NA  
              levels.curr[1:(nbasis.curr[j]+1),1,,j] <- levels.curr[c(setdiff(1:(nbasis.curr[j]+1), tokill), tokill),1,,j]
              
            } else {
              knots.curr[tokill,,j] <- NA
              knots.curr[1:(nbasis.curr[j]+1),,j] <- knots.curr[c(setdiff(1:(nbasis.curr[j]+1),tokill),tokill),,j]
              signs.curr[tokill,,j] <- NA
              signs.curr[1:(nbasis.curr[j]+1),,j] <- signs.curr[c(setdiff(1:(nbasis.curr[j]+1),tokill),tokill),,j]
              vars.curr[tokill,,j] <- NA
              vars.curr[1:(nbasis.curr[j]+1),,j] <- vars.curr[c(setdiff(1:(nbasis.curr[j]+1),tokill),tokill),,j]
              levels.curr[tokill,,,j] <- NA  
              levels.curr[1:(nbasis.curr[j]+1),,,j] <- levels.curr[c(setdiff(1:(nbasis.curr[j]+1), tokill), tokill),,,j]
            }
            
            # update probabilities if weighted (before deleting)
            if(birth.prop != "uniform"){
              p.int <- p.int.temp
              p.vars <- p.vars.temp
            }
          }
          accept <- TRUE
        }
        
      } else {
        valid.proposal.temp <- FALSE
        count.temp <- 0
        tochange<-sample(nbasis.curr[j],1) # which basis function we will change
        
        while(!valid.proposal.temp & count.temp < max.prop){
          count.temp <- count.temp + 1
          if(max.int==1){
            tochange2 <- 1
            knots.cand<-knots.curr[tochange,tochange2,j] # copy
            vars.to.change <- vars.curr[tochange,tochange2,j]  # variable (single!) in change basis function
            levels.cand <- levels.curr[tochange,1,,j]
            levels.cand <- matrix(levels.cand, nrow=1)
            signs.cand<-signs.curr[tochange,tochange2,j] # copy
          } else {
            tochange2<-sample(nint.curr[tochange,j],1) # which element in the basis function tensor product we will change
            knots.cand<-knots.curr[tochange,1:nint.curr[tochange,j],j] # copy
            vars.to.change <- vars.curr[tochange,tochange2,j]  # variable (single!) in change basis function
            levels.cand <- levels.curr[tochange,,,j]
            levels.cand <- levels.cand[1:nint.curr[tochange,j],,drop=FALSE]
            signs.cand<-signs.curr[tochange,1:nint.curr[tochange,j],j] # copy
          }
          
          # Sample signs for all cases
          signs.cand[tochange2] <- sample(c(-1,1), 1)
          
          # knots.cand[tochange2]<-runif(1) # change one element
          if(ordinal[vars.to.change] & knot.prop == "discrete" ){
            knots.cand[tochange2]<- sample(vars.list[[vars.to.change]], 1)
            # ll.prop.knot <- 0  # proposal probabilities cancel in change
          } else if (ordinal[vars.to.change] & knot.prop != "discrete"){
            knots.cand[tochange2]<-runif(1) # change one element
            
          } else if (!ordinal[vars.to.change]){
            # p.old <- sum(!is.na(levels.cand[tochange2,]))
            # if(max.int==1){
            #   levels.cand[tochange2,] <- NA  # overwrite levels we will change
            # } else {
            #   levels.cand[tochange2,] <- NA  # overwrite levels we will change
            # }
            # 
            # cat.poss <- vars.list[[vars.to.change]]  #sort(unique(tX[tochange, ]))
            # if(length(cat.poss) <= 3){
            #   prob.temp <- rep(1, length(cat.poss)-1)
            # } else {
            #   prob.temp <- 1/sqrt(choose(length(cat.poss), 1:(length(cat.poss)-1)))   # wtd prob SQRT
            # }
            # u <- sample(1:(length(cat.poss)-1), 1, prob=prob.temp)  # at least 1 and not all
            cat.poss <- vars.list[[vars.to.change]]  #sort(unique(tX[tochange, ]))
            u = sum(!is.na(levels.cand[tochange2,]))
            levels.temp <- sample(cat.poss, u, replace=FALSE)
            levels.cand[tochange2,1:u] <- levels.temp
          }
          
          if(max.int==1){
            vars.basis <- vars.curr[tochange,1,j]
          } else {
            vars.basis <- vars.curr[tochange,1:nint.curr[tochange,j],j]
          }
          basis<-makeBasis.2(signs.cand,vars.basis,knots.cand,Xt,ordinal=ordinal,
                             levels=levels.cand[1:nint.curr[tochange,j],,drop=FALSE]) # make the new basis function
          X.cand<-X.curr[[j]]
          X.cand[,tochange+1]<-basis # replace with our new basis function (+1 for intercept)
          # Vinv.cand<-crossprod(X.cand)+diag(nbasis+1)/tau2
          Vinv.cand <- Vinv.curr[[j]]
          x.new <- c(crossprod(basis,X.cand))
          Vinv.cand[tochange+1,] <- x.new
          Vinv.cand[,tochange+1] <- x.new
          Vinv.cand[tochange+1,tochange+1] <- sum(basis^2) + 1/tau2.curr[j]
          chol.cand <- chol.default(Vinv.cand)
          diag.cand <- diag(chol.cand)
          cond.cand <- (min(diag.cand)/max(diag.cand))^2  # approx condition number
          
          nonzero <- sum(abs(basis) > eps) >= min.nonzero
          nonzero2 <- cond.cand > eps   # computationally nonsingular
          
          valid.proposal.temp <- nonzero & nonzero2
        }
        
        if(valid.proposal.temp){  
          log.det.cand <- sum(log(diag.cand))   # 1/2 log determinant of Vinv
          Xy <- crossprod(X.cand,z[,j])
          bhat.cand <- chol2inv(chol.cand)%*%Xy
          d.cand <- sum(Xy*bhat.cand)/2
          
          llik.alpha <- -log.det.cand + log.det.curr[[j]] + d.cand - d.curr[[j]]
          alpha <- llik.alpha    # the proposals and priors all cancel since there is no dimension change
        } else {
          alpha <-  -Inf
        }
        
        
        if(is.na(alpha)){
          cat("NA alpha, change \n")
          cat("max.int=", max.int, "\n")
          cat("i=", i, "\n")
          cat("j=", j, "\n")
          cat("range z=", range(z[,j]), "\n")
          cat("nbasis.curr[j]=", nbasis.curr[j], "\n")
          cat("lam.curr[j]=", lam.curr[j], "\n")
          cat("d.cand=", d.cand, "\n")
          cat("d.curr=", d.curr[[j]], "\n")
          cat("tau2.curr=", tau2.curr[j], "\n")
          cat("log.det.cand=", log.det.cand, "\n")
          cat("log.det.curr=", log.det.curr[[j]], "\n")
          cat("llam=", llam, "\n")
          cat("ll.birth.vars=", ll.birth.vars, "\n")
          cat("ll.birth.int=", ll.birth.int, "\n")
        }
        
        
        if(log(runif(1))<alpha){ # accept, update
          X.curr[[j]]<-X.cand
          # bhat<-bhat.cand
          Vinv.curr[[j]]<-Vinv.cand
          d.curr[[j]]<-d.cand
          log.det.curr[[j]] <- log.det.cand
          if(max.int == 1){
            knots.curr[tochange,1,j]<-knots.cand
            signs.curr[tochange,1,j]<-signs.cand
            levels.curr[tochange,1,,j]<-levels.cand 
          } else {
            knots.curr[tochange,1:nint.curr[tochange,j],j]<-knots.cand
            signs.curr[tochange,1:nint.curr[tochange,j],j]<-signs.cand
            levels.curr[tochange,1:nint.curr[tochange,j],,j]<-levels.cand
          }
          accept <- TRUE
        }
        
      }
      
      # Gibbs steps
      # S<-solve(crossprod(X.curr)/s2+diag(nbasis+1)/tau2) # covariance matrix for beta update
      # S <- solve(Vinv.curr[[j]])
      # S2 <- chol.default(S)
      S2 <- my_symm_square_root(Vinv.curr[[j]], invert=TRUE)
      mean.temp <- solve(Vinv.curr[[j]],crossprod(X.curr[[j]],z[,j]))
      beta.curr[1:(nbasis.curr[j]+1),j]<- mean.temp + crossprod(S2, rnorm(length(mean.temp))) # update beta
      
      res <- sample.lambda.mh(lam.curr[j], nbasis.curr[j], 
                              max.basis, h1, h2)
      lam.curr[j] <- res$lambda  
      accept.lambda <- res$accept
      
      # sample tau2
      ssb <- sum(beta.curr[1:(nbasis.curr[j]+1),j]^2)/2
      p0 <- nbasis.curr[j]+1
      tau2.curr[j] <- 1/rgamma(1, p0/2 + 1, ssb + tau2.0)
      
      if(accept){
        count[j,move.type] <- count[j,move.type] + 1
        # mu[,j] <- X.curr[[j]]%*%beta.curr[1:ncol(X.curr[[j]]),j]
      }
      
      valid.prop.curr <- valid.proposal.temp
      n.prop.curr <- count.temp
    }
    
    # # update lambda
    # # res <- sample.lambda.mh.multiple(lam.curr, nvars.curr, 
    # #                                  max.vars, h1, h2)
    # res <- sample.lambda.mh.multiple(lam.curr, nbasis.curr, 
    #                                  max.basis, h1, h2)
    # lam.curr <- res$lambda  
    # accept.lambda <- res$accept
    
    # Writeout
    if(verbose & i %% print.interval == 0){
      cat("MCMC iteration ", i, "\n")
      cat("Number of basis functions:\n")
      print(nbasis.curr)
      cat("Lambda est:", round(lam.curr, 3), "\n")
      cat("Proposal acceptance rate:", round(sum(count) / i / d, 3), "\n")
      cat("Lambda acceptance rate:",round(accept.lambda, 3), "\n")
      cat("######################\n")
    }
    
    
    
    # Save if i is large enough
    if(i > nburn & (i %% nthin==0)){
      jj <- which(i.save == i)
      nint[jj,,] <- nint.curr
      beta[jj,,] <- beta.curr
      
      knots[jj,,,] <- knots.curr
      signs[jj,,,] <- signs.curr
      vars[jj,,,] <- vars.curr
      
      levels[jj,,,,] <- levels.curr
      
      nbasis[jj,] <- nbasis.curr
      nvars[jj,] <- nvars.curr
      lam[jj,] <- lam.curr
      tau2[jj,] <- tau2.curr
      
      valid.prop[jj] <- valid.prop.curr
      n.prop[jj] <- n.prop.curr
      
      if(writeout & is.character(writedir)){
        write.table(matrix(c(nint.curr, nrow=1)),
                    file=file.path(writedir, "nint.txt"), append=TRUE, col.names = FALSE)
        write.table(matrix(c(beta.curr, nrow=1)),
                    file=file.path(writedir, "beta.txt"), append=TRUE, col.names = FALSE)
        write.table(matrix(c(knots.curr, nrow=1)),
                    file=file.path(writedir, "knots.txt"), append=TRUE, col.names = FALSE)
        write.table(matrix(c(signs.curr, nrow=1)),
                    file=file.path(writedir, "signs.txt"), append=TRUE, col.names = FALSE)
        write.table(matrix(c(levels.curr, nrow=1)),
                    file=file.path(writedir, "levels.txt"), append=TRUE, col.names = FALSE)
        write.table(matrix(c(vars.curr, nrow=1)),
                    file=file.path(writedir, "vars.txt"), append=TRUE, col.names = FALSE)
        save(z,X.curr, Vinv.curr, d.curr, log.det.curr,
             file=file.path(writedir, "state_variables.rda"))
      }
      
      
      
      if(res$accept){
        basis.samples[jj,j,res$move.type] <- 1
      }
    }
  }
  
  return(list(X=X.curr,b=0,count=count,knots=knots,signs=signs,vars=vars,
              nint=nint,nbasis=nbasis,beta=beta,lam=lam,
              accept.lambda=accept.lambda,basis.samples=basis.samples,
              p.int=p.int, p.vars=p.vars, ordinal=ordinal, levels=levels, 
              tau2=tau2, classes=classes,
              nburn=nburn, nmcmc=nmcmc, nthin=nthin,
              V=Vinv.curr, d=d.curr, log.det=log.det.curr,
              z=z, max.int=max.int,
              valid.prop=valid.prop, n.prop=n.prop,
              nvars=nvars, h1=h1, h2=h2, 
              ranges=ranges.X, mins=mins.X
              # inf.flags=inf.flags
  ))
}





# Initialize model with single interaction fit
init.mod.1 <- function(X, y, 
                             tau2=1e3, 
                             max.basis = ncol(X)*2,
                             h1=4, h2=4*(length(unique(y))-1)/length(y),
                             nmcmc=5e3, nthin=10, nburn=3e3,
                             birth.prop="uniform", 
                             ordinal=NULL)
{
  mod1 <-  fit.cbass(X, 
                     y,
                     ordinal=ordinal,
                     max.int=1,
                     init1=FALSE,
                     nmcmc=nmcmc, 
                     nburn=nburn, 
                     nthin=nthin,
                     tau2=tau2,
                     max.basis = max.basis,
                     h1=h1, h2=h2)
  
  # names(mod)
  
  nn <- nrow(mod1$nbasis)
  return(list(z=mod1$z, 
              X=mod1$X,
              Vinv=mod1$V,
              d=mod1$d,
              log.det=mod1$log.det,
              nbasis=mod1$nbasis[nn,],
              nvars=mod1$nvars[nn,],
              nint=mod1$nint[nn,,],
              knots=mod1$knots[nn,,,], 
              signs=mod1$signs[nn,,,], 
              vars=mod1$vars[nn,,,],
              levels=mod1$levels[nn,,,,],
              beta=mod1$beta[nn,,],
              lam=mod1$lam[nn,],
              tau2=mod1$tau2[nn,],
              # h1=mean(tail(mod1$lam, round(nn/2)))
              h1=mean(mod1$nbasis[tail(1:nn, round(nn/2)),])        
  ))
}




