library(dirichletprocess)
PriorDraw.mvnormal2.my<-function (mdObj, n = 1){
  priorParameters <- mdObj$priorParameters
  
  sig <- array(sapply(seq_len(n), function(x) solve(rWishart(1,
                                                             priorParameters$nu0, solve(priorParameters$phi0))[, ,
                                                                                                               1])),dim=c(2,2,n))
  mu <- array(t(mvtnorm::rmvnorm(n,
                                 priorParameters$mu0, priorParameters$sigma0)),dim=c(1,2,n))
  
  theta <- list(mu = mu, sig = sig)
  return(theta)
}
#
ClusterParameterUpdate.nonconjugate.my<-function (dpObj){
  y <- dpObj$data
  numLabels <- dpObj$numberClusters
  clusterLabels <- dpObj$clusterLabels
  clusterParams <- dpObj$clusterParameters
  mdobj <- dpObj$mixingDistribution
  mhDraws <- dpObj$mhDraws
  accept_ratio <- numeric(numLabels)
  for (i in 1:numLabels) {
    pts <- y[which(clusterLabels == i), , drop = FALSE]
    
    parameter_samples <- PosteriorDraw(mdobj, pts, mhDraws)
    for (j in seq_along(clusterParams)) {
      clusterParams[[j]][, , i] <- parameter_samples[[j]][,
                                                          , mhDraws]
    }
    accept_ratio[i] <- length(unique(parameter_samples[[1]]))/mhDraws
  }
  dpObj$clusterParameters <- clusterParams
  return(dpObj)
}

DirichletProcessGaussian.my<-function (y, g0Priors = c(0, 1, 1, 1), 
                                       alphaPriors = c(10, 0.5)){
  mdobj <- GaussianMixtureCreate(g0Priors)
  dpobj <- DirichletProcessCreate(y, mdobj, alphaPriors)
  dpobj <- Initialise(dpobj)
  return(dpobj)
}

Initialise.nonconjugate.my<-function (dpObj, posterior = TRUE, 
                                      updatePrior = FALSE, m = 3, verbose = TRUE){
  dpObj$clusterLabels <- rep(1, dpObj$n)
  dpObj$numberClusters <- 1
  dpObj$pointsPerCluster <- dpObj$n
  if (posterior) {
    post_draws <- PosteriorDraw(dpObj$mixingDistribution, 
                                dpObj$data, 1000)
    if (verbose) 
      cat(paste("Accept Ratio: ", length(unique(c(post_draws[[1]])))/1000, 
                "\n"))
    dpObj$clusterParameters <- list(post_draws[[1]][, , 1000, 
                                                    drop = FALSE], post_draws[[2]][, , 1000, drop = FALSE])
  }
  else {
    dpObj$clusterParameters <- PriorDraw(dpObj$mixingDistribution, 
                                         1)
  }
  dpObj$nu <- post_draws[[3]]
  dpObj$phi <- post_draws[[4]]
  dpObj$m <- m
  return(dpObj)
}

GaussianMixtureCreate.my<-function (priorParameters = c(0, 1e3, 1, 1))
{
  mdobj <- MixingDistribution("normal", priorParameters,
                              "nonconjugate")
  return(mdobj)
}

PriorDraw.normal.my<-function (mdObj, n = 1){
  priorParameters <- mdObj$priorParameters
  mu.0 <- priorParameters[1]
  tau2.0 <- priorParameters[2]
  nu.0 <- priorParameters[3]
  sigma2.0 <- priorParameters[4]
  
  sigma2.samp <- 1/rgamma(n, nu.0/2, sigma2.0 * nu.0/2)
  mu.samp <- rnorm(n, mu.0, sqrt(tau2.0))
  
  theta <- list(array(mu.samp, dim = c(1, 1, n)), array(sqrt(sigma2.samp),
                                                        dim = c(1, 1, n)))
  return(theta)
}

PosteriorDraw.normal.my<-function (mdObj, x, n = 1, ...){
  n.x <- length(x)
  mu.0 <- mdObj$priorParameters[1]
  tau2.0 <- mdObj$priorParameters[2]
  nu.0 <- mdObj$priorParameters[3]
  sigma2.0 <- mdObj$priorParameters[4]
  
  muSamples <- array(dim = c(1,1,n))
  sigSamples <- array(dim = c(1,1,n))
  
  mean.x <- mean(x)
  
  muSamp <- 0
  precSamp <- 1/1e6
  for (i in seq_len(n)) {
    mu.n <- (mu.0/tau2.0+n.x*mean.x*precSamp)/(1/tau2.0+n.x*precSamp)
    tau2.n <- 1/(1/tau2.0+n.x*precSamp)
    muSamp <- rnorm(1,mean = mu.n,sqrt(tau2.n))
    muSamples[,,i] <- muSamp
    
    nu.n <- nu.0+n.x
    sigma2.n <- (nu.0*sigma2.0+sum((x-muSamp)^2))/nu.n
    precSamp <- rgamma(1,nu.n/2,nu.n*sigma2.n/2)
    sigSamples[,,i] <- sqrt(1/precSamp)
    
  }
  return(list(mu = muSamples, sig = sigSamples, nu = nu.n, sigma = sigma2.n))
}

PosteriorDraw.mvnormal2.my <- function (mdObj, x, n = 1, ...){
  if (!is.matrix(x)) {
    x <- matrix(x, ncol = length(x))
  }
  phi0 <- mdObj$priorParameters$phi0
  mu0 <- mdObj$priorParameters$mu0
  sigma0 <- mdObj$priorParameters$sigma0
  muSamples <- array(dim = c(dim(mu0), n))
  sigSamples <- array(dim = c(dim(phi0), n))
  muSamp <- matrix(rep_len(0, ncol(mu0)), ncol = ncol(mu0))
  for (i in seq_len(n)) {
    nuN <- nrow(x) + mdObj$priorParameters$nu0
    phiN <- phi0 + 
      t(x-matrix(1,nrow = nrow(x))%*%muSamp)%*%(x-matrix(1,nrow = nrow(x))%*%muSamp)
    sigSamp <- solve(rWishart(1, nuN, solve(phiN))[, , 1])
    sigN <- solve(solve(sigma0) + nrow(x) * solve(sigSamp))
    muN <- sigN %*% (nrow(x) * solve(sigSamp) %*% colMeans(x) + 
                       solve(sigma0) %*% c(mu0))
    muSamp <- mvtnorm::rmvnorm(1, muN, sigN)
    muSamples[, , i] <- muSamp
    sigSamples[, , i] <- sigSamp
  }
  return(list(mu = muSamples, sig = sigSamples, nu = nuN, phi = phiN))
}

UpdateAlpha.default.my<-function (dpobj){
  newAlpha <- dirichletprocess:::update_concentration(dpobj$alpha, dpobj$n, dpobj$numberClusters,
                                                      dpobj$alphaPriorParameters)
  dpobj$alpha <- newAlpha$alpha
  dpobj$alphapost <- c(newAlpha$postpar1,newAlpha$postpar2)
  return(dpobj)
}

update_concentration.my<-function (oldParam, n, nParams, priorParameters){
  x <- rbeta(1, oldParam + 1, n)
  pi1 <- priorParameters[1] + nParams - 1
  pi2 <- n * (priorParameters[2] - log(x))
  pi1 <- pi1/(pi1 + pi2)
  postParams1 <- priorParameters[1] + nParams
  postParams2 <- priorParameters[2] - log(x)
  if (runif(1) > pi1) {
    postParams1 <- postParams1 - 1
  }
  new_alpha <- rgamma(1, postParams1, postParams2)
  return(list(alpha = new_alpha, postpar1 =postParams1,postpar2 =postParams2))
}


assignInNamespace("PriorDraw.mvnormal2",PriorDraw.mvnormal2.my,"dirichletprocess")
assignInNamespace("ClusterParameterUpdate.nonconjugate",ClusterParameterUpdate.nonconjugate.my,
                  "dirichletprocess")

assignInNamespace("DirichletProcessGaussian",DirichletProcessGaussian.my,
                  "dirichletprocess")

assignInNamespace("GaussianMixtureCreate",GaussianMixtureCreate.my,
                  "dirichletprocess")
assignInNamespace("PriorDraw.normal",PriorDraw.normal.my,
                  "dirichletprocess")

assignInNamespace("Initialise.nonconjugate",Initialise.nonconjugate.my,
                  "dirichletprocess")

assignInNamespace("PosteriorDraw.mvnormal2",PosteriorDraw.mvnormal2.my,
                  "dirichletprocess")
assignInNamespace("PosteriorDraw.normal",PosteriorDraw.normal.my,
                  "dirichletprocess")

assignInNamespace("update_concentration",update_concentration.my,
                  "dirichletprocess")
assignInNamespace("UpdateAlpha.default",UpdateAlpha.default.my,
                  "dirichletprocess")
