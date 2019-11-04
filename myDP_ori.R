rm(list=ls())
library(mvtnorm)
library(LaplacesDemon)
library(mixtools)
# reinstall dirichletprocess package if it is needed to be used.
# Initialization:
FlatMixtureCreate<-function(priorParameters = list(beta=array(c(4000,-7000),
                                                              dim = c(2,1)),
                                                   nu0=2,sigma02=100^2,
                                                   D=matrix(c(4*1e6,0,0,1e8),nrow = 2,ncol = 2)))
{
  mdobj <- MixingDistribution_my("Flat", priorParameters, 
                                 "nonconjugate",
                                 mhStepSize = c(100, 5))
  return(mdobj)  
}


MixingDistribution_my<-function (distribution, priorParameters, conjugate, mhStepSize = NULL, 
          hyperPriorParameters = NULL) 
{
  mdObj <- list(distribution = distribution, priorParameters = priorParameters, 
                conjugate = conjugate, mhStepSize = mhStepSize, hyperPriorParameters = hyperPriorParameters)
  class(mdObj) <- append(class(mdObj), c(distribution, conjugate))
  return(mdObj)
}

DirichletProcessCreate_my<-function (x, y, mdObject=FlatMixtureCreate(), 
                                     alphaPriorParameters = c(10, 0.5), mhDraws = 30) 
{
  # y is response, x is data
  if (!is.matrix(y)) {
    y <- matrix(y, ncol = 1)
  }
  dpObj <- list(data = x, response = y, mixingDistribution = mdObject, n = dim(x)[1], 
                alphaPriorParameters = alphaPriorParameters, alpha = rgamma(1, 
                                                                            alphaPriorParameters[1], alphaPriorParameters[2]), 
                mhDraws = mhDraws)
  class(dpObj) <- append(class(dpObj), c("dirichletprocess", 
                                         class(mdObject)[-1]))
  return(dpObj)
}

# PriorDraw
PriorDraw<-function (dpObj, n=1){
  UseMethod("PriorDraw", dpObj)
}

PriorDraw.OLS <- function(dpObj, n=1){
  n.obs <- dpObj$n
  mdObj <- dpObj$mixingDistribution
  priorParameters <- mdObj$priorParameters
  
  x.mat <- dpObj$data
  x.mat <- solve(t(x.mat)%*%x.mat)
  
  beta.sig <- priorParameters$beta.sig
  
  beta.sig <- array(rgamma(n, shape = beta.sig[1], scale = beta.sig[2]),dim = c(1,1,n))
  #using g-prior here, initialized as OLS result:
  n.per.cluster<-ceiling(n.obs/n)
  beta <- array(mvtnorm:::rmvnorm(n,mean = priorParameters$beta,
                        sigma = beta.sig[,,n]^2*n.per.cluster*x.mat),dim = c(2,1,n))
  
  theta <- list(beta, beta.sig)
  return(theta)
}

PriorDraw.Wishart <- function(dpObj, n=1){
  n.obs <- dpObj$n
  mdObj <- dpObj$mixingDistribution
  priorParameters <- mdObj$priorParameters
  
  beta.sig <- priorParameters$beta.sig
  
  beta.sig <- array(rgamma(n, shape = beta.sig[1], scale = beta.sig[2]),dim = c(1,1,n))
  
  D <- priorParameters$D
  
  D <- rWishart(n,df = 2,Sigma = D)/2
  
  beta <- array(unlist(lapply(seq_len(n), function(i){mvtnorm:::rmvnorm(1,mean = priorParameters$beta,
                                                 sigma = D[,,i])})),dim = c(2,1,n))
  
  theta <- list(beta, beta.sig, D)
  return(theta)
}

PriorDraw.Flat <- function(dpObj, n=1){
  n.obs <- dpObj$n
  mdObj <- dpObj$mixingDistribution
  priorParameters <- mdObj$priorParameters
  
  nu.0 <- priorParameters$nu0
  sigma0.2 <- priorParameters$sigma02
  
  beta.sig <- array(sqrt(1/rgamma(n, nu.0/2,sigma0.2*nu.0/2)),dim = c(1,1,n))
  
  D <- priorParameters$D
  
  beta <- array(mvtnorm:::rmvnorm(n,mean = priorParameters$beta,
                        sigma = D),dim = c(2,1,n))
  
  theta <- list(beta, beta.sig)
  return(theta)
}

Initialise_my <- function(dpObj, posterior = TRUE, m=3, verbose=TRUE){
  UseMethod("Initialise_my", dpObj)
}

# Maybe it's conjugate
Initialise_my.nonconjugate <- function(dpObj, m = 100, verbose = TRUE) {
  dpObj$clusterParameters <- PriorDraw(dpObj, m+1)
  
  pstick <- StickBreaking_my(alpha = dpObj$alpha,N = m)
  pend <- ((1-sum(pstick))>=0)*(1-sum(pstick))
  pstick <- c(pstick,pend)
  
  # labelling:
  label <- c()
  for(j in seq_len(dpObj$n)){
    init.llh <- sapply(seq_len(m+1), 
                       function(i) return(dnorm(dpObj$response[j],
                                                dpObj$data[j,]%*%dpObj$clusterParameters[[1]][,,i],
                                                sqrt(dpObj$clusterParameters[[2]][,,i]))))
    probs <- init.llh*pstick
    if(all(probs==0)){
      probs <- rep(1,m+1)
    }
    label[j] <- sample.int(m+1, 1, prob = probs)
  }
  table.label <- table(label)
  table.length <- length(table.label)
  clusterlabel <- as.numeric(names(table.label))
  names(table.label) <- NULL
  for(j in seq_len(table.length)){
    indj <- which(label==clusterlabel[j])
    label[indj] <- j
  }
  
  pointspercluster <- table.label
  dpObj$clusterParameters[[1]] <- array(dpObj$clusterParameters[[1]][,,clusterlabel],
                                        dim=c(2,1,table.length))
  dpObj$clusterParameters[[2]] <- array(dpObj$clusterParameters[[2]][,,clusterlabel],
                                        dim=c(1,1,table.length))
  dpObj$pointsPerCluster <- pointspercluster
  dpObj$numberClusters <- table.length
  dpObj$clusterLabels <- label
  
  dpObj$m <- 3
  
  return(dpObj)
}

# Fitting:
Fit_my <- function(dpObj, its, alphaInterval= 1,updatePrior = FALSE, progressBar=TRUE) UseMethod("Fit_my", dpObj)

Fit_my.default <- function(dpObj, its, alphaInterval= 1, updatePrior = FALSE, progressBar=TRUE) {
  
  if (progressBar){
    pb <- txtProgressBar(min=0, max=its, width=50, char="-", style=3)
  }
  
  alphaChain <- numeric(its)
  likelihoodChain <- numeric(its)
  weightsChain <- vector("list", length = its)
  clusterParametersChain <- vector("list", length = its)
  priorParametersChain <- vector("list", length = its)
  labelsChain <- vector("list", length = its)
  numLabelChain <- vector("numeric",length = its)
  pointsPerClusterChain <- vector("list",length = its)
  
  for (i in seq_len(its)) {
    
    dpObj <- ClusterParameterUpdate_my(dpObj)
    dpObj <- ClusterComponentUpdate_my(dpObj)
    if(i %% alphaInterval ==0){
      dpObj <- UpdateAlpha_my(dpObj)
    }
    
    alphaChain[i] <- dpObj$alpha
    weightsChain[[i]] <- dpObj$pointsPerCluster / dpObj$n
    clusterParametersChain[[i]] <- dpObj$clusterParameters
    priorParametersChain[[i]] <- dpObj$mixingDistribution$priorParameters
    labelsChain[[i]] <- dpObj$clusterLabels
    numLabelChain[i] <- dpObj$numberClusters
    pointsPerClusterChain[[i]] <- dpObj$pointsPerCluster
    
    likelihoodChain[i] <- sum(log(LikelihoodDP_my(dpObj)))
    
    # if (updatePrior) {
    #   dpObj$mixingDistribution <- PriorParametersUpdate(dpObj$mixingDistribution,
    #                                                     dpObj$clusterParameters)
    # }
    if (progressBar){
      setTxtProgressBar(pb, i)
    }
  }
  
  dpObj$weights <- dpObj$pointsPerCluster / dpObj$n
  dpObj$alphaChain <- alphaChain
  dpObj$likelihoodChain <- likelihoodChain
  dpObj$weightsChain <- weightsChain
  dpObj$clusterParametersChain <- clusterParametersChain
  dpObj$priorParametersChain <- priorParametersChain
  dpObj$labelsChain <- labelsChain
  dpObj$numLabelChain <- numLabelChain
  dpObj$pointsPerClusterChain <- pointsPerClusterChain
  
  if (progressBar) {
    close(pb)
  }
  return(dpObj)
}

# Likelihood:
LikelihoodDP_my <- function(dpobj){
  
  allParameters <- lapply(seq_along(dpobj$clusterParameters),
                          function(i) dpobj$clusterParameters[[i]][,,dpobj$clusterLabels, drop=FALSE])
  
  likelihoodValues <- Likelihood(dpobj, allParameters)
  
  return(likelihoodValues)
  
}

Likelihood<-function (dpObj, theta, idp=NaN){
  UseMethod("Likelihood", dpObj) 
}

Likelihood.default <- function(dpobj, theta, ifp=NaN){
  # No matter Wishart model or OLS model is used, likelihood will not change
  # Likelihood depends only on the linear regression format.
  # That is, only beta and beta.sig will be taken into account.
    
  numLabels <- dpobj$numberClusters
  clusterLabels <- dpobj$clusterLabels
  clusterParams <- dpobj$clusterParameters

  if(any(is.nan(ifp))){# compute likelihood for each cluster
    x.mat <- dpobj$data
    response <- dpobj$response
    
    mean <- c()
    sigma <- c()
    
    for(i in 1:numLabels){
      ind.i <- which(clusterLabels==i)
      x.mat <- dpobj$data[ind.i,]
      
      if(length(ind.i)==1){beta <- theta[[1]][,,ind.i]}else{
        beta <- theta[[1]][,,ind.i][,1] # take one of them  
      }
      sigma.cluster <- theta[[2]][,,ind.i]
      
      mean[ind.i]<-x.mat%*%beta
      sigma[ind.i]<-sigma.cluster
      
      y <- dnorm(response,mean,sigma)
    }
  }else{# compute the likelihood for a given observation in each cluster
    x.mat <- dpobj$data[ifp,]
    response <- dpobj$response[ifp]
    
    beta <- theta[[1]][,,,drop=TRUE]
    sigma.cluster <- theta[[2]][,,,drop=TRUE]
    
    mean <- x.mat%*%beta
    y <- dnorm(response, mean, sigma.cluster)
  }  
  
  return(y)
}

# Label Update:
ClusterComponentUpdate_my<-function (dpObj){
  UseMethod("ClusterComponentUpdate_my", dpObj)
}

ClusterComponentUpdate_my.nonconjugate <- function(dpObj) {
  
  n <- dpObj$n
  alpha <- dpObj$alpha
  
  clusterLabels <- dpObj$clusterLabels
  clusterParams <- dpObj$clusterParameters
  numLabels <- dpObj$numberClusters
  
  mdObj <- dpObj$mixingDistribution
  m <- dpObj$m
  
  pointsPerCluster <- dpObj$pointsPerCluster
  
  aux <- vector("list", length(clusterParams))
  
  for (i in seq_len(n)) {
    
    currentLabel <- clusterLabels[i]
    
    pointsPerCluster[currentLabel] <- pointsPerCluster[currentLabel] - 1
    
    if (pointsPerCluster[currentLabel] == 0) {
      
      priorDraws <- PriorDraw(dpObj, m - 1)
      
      for (j in seq_along(priorDraws)) {
        aux[[j]] <- array(c(clusterParams[[j]][, , currentLabel], priorDraws[[j]]),
                          dim = c(dim(priorDraws[[j]])[1:2], m))
      }
    } else {
      aux <- PriorDraw(dpObj, m)
    }
    
    probs <- c(
      pointsPerCluster * Likelihood(dpObj, clusterParams, i)[1,],
      (alpha/m) * Likelihood(dpObj, aux, i))
    
    if (any(is.nan(probs))) {
      probs[is.nan(probs)] <- 0
    }
    
    
    probs[is.na(probs)] <- 0
    
    
    if (any(is.infinite(probs))) {
      probs[is.infinite(probs)] <- 1
      probs[-is.infinite(probs)] <- 0
    }
    
    if (all(probs == 0)) {
      probs <- rep_len(1, length(probs))
    }
    newLabel <- sample.int(numLabels + m, 1, prob = probs)
    
    dpObj$pointsPerCluster <- pointsPerCluster
      
    dpObj <- ClusterLabelChange_my(dpObj, i, newLabel, currentLabel, aux)
    
    pointsPerCluster <- dpObj$pointsPerCluster
    clusterLabels <- dpObj$clusterLabels
    clusterParams <- dpObj$clusterParameters
    numLabels <- dpObj$numberClusters
    
  }
  
  dpObj$pointsPerCluster <- pointsPerCluster
  dpObj$clusterLabels <- clusterLabels
  dpObj$clusterParameters <- clusterParams
  dpObj$numberClusters <- numLabels
  return(dpObj)
}

ClusterLabelChange_my <- function(dpObj, x, newLabel, currentLabel, ...){
  UseMethod("ClusterLabelChange_my", dpObj)
}

ClusterLabelChange_my.nonconjugate <- function(dpObj, i, newLabel, currentLabel, aux) {
  
  pointsPerCluster <- dpObj$pointsPerCluster
  clusterLabels <- dpObj$clusterLabels
  clusterParams <- dpObj$clusterParameters
  numLabels <- dpObj$numberClusters
  # mdObj <- dpObj$mixingDistribution
  
  if (newLabel <= numLabels) {
    pointsPerCluster[newLabel] <- pointsPerCluster[newLabel] + 1
    clusterLabels[i] <- newLabel
    
    if (pointsPerCluster[currentLabel] == 0) {
      # print('B') Removing the Empty Cluster ###
      numLabels <- numLabels - 1
      pointsPerCluster <- pointsPerCluster[-currentLabel]
      # clusterParams <- clusterParams[-currentLabel, ,drop=FALSE]
      clusterParams <- lapply(clusterParams, function(x) x[, , -currentLabel,
                                                           drop = FALSE])
      
      inds <- clusterLabels > currentLabel
      clusterLabels[inds] <- clusterLabels[inds] - 1
    }
  } else {
    
    if (pointsPerCluster[currentLabel] == 0) {
      # print('C') clusterParams[currentLabel, ] = aux[newLabel-numLabels, ]
      
      for (j in seq_along(clusterParams)) {
        clusterParams[[j]][, , currentLabel] <- aux[[j]][, , newLabel - numLabels]
      }
      pointsPerCluster[currentLabel] <- pointsPerCluster[currentLabel] + 1
      
    } else {
      # print('D')
      clusterLabels[i] <- numLabels + 1
      pointsPerCluster <- c(pointsPerCluster, 1)
      # clusterParams = rbind(clusterParams, aux[newLabel-numLabels, ])
      
      for (j in seq_along(clusterParams)) {
        clusterParams[[j]] <- array(c(clusterParams[[j]],
                                      aux[[j]][, , newLabel - numLabels]),
                                    dim = c(dim(clusterParams[[j]])[1:2],
                                            dim(clusterParams[[j]])[3] + 1))
      }
      
      numLabels <- numLabels + 1
    }
  }
  
  dpObj$pointsPerCluster <- pointsPerCluster
  dpObj$clusterLabels <- clusterLabels
  dpObj$clusterParameters <- clusterParams
  dpObj$numberClusters <- numLabels
  return(dpObj)
}

# Parameter Update:
ClusterParameterUpdate_my<-function (dpObj){
  UseMethod("ClusterParameterUpdate_my", dpObj)
}

ClusterParameterUpdate_my.nonconjugate <- function(dpObj) {
  
  numLabels <- dpObj$numberClusters
  
  clusterLabels <- dpObj$clusterLabels
  clusterParams <- dpObj$clusterParameters
  
  mdobj <- dpObj$mixingDistribution
  mhDraws <- dpObj$mhDraws
  
  accept_ratio <- numeric(numLabels)
  
  start_pos <- PriorDraw(dpObj)
  
  for (i in 1:numLabels) {
    ind <- which(clusterLabels == i)
    # annotation is for mh:
    
    # for (j in seq_along(clusterParams)) {
    #   start_pos[[j]] <- clusterParams[[j]][, , i, drop = FALSE]
    # }
    
    #parameter_samples <- PosteriorDraw_my(dpObj, ind, start_pos = start_pos, n = mhDraws)
    
    parameter_samples <- PosteriorDraw_my(dpObj, ind, n = mhDraws)
    
    for (j in seq_along(clusterParams)) {
      clusterParams[[j]][, , i] <- parameter_samples[[j]][, , mhDraws]
    }
    
    
    accept_ratio[i] <- length(unique(parameter_samples[[1]]))/mhDraws
  }
  dpObj$clusterParameters <- clusterParams
  return(dpObj)
}

PosteriorDraw_my<-function (dpObj, ind, start_pos, n = 1, ...){
  UseMethod("PosteriorDraw_my", dpObj)
}

PosteriorDraw_my.nonconjugate<-function (dpObj, ind, n = 1, ...) {
  x <- dpObj$data[ind,,drop=FALSE]
  y <- dpObj$response[ind]
  n.x <- dim(x)[1]
  mu.0 <- dpObj$mixingDistribution$priorParameters$beta
  tau2.0 <- dpObj$mixingDistribution$priorParameters$D[1,1]
  nu.0 <- dpObj$mixingDistribution$priorParameters$nu0
  sigma2.0 <- dpObj$mixingDistribution$priorParameters$sigma02

  muSamples <- array(dim = c(2,1,n))
  sigSamples <- array(dim = c(1,1,n))

  muSamp <- rep(0,2)
  sig2Samp <- 1/rgamma(1,nu.0/2,nu.0*sigma2.0/2)
  for (i in seq_len(n)) {
    tau2.n <- solve(t(x)%*%x/sig2Samp+diag(2)/tau2.0)
    mu.n <- tau2.n%*%(t(x)%*%y/sig2Samp+mu.0/tau2.0)
    muSamp <- mvtnorm::rmvnorm(1,mean = mu.n,sigma = tau2.n)
    muSamples[,,i] <- muSamp

    nu.n <- nu.0+n.x
    sigma2.n <- (nu.0*sigma2.0+sum((t(y)-muSamp%*%t(x))^2))/nu.n
    sig2Samp <- 1/rgamma(1,nu.n/2,nu.n*sigma2.n/2)
    sigSamples[,,i] <- sqrt(sig2Samp)
  }
  return(list(muSamples,sigSamples))
}

# Metropolis:
MetropolisHastings_my<-function (dpObj, ind, start_pos, no_draws = 100) {
  UseMethod("MetropolisHastings_my", dpObj)
}

MetropolisHastings_my.default<-function (dpObj, ind, start_pos, no_draws = 100) {
  parameter_samples <- vector("list", length(start_pos))
  for (i in seq_along(start_pos)) {
    parameter_samples[[i]] <- array(dim = c(dim(start_pos[[i]])[1:2], 
                                            no_draws))
    parameter_samples[[i]][, , 1] <- start_pos[[i]][, , 1]
  }
  accept_count <- 0
  old_param <- start_pos
  old_prior <- log(PriorDensity_my(dpObj, old_param))
  
  old_all_param <- lapply(seq_along(old_param),
                          function(i) old_param[[i]][,,rep(1,dpObj$n), 
                                                    drop=FALSE])
  old_Likelihood <- sum(log(Likelihood(dpObj,old_all_param,ind)))
  
  for (i in seq_len(no_draws - 1)) {
    prop_param <- MhParameterProposal_my(dpObj, 
                                      old_param)
    # 
    new_prior <- log(PriorDensity_my(dpObj, prop_param))
    
    prop_all_param <- lapply(seq_along(prop_param),
                             function(i) prop_param[[i]][,,rep(1,dpObj$n), 
                                                        drop=FALSE])
    new_Likelihood <- sum(log(Likelihood(dpObj, 
                                         prop_all_param, ind)))
    accept_prob <- min(1, exp(new_prior + new_Likelihood - 
                                old_prior - old_Likelihood))
    if (is.na(accept_prob) | !length(accept_prob)) {
      accept_prob <- 0
    }
    if (runif(1) < accept_prob) {
      accept_count <- accept_count + 1
      sampled_param <- prop_param
      old_Likelihood <- new_Likelihood
      old_prior <- new_prior
    }
    else {
      sampled_param <- old_param
    }
    old_param <- sampled_param
    for (j in seq_along(start_pos)) {
      parameter_samples[[j]][, , i + 1] <- sampled_param[[j]]
    }
  }
  accept_ratio <- accept_count/no_draws
  return(list(parameter_samples = parameter_samples, accept_ratio = accept_ratio))
}

MhParameterProposal_my<-function (dpObj, old_params) 
{
  UseMethod("MhParameterProposal_my", dpObj)
}

MhParameterProposal_my.OLS <- function(dpObj, old_params){
  mdobj <- dpObj$mixingDistribution
  mhStepSize <- mdobj$mhStepSize
  new_params <- old_params
  
  # update beta:
  new_params[[1]][,,1] <- old_params[[1]][,,1]+
    mhStepSize[1]*mvtnorm:::rmvnorm(1,mean = rep(0,2),sigma = diag(rep(1,2)))
  # update beta.sig:
  new_params[[2]][,,1] <- old_params[[2]][,,1]+mhStepSize[2]*rnorm(1)
  if(new_params[[2]][,,1]<0){
    new_params[[2]][,,1]<-old_params[[2]][,,1]
  }

  return(new_params)
}

MhParameterProposal_my.Wishart <- function(dpObj, old_params){
  mdobj <- dpObj$mixingDistribution
  mhStepSize <- mdobj$mhStepSize
  new_params <- old_params
  
  # update beta:
  new_params[[1]][,,1] <- old_params[[1]][,,1]+
    mhStepSize[1]*mvtnorm:::rmvnorm(1,mean = rep(0,2),sigma = diag(rep(1,2)))
  
  # update beta.sig:
  new_params[[2]][,,1] <- old_params[[2]][,,1]+mhStepSize[2]*rnorm(1)
  if(new_params[[2]][,,1]<0){
    new_params[[2]][,,1]<-old_params[[2]][,,1]
  }
  
  # update D:
  new_params[[3]][,,1] <- old_params[[3]][,,1]+mhStepSize[3]*rWishart(1,
                                                                      df = 2,
                                                                      Sigma = diag(rep(1,2)))[,,1]/2
  
  return(new_params)
}

MhParameterProposal_my.Flat <- function(dpObj, old_params){
  mdobj <- dpObj$mixingDistribution
  mhStepSize <- mdobj$mhStepSize
  new_params <- old_params
  
  # update beta:
  new_params[[1]][,,1] <- old_params[[1]][,,1]+
    mhStepSize[1]*mvtnorm:::rmvnorm(1,mean = rep(0,2),sigma = diag(rep(1,2)))
  
  # update beta.sig:
  new_params[[2]][,,1] <- old_params[[2]][,,1]+mhStepSize[2]*rnorm(1)
  if(new_params[[2]][,,1]<0){
    new_params[[2]][,,1]<-old_params[[2]][,,1]
  }
  return(new_params)
}

# Update Alpha:
UpdateAlpha_my<-function(dpobj){
  UseMethod("UpdateAlpha_my", dpobj)
}

UpdateAlpha_my.default<-function (dpobj){
  newAlpha <- update_concentration_my(dpobj$alpha, dpobj$n, dpobj$numberClusters, 
                                   dpobj$alphaPriorParameters)
  dpobj$alpha <- newAlpha
  return(dpobj)
}

update_concentration_my<-function (oldParam, n, nParams, priorParameters){
  # x <- rbeta(1, oldParam + 1, n)
  # pi1 <- priorParameters[1] + nParams - 1
  # pi2 <- n * (priorParameters[2] - log(x))
  # pi1 <- pi1/(pi1 + pi2)
  x <- oldParam/(oldParam+n)
  postParams1 <- priorParameters[1] + nParams
  postParams2 <- priorParameters[2] - log(x)
  # if (runif(1) > pi1) {
  #   postParams1 <- postParams1 - 1
  # }
  new_alpha <- rgamma(1, postParams1, postParams2)
  return(new_alpha)
}

# Prior Density:
PriorDensity_my<-function (dpObj, theta){
  UseMethod("PriorDensity_my", dpObj)  
}

PriorDensity_my.Flat <- function(dpobj, theta){
  mdobj <- dpobj$mixingDistribution
  priorParameters <- mdobj$priorParameters
  
  beta <- priorParameters$beta
  beta.sig <- priorParameters$beta.sig
  D <- priorParameters$D
  
  sigma.density <- as.numeric(dgamma(theta[[2]][,,,drop=TRUE], 
                                     shape = beta.sig[1], 
                                     scale = beta.sig[2]))
  
  beta.density <- mvtnorm:::dmvnorm(x = theta[[1]][,,,drop=TRUE],
                          mean = priorParameters$beta,
                          sigma = D)
  
  prior.density <- sigma.density*beta.density
  return(prior.density)
}
# posterior inference:
post.region <- function(dpObj, num.ind, cluster.ind){
  mu.0 <- dpObj$mixingDistribution$priorParameters$beta
  tau2.0 <- dpObj$mixingDistribution$priorParameters$D[1,1]
  ind <- which(dpObj$labelsChain[[num.ind]]==cluster.ind)
  x.i <- dpObj$data[ind,,drop=FALSE]
  y.i <- dpObj$response[ind,]
  sigma2.n <- dpObj$clusterParametersChain[[num.ind]][[2]][,,cluster.ind]^2
  
  SIG2.n <- solve(t(x.i)%*%x.i/sigma2.n+diag(2)/tau2.0)
  MU.n <- SIG2.n%*%(t(x.i)%*%y.i/sigma2.n+mu.0/tau2.0)
  return(list(MU.n,SIG2.n))
}

library(plotrix)
# plot dp:
plot_my<-function(dpObj, centroid, num.ind, post=FALSE){
  UseMethod("plot_my", dpObj)  
}

plot_my.default<-function(dpObj, centroid, num.ind, post=FALSE){
  cluster_num <- dpObj$numLabelChain[[num.ind]]
  pointpercluster <- dpObj$pointsPerClusterChain[[num.ind]]
  if(post){
    MU.mat <- matrix(0,nrow = 2,ncol = cluster_num)
    SIG.list <- list()
    for(i in 1:cluster_num){
      MU.mat[,i] <- post.region(dpObj, num.ind, i)[[1]]
      SIG.list[[i]] <- post.region(dpObj, num.ind, i)[[2]]
    }
    xbound <- c(0,max(max(MU.mat[1,])+500,centroid[1,1][[1]]+500,max(dpObj$response)))
    ybound <- c(min(min(MU.mat[2,])-500,centroid[1,2][[1]]-500),0)
    plot(centroid[,1],centroid[,2],col="blue",
         lwd=5,
         pch=4,
         xlab = "x direction",
         ylab = "distance from wall",
         ylim = ybound,
         xlim = xbound)
    for(i in 1:cluster_num){
      points(x = MU.mat[1,i], y = MU.mat[2,i],col="red",lwd=0.1)
      textbox(c(MU.mat[1,i]-20,MU.mat[1,i]+20),
              MU.mat[2,i]-30,pointpercluster[i],box = FALSE)
      
      ellipse(mu = MU.mat[,i],
              sigma = SIG.list[[i]],alpha = 0.05,npoints = 500,
              draw = TRUE,newplot = FALSE,col="red",lwd=1e-10)    
    }     
  }else{
    cluster_params <- dpObj$clusterParametersChain[[num.ind]]
    label <- dpObj$labelsChain[[num.ind]]
    pointPerCluster<-dpObj$pointsPerClusterChain[[num.ind]]
    x.dat <- dpObj$data
    params <- cluster_params[[1]][,,,drop=TRUE]
    
    xbound <- c(0,max(max(params[1,])+500,centroid[1,1][[1]]+500,max(dpObj$response)))
    ybound <- c(min(min(params[2,])-500,centroid[1,2][[1]]-500),0)
    
    plot(centroid[,1],centroid[,2],col="blue",
         lwd=5,
         pch=4,
         xlab = "x direction",
         ylab = "distance from wall",
         ylim = ybound,
         xlim = xbound)
    for(i in 1:cluster_num){
      ind <- which(label==i)
      X.mat<-x.dat[ind,]
      try.solve <- try(solve(t(X.mat)%*%X.mat),silent = TRUE)
      points(x = cluster_params[[1]][,,i][1],y=cluster_params[[1]][,,i][2],col="red",lwd=0.1)
      if(pointPerCluster[i]!=1 & class(try.solve)!="try-error"){
        ellipse(mu =cluster_params[[1]][,,i],
                sigma = cluster_params[[2]][,,i]^2*solve(t(X.mat)%*%X.mat),alpha = 0.05,npoints = 500,
                draw = TRUE,newplot = FALSE,col="red",lwd=1e-10)    
      }
    }    
  }
}

# posterior inference:
piDirichlet_my<-function(betas){
  pis <- numeric(length(betas))
  pis[1] <- betas[1]
  for (i in 2:length(betas)) {
    pis[i] <- betas[i] * prod(1 - betas[1:(i - 1)])
  }
  return(pis)
}

StickBreaking_my<-function(alpha, N){
  betas <- rbeta(N, 1, alpha)
  pis <- piDirichlet_my(betas)
  return(pis)
}

PosteriorClusters<-function (dpobj, ind){
  UseMethod("PosteriorClusters", dpobj) 
}

PosteriorClusters.dirichletprocess<-function (dpobj, ind){
  if (!missing(ind)) {
    pointsPerCluster <- dpobj$weightsChain[[ind]] * dpobj$n
    alpha <- dpobj$alphaChain[ind]
    clusterParams <- dpobj$clusterParametersChain[[ind]]
  }
  else {
    pointsPerCluster <- dpobj$pointsPerCluster
    alpha <- dpobj$alpha
    clusterParams <- dpobj$clusterParameters
  }
  numLabels <- length(pointsPerCluster)
  mdobj <- dpobj$mixingDistribution
  dirichlet_draws <- gtools::rdirichlet(1, c(pointsPerCluster, 
                                             alpha))
  numBreaks <- ceiling(alpha + numLabels) * 20 + 5
  sticks <- StickBreaking_my(alpha + numLabels, numBreaks)
  sticks <- sticks * dirichlet_draws[numLabels + 1]
  sticks <- c(dirichlet_draws[-(numLabels + 1)], sticks)
  n_smps <- numBreaks + numLabels
  PriorDraws <- PriorDraw(dpobj, numBreaks)
  postParams <- list()
  for (i in seq_along(clusterParams)) {
    postParams[[i]] <- array(c(clusterParams[[i]], PriorDraws[[i]]), 
                             dim = c(dim(PriorDraws[[i]])[1:2], numBreaks + numLabels))
  }
  returnList <- list(weights = sticks, params = postParams)
  return(returnList)
}
# Posterior density sampling:
postSampling <- function(n, weights, params, dpObj){
  probs <- unlist(weights)
  beta.mat <- matrix(0,nrow = n,ncol = 2)
  X.mat.list <- list()
  for(i in 1: dpObj$numberClusters){
    X.mat.i <- dpObj$data[which(dpObj$clusterLabels==i),]
    X.mat.list[[i]] <- solve(t(X.mat.i)%*%X.mat.i)
  }
  for(i in seq_len(n)){
    ind.i <- which(rmultinom(1,1,probs)==1)
    if(ind.i <= dpObj$numberClusters){
      cov.inv <- X.mat.list[[ind.i]]*params[[2]][,,ind.i]^2
    }else{
      cov.inv <- matrix(c(1e8,0,0,1e8),nrow = 2,ncol = 2)
    }
    beta <- mvtnorm:::rmvnorm(1,mean = params[[1]][,,ind.i],sigma = cov.inv)
    beta.mat[i,] <- beta
  }
}

# splitting event simulation:
# beta1 <- c(3000,-8000)
# beta2 <- c(4000,-7000)
# beta3 <- c(5000,-8000)
# beta0 <- c(4000,-9000)
# val1 <- (beta1[1]-beta0[1])/(beta1[2]-beta0[2])
# val2 <- (beta2[1]-beta0[1])/(beta2[2]-beta0[2])
# val3 <- (beta3[1]-beta0[1])/(beta3[2]-beta0[2])
# 
# x1 <- rnorm(100,val1,0.5);x1.mat <- model.matrix(~x1)
# x2 <- rnorm(100,val2,0.5);x2.mat <- model.matrix(~x2)
# x3 <- rnorm(100,val3,0.5);x3.mat <- model.matrix(~x3)
# 
# y1 <- x1.mat%*%beta1 + rnorm(100,mean = 0,sd = 50)
# y2 <- x2.mat%*%beta2 + rnorm(100,mean = 0,sd = 50)
# y3 <- x3.mat%*%beta3 + rnorm(100,mean = 0,sd = 50)
# 
# x <- model.matrix(~c(x1,x2,x3))
# y <- c(y1,y2,y3)
# 
# lm1 <- lm(y~x[,2])
# betain <- lm1$coefficient
# 
# dpObj <- DirichletProcessCreate_my(x,y,
#                                    mdObject = FlatMixtureCreate(priorParameters = list(beta=array(betain,
#                                                                                 dim = c(2,1)),
#                                                                      nu0=2,sigma02=100^2,
#                                                                      D=matrix(c(1e8,0,0,1e8),nrow = 2,ncol = 2))),
#                                    alphaPriorParameters = c(50,1))
# dpObj <- Initialise_my(dpObj)
# dpObj <- Fit_my(dpObj,its = 10000)
# 


