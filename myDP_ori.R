rm(list=ls())
library(mvtnorm)
library(LaplacesDemon)
library(mixtools)
registerDoParallel(cores=4)
# Initialization:
OLSMixtureCreate<-function (priorParameters = list(beta=array(c(5000,-3000),
                                                              dim = c(2,1)),beta.sig=c(1e2,1e1))) 
{
  mdobj <- MixingDistribution_my("OLS", priorParameters, 
                              "nonconjugate",
                              mhStepSize = c(100, 5))
  return(mdobj)
}

WishartMixtureCreate<-function(priorParameters = list(beta=array(c(5000,-3000),
                                                      dim = c(2,1)),
                                                      beta.sig=c(1e2,1e1),
                                                      D=matrix(c(1e8,0,0,1e8),nrow = 2,ncol = 2)))
{
  mdobj <- MixingDistribution_my("Wishart", priorParameters, 
                                 "nonconjugate",
                                 mhStepSize = c(100, 5, 1e4))
  return(mdobj)  
}

FlatMixtureCreate<-function(priorParameters = list(beta=array(c(5000,-3000),
                                                              dim = c(2,1)),
                                                   beta.sig=c(1e2,1e1),
                                                   D=matrix(c(1e8,0,0,1e8),nrow = 2,ncol = 2)))
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

DirichletProcessCreate_my<-function (x, y, mdObject=OLSMixtureCreate(), 
                                     alphaPriorParameters = c(1, 1), mhDraws = 250) 
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
  
  beta.sig <- priorParameters$beta.sig
  
  beta.sig <- array(rgamma(n, shape = beta.sig[1], scale = beta.sig[2]),dim = c(1,1,n))
  
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
Initialise_my.nonconjugate <- function(dpObj, posterior = TRUE, m = 3, verbose = TRUE) {
  
  # dpObj$clusterLabels <- 1:dpObj$n dpObj$numberClusters <- dpObj$n
  # dpObj$pointsPerCluster <- rep(1, dpObj$n) dpObj$clusterParameters <-
  # PosteriorDraw(dpObj$MixingDistribution, dpObj$data, dpObj$n)
  dpObj$clusterLabels <- rep(1, dpObj$n)
  dpObj$numberClusters <- 1
  dpObj$pointsPerCluster <- dpObj$n
  
  if (posterior) {
    post_draws <- PosteriorDraw(dpObj$mixingDistribution, dpObj$data, 1000)
    
    if (verbose)
      cat(paste("Accept Ratio: ",
                length(unique(c(post_draws[[1]])))/1000,
                "\n"))
    
    dpObj$clusterParameters <- list(post_draws[[1]][, , 1000, drop = FALSE],
                                    post_draws[[2]][, , 1000, drop = FALSE])
  } else {
    dpObj$clusterParameters <- PriorDraw(dpObj, 1)
  }
  
  dpObj$m <- m
  
  return(dpObj)
}

# Fitting:
Fit_my <- function(dpObj, its, updatePrior = FALSE, progressBar=TRUE) UseMethod("Fit_my", dpObj)

Fit_my.default <- function(dpObj, its, updatePrior = FALSE, progressBar=TRUE) {
  
  if (progressBar){
    pb <- txtProgressBar(min=0, max=its, width=50, char="-", style=3)
  }
  
  alphaChain <- numeric(its)
  likelihoodChain <- numeric(its)
  weightsChain <- vector("list", length = its)
  clusterParametersChain <- vector("list", length = its)
  priorParametersChain <- vector("list", length = its)
  labelsChain <- vector("list", length = its)
  
  for (i in seq_len(its)) {
    
    alphaChain[i] <- dpObj$alpha
    weightsChain[[i]] <- dpObj$pointsPerCluster / dpObj$n
    clusterParametersChain[[i]] <- dpObj$clusterParameters
    priorParametersChain[[i]] <- dpObj$mixingDistribution$priorParameters
    labelsChain[[i]] <- dpObj$clusterLabels
    
    
    likelihoodChain[i] <- sum(log(LikelihoodDP_my(dpObj)))
    
    dpObj <- ClusterComponentUpdate_my(dpObj)
    dpObj <- ClusterParameterUpdate_my(dpObj)
    dpObj <- UpdateAlpha_my(dpObj)
    
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
      pointsPerCluster * Likelihood(dpObj, clusterParams, i),
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
    
    for (j in seq_along(clusterParams)) {
      start_pos[[j]] <- clusterParams[[j]][, , i, drop = FALSE]
    }
    
    parameter_samples <- PosteriorDraw_my(dpObj, ind, start_pos = start_pos, n = mhDraws)
    
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

PosteriorDraw_my.nonconjugate<-function (dpObj, ind, start_pos, n = 1, ...) {
  # if (missing(...)) {
  #   start_pos <- PenalisedLikelihood(mdObj, x)
  # }
  # else {
  #   start_pos <- list(...)$start_pos
  # }
  mh_result <- MetropolisHastings_my(dpObj, ind, start_pos, no_draws = n)
  theta <- vector("list", length(mh_result))
  for (i in seq_along(mh_result$parameter_samples)) {
    theta[[i]] <- array(mh_result$parameter_samples[[i]], 
                        dim = c(dim(mh_result$parameter_sample[[i]])[1:2], 
                                n))
  }
  return(theta)
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
  return(new_alpha)
}

# Prior Density:
PriorDensity_my<-function (dpObj, theta){
  UseMethod("PriorDensity_my", dpObj)  
}

PriorDensity_my.OLS <- function(dpobj, theta){
  x.mat <- dpobj$data
  x.mat <- solve(t(x.mat)%*%x.mat)
  
  mdobj <- dpobj$mixingDistribution
  priorParameters <- mdobj$priorParameters
  beta.sig <- priorParameters$beta.sig
  
  sigma.density <- as.numeric(dgamma(theta[[2]][,,,drop=TRUE], 
                                     shape = beta.sig[1], 
                                     scale = beta.sig[2]))
  # g prior:
  beta.density <- mvtnorm:::dmvnorm(x = theta[[1]][,,,drop=TRUE],
                          mean = priorParameters$beta,
                          sigma = theta[[2]][,,,drop=TRUE]^2*dpobj$n*x.mat)
  prior.density <- sigma.density*beta.density  
  return(prior.density)
}

PriorDensity_my.Wishart <- function(dpobj, theta){
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
                          sigma = theta[[3]][,,,drop=TRUE])
  
  D.density <- dwishart(Omega = theta[[3]][,,,drop=TRUE]*2,nu = 2, S= D)
  
  prior.density <- sigma.density*beta.density*D.density
  return(prior.density)
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

# plot dp:
plot_my<-function(dpObj, centroid){
  UseMethod("plot_my", dpObj)  
}

plot_my.default<-function(dpObj, centroid){
  cluster_params <- dpObj$clusterParameters
  cluster_num <- dpObj$numberClusters
  label <- dpObj$clusterLabels
  pointPerCluster<-dpObj$pointsPerCluster
  x.dat <- dpObj$data
  params <- dpObj$clusterParameters[[1]][,,,drop=TRUE]
  
  xbound <- c(min(min(params[1,])-500,centroid[1,1][[1]]-500),max(max(params[1,])+500,centroid[1,1][[1]]+500))
  ybound <- c(min(min(params[2,])-500,centroid[1,2][[1]]-500),max(max(params[2,])+500,centroid[1,2][[1]]+500))
  
  plot(centroid[,1],centroid[,2],col="blue",
       lwd=0.1,
       xlab = "x direction",
       ylab = "distance from wall",
       ylim = ybound,
       xlim = xbound)
  
  for(i in 1:cluster_num){
    ind <- which(label==i)
    if(pointPerCluster[i]!=1){
      X.mat<-x.dat[ind,]
      points(x = cluster_params[[1]][,,i][1],y=cluster_params[[1]][,,i][2],col="red",lwd=0.1)
      ellipse(mu =cluster_params[[1]][,,i],
              sigma = cluster_params[[2]][,,i]^2*solve(t(X.mat)%*%X.mat),alpha = 0.05,npoints = 500,
              draw = TRUE,newplot = FALSE,col="red",lwd=1e-10)    
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
# # ## Debugging:
# centroid<-matrix(c(4000, 5000, 6000,
#                    -6000, -2000, -5000),nrow = 3,ncol = 2)
# x1<-scale((1:100));x2<-scale((1:100));x3<-scale((1:100))
# x<- model.matrix(~c(x1,x2,x3))
# y<-c(model.matrix(~x1)%*%c(6000,-5000),
#      model.matrix(~x2)%*%c(4000,-6000),
#      model.matrix(~x3)%*%c(5000,-2000))
# 
# # # well separated
# # # sd<-c(rep(100,100),rep(100,100),rep(100,100))
# # # y<-rnorm(300,mean=y,sd=sd)
# # # 
# # # overlapped
# sd<-c(rep(500,100),rep(500,100),rep(500,100))
# y<-rnorm(300,mean=y,sd=sd)
# 
# # # # extremely overlapped
# # # sd<-c(rep(1000,100),rep(1000,100),rep(1000,100))
# # # y<-rnorm(300,mean=y,sd=sd)
# # # 
# dp1<-DirichletProcessCreate_my(x,y,mdObject = FlatMixtureCreate(),alphaPriorParameters = c(10,0.5))
# dp1<-Initialise_my(dp1,FALSE)
# dp1<- Fit_my(dp1,50)
# post.dp1<-PosteriorClusters(dp1)
# # 
# # # posterior plot:
# library(plotly)
# x_vec <- c(seq(3000,7000,20))
# z_vec <- c(seq(-7000,-1000,20))
# mean.list <- post.dp1$params[[1]]
# sd.list <- post.dp1$params[[2]]
# weight.list <- post.dp1$weights
# post.ind <- dp1$clusterLabels
# dat.mat <- dp1$data
# post.density <- function(x){
#   numL <- length(unique(post.ind))
#   n <- length(weight.list)
#   sum0<-0
#   for(i in 1:n){
#     x.mat <- matrix(c(1e8,0,0,1e8),nrow = 2,ncol = 2)
#     if(i <= 3){
#       ind.i <- which(post.ind==i)
#       x.mat <- dat.mat[ind.i,]
#       x.mat <- solve(t(x.mat)%*%x.mat)*sd.list[i]^2
#     }
#     sum0 <- sum0+mvtnorm::dmvnorm(x,mean = mean.list[,,i],sigma = x.mat)
#   }
#   return(sum0)
# }
# #post.density(c(4000,-5000),mean.list,sd.list,weight.list,post.ind,dat.mat)
# val <- matrix(0,301,201)
# for(i in 1:201){
#   for(j in 1:301){
#     val[j,i] <- post.density(c(x_vec[i],z_vec[j]))
#   }
# }
# val.95contour <- matrix(quantile(val,0.05),301,201)
# p <- plot_ly(x=~x_vec,y=~z_vec,z=~val,type = "surface")
# p
# 
# val.max <- max(val)
# 
# reject.sample.2d <- function(n,pdf,maxval,xlim,ylim)
# {
#   smpl <- data.frame(x=numeric(n),y=numeric(n))
#   i <- 0
#   while (i<n)
#   {
#     xval <- runif(1,xlim[1],xlim[2])
#     yval <- runif(1,ylim[1],ylim[2])
#     if (runif(1)<pdf(c(xval,yval))/maxval)
#     {
#       i <- i+1
#       smpl[i,] <- c(xval,yval)
#     }
#   }
#   return(smpl)
# }
# 
# source("C:\\Users\\tiany\\Documents\\HPDregionplot.R")
# library(coda)
# res <- reject.sample.2d(1e3,post.density,val.max+0.1e-5,c(3000,7000),c(-7000,-1000))
# HPDregionplot(res,lims = c(3000,7000,-7000,-1000),prob = 0.9)

# data <- list(
#   x = x_vec,
#   y = z_vec,
#   z = val,
#   type = "surface")
# 
# axx <- list(
#   nticks = 4,
#   range = c(3e3,7e3)
# )
# 
# axy <- list(
#   nticks = 4,
#   range = c(-7e3,-1e3)
# )
# 
# axz <- list(
#   nticks = 4,
#   range = c(0,1e-3)
# )

# plot_my(dp1,centroid)
# random.var.list<-list()
# for(i in 1:5){
#   blood.stain <- c()
#   ind<-dp1$clusterLabels
#   sigmaI <- diag(sapply(ind, function(i) return(dp1$clusterParameters[[2]][,,i])))
#   blood.stain$int <- solve(sigmaI)%*%x[,1]
#   blood.stain$x <- solve(sigmaI)%*%x[,2];blood.stain$y <- solve(sigmaI)%*%y
#   blood.stain$id <- dp1$clusterLabels;
#   
#   fit.lme <- lme(y~ -1+int+x,method = "REML",random = reStruct(~-1+int+x|id,pdClass = "pdSymm"),
#                  data = blood.stain)
#   beta.fix <- fixed.effects(fit.lme)
#   v <- VarCorr(fit.lme)
#   D.var <- as.numeric(v[,"Variance"])[1:2]
#   D.corr <- as.numeric(v[,"Corr"][2])
#   D.cov <- sqrt(D.var[1]*D.var[2])*D.corr
#   D <- matrix(c(D.var[1],D.cov,D.cov,D.var[2]),nrow = 2)
#   dp1$mixingDistribution$priorParameters$D <- D
#   dp1$mixingDistribution$priorParameters$beta <- matrix(beta.fix,nrow = 2)
#   dp1 <- Fit_my(dp1,50)
#   plot_my(dp1,centroid)
# }
# plot_my(dp1,centroid)
# points(x = beta.fix[1],y=beta.fix[2],col="green",lwd=0.1)
# 
# # plot the posterior variance region:
# n.cluster <- dp1$numberClusters
# x.mat <- dp1$data
# for(i in 1:n.cluster){
#   n.i <- dp1$pointsPerCluster[i]
#   ind.i <- which(dp1$clusterLabels==i)
#   z <- x.mat[ind.i,]
#   R.i <- diag(rep(dp1$clusterParameters[[2]][,,i]^2,n.i))
#   sigma.s <- D-D%*%t(z)%*%solve(z%*%D%*%t(z)+R.i)%*%z%*%D
#   random.mat.list[[i]] <- sigma.s
# }
# for(i in 1:n.cluster){
#   ellipse(mu =dp1$clusterParameters[[1]][,,i],
#           sigma = random.mat.list[[i]],alpha = 0.05,npoints = 500,
#           draw = TRUE,newplot = FALSE,col="green",lwd=1e-10)  
# }

# wishart plot:
# plot(c(6000,4000,5000),c(-2000,-3000,-4000),col="blue",
#      lwd=0.1,
#      xlab = "x direction",
#      ylab = "distance from wall")
# for(i in 1:cluster_num){
#   if(pointPerCluster[i]!=1){
#     X.mat<-x[ind,]
#     points(x = cluster_params[[1]][,,i][1],y=cluster_params[[1]][,,i][2],col="red",lwd=0.1)
#     ellipse(mu =cluster_params[[1]][,,i],
#             sigma = cluster_params[[3]][,,i],alpha = 0.05,npoints = 500,
#             draw = TRUE,newplot = FALSE,col="red",lwd=1e-10)    
#   }
# }
