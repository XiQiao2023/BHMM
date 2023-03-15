####################################################################
############# Univariate Mediation with Mediation package ##########
####################################################################
#------------------ Parallel Running ----------------#
library(foreach)
library(parallel)
library(doParallel)
library(coda)
library(mediation)
library(lme4)
library(mvnfast)

## Working Cores
my.cluster = makeCluster(4, type = "PSOCK")
my.cluster 
registerDoParallel(cl = my.cluster)
getDoParRegistered()
getDoParWorkers()

detach_package <- function(pkg, character.only = FALSE)
{
  if(!character.only)
  {
    pkg <- deparse(substitute(pkg))
  }
  search_item <- paste("package", pkg, sep = ":")
  while(search_item %in% search())
  {
    detach(search_item, unload = TRUE, character.only = TRUE)
  }
}

detach_package(lmerTest)

multi.p.data = function(n1,n2,c,
                        Alpha0, Salpha0, alphaT, alphaX, sigmaM,
                        Beta0, Sbeta0, betaT, betaM, betaX,sigmaY)
{
  ## Sample Size
  N = n1 * n2 
  
  ## Group ID
  Group = matrix(rep(1:n2, n1), N, 1,byrow = TRUE) 
  
  ## Exposure
  T = matrix(rep(0:1, each = n2), N, 1,byrow = TRUE) 
  
  ## c continuous covariates 
  X = matrix(rnorm(N * c, 0, 1), N, c, byrow = TRUE)
  
  ## Random Intercept
  ## Uncorrelated Mediators
  ##alpha0 = rmvn(n2, Alpha0, diag(Salpha0,p))
  ## Correlated Mediators
  corr0 = diag(p) + 0.01 * (matrix(1, p, p) - diag(p))
  alpha0 = rmvn(n2, Alpha0, diag(sqrt(Salpha0), p)%*% corr0 %*%diag(sqrt(Salpha0), p))
  
  beta0 = rnorm(n2, Beta0, sqrt(Sbeta0))
  
  ## data generation 
  mean.M = matrix(c(rep(0, N * p)), N, p, byrow = TRUE)
  M = matrix(c(rep(0, N * p)), N, p, byrow = TRUE)
  
  mean.Y =  matrix(rep(0, N), N, 1 ,byrow = TRUE)
  Y =  matrix(rep(0, N), N, 1 ,byrow = TRUE)
  
  for (i in 1:N) {
    mean.M[i, ] = alpha0[Group[i],] + T[i] * alphaT + X[i,] %*% t(alphaX)
    ## Uncorrelated Mediators
    ##M[i, ] = rmvn(1, mean.M[i,], diag(sigmaM, p)) 
    ## Correlated Mediators
    corrM = diag(p) + 0.01* (matrix(1, p, p) - diag(p))
    M[i, ] = rmvn(1, mean.M[i,], diag(sqrt(sigmaM), p)%*% corrM %*%diag(sqrt(sigmaM), p))
    
    mean.Y[i] = beta0[Group[i]] + T[i] * betaT + M[i,] %*% betaM + X[i,] %*% betaX 
    Y[i] = rnorm(1, mean.Y[i], sqrt(sigmaY))
  }
  
  ## return data
  colnames(M)= paste0("cg",1:p)
  data = list(Group,T,X,M,Y)
  names(data) = c("Group","T","X","M","Y")
  return(data)
}

SingleMediation = function(M, nsim = 2003){
  #M = scale(M, center = TRUE, scale = TRUE)
  fit.mediator= lmer(M ~ T + X + (1|Group), REML = FALSE)
  
  #Y = scale(Y, center = TRUE, scale = TRUE)
  fit.outcome = lmer(Y ~ T + M + X + (1|Group), REML = FALSE)
  
  results = mediate(fit.mediator,fit.outcome,treat = "T", mediator = "M", sims = nsim)
  
  ACME = results$d0
  ACME_pvalue = results$d0.p
  ADE = results$z0
  ADE_pvalue = results$z0.p
  proportion_mediated = results$n0
  
  out = c(ACME, ACME_pvalue, ADE, ADE_pvalue, proportion_mediated)
  names(out) = c("ACME", "ACME_pvalue", "ADE", "ADE_pvalue", "proportion_mediated")
  return(out)
}

replicate = 4
p = 500 ## p >= 10
c = 2
n1 = 2
n2 = 100

results = 
  foreach (r = 1:replicate,
           .packages = c("MASS","mediation","lme4")) %dopar% {
             tmp = multi.p.data(n1 = n1, n2 = n2, c = c,
                                alphaT = c(rep(1, 0.2 * p), rep(0, 0.8 * p)),
                                betaM =  c(rep(1, 0.1 * p), rep(0, 0.1 * p), rep(1, 0.1 * p), rep(0, 0.7 * p)),
                                Alpha0 = rep(1,p), Salpha0 = 0.01,
                                alphaX = matrix(rep(1,p*c),p,c,byrow = TRUE),sigmaM = 0.01,
                                Beta0 = 1,Sbeta0 = 0.01,betaT = 1,
                                betaX = rep(1,c),sigmaY = 0.01)
             
             Group = tmp$Group
             X = tmp$X
             T = tmp$T
             M = tmp$M
             Y = tmp$Y
             n = nrow(X)
             n2 = length(unique(Group))
             p = ncol(M)
             c = ncol(X)
             
             apply(M, 2, SingleMediation)
           }


save(results,file = paste0("simulation/Univariate/M", p, "A", 0.1 * p, "N", n1 * n2,"C",c,".Rdata"))

