#------------------ Parallel Running ----------------#
library(foreach)
library(parallel)
library(doParallel)
library(coda)
library(mvnfast)
library(MCMCpack)


## Working Cores
my.cluster = makeCluster(4, type = "PSOCK")
registerDoParallel(cl = my.cluster)
getDoParRegistered()
getDoParWorkers()



source("Gibbs Samling with Nested Model.R")
#source("Gibbs Samling with Nested Model Cor.R")

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
  #alpha0 = rmvn(n2, Alpha0, diag(Salpha0,p))
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
    corrM = diag(p) + 0.01 * (matrix(1, p, p) - diag(p))
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


replicate = 4
p = 10 ## p >= 10
c = 2
n1 = 2
n2 = 100



results = 
  foreach (r = 1:replicate,
           .packages = c("Rlab","invgamma","MASS","LaplacesDemon","mvnfast","MCMCpack")) %dopar% {
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
             
             pos = Gibbs.Posterior(iteration  = 100, burn = 50, chain = 1)
           }

save(results,file = paste0("simulation/Gibbs/M", p, "A", 0.1 * p, "N", n1 * n2,"C",c,".Rdata"))































#############################################
## Threshold Determination at Certain FDR ###
#### For simulate all CORRECTED PIP##########
#############################################

load("simulation/Gibbs/M200A20N200C2.Rdata")

pos = results[[1]][[1]]
p=100

getPIP = function(pos){
  PIPa = which(colMeans(pos$` PIPalphaT`)>0.5)
  PIP = colMeans(pos$` PIPbetaM`)[PIPa]
  names(PIP) = paste0("cg",PIPa)
  return(PIP)
}
allPIP = lapply(results, sapply, getPIP)

FDRAdjust = function(PIP){
  
  allm = paste0("cg",1:p)
  trueactive = c(rep(1, p * 0.1), rep(0, 0.9 * p))
  
  FDRTable = as.data.frame(cbind(PIP,rep(0,length(PIP)),rep(0,length(PIP)),rep(0,length(PIP)),rep(0,length(PIP))))
  
  colnames(FDRTable) = c("PIPbetaM","FDR","FPR","TPR","TNR")
  
  for (i in 1:length(PIP)) {
    threshhold = FDRTable$PIPbetaM[i]
    postactiveCpG = rownames(FDRTable)[which(FDRTable$PIPbetaM >= threshhold)]
    postactive = ifelse(allm %in% postactiveCpG, 1, 0)
    
    False.Positive = ifelse(trueactive == 0 & postactive == 1, 1, 0)
    True.Postive   = ifelse(trueactive == 1 & postactive == 1, 1, 0)
    True.Negative  = ifelse(trueactive == 0 & postactive == 0, 1, 0)
    
    FDRTable[i,"FDR"] = sum(False.Positive) / sum(postactive)
    FDRTable[i,"FPR"] = sum(False.Positive) / (p - sum(trueactive))
    FDRTable[i,"TPR"] = sum(True.Postive) / sum(trueactive)
    FDRTable[i,"TNR"] = sum(True.Negative) / (p - sum(trueactive))
  }
  
  FDRcontrol = FDRTable[FDRTable[,"FPR"] <= FPR,]
  finalthresh = min(FDRcontrol[,"PIPbetaM"])
  Decision = FDRcontrol[FDRcontrol[,"PIPbetaM"] == finalthresh,][1,]
  
  return(Decision)
}

p=500
FPR = 0.1
rate = lapply(allPIP, FDRAdjust)
dfrate = as.data.frame(matrix(unlist(rate),length(rate),5,byrow=T))
colnames(dfrate) = c("Threshhold","FDR","FPR","TPR","TNR")
colMeans(na.omit(dfrate))


#############################################
## Threshold Determination at Certain FDR ###
##### For Nested Model corrected PIP ########
#############################################

load("simulation/Gibbs/M200A20N200C2.Rdata")

getPIP = function(pos){
  PIP = pos$PIPbetaM
  return(PIP)
}

allPIP = lapply(results, sapply, getPIP)

FDRAdjust = function(PIP){
  
  allm = paste0("cg",1:p)
  trueactive = c(rep(1, p * 0.1), rep(0, 0.9 * p))
  
  FDRTable = as.data.frame(cbind(PIP,rep(0,length(PIP)),rep(0,length(PIP)),rep(0,length(PIP)),rep(0,length(PIP))))
  
  colnames(FDRTable) = c("PIPbetaM","FDR","FPR","TPR","TNR")
  
  for (i in 1:length(PIP)) {
    threshhold = FDRTable$PIPbetaM[i]
    postactiveCpG = rownames(FDRTable)[which(FDRTable$PIPbetaM >= threshhold)]
    postactive = ifelse(allm %in% postactiveCpG, 1, 0)
    
    False.Positive = ifelse(trueactive == 0 & postactive == 1, 1, 0)
    True.Postive   = ifelse(trueactive == 1 & postactive == 1, 1, 0)
    True.Negative  = ifelse(trueactive == 0 & postactive == 0, 1, 0)
    
    FDRTable[i,"FDR"] = sum(False.Positive) / sum(postactive)
    FDRTable[i,"FPR"] = sum(False.Positive) / (p - sum(trueactive))
    FDRTable[i,"TPR"] = sum(True.Postive) / sum(trueactive)
    FDRTable[i,"TNR"] = sum(True.Negative) / (p - sum(trueactive))
  }
  
  FDRcontrol = FDRTable[FDRTable[,"FPR"] <= FPR,]
  finalthresh = min(FDRcontrol[,"PIPbetaM"])
  Decision = FDRcontrol[FDRcontrol[,"PIPbetaM"] == finalthresh,][1,]
  
  return(Decision)
}

p=100
FPR = 0.1
rate = lapply(allPIP, FDRAdjust)
dfrate = as.data.frame(matrix(unlist(rate),length(rate),5,byrow=T))
colnames(dfrate) = c("Threshhold","FDR","FPR","TPR","TNR")
colMeans(na.omit(dfrate))

