library(Rcpp)
library(RcppArmadillo)
library(mvnfast)
library(MCMCpack)

sourceCpp("Rcpp.cpp")

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
  ## Specify Correlation between of random intercept
  corr0 = diag(p) + 0 * (matrix(1, p, p) - diag(p))
  alpha0 = mvrnorm(n2, Alpha0, diag(sqrt(Salpha0), p)%*% corr0 %*%diag(sqrt(Salpha0), p))
  
  beta0 = rnorm(n2, Beta0, sqrt(Sbeta0))
  
  ## data generation 
  mean.M = matrix(c(rep(0, N * p)), N, p, byrow = TRUE)
  M = matrix(c(rep(0, N * p)), N, p, byrow = TRUE)
  
  mean.Y =  matrix(rep(0, N), N, 1 ,byrow = TRUE)
  Y =  matrix(rep(0, N), N, 1 ,byrow = TRUE)
  
  for (i in 1:N) {
    mean.M[i, ] = alpha0[Group[i],] + T[i] * alphaT + X[i,] %*% t(alphaX)
    ## Specify Correlation between Mediators
    corrM = diag(p) + 0 * (matrix(1, p, p) - diag(p))
    M[i, ] = mvrnorm(1, mean.M[i,], diag(sqrt(sigmaM), p)%*% corrM %*%diag(sqrt(sigmaM), p))
  }  
  
  M = apply(M, 2, function(x){scale(x,center = TRUE, scale = TRUE)})
  
  for (i in 1:N) {
    mean.Y[i] = beta0[Group[i]] + T[i] * betaT + M[i,] %*% betaM + X[i,] %*% betaX 
    Y[i] = rnorm(1, mean.Y[i], sqrt(sigmaY))
  }
  
  ## return data
  colnames(M)= paste0("cg",1:p)
  data = list(Group,T,X,M,Y)
  names(data) = c("Group","T","X","M","Y")
  return(data)
}


p = 10 ## p >= 10
c = 2
n1 = 2
n2 = 100

tmp = multi.p.data(n1 = n1, n2 = n2, c = c,
                   alphaT = c(rep(1, 0.2 * p), rep(0, 0.8 * p)),
                   betaM =  c(rep(1, 0.1 * p), rep(0, 0.1 * p), rep(1, 0.1 * p), rep(0, 0.7 * p)),
                   Alpha0 = rep(1,p), Salpha0 = 0.01,
                   alphaX = matrix(rep(1,p*c),p,c,byrow = TRUE),sigmaM = 0.01,
                   Beta0 = 1,Sbeta0 = 0.01,betaT = 1, betaX = rep(1,c),sigmaY = 0.01)

Group = tmp$Group
X = tmp$X
T = tmp$T
M = tmp$M
Y = tmp$Y
n = nrow(X)
n2 = length(unique(Group))
p = ncol(M)
c = ncol(X)
             
pos = Gibbs_Posterior(Group,M,T,X,Y,iteration = 80000, burn = 40000,chain = 1,p,n2,n,c)

save(pos,file = paste0("simulation1/M", p, "A", 0.1 * p, "N", n1 * n2,"C",c,"_", Sys.getenv("SLURM_ARRAY_TASK_ID"),".Rdata"))


