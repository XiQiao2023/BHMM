############################################################################################
##################### Multivariate multilevel mediation analysis ###########################
############################# Bayesian Variable Selection ##################################
#################################### Corrected PIP #########################################
############################################################################################

Gibbs.Posterior = function(iteration, burn, chain){
  
  posterior = vector(mode='list', length = chain)
  
  for (l in 1:chain) {
    
    iter = iteration + burn
    
    ###########################################################################
    ########################## Mediator Model #################################
    ###########################################################################
    
    #-------- Mediator Model Parameters -------#
    #pos.Alpha0   = matrix(rep(0, p * iter), iter, p, byrow = T)
    #pos.Salpha0  = array(rep(diag(p), iter), dim = c(p, p, iter))
    
    #pos.PalphaT  = matrix(rep(0, 1 * iter), iter, 1, byrow = T)
    PIP.alphaT   = matrix(rep(1, p * iter), iter, p, byrow = T)
    #pos.IalphaT  = matrix(rep(1, p * iter), iter, p, byrow = T)
    pos.alphaT   = matrix(rep(0, p * iter), iter, p, byrow = T)
    #pos.SalphaT1 = matrix(rep(1, 1 * iter), iter, 1, byrow = T)
    #pos.SalphaT0 = matrix(rep(1, 1 * iter), iter, 1, byrow = T)
    
    #pos.sigmaM   = array(rep(diag(p), iter), dim = c(p, p, iter))
    
    ## Mediator Model Initial Value
    alpha0 = matrix(rep(0, p * n2), n2, p, byrow = T)
    alphaX = matrix(rep(0, p * c) , p, c,  byrow = T)
    for (m in 1:p) {
      alpha0[,m] = rep(0, n2)
      alphaX[m,] = rep(0, c)
    }
    Alpha0 = rep(0, p)
    Salpha0 = diag(p)
    Mintercept = matrix(rep(0, p * n), n, p, byrow = T)
    
    PalphaT = 0
    PIPalphaT = rep(0, p)
    IalphaT = rep(0, p)
    alphaT = rep(0, p)
    SalphaT1 = 1
    SalphaT0 = 1
    
    sigmaM = diag(p)
    
    ## Mediator Model Sampling Algorithm
    for (i in 1:iter) {
      
      ## sampling alpha0
      for (j in 1:n2) {
        alpha0j = ifelse(Group == j, 1, 0)
        va0 = chol2inv(chol(chol2inv(chol(Salpha0)) + sum(alpha0j ^ 2) * chol2inv(chol(sigmaM))))
        sa0 = t(M - T%*%alphaT- X%*%t(alphaX)) %*% alpha0j
        ea0 = va0 %*% (chol2inv(chol(Salpha0)) %*% Alpha0 + chol2inv(chol(sigmaM)) %*% sa0) 
        alpha0[j,] = rmvn(1, as.vector(ea0), va0)
      }
      
      for (s in 1:n) {Mintercept[s,] = alpha0[Group[s],]}
      
      ## sampling Alpha0
      vA0 = Salpha0 / n2 
      eA0 = vA0  %*% chol2inv(chol(Salpha0)) %*% t(alpha0) %*% (rep(1,n2))
      Alpha0 = as.vector(rmvn(1, as.vector(eA0), vA0))
      
      ## sampling Salpha0
      ln = as.vector(rep(1,n2))
      ssalpha0 = t(alpha0 - ln %*% t(Alpha0)) %*% (alpha0 - ln %*% t(Alpha0))
      Salpha0 = riwish(p + 2 + n2, diag(p) + ssalpha0)
      
      ## sampling Indicator of alphaT
      for (m in 1:p) {
        a = dnorm(alphaT[m], 0, sqrt(SalphaT1)) * PalphaT
        b = dnorm(alphaT[m], 0, sqrt(SalphaT0)) * (1 - PalphaT)
        PIPalphaT[m] =  a / (a+b)
        IalphaT[m] = rbinom(1, 1, a / (a+b))
      }
      
      ## sampling alphaT
      vaT0 = diag(IalphaT) * SalphaT1 + diag(1 - IalphaT) * SalphaT0
      vaT = chol2inv(chol((chol2inv(chol(vaT0))) + sum(T^2) * chol2inv(chol(sigmaM))))
      saT = t(M - Mintercept - X %*% t(alphaX)) %*% T
      eaT = vaT %*% chol2inv(chol(sigmaM)) %*% saT 
      alphaT = as.vector(rmvn(1, as.vector(eaT), vaT))
      
      ## sampling alphaX
      for (e in 1:c) {
        vaX = sigmaM / sum(X[,e] ^ 2) 
        saX = t(M - Mintercept - T%*%alphaT - X%*%t(alphaX) + X[,e]%*%(t(alphaX[,e]))) %*% X[,e]
        eaX = vaX %*% chol2inv(chol(sigmaM)) %*% saX 
        alphaX[,e] = rmvn(1, as.vector(eaX), vaX)
      }
    
    ## sampling PalphaT
    PalphaT = rbeta(1, 1 + sum(IalphaT), 1 + p - sum(IalphaT))
    
    ## sampling SalphaT1 and SalphaT0
    SalphaT1 = rinvgamma(1, sum(IalphaT) / 2 + 1, sum((IalphaT * alphaT) ^ 2) / 2 + 1)
    SalphaT0 = rinvgamma(1, (p - sum(IalphaT)) / 2 + 1, sum(((1-IalphaT) * alphaT) ^ 2) / 2 + 10E-4)
    
    ## sampling sigmaM
    MSSR = t(M - Mintercept - T%*%t(alphaT) - X%*%t(alphaX)) %*% (M - Mintercept - T%*%t(alphaT) - X%*%t(alphaX))
    sigmaM = riwish(p + 2 + n, diag(p) + MSSR)
    
    ## Updated Value
    PIP.alphaT[i,] = PIPalphaT
    pos.alphaT[i,] = alphaT

    gc()
  }
  
    ###########################################################################
    ##################### Mediator Model selection ############################
    ###########################################################################
    PIPalphaT = colMeans(PIP.alphaT[-(1:burn),])
    names(PIPalphaT) = colnames(M)
    
    seleted1 = colnames(M)[which(PIPalphaT > 0.5)]
    M = M[,seleted1]
    p = ncol(M)
    
    ###########################################################################
    ########################## Outcome Model ##################################
    ###########################################################################
    
    #------- Outcome Model Parameters -------#
    #pos.Sbeta0  = matrix(rep(1, 1 * iter), iter, 1, byrow = TRUE)
    
    #pos.betaT   = matrix(rep(0, 1 * iter), iter, 1, byrow = TRUE)
    
    #pos.PbetaM  = matrix(rep(0, 1 * iter), iter, 1, byrow = TRUE)
    PIP.betaM   = matrix(rep(1, p * iter), iter, p, byrow = T)
    #pos.IbetaM  = matrix(rep(1, p * iter), iter, p, byrow = TRUE)
    pos.betaM   = matrix(rep(0, p * iter), iter, p, byrow = TRUE)
    #pos.SbetaM1 = matrix(rep(1, 1 * iter), iter, 1, byrow = TRUE)
    #pos.SbetaM0 = matrix(rep(1, 1 * iter), iter, 1, byrow = TRUE)
    
    #pos.sigmaY  = matrix(rep(1, 1 * iter), iter, 1, byrow = TRUE)
    
    
    ## Outcome Model Initial Value
    beta0 = rep(0, n2)
    Beta0 = 0
    Sbeta0 = 1
    Yintercept = rep(0,n)
    
    betaT = 0
    PbetaM = 0
    PIPbetaM = rep(0, p)
    IbetaM = rep(0, p)
    betaM = rep(0, p)
    SbetaM0 = 1
    SbetaM1 = 1
    
    betaX = rep(0, c)
    
    sigmaY = 1
    
    ## Outcome Model Sampling Algorithm
    for (i in 1:iter) {
      
      ## sampling beta0
      for (j in 1:n2) {
        beta0j = ifelse(Group == j, 1, 0)
        vb0 = (1 / Sbeta0 + sum(beta0j ^ 2) / sigmaY) ^ (-1)
        sb0 = sum(beta0j * (Y - T*betaT - M%*%betaM - X%*%betaX))
        eb0 = vb0 * (Beta0  / Sbeta0 + sb0 / sigmaY) 
        beta0[j] = rnorm(1, mean = eb0, sd = sqrt(vb0))
      }
      
      for (s in 1:n) {Yintercept[s] = beta0[Group[s]]}
      
      ## sampling Beta0
      vB0 = (10E-4 + n2 / Sbeta0) ^ (-1)
      eB0 = vB0 * n2 * mean(beta0)/Sbeta0
      Beta0 = rnorm(1, mean = eB0, sd = sqrt(vB0))
      
      ## sampling Sbeta0
      ssbeta0 = sum((beta0 - Beta0) ^ 2)
      Sbeta0 = rinvgamma(1, n2/2 + 10E-4, ssbeta0/2 + 10E-4)
      
      ## sampling betaT
      vbT = (10E-4 + sum(T^2) / sigmaY) ^ (-1)
      sbT = sum(T * (Y - Yintercept - M%*%betaM - X%*%betaX))
      ebT = vbT * sbT / sigmaY
      betaT = rnorm(1, mean = ebT, sd = sqrt(vbT))
      
      for (m in 1:p) {
        ## sampling Indicator of betaM
        a = dnorm(betaM[m], 0, sqrt(SbetaM1)) * PbetaM
        b = dnorm(betaM[m], 0, sqrt(SbetaM0)) * (1 - PbetaM)
        PIPbetaM[m] = a / (a + b)
        IbetaM[m] = rbinom(1, 1, a / (a + b))
        
        ## sampling betaM
        vbM0 = IbetaM[m] * SbetaM1 + (1 - IbetaM[m]) * SbetaM0
        vbM = (1 / vbM0 + sum(M[,m] ^ 2)/ sigmaY) ^ (-1)
        sbM = sum(M[,m] * (Y - Yintercept - T*betaT - X%*%betaX - M%*%betaM + M[,m]*betaM[m]))
        ebM = vbM * sbM / sigmaY
        betaM[m] = rnorm(1, mean = ebM, sd = sqrt(vbM))
      }
      
      ## sampling pi of betaM
      PbetaM = rbeta(1, 1 + sum(IbetaM), 1 + p - sum(IbetaM))
      
      ## sampling SbetaM1 and SbetaM0
      SbetaM1 = rinvgamma(1, sum(IbetaM) / 2 + 1, sum((IbetaM * betaM) ^ 2) / 2  + 1)
      SbetaM0 = rinvgamma(1, (p - sum(IbetaM)) / 2 + 1, sum(((1 - IbetaM) * betaM) ^ 2) / 2  + 10E-4)
      
      ## sampling betaX
      for (e in 1:c) {
        vbX = (10E-4 + sum(X[,e] ^ 2) / sigmaY) ^ (-1)
        sbX = sum(X[,e] * (Y - Yintercept - T*betaT - M%*%betaM - X%*%betaX + X[,e]*betaX[e]))
        ebX = vbX * sbX / sigmaY
        betaX[e] = rnorm(1, mean = ebX, sd = sqrt(vbX))
      }
      
      ## sampling sigmaY
      YSSR = sum((Y - Yintercept - T * betaT - M %*% betaM - X %*% betaX)^2)
      sigmaY = rinvgamma(1, n/2 + 10E-4, YSSR/2 + 10E-4)
      
      ## Outcome Model Updated Value
      PIP.betaM[i,] =  PIPbetaM 
      pos.betaM[i,] = betaM
      
      gc()

    }
    
    ###########################################################################
    ##################### Outcome Model selection ############################
    ###########################################################################
    PIPbetaM = colMeans(PIP.betaM[-(1:burn),])
    names(PIPbetaM) = colnames(M)
   
    
    ###########################################################################
    ############################### Output ####################################
    ###########################################################################
    
    posterior.sample = list(PIPbetaM, pos.betaM[-(1:burn),], 
                            PIPalphaT,pos.alphaT[-(1:burn),which(PIPalphaT > 0.5)])
    
    names(posterior.sample) = c("PIPbetaM", "betaM", "PIPalphaT","alphaT")
    
    posterior[[l]] = posterior.sample
  }
  
  return(posterior)
}




