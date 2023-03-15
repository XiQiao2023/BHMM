############################################################################################
######################## Multivariate multilevel mediation analysis ########################
############################## Bayesian Variable Selection #################################
################################## Correlated mediators ####################################
###################################### Simulate All ########################################
############################################################################################

Gibbs.Posterior = function(iteration, burn, chain){
  
  posterior = vector(mode='list', length = chain)
  
  for (l in 1:chain) {
    
    iter = iteration + burn
    
    #------- Outcome Model Parameters -------#
    pos.Beta0   = matrix(rep(0, 1 * iter), iter, 1, byrow = T)
    pos.Sbeta0  = matrix(rep(1, 1 * iter), iter, 1, byrow = T)
    
    pos.betaT   = matrix(rep(0, 1 * iter), iter, 1, byrow = T)
    
    pos.PbetaM  = matrix(rep(0, 1 * iter), iter, 1, byrow = T)
    PIP.betaM   = matrix(rep(1, p * iter), iter, p, byrow = T)
    pos.IbetaM  = matrix(rep(1, p * iter), iter, p, byrow = T)
    pos.betaM   = matrix(rep(0, p * iter), iter, p, byrow = T)
    pos.SbetaM1 = matrix(rep(1, 1 * iter), iter, 1, byrow = T)
    pos.SbetaM0 = matrix(rep(1, 1 * iter), iter, 1, byrow = T)
    
    pos.sigmaY  = matrix(rep(1, 1 * iter), iter, 1, byrow = T)
    
    #-------- Mediator Model Parameters -------#
    
    pos.Alpha0   = matrix(rep(0, p * iter), iter, p, byrow = T)
    pos.Salpha0  = matrix(rep(1, 1 * iter), iter, 1, byrow = T)
    
    pos.PalphaT  = matrix(rep(0, 1 * iter), iter, 1, byrow = T)
    PIP.alphaT   = matrix(rep(1, p * iter), iter, p, byrow = T)
    pos.IalphaT  = matrix(rep(1, p * iter), iter, p, byrow = T)
    pos.alphaT   = matrix(rep(0, p * iter), iter, p, byrow = T)
    pos.SalphaT1 = matrix(rep(1, 1 * iter), iter, 1, byrow = T)
    pos.SalphaT0 = matrix(rep(1, 1 * iter), iter, 1, byrow = T)
    
    pos.sigmaM   = matrix(rep(1, 1 * iter), iter, 1, byrow = T)
    
    ###########################################################################
    ########################## Gibbs sampling #################################
    ###########################################################################
    
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
    
    ## Mediator Model Initial Value
    alpha0 = matrix(rep(0, p * n2), p, n2, byrow = T)
    alphaX = matrix(rep(0, p * c) , p, c,  byrow = T)
    for (m in 1:p) {
      alpha0[m,] = rep(0, n2)
      alphaX[m,] = rep(0, c)
    }
    Alpha0 = rep(0, p)
    Salpha0 = 1
    Mintercept = matrix(rep(0, p * n), n, p, byrow = T)
    
    PalphaT = 0
    PIPalphaT = rep(0, p)
    IalphaT = rep(0, p)
    alphaT = rep(0, p)
    SalphaT1 = 1
    SalphaT0 = 1
    
    sigmaM = 1
    
    for (i in 1:iter) {
      
      ## Outcome Model Sampling Algorithm
      
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
        IbetaM[m] = rbinom(1, 1, PIPbetaM[m])
        
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
      
      ## Mediator Model Sampling Algorithm
      
      for (m in 1:p) {
        
        ## sampling alpha0
        for (j in 1:n2) {
          alpha0j = ifelse(Group == j, 1, 0)
          va0 = (1 / Salpha0 + sum(alpha0j ^ 2) / sigmaM) ^ (-1)
          sa0 = sum(alpha0j * (M[,m] - T*alphaT[m] - X%*%alphaX[m,]))
          ea0 = va0 * (Alpha0[m] / Salpha0 + sa0 / sigmaM) 
          alpha0[m,j] = rnorm(1, mean = ea0, sd = sqrt(va0))
        }
        
        for (s in 1:n) {Mintercept[s, m] = alpha0[m, Group[s]]}
        
        ## sampling Alpha0
        vA0 = (10E-4 + n2 / Salpha0) ^ (-1)
        eA0 = vA0 * n2 * mean(alpha0[m, ]) / Salpha0
        Alpha0[m] = rnorm(1, mean = eA0, sd = sqrt(vA0))
        
        ## sampling Indicator of alphaT
        a = dnorm(alphaT[m], 0, sqrt(SalphaT1)) * PalphaT
        b = dnorm(alphaT[m], 0, sqrt(SalphaT0)) * (1 - PalphaT)
        PIPalphaT[m] =  a / (a+b)
        IalphaT[m] = rbinom(1, 1, PIPalphaT[m])
        
        ## sampling alphaT
        vaT0 = IalphaT[m] * SalphaT1 + (1 - IalphaT[m]) * SalphaT0
        vaT = (1 / vaT0 + sum(T^2) / sigmaM) ^ (-1)
        saT = sum(T * (M[,m] - Mintercept[,m] - X %*% alphaX[m,]))
        eaT = vaT * saT / sigmaM
        alphaT[m] = rnorm(1, mean = eaT, sd = sqrt(vaT))
        
        ## sampling alphaX
        for (e in 1:c) {
          vaX = (10E-4 + sum(X[,e] ^ 2) / sigmaM) ^ (-1)
          saX = sum(X[,e] * (M[,m] - Mintercept[,m] - T*alphaT[m] - X%*%alphaX[m,] + X[,e]*alphaX[m,e]))
          eaX = vaX * saX /sigmaM
          alphaX[m,e] = rnorm(1, mean = eaX, sd = sqrt(vaX))
        }
      }
      
      ## sampling Salpha0
      Alpha0Matrix = matrix(rep(Alpha0[m],n2), p, n2,byrow = FALSE)
      ssalpha0 = sum((alpha0 - Alpha0Matrix) ^ 2)
      Salpha0 = rinvgamma(1, n2 * p /2 + 10E-4, ssalpha0/2 + 10E-4)
      
      ## sampling PalphaT
      PalphaT = rbeta(1, 1 + sum(IalphaT), 1 + p - sum(IalphaT))
      
      ## sampling SalphaT1 and SalphaT0
      SalphaT1 = rinvgamma(1, sum(IalphaT) / 2 + 1, sum((IalphaT * alphaT) ^ 2) / 2 + 1)
      SalphaT0 = rinvgamma(1, (p - sum(IalphaT)) / 2 + 1, sum(((1-IalphaT) * alphaT) ^ 2) / 2 + 10E-4)
      
      ## sampling sigmaM
      MSSR = sum((M - Mintercept - T%*%t(alphaT) - X%*%t(alphaX))^2)
      sigmaM = rinvgamma(1, n * p / 2 + 10E-4, MSSR/2 + 10E-4)
      
      ## Outcome Model Updated Value
      pos.Beta0[i,] = Beta0
      pos.Sbeta0[i,] = Sbeta0
      
      pos.betaT[i,] = betaT
      
      pos.PbetaM[i,] = PbetaM
      PIP.betaM[i,] =  PIPbetaM 
      pos.IbetaM[i,] =  IbetaM
      pos.betaM[i,] = betaM
      pos.SbetaM1[i,] = SbetaM1
      pos.SbetaM0[i,] = SbetaM0
      
      pos.sigmaY[i,] = sigmaY
      
      ## Mediator Model Previous Value
      pos.Alpha0[i,] = Alpha0
      pos.Salpha0[i,] = Salpha0
      
      pos.PalphaT[i,] = PalphaT
      PIP.alphaT[i,] = PIPalphaT
      pos.IalphaT[i,] = IalphaT
      pos.alphaT[i,] = alphaT
      pos.SalphaT1[i,] = SalphaT1
      pos.SalphaT0[i,] = SalphaT0
      
      pos.sigmaM[i,] = sigmaM 
    }
    
    ###########################################################################
    ############################### Selection #################################
    ###########################################################################

    posterior.sample = list(pos.betaT[-(1:burn),],  pos.sigmaY[-(1:burn),],
                            pos.PbetaM[-(1:burn),], PIP.betaM[-(1:burn),], pos.betaM[-(1:burn),], 
                            pos.SbetaM1[-(1:burn),], pos.SbetaM0[-(1:burn),],
                            pos.sigmaM[-(1:burn),],
                            pos.PalphaT[-(1:burn),], PIP.alphaT[-(1:burn),], pos.alphaT[-(1:burn),],
                            pos.SalphaT1[-(1:burn),], pos.SalphaT0[-(1:burn),])
    
    names(posterior.sample) = c("betaT",  "sigmaY",
                                "PbetaM", " PIPbetaM", "betaM",   
                                "SbetaM1", "SbetaM0",
                                "sigmaM",
                                "PalphaT"," PIPalphaT","alphaT",  
                                "SalphaT1","SalphaT0")
    
    posterior[[l]] = posterior.sample
  }

  return(posterior)
}




