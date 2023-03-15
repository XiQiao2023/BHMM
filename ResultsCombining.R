setwd("/Users/xiqiao/Desktop/Server/simulation3")
p=500

data_dir <- getwd()

# get a list of all the .RData files in the directory
data_files <- list.files(data_dir, full.names = TRUE)

# create an empty list to store the loaded data
tmp <- list()

# loop through each .RData file and load its contents into the list
for (file in data_files) {
  # load the object into the list with the same name as the file
  tmp[[file]] <- get(load(file))
}


n1=2
n2=100
c=2


save(results,file = paste0("M", p, "A", 0.1 * p, "N", n1 * n2,"C",c,".Rdata"))



## Organize Rcpp Results
convertResults = function(pos){
  pos$PIPbetaM = setNames(as.vector(pos$PIPbetaM)[["PIPbetaM"]],row.names(pos$PIPbetaM))
  pos$PIPalphaT = setNames(as.vector(pos$PIPalphaT)[["PIPalphaT"]],row.names(pos$PIPalphaT))
  names(pos) = c("PIPbetaM", "betaM", "PIPalphaT","alphaT")
  return(pos)
}


results =  list()
for (i in 1:100) {
  results[[i]] = lapply(tmp[[i]], function(list) convertResults(list))
}

save(results,file = paste0("M", p, "A", 0.1 * p, "N", n1 * n2,"C",c,".Rdata"))





#############################################
## Threshold Determination at Certain FDR ###
##### For Nested Model corrected PIP ########
#############################################

getPIP = function(pos){
  PIP = pos$PIPbetaM
  return(PIP)
}


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



allPIP = lapply(results, sapply, getPIP)

p=500
FPR = 0.01
rate = lapply(allPIP, FDRAdjust)
dfrate = as.data.frame(matrix(unlist(rate),length(rate),5,byrow=T))
colnames(dfrate) = c("Threshhold","FDR","FPR","TPR","TNR")
colMeans(na.omit(dfrate))




#############################################
## Threshold Determination at Certain FDR ###
########### For Uni-variate Model ############
#############################################

FDRAdjust = function(result){
  
  Pvalues=result["ACME_pvalue",]
  allm = paste0("cg",1:p)
  trueactive = c(rep(1, p * 0.1), rep(0, 0.9 * p))
  
  FDRTable = as.data.frame(cbind(Pvalues,rep(0,length(Pvalues)),rep(0,length(Pvalues)),rep(0,length(Pvalues)),rep(0,length(Pvalues))))
  
  colnames(FDRTable) = c("P_thresh","FDR","FPR","TPR","TNR")
  
  for (i in 1:length(Pvalues)) {
    threshhold = FDRTable$P_thresh[i]
    postactiveCpG = rownames(FDRTable)[which(FDRTable$P_thresh <= threshhold)]
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
  finalthresh = max(FDRcontrol[,"P_thresh"])
  Decision = FDRcontrol[FDRcontrol[,"P_thresh"] == finalthresh,][1,]
  
  return(Decision)
}


p=1000
FPR = 0.01
rate = lapply(results, FDRAdjust)
dfrate = as.data.frame(matrix(unlist(rate),length(rate),5,byrow=T))
colnames(dfrate) = c("Threshhold","FDR","FPR","TPR","TNR")
colMeans(na.omit(dfrate))


