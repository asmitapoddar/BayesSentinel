percent <- function(list1,list2){
  diff = list1 == list2
  sum(diff)/length(diff)*100
}


# separera la data en 2 parties
testTrain <- function(data,pTrain){
  n = length(data[[1]])
  sampl = runif(n)
  test <- lapply(data[[3]], function(mat,list){mat[list,]},list=which(sampl<pTrain))
  train <- lapply(data[[3]], function(mat,list){mat[list,]},list=which(sampl>=pTrain))
  test <- list(data[[1]][which(sampl<pTrain)],data[[2]],test,data[[4]][which(sampl<pTrain)])
  train <- list(data[[1]][which(sampl>=pTrain)],data[[2]],train,data[[4]][which(sampl>=pTrain)])
  list(test,train)
}


# pourcent predict a partir d'une data
fitPredPrecent <- function(data,modelname = "full", spectra = "diag", time = "diag", kerneltypeSpectra = "exponential", kerneltypeTime = "exponential", h = 10,lambda=0.3,model="gaussian",pTrain=0.1){
  l = testTrain(data,pTrain)
  percent(predict(l[[2]],fit(l[[1]],modelname,spectra,time,kerneltypeSpectra,kerneltypeTime,h),lambda,model),l[[2]][[1]])
}


#### inversion matrice ###
inversion <- function(mat){
  chol2inv(chol(mat))
}




## tous les test possibles :
tablesModel = function(data,lkernelS=list(),lkernelT=list(),lhS=list(),lhT=list(),listLambdaS=seq(from=0,to=10,by=0.1),listLambdaT=seq(from=0,to=10,by=0.1)){
  l = testTrain(data,0.1)
  nr = 2 + length(lkernelS)
  nc = 2 + length(lkernelT)
  cov <- fit(l[[1]])
  cv <- bestLambda(data = l[[2]],fittedCov = cov,listLambdaS = listLambdaS,listLambdaT = listLambdaT,model = "gaussian")
  print(cv$lambdaS)
  fuulG <- cv$percent
  pG = matrix(ncol = nc, nrow = nr)
  for(i in 1:nr){
    if(i==1){
      spectra = "diag"
    }
    if(i==2){
      spectra = "unknown"
    }
    kerneltypeSpectra = ""
    h=0
    if(i>2){
      spectra = "kernel"
      kerneltypeSpectra = lkernelS[[i-2]]
      h = lhS[[i-2]]
    }

    for(j in 1:nc){
      if(j==1){
        time = "diag"
      }
      if(j==2){
        time = "unknown"
      }
      kerneltypeTime = ""
      if(j>2){
        time = "kernel"
        kerneltypeTime = lkernelT[[j-2]]
        h = lhT[[j-2]]
      }
      print(c(spectra,time))
      cov <- fit(l[[1]],"parsimonious",spectra,time,kerneltypeSpectra,kerneltypeTime,h,s,lambdaS,lambdaT)
      cv <- bestLambda(data = l[[2]],fittedCov = cov,listLambdaS = listLambdaS,listLambdaT = listLambdaT,model = "gaussian")
      print(c(cv$lambdaS,cv$lambdaT))
      pG[i,j] <- cv$percent
    }
  }
  rownames(pG) <- c("diag","unknown",as.character(lkernelS))
  colnames(pG) <- c("diag","unknown",as.character(lkernelT))

  # pour fisher
  cov <- fit(l[[1]])
  cv <- bestLambda(data = l[[2]],fittedCov = cov,listLambdaS = listLambdaS,listLambdaT = listLambdaT,model = "fisher")
  print(cv$lambdaS)
  fuulF <- cv$percent
  pF = matrix(ncol = nc, nrow = nr)
  for(i in 1:nr){
    if(i==1){
      spectra = "diag"
    }
    if(i==2){
      spectra = "unknown"
    }
    kerneltypeSpectra = ""
    h=0
    if(i>2){
      spectra = "kernel"
      kerneltypeSpectra = lkernelS[[i-2]]
      h = lhS[[i-2]]
    }

    for(j in 1:nc){
      if(j==1){
        time = "diag"
      }
      if(j==2){
        time = "unknown"
      }
      kerneltypeTime = ""
      if(j>2){
        time = "kernel"
        kerneltypeTime = lkernelT[[j-2]]
        h = lhT[[j-2]]
      }
      print(c(spectra,time))
      cov <- fit(l[[1]],"parsimonious",spectra,time,kerneltypeSpectra,kerneltypeTime,h)
      cv <- bestLambda(data = l[[2]],fittedCov = cov,listLambdaS = listLambdaS,listLambdaT = listLambdaT,model = "fisher")
      print(c(cv$lambdaS,cv$lambdaT))
      pF[i,j] <- cv$percent
    }
  }
  rownames(pF) <- c("diag","unknown",as.character(lkernelS))
  colnames(pF) <- c("diag","unknown",as.character(lkernelT))

  list(fullG=fullG,pG=pG,fullF=fullF,pF=pF)
}


