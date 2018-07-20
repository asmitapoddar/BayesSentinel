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
tablesModel = function(data,lkernelS=list(),lkernelT=list(),lhS=list(),lhT=list(),listLambdaS=seq(from=0,to=10,by=0.1),listLambdaT=seq(from=0,to=10,by=0.1),listModel=list("gaussian")){
  l = testTrain(data,0.8)
  nr = 2 + length(lkernelS)
  nc = 2 + length(lkernelT)
  res = list()
  for(k in 1:length(listModel))
  {
    cov <- new("fitSpectra",l[[1]])
    cov <- fit(cov)
    pred <- new("predictClass",m=l[[2]],fittedCov=cov,model=listModel[[k]],validation=TRUE)
    pred <- predict(pred)
    full <- pred@accuracy
    p = matrix(ncol = nc, nrow = nr)
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
        cov <- new("fitSpectra",m=l[[1]],modelname="parsimonious",spectra=spectra,time=time,
                   kerneltypeSpectra=kerneltypeSpectra,kerneltypeTime=kerneltypeTime,h=h,s=s,
                   model=listModel[[k]],lambdaS=lambdaS,lambdaT=lambdaT,validation=TRUE)
        cov <- fit(cov)
        pred <- new("predictClass",m=l[[2]],fittedCov=cov,model=listModel[[k]],validation=TRUE)
        pred <- predict(pred)
        p[i,j] <- pred@accuracy
      }
    }

    rownames(p) <- c("diag","unknown",as.character(lkernelS))
    colnames(p) <- c("diag","unknown",as.character(lkernelT))
    noms = paste(c(full,p),"_",listModel[[k]],sep="")
    res = c(res,list(noms[1]=full,noms[2]=p))
  }

  res
}


