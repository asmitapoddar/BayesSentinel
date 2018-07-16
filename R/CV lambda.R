bestLambda = function(data,fittedCov,listLambdaS=c(0),listLambdaT = c(0),model = "gaussian"){
  if(fittedCov$modelname == "full"){
    listLambdaS = unique(c(listLambdaS , listLambdaT))
    listLambdaT = c(0)
  }
  else{
    if(fittedCov$spectra != "unknown"){
      listLambdaS = c(0)
    }
    if(fittedCov$time != "unknown"){
      listLambdaT = c(0)
    }
  }
  p = lapply(X = listLambdaS,FUN = predict2, data=data,fittedCov=fittedCov,listLambdaT=listLambdaT,model=model)
  perc = vapply(p, percent2,list2=data[[1]],FUN.VALUE = vector('double',length = length(listLambdaT)))

  if(length(listLambdaS)==1 | length(listLambdaS)==1){
    plot(x=listLambda,y=perc,type = 'l')
    title(paste(fittedCov$modelname,fittedCov$spectra,fittedCov$time))
  }

  list(lambdaS = listLambdaS[matxMax(perc)[1]],lambdaT = listLambdaT[matxMax(perc)[2]],predict = p[[matxMax(perc)[1]]][[matxMax(perc)[2]]],percent = max(perc))
}



predict2 = function(data,fittedCov,lambdaS,listLambdaT,model){
  lapply(listLambdaT,predict,m=data,fittedCov=fittedCov,lambdaS=lambdaS, model=model)
}

percent2 = function(Llist,list2){
  vapply(Llist,percent,list2 = list2,FUN.VALUE = vector('double',1))
}

matxMax <- function(mat)
{
  m = which.max(mat)
  colmn <- (m-1) %/% nrow(mat) + 1
  row <- (m-1) %% nrow(mat) + 1
  c(row, colmn)
}





bestFitLambda = function(objFit,listS,listT,model){
  perc = vapply(listS, function(objFit,lambdaS,listT,model){vapply(listT,fitChangLambda,objFit = objFit,lambdaS = lambdaS,model =model,FUN.VALUE = vector('double',length(1)))},objFit=objFit,listT=listT,model=model,FUN.VALUE = vector('double',length = length(listT)))
  list(lambdaS = listS[[matxMax(perc)[1]]], lambdaT = listS[[matxMax(perc)[2]]])
}

fitChangLambda = function(objFit,lambdaS,lambdaT,model)
{
  objFit@lambdaS = lambdaS
  objFit@lambdaT = lambdaT
  p = new("predictClass",objFit@m,fit(objFit),lambdaS = lambdaS, lambdaT = lambdaT, model= model)
  predict(p)@accuracy
}
