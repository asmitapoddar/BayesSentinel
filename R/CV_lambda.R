#-----------------------------------------------------------------------
#' Find best fitted lambda
#'
#' Find best fitted lambda
#'
#' @param objPred prediction object
#' @param listLambdaS list of lambda for spectra
#' @param listLambdaT list of lambda for time
#'
#' @return list of predicted lambda for spectra and time
#'
#' @author Asmita Poddar & Florent Latimier


bestPredLambda = function(objPred)
{
  objPred@validation = FALSE
  if(objPred@fittedCov$modelname == "full"){
    listLambdaS = unique(c(objPred@listLambdaS , objPred@listLambdaT))
    listLambdaT = c(0)
  }
  else
  {
    listLambdaS = objPred@listLambdaS
    listLambdaT = objPred@listLambdaT

    if(objPred@fittedCov$spectra != "unknown"){
      listLambdaS = c(0)
    }
    if(objPred@fittedCov$time != "unknown"){
      listLambdaT = c(0)
    }
    if(objPred@fittedCov$spectra == "unknown" && objPred@fittedCov$time == "unknown")
    {
      p = predict(objPred)
      return(list(lambdaS = objPred@lambdaS, lambdaT = objPred@lambdaT, predicted = p@predicted_labels, percent = p@accuracy))
    }
  }

  predict2 = function(objPred,lambdaS,listLambdaT)
  {
    objPred@lambdaS = lambdaS
    lapply(listLambdaT,
           function(objPred,lambdaT)
           { objPred@lambdaT = lambdaT
           predict(objPred)
           }
           ,objPred = objPred)
  }
  p = lapply(listLambdaS,predict2, objPred = objPred,listLambdaT=listLambdaT)
  print(length(p))
  print(length(p[[1]]))
  perc = vapply(p,function(list){vapply(list,function(pred){pred@accuracy},FUN.VALUE = vector('double',length = 1))},FUN.VALUE = vector('double',length = length(p[[1]])))
  print(perc)


  lambda = c(0)
  if(length(listLambdaS)==1 && length(listLambdaT)!=1)
  {
    print(which.max(perc))
    lambda = listLambdaT
    plot(x=lambda,y=perc,type = 'l')
    title(paste(objPred@fittedCov$modelname,objPred@fittedCov$spectra,objPred@fittedCov$time))
    l = list(lambdaS = lambda[which.max(perc)],lambdaT = lambda[which.max(perc)],predicted = p[[1]][[which.max(perc)]]@predicted_labels,percent = max(perc))
  }
  if(length(listLambdaS)!=1 && length(listLambdaT)==1)
  {
    lambda = listLambdaS
    plot(x=lambda,y=perc,type = 'l')
    title(paste(objPred@fittedCov$modelname,objPred@fittedCov$spectra,objPred@fittedCov$time))
    l = list(lambdaS = lambda[which.max(perc)],lambdaT = lambda[which.max(perc)],predicted = p[[which.max(perc)]][[1]]@predicted_labels,percent = max(perc))
  }
  if(length(listLambdaS)==1 && length(listLambdaT)==1)
  {
    l = list(lambdaS = 0,lambdaT = 0,predicted = p[[1]][[1]]@predicted_labels,percent = perc)
  }

  l
}

#-----------------------------------------------------------------------
#' Find row and column of maximum element
#' Find row and column of maximum element
#'
#' @param mat Matrix
#'
#' @return Vector containing the row and column of maximum element
#'
#' @author Asmita Poddar & Florent Latimier


matxMax <- function(mat)
{
  m = which.max(mat)
  colmn <- (m-1) %/% nrow(mat) + 1
  row <- (m-1) %% nrow(mat) + 1
  c(row, colmn)
}

#-----------------------------------------------------------------------
#' Find best fitted lambda
#'
#' Find best fitted lambda
#'
#' @param objFit prediction object
#' @param listS list of lambda for spectra
#' @param listT list of lambda for time
#'
#' @return list of predicted lambda for spectra and time
#'
#' @author Asmita Poddar & Florent Latimier
#'


bestFitLambda = function(objFit)
{
  objFit@validation = FALSE
  fitChangLambda = function(objFit, lambdaS, lambdaT, model)
  {
    objFit@lambdaS = lambdaS
    objFit@lambdaT = lambdaT
    p = new("predictClass", objFit@m,fit(objFit), lambdaS = lambdaS, lambdaT = lambdaT
            , model= model)
    predict(p)@accuracy
  }

  model = objFit@model
  perc = vapply(objFit@listLambdaS, function(objFit, lambdaS, model)
    {vapply(objFit@listLambdaT, fitChangLambda, objFit = objFit, lambdaS = objFit@lambdaS, model =model
            ,FUN.VALUE = vector('double',length(1)))}
    , objFit=objFit, model=model
    , FUN.VALUE = vector('double',length = length(objFit@listLambdaT)))

  list(lambdaS = objFit@listLambdaS[[matxMax(perc)[2]]], lambdaT = objFit@listLambdaT[[matxMax(perc)[1]]])
}


