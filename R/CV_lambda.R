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
#' @examples
#' bestlambda = bestPredLambda(o)
#'
#' @author Asmita Poddar & Florent Latimier


bestPredLambda = function(objPred,listLambdaS=c(0),listLambdaT = c(0))
{
  objPred@validation = FALSE
  if(objPred@fittedCov$modelname == "full"){
    listLambdaS = unique(c(listLambdaS , listLambdaT))
    listLambdaT = c(0)
  }
  else
  {
    if(objPred@fittedCov$spectra != "unknown"){
      listLambdaS = c(0)
    }
    if(objPred@fittedCov$time != "unknown"){
      listLambdaT = c(0)
    }
    if(objPred@fittedCov$spectra == "unknown" && objPred@fittedCov$time == "unknown")
    {
      p = predict(objPred)
      return(list(lambdaS = objPred@lamdbaS, lambdaT = objPred@lambdaT, predicted = p@perdict, percent = p@accuracy))
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
  perc = vapply(p,function(list){vapply(list,function(pred){pred@accuracy},FUN.VALUE = vector('double',length = 1))},FUN.VALUE = vector('double',length = length(p[[1]])))


  if(length(listLambdaS)==1 | length(listLambdaT)==1)
  {
    plot(x=listLambdaS,y=perc,type = 'l')
    title(paste(objPred@fittedCov$modelname,objPred@fittedCov$spectra,objPred@fittedCov$time))
    l = list(lambdaS = listLambdaS[which.max(perc)],lambdaT = listLambdaS[which.max(perc)],predicted = p[[which.max(perc)]][[1]]@predicted_labels,percent = max(perc))
  }
  else
  {
    l = list(lambdaS = listLambdaS[matxMax(perc)[2]],lambdaT = listLambdaT[matxMax(perc)[1]],predicted = p[[matxMax(perc)[2]]][[matxMax(perc)[1]]]@predicted_labels,percent = max(perc))
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
#' @examples
#' rowcol = matxMax(o)
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
#' @param objPred prediction object
#' @param listS list of lambda for spectra
#' @param listT list of lambda for time
#'
#' @return list of predicted lambda for spectra and time
#'
#' @examples
#' bestlambda = bestFitLambda(o,listS,listT)
#'
#' @author Asmita Poddar & Florent Latimier
#'
bestFitLambda = function(objFit, listS, listT)
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
  perc = vapply(listS, function(objFit, lambdaS, listT, model)
    {vapply(listT, fitChangLambda, objFit = objFit, lambdaS = lambdaS, model =model
            ,FUN.VALUE = vector('double',length(1)))}
    , objFit=objFit, listT=listT, model=model
    , FUN.VALUE = vector('double',length = length(listT)))

  list(lambdaS = listS[[matxMax(perc)[2]]], lambdaT = listS[[matxMax(perc)[1]]])
}


