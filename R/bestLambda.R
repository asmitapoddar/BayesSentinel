#-----------------------------------------------------------------------
#' Best regularisation for prediction
#'
#' Find best fitted lambda
#'
#' @param objPred prediction object
#'
#' @return list of predicted lambda for row and column and the prediction with accuracy
#'
#' @author Asmita Poddar & Florent Latimier


bestPredLambda = function(objPred)
{
  objPred@validation = FALSE
  if(objPred@fittedCov$modelname == "full"){
    listLambdaR = unique(c(objPred@listLambdaR , objPred@listLambdaC))
    listLambdaC = c(0)
  }
  else
  {
    listLambdaR = objPred@listLambdaR
    listLambdaC = objPred@listLambdaC

    if(objPred@fittedCov$row != "unknown"){
      listLambdaR = c(0)
    }
    if(objPred@fittedCov$column != "unknown"){
      listLambdaC = c(0)
    }
    if(objPred@fittedCov$row == "unknown" && objPred@fittedCov$column == "unknown")
    {
      p = predict(objPred)
      return(list(lambdaR = objPred@lambdaR, lambdaC = objPred@lambdaC, predicted = p@predicted_labels, percent = p@accuracy))
    }
  }

  predict2 = function(objPred,lambdaR,listLambdaC)
  {
    objPred@lambdaR = lambdaR
    lapply(listLambdaC,
           function(objPred,lambdaC)
           { objPred@lambdaC = lambdaC
           predict(objPred)
           }
           ,objPred = objPred)
  }
  p = lapply(listLambdaR,predict2, objPred = objPred,listLambdaC=listLambdaC)
  print(length(p))
  print(length(p[[1]]))
  perc = vapply(p,function(list){vapply(list,function(pred){pred@accuracy},FUN.VALUE = vector('double',length = 1))},FUN.VALUE = vector('double',length = length(p[[1]])))
  print(perc)


  lambda = c(0)
  if(length(listLambdaR)==1 && length(listLambdaC)!=1)
  {
    print(which.max(perc))
    lambda = listLambdaC
    plot(x=lambda,y=perc,type = 'l')
    title(paste(objPred@fittedCov$modelname,objPred@fittedCov$row,objPred@fittedCov$column))
    l = list(lambdaR = lambda[which.max(perc)],lambdaC = lambda[which.max(perc)],predicted = p[[1]][[which.max(perc)]]@predicted_labels,percent = max(perc))
  }
  if(length(listLambdaR)!=1 && length(listLambdaC)==1)
  {
    lambda = listLambdaR
    plot(x=lambda,y=perc,type = 'l')
    title(paste(objPred@fittedCov$modelname,objPred@fittedCov$row,objPred@fittedCov$column))
    l = list(lambdaR = lambda[which.max(perc)],lambdaC = lambda[which.max(perc)],predicted = p[[which.max(perc)]][[1]]@predicted_labels,percent = max(perc))
  }
  if(length(listLambdaR)==1 && length(listLambdaC)==1)
  {
    l = list(lambdaR = 0,lambdaC = 0,predicted = p[[1]][[1]]@predicted_labels,percent = perc)
  }

  l
}



#-----------------------------------------------------------------------
#' matxMax
#'
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
#' bestFitLambda
#'
#' Find best fitted lambda for row and column
#'
#' @param objFit fit object
#'
#' @return list of covariance matrix with the best lambda for row and column
#'
#' @author Asmita Poddar & Florent Latimier
#'
bestFitLambda = function(objFit)
{
  objFit@validation = FALSE
  fitChangLambda = function(objFit, lambdaR, lambdaC, model)
  {
    objFit@lambdaR = lambdaR
    objFit@lambdaC = lambdaC
    p = new("predictClass", objFit@m,fit(objFit), lambdaR = lambdaR, lambdaC = lambdaC
            , model= model)
    predict(p)@accuracy
  }

  model = objFit@model
  perc = vapply(objFit@listLambdaR, function(objFit, lambdaR, model)
  {vapply(objFit@listLambdaC, fitChangLambda, objFit = objFit, lambdaR = objFit@lambdaR, model =model
          ,FUN.VALUE = vector('double',length(1)))}
  , objFit=objFit, model=model
  , FUN.VALUE = vector('double',length = length(objFit@listLambdaC)))

  list(lambdaR = objFit@listLambdaR[[matxMax(perc)[2]]], lambdaC = objFit@listLambdaC[[matxMax(perc)[1]]])
}


