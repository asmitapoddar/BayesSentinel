#-----------------------------------------------------------------------
#' Calculate accuracy of prediction
#'
#' Take two lists containing the predicted and actual labels of the data
#'
#' @param list1 list containing predicted labels
#' @param list2 list containing actual ground truth labels
#'
#' @examples
#' p = percent(p, m$labels)
#'
#' @return A floating value of the accuracy percentage
#'
#' @author Asmita Poddar & Florent Latimier
#'
#'

percent <- function(list1,list2){
  diff = list1 == list2
  sum(diff)/length(diff)*100
}

#-----------------------------------------------------------------------
#' Split data into training and test set
#'
#' Split data into training and test set depending on percentage set by the user
#'
#' @param data spectroscopic data
#' @param pTrain percentage of data to be used for training
#'
#' @examples
#' t = testTrain(data, 0.2)
#'
#' @return A list with training and test data
#'
#' @author Asmita Poddar & Florent Latimier
#'
#'

testTrain <- function(data,pTrain)
{
  n = length(data[[1]])
  sampl = runif(n)
  test <- lapply(data[[3]], function(mat,list){mat[list,]},list=which(sampl<pTrain))
  train <- lapply(data[[3]], function(mat,list){mat[list,]},list=which(sampl>=pTrain))
  test <- list(data[[1]][which(sampl<pTrain)],data[[2]],data[[3]],data[[4]],data[[5]],
               test,data[[7]][which(sampl<pTrain)])
  train <- list(data[[1]][which(sampl>=pTrain)],data[[2]],data[[3]],data[[4]],data[[5]],
                train,data[[7]][which(sampl>=pTrain)])
  list(test,train)
}

#-----------------------------------------------------------------------
#' Predict the label classes of the data
#'
#' Simulate one or more Gaussian spectra at regularly sampling time
#'
#' @param data spectroscopic data
#' @param modelname fitted covariance matrix for the data
#' @param spectra type of spectra. Available models are "diag", "unknown" and "kernel".
#' Default is "diag".
#' @param time type of time. Available models are "diag", "unknown" and "kernel".
#' Default is "diag".
#' @param kerneltypeSpectra kernel to be used for covariance matrix of spectra
#' Available kernels are "epanechnikov", "gaussian", "exponential", "uniform",
#' "quadratic", "circular", "triangular", "rational quadratic", "inverse multiquadratic".
#' Default is "exponential".
#' @param kerneltypeTime kernel to be used for covariance matrix of time
#' Available kernels are "epanechnikov", "gaussian", "exponential", "uniform",
#' "quadratic", "circular", "triangular", "rational quadratic", "inverse multiquadratic".
#' Default is "exponential".
#' @param h used for kernel calculation
#' @param lambda parameter for regularisation
#' @param model type of model to be used for prediction of labels
#' Available models are "gaussian", "fisher". Default is "gaussian".
#' @param pTrain percentage of data to be used for training.
#'
#' @examples
#' p = fitPredPrecent(m)
#'
#' @return the percentage of accuracy
#'
#' @author Asmita Poddar & Florent Latimier
#'

fitPredPrecent <- function(data, modelname = "full", spectra = "diag", time = "diag"
                           , kerneltypeSpectra = "exponential", kerneltypeTime = "exponential"
                           , h = 10, lambda=0.3, model="gaussian", pTrain=0.1
                           )
{
  l = testTrain(data,pTrain)
  percent(predict(l[[2]],fit(l[[1]],modelname,spectra,time,
                             kerneltypeSpectra,kerneltypeTime,h),lambda,model),l[[2]][[1]])
}

#-----------------------------------------------------------------------
#' Invert a matrix
#'
#' Invert a matrix
#'
#' @param mat matrix to be inverted
#'
#' @examples
#' i = inversion(mat)
#'
#' @return Inverted matrix
#'
#' @author Asmita Poddar & Florent Latimier
#'

inversion <- function(mat){
  cma <- chol(mat)
  chol2inv(cma)
}

