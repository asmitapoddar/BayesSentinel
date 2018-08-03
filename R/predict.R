
#-----------------------------------------------------------------------
#' Prediction
#'
#' Predict the label classes of the data with 3 dimensions
#'
#' @slot m the 3 dimentional data
#' @slot fittedCov fitted covariance matrix for the data
#' @slot lambdaR parameter for regularisation of row
#' @slot lambdaC parameter for regularisation of column
#' @slot model type of model to be used for prediction of labels
#' Available models are "gaussian", "tstudent". Default is "gaussian".
#' @slot validation logical to optimize the lambda.
#' @slot listLambdaR list of parameter for regularisation of row when do validation
#' @slot listLambdaC list of parameter for regularisation of column when do validation
#' @slot predicted_labels predicted class labels
#' @slot accuracy accracy of prediction
#'
#' @return A list with the row
#'
#' @author Asmita Poddar & Florent Latimier
#'
#' @name predictClass
#' @aliases predictClass-class
#' @rdname predictClass-class
#'

setClass(
  Class="predictClass",
  representation( m                     = "list"
                  , fittedCov           = "list"
                  , lambdaR             = "numeric"
                  , lambdaC             = "numeric"
                  , model               = "character"
                  , validation          = "logical"
                  , listLambdaR         = "numeric"
                  , listLambdaC         = "numeric"
                  , predicted_labels    = "integer"
                  , accuracy            = "numeric"
  ),
  prototype( m                   = list(0)
             , fittedCov         = list(0)
             , lambdaR           = 0.3
             , lambdaC           = 0.3
             , listLambdaR       = seq(from=0,to=1,by=0.1)
             , listLambdaC       = seq(from=0,to=1,by=0.1)
             , model = "gaussian"
             , validation        = FALSE
  ),
  # validity function
  validity = function(object)
  {
    #if (length(object@m)!=7)
    # stop("Enter correct format of data to be predicted.")
    #if (length(object@fittedCov)!=7)
    # stop("Enter correct format of covariance matrix to be predicted.")
    #if ( round(object@lambda) != object@lambda)
    # stop("lambda must be an integer.")
    if (object@model != "gaussian" && object@model !="tstudent")
    { stop("model must be either \"gaussian\", \"tstudent\".")}
    if (object@validation != TRUE && object@validation != FALSE)
    { stop("validation must be logical.")}
    return(TRUE)
  }
)

#' Method num.
#'
#' @name predict
#' @rdname predict-method
#' @exportMethod predict

setGeneric("predict",
           def=function(Object)
           {
             standardGeneric("predict")
           }
)

#' Method num.
#'
#' @param Object object to be input
#'
#' @rdname predict-method
#' @aliases predict

setMethod(
  f = "predict",
  signature = "predictClass",
  definition=function(Object)
  {
    if(Object@validation)
    {
      res = bestPredLambda(Object)
      Object@lambdaR = res$lambdaR
      Object@lambdaC = res$lambdaC
      Object@predicted_labels = res$predicted
      Object@accuracy = res$percent
    }

    else
    {
      mvnorm = function(data,reg,mean,weight, X)
      {
        weight * dmvnorm(X,mean,reg,log=TRUE)
      }

      nbLabel = length(unique(Object@m[[1]]))
      nbSample = length(Object@m[[1]])

      p = matrix(0, nbLabel, nbSample)
      weight  = Object@fittedCov$weight
      mean = Object@fittedCov$mean
      reg = regularisation(Object@fittedCov, Object@lambdaR, Object@lambdaC)

      powerLabelG <- function(data,inv,mean,weight,X)
      {
        X = X - mean
        power = rowSums((X %*% inv) * X)
        log(weight) +  log(sqrt(abs(det(inv)))) + (-power/2)
      }

      powerLabelF <-function(data, inv, mean, weight, X, nbRow, nbSample)
      {
        X = X - mean
        power = rowSums((X %*% inv) * X)
        log(weight) - 1/2*log(abs(det(inv))) -(nbRow*nbSample+3)/2*log(1+power)

      }

      power <- function(data,invers,mean,weight)
      {
        nbLabel = length( unique( Object@m[[1]] ) )
        nbSample = length( Object@m[[1]] )
        p = matrix(0, nbLabel, nbSample)
        if(Object@model == "gaussian"){powerL = powerLabelG}
        if(Object@model == "tstudent"){powerL = powerLabelF}
        if(Object@model == "mvnorm"){powerL = mvnorm}
        X = do.call('cbind',Object@m[[3]])
        for(i in 1:nbLabel)
        {
          if (Object@model == "tstudent")
            p[i,] = powerL(data,invers[[i]],mean[[i]],weight[[i]], X
                           , Object@m$nbRow, Object@m$nbSample)
          else
            p[i,] = powerL(data,invers[[i]],mean[[i]],weight[[i]], X)
        }
        p
      }

      Object@predicted_labels = max.col(t(power(Object@m,reg,mean,weight)))
      Object@accuracy = percent(Object@m[[1]], Object@predicted_labels)

    }

    Object
  }
)


#-----------------------------------------------------------------------
#' Initialize an instance of a predictClass S4 class.
#'
#' Initialization method of the predictClass class.
#'
#' @param .Object object of class predictClass
#' @param m the 3 dimentional data
#' @param fittedCov fitted covariance matrix for the data
#' @param lambdaR parameter for regularisation of row
#' @param lambdaC parameter for regularisation of column
#' @param model type of model to be used for prediction of labels
#' Available models are "gaussian", "tstudent". Default is "gaussian".
#' @param validation logical to optimize the lambda.
#' @param predicted_labels predicted class labels
#' @param accuracy accracy of prediction
#'
#' @name initialize
#' @rdname initialize-method
#' @keywords internal
#'
setMethod(
  "initialize",
  "predictClass",
  function(.Object, m = list(0), fittedCov = list(0), lambdaR = 0.3, lambdaC = 0.3, model = "gaussian"
           , validation = FALSE, listLambdaR = seq(from=0,to=1,by=0.1), listLambdaC = seq(from=0,to=1,by=0.1))
  { .Object@m = m
  .Object@fittedCov = fittedCov
  .Object@lambdaR = lambdaR
  .Object@lambdaC = lambdaC
  .Object@model = model
  .Object@validation = validation
  .Object@listLambdaR = listLambdaR
  .Object@listLambdaC = listLambdaC
  return(.Object)
  }
)


#' Prediction of labels
#'
#' Predict the label classes of the data according to the fitted arguments.
#'
#' @param m the 3 dimentional data
#' @param fittedCov fitted covariance matrix for the data
#' @param lambdaR parameter for regularisation of row
#' @param lambdaC parameter for regularisation of column
#' @param model type of model to be used for prediction of labels
#' Available models are "gaussian", "tstudent". Default is "gaussian".
#' @param validation logical to optimize the lambda.
#' @param listLambdaR list of parameter for regularisation of row when do validation
#' @param listLambdaC list of parameter for regularisation of column when do validation
#' @param predicted_labels predicted class labels
#' @param accuracy accracy of prediction
#'
#' @return A list with the row
#'
#' @name predictData
#' @export predictData

predictData <- function(...)
{
  o = new("predictClass", ...)
  predict(o)
}
