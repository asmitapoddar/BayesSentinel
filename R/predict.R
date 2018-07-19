
#-----------------------------------------------------------------------
#' Predict the label classes of the data
#'
#' Simulate one or more Gaussian spectra at regularly sampling time
#'
#' @slot m spectroscopic data
#' @slot fittedCov fitted covariance matrix for the data
#' @slot lambdaS parameter for regularisation of spectra
#' @slot lambdaT parameter for regularisation of time
#' @slot model type of model to be used for prediction of labels
#' Available models are "gaussian", "fisher". Default is "gaussian".
#' @slot validation logical to optimize the lambda.
#' @slot predicted_labels predicted class labels
#' @slot accuracy accracy of prediction
#'
#' @examples
#' p = predict(m, fittedCov)
#'
#' @return A list with the spectra
#'
#' @author Asmita Poddar & Florent Latimier
#'

setClass(
  Class="predictClass",
  representation( m                     = "list"
                  , fittedCov           = "list"
                  , lambdaS             = "numeric"
                  , lambdaT             = "numeric"
                  , model               = "character"
                  , validation          = "logical"
                  , predicted_labels    = "numeric"
                  , accuracy            = "numeric"
  ),
  prototype( m                   = list(0)
             , fittedCov         = list(0)
             , lambdaS           = 0.3
             , lambdaT           = 0.3
             , model             = "gaussian"
             , validation        = FALSE
             , predicted_labels  = 0
             , accuracy          = 0
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
    if (object@model != "gaussian" && object@model !="fisher")
    { stop("model must be either \"gaussian\", \"fisher\".")}
    if (object@validation != TRUE && object@validation != FALSE)
    { stop("validation must be logical.")}
    return(TRUE)
  }
)

setGeneric("predict",
           def=function(Object)
           {
             standardGeneric("predict")
           }
)

setMethod(
  f = "predict",
  signature = "predictClass",
  definition=function(Object)
  {
    if(Object@validation)
    {
      res = bestPredLambda(Object,listLambdaS = seq(from=0,to=10,by=0.1),listLambdaT = seq(from=0,to=10,by=0.1))
      Object@lambdaS = res$lambdaS
      Object@lambdaT = res$lambdaT
      Object@predicted_labels = res$perdict
      Object@accuracy = res$percent
    }

    else
    {
      mvnorm = function(data,reg,mean,weight, X)
      {
        weight * dmvnorm(X,mean,reg,log=TRUE)
      }

      nbLabel = length(unique(Object@m[[1]]))
      nbPixel = length(Object@m[[1]])

      p = matrix(0, nbLabel, nbPixel)
      weight  = Object@fittedCov$weight
      mean = Object@fittedCov$mean
      reg = regularisation(Object@fittedCov, Object@lambdaS, Object@lambdaT)

      powerLabelG <- function(data,inv,mean,weight,X)
      {
        X = X - mean
        power = rowSums((X %*% inv) * X)
        log(weight) +  log(sqrt(abs(det(inv)))) + (-power/2)
      }

      powerLabelF <-function(data,inv,mean,weight,X)
      {
        X = X - mean
        power = rowSums((X %*% inv) * X)
        log(weight) - 1/2*log(abs(det(inv))) -165*log(1+power)

      }

      power <- function(data,invers,mean,weight)
      {
        nbLabel = length( unique( Object@m[[1]] ) )
        nbPixel = length( Object@m[[1]] )
        p = matrix(0, nbLabel, nbPixel)
        if(Object@model == "gaussian"){powerL = powerLabelG}
        if(Object@model == "fisher"){powerL = powerLabelF}
        if(Object@model == "mvnorm"){powerL = mvnorm}
        X = do.call('cbind',Object@m[[3]])
        for(i in 1:nbLabel)
          p[i,] = powerL(data,invers[[i]],mean[[i]],weight[[i]], X)
        p
      }

      Object@predicted_labels = max.col(t(power(Object@m,reg,mean,weight)))
      Object@accuracy = percent(Object@m$labels, Object@predicted_labels)

    }

    Object
  }
)

setMethod(
  "initialize",
  "predictClass",
  function(.Object, m = list(0), fittedCov = list(0), lambdaS = 0.3, lambdaT = 0.3
           , model = "gaussian", validation = FALSE)
  { .Object@m = m
  .Object@fittedCov = fittedCov
  .Object@lambdaS = lambdaS
  .Object@lambdaT = lambdaT
  .Object@model = model
  .Object@validation = validation
  return(.Object)
  }
)
