#-----------------------------------------------------------------------
#' Create a list with covariance matrices of a data with matrix observation.
#'
#' Return the covariance matrices
#'
#' @slot m the 3 dimentional data
#' @slot modelname name of model to be used for calculating the covariance matrix. Available models are
#' "full", "parsimonious". Default is "full".
#' @slot rows type of rows. Available models are "diag", "unknown" and "kernel".
#' Default is "diag".
#' @slot column type of column. Available models are "diag", "unknown" and "kernel".
#' Default is "diag".
#' @slot kerneltypeRow kernel to be used for covariance matrix of rows
#' Available kernels are "epanechnikov", "gaussian", "exponential", "uniform",
#' "quadratic", "circular", "triangular", "rational quadratic", "inverse multiquadratic".
#' Default is "exponential".
#' @slot kerneltypeCol kernel to be used for covariance matrix of column
#' Available kernels are "epanechnikov", "gaussian", "exponential", "uniform",
#' "quadratic", "circular", "triangular", "rational quadratic", "inverse multiquadratic".
#' Default is "exponential".
#' @slot h used for kernel calculation
#' @slot s correction limit paramater for flip flop algorithm
#' @slot lambdaR regularisation for rows for flip flop algorithm
#' @slot lambdaC regularisation for rows for flip flop algorithm
#' @slot validation to optimize lambda in case of th model is : M = parsimonious, S=unknown, T=unknow
#' @slot listLambdaR list of lambdaR used in prediction in case validation is TRUE
#' @slot listLambdaC list of lambdaC used in prediction in case validation is TRUE
#' @slot model use in prediction in case of validation is TRUE
#' @slot covMat returning the covariance matrx
#'
#' @name fitData-class
#' @aliases fitData-class
#' @rdname fitData-class
#'
#' @author Asmita Poddar & Florent Latimier
#'
setClass(
  Class="fitData",
  representation( m                     = "list"
                  , modelname           = "character"
                  , rows             = "character"
                  , column                = "character"
                  , kerneltypeRow   = "character"
                  , kerneltypeCol      = "character"
                  , h                   = "numeric"
                  , s                   = "numeric"
                  , lambdaR             = "numeric"
                  , lambdaC             = "numeric"
                  , validation          = "logical"
                  , listLambdaR         = "numeric"
                  , listLambdaC         = "numeric"
                  , model               = "character"
                  , covMat              = "list"
  ),
  prototype( m                   = list(0)
             , modelname         = "full"
             , rows           = "diag"
             , column              = "diag"
             , kerneltypeRow = "exponential"
             , kerneltypeCol    = "exponential"
             , h                 = 10
             , s                 = 0.01
             , lambdaR           = 0.3
             , lambdaC           = 0.3
             , validation        = FALSE
             , listLambdaR       = seq(from=0.1,to=0.3,by=0.1)
             , listLambdaC       = seq(from=0.1,to=0.3,by=0.1)
             , model             = "gaussian"
  ),
  # validity function
  validity = function(object)
  {
    if (object@modelname != "full" && object@modelname != "parsimonious")
    { stop("modelname must be either \"full\", \"parsimonious\".")}
    if (object@rows != "diag" && object@rows != "unknown" && object@rows != "kernel")
    { stop("rows must be either \"diag\", \"unknown\", \"kernel\".")}
    if (object@column != "diag" && object@column != "unknown" && object@column != "kernel")
    { stop("column must be either \"diag\", \"unknown\", \"kernel\".")}
    if (object@kerneltypeRow != "epanechnikov" && object@kerneltypeRow !="gaussian"
        && object@kerneltypeRow != "exponential" && object@kerneltypeRow !="uniform"
        &&object@kerneltypeRow != "quadratic" && object@kerneltypeRow != "circular"
        &&object@kerneltypeRow !="triangular" && object@kerneltypeRow !="rational quadratic"
        &&object@kerneltypeRow !="inverse multiquadratic")
    { stop("wrong kerneltypeRow entered. ")}
    if (object@kerneltypeCol != "epanechnikov" && object@kerneltypeCol !="gaussian"
        && object@kerneltypeCol != "exponential" && object@kerneltypeCol !="uniform"
        &&object@kerneltypeCol != "quadratic" && object@kerneltypeCol != "circular"
        &&object@kerneltypeCol !="triangular" && object@kerneltypeCol !="rational quadratic"
        &&object@kerneltypeCol !="inverse multiquadratic")
    { stop("wrong kerneltypeCol entered. ")}
    if (round(object@h) != object@h)
    { stop("h must be an integer.")}
    if (object@lambdaR < 0 | object@lambdaC < 0 )
    { stop("lambdaR and lambdaC must be positive.")}
    if (object@validation != TRUE && object@validation != FALSE )
    { stop("validation is a logical argument.")}
    if(object@validation)
    { if(object@model != "gaussian" && object@model != "fisher")
    { stop("With validation, the model must be either \"gaussian\", \"fisher\".")}
    }
    return(TRUE)
  }
)

#' Method num.
#'
#' @name fit
#' @rdname fit-method
#' @exportMethod fit

setGeneric("fit",
           def=function(Object)
           {
             standardGeneric("fit")
           }
)

#' #' Method num.
#'
#' @param Object object to be input
#'
#' @rdname fit-method
#' @aliases fit

setMethod(
  f = "fit",
  signature = "fitData",
  definition=function(Object)
  {

    #return a list of the list of integer elements
    listof = function(list, int)
    { lapply(list, function(l,int){l[[int]]}, int = int)
    }

    weight = lapply(levels(factor(Object@m[[1]])), function(k,data){length(data[which(data==k)])/length(data)}
                    , data = Object@m[[1]])


    mean = meanData(Object@m)

    if (Object@modelname=="full")
    {
      covmat = full(Object@m)
    }
    if (Object@modelname=="parsimonious")
    {
      if (Object@rows == "diag")
      {
        covRow = diagRow(Object@m)
      }
      if (Object@rows == "unknown")
      {
        if(Object@column=="unknown")
        {
          if(Object@validation)
          {
            lambda = bestFitLambda(Object)
            Object@lambdaR = lambda[[1]]
            Object@lambdaC = lambda[[2]]
          }
          ff = flipflop(data = Object@m, lmean = mean, s=Object@s, lambdaR = Object@lambdaR,lambdaC = Object@lambdaC)
          covRow = listof(ff,1)
          covCol = listof(ff,2)
        }
        else
          covRow = unknownRow(Object@m)
      }
      if (Object@rows == "kernel")
      {
        covRow = kernelRow(Object@m, Object@kerneltypeRow, Object@h)
      }

      if (Object@column == "diag")
      {
        covCol = diagCol(Object@m)
      }
      if (Object@column == "unknown")
      {
        if(Object@rows!="unknown")
          covCol = unknownCol(Object@m)
      }
      if (Object@column == "kernel")
      {
        covCol = kernelCol(Object@m, Object@kerneltypeCol, Object@h)
      }
    }

    if (Object@modelname=="full")
    {
      Object@covMat = list (covR = covmat, covC = 1, modelname = "full", rows = Object@rows
                            , column = Object@column, weight=weight, mean=mean)
    }
    else
    {
      Object@covMat = list (covR = covRow, covC = covCol, modelname = Object@modelname
                            , rows = Object@rows, column = Object@column, weight=weight, mean=mean)
    }

    Object@covMat
  }
)

#-----------------------------------------------------------------------
#' Initialize an instance of a fitData S4 class.
#'
#' Initialization method of the fitData class.
#'
#' @param .Object object of class fitData
#' @param m the 3 dimentional data
#' @param modelname name of model to be used for calculating the covariance matrix. Available models are
#' "full", "parsimonious". Default is "full".
#' @param rows type of rows. Available models are "diag", "unknown" and "kernel".
#' Default is "diag".
#' @param column type of column. Available models are "diag", "unknown" and "kernel".
#' Default is "diag".
#' @param kerneltypeRow kernel to be used for covariance matrix of rows
#' Available kernels are "epanechnikov", "gaussian", "exponential", "uniform",
#' "quadratic", "circular", "triangular", "rational quadratic", "inverse multiquadratic".
#' Default is "exponential".
#' @param kerneltypeCol kernel to be used for covariance matrix of column
#' Available kernels are "epanechnikov", "gaussian", "exponential", "uniform",
#' "quadratic", "circular", "triangular", "rational quadratic", "inverse multiquadratic".
#' Default is "exponential".
#' @param h used for kernel calculation
#' @param s regularisation paramater for flip flop algorithm
#' @param lambdaR regularisation for rows for flip flop algorithm
#' @param lambdaC regularisation for rows for flip flop algorithm
#' @param validation to optimize lambda or not
#' @param listLambdaR list of lambdaR used in prediction in case validation is TRUE
#' @param listLambdaC list of lambdaC used in prediction in case validation is TRUE
#' @param model use in prediction in case of validation is TRUE
#'
#' @name initialize
#' @rdname initialize-method
#' @keywords internal
#'

setMethod(
  "initialize",
  "fitData",
  function(.Object, m = list(0), modelname = "full", rows = "diag", column = "diag"
           , kerneltypeRow = "exponential", kerneltypeCol    = "exponential"
           , h = 10, s = 0.01, lambdaR = 0.3, lambdaC = 0.3, validation = FALSE
           , listLambdaR = seq(from = 0.1, to=0.3, by=0.1)
           , listLambdaC = seq(from=0.1, to=0.3, by=0.1), model = "gaussian")
  { .Object@m = m
  .Object@modelname = modelname
  .Object@rows = rows
  .Object@column = column
  .Object@kerneltypeRow = kerneltypeRow
  .Object@kerneltypeCol = kerneltypeCol
  .Object@h = h
  .Object@s = s
  .Object@lambdaR = lambdaR
  .Object@lambdaC = lambdaC
  .Object@validation = validation
  .Object@listLambdaR = listLambdaR
  .Object@listLambdaC = listLambdaC
  .Object@model = model
  return(.Object)
  }
)



#' fitDataMatrix
#'
#' Fit the mean and the covariance matrix according to the data
#'
#' @param .Object object of class fitData
#' @param m the 3 dimentional data
#' @param modelname name of model to be used for calculating the covariance matrix. Available models are
#' "full", "parsimonious". Default is "full".
#' @param rows type of rows. Available models are "diag", "unknown" and "kernel".
#' Default is "diag".
#' @param column type of column. Available models are "diag", "unknown" and "kernel".
#' Default is "diag".
#' @param kerneltypeRow kernel to be used for covariance matrix of rows
#' Available kernels are "epanechnikov", "gaussian", "exponential", "uniform",
#' "quadratic", "circular", "triangular", "rational quadratic", "inverse multiquadratic".
#' Default is "exponential".
#' @param kerneltypeCol kernel to be used for covariance matrix of column
#' Available kernels are "epanechnikov", "gaussian", "exponential", "uniform",
#' "quadratic", "circular", "triangular", "rational quadratic", "inverse multiquadratic".
#' Default is "exponential".
#' @param h used for kernel calculation
#' @param s regularisation paramater for flip flop algorithm
#' @param lambdaR regularisation for rows for flip flop algorithm
#' @param lambdaC regularisation for rows for flip flop algorithm
#' @param validation to optimize lambda or not
#' @param listLambdaR list of lambdaR used in prediction in case validation is TRUE
#' @param listLambdaC list of lambdaC used in prediction in case validation is TRUE
#' @param model use in prediction in case of validation is TRUE
#'
#' @return a list with covariance matrices for each clusther, of a data with matrix observation
#'
#' @name fitDataMatrix
#' @export fitDataMatrix
#'
fitDataMatrix <- function(m             = list(0)
                    , modelname         = "full"
                    , rows              = "diag"
                    , column            = "diag"
                    , kerneltypeRow     = "exponential"
                    , kerneltypeCol     = "exponential"
                    , h                 = 10
                    , s                 = 0.01
                    , lambdaR           = 0.3
                    , lambdaC           = 0.3
                    , validation        = FALSE
                    , listLambdaR       = seq(from=0.1,to=0.3,by=0.1)
                    , listLambdaC       = seq(from=0.1,to=0.3,by=0.1)
                    , model             = "gaussian")
{
  o = new("fitData", m = m, modelname = modelname, rows = rows, column = column, kerneltypeRow = kerneltypeRow
          , kerneltypeCol    = kerneltypeCol, h = h, s = s, lambdaR = lambdaR, lambdaC = lambdaC, validation = validation
          , listLambdaR  = listLambdaR, listLambdaC = listLambdaC, model = model)
  fit(o)
}
