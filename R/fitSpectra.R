#-----------------------------------------------------------------------
#' Create a list with covariance matrices of the spectra and times
#'
#' Return the covariance matrices
#'
#' @slot m spectroscopic data
#' @slot modelname name of model to be used for calculating the covariance matrix. Available models are
#' "full", "parsimonious". Default is "full".
#' @slot spectra type of spectra. Available models are "diag", "unknown" and "kernel".
#' Default is "diag".
#' @slot time type of time. Available models are "diag", "unknown" and "kernel".
#' Default is "diag".
#' @slot kerneltypeSpectra kernel to be used for covariance matrix of spectra
#' Available kernels are "epanechnikov", "gaussian", "exponential", "uniform",
#' "quadratic", "circular", "triangular", "rational quadratic", "inverse multiquadratic".
#' Default is "exponential".
#' @slot kerneltypeTime kernel to be used for covariance matrix of time
#' Available kernels are "epanechnikov", "gaussian", "exponential", "uniform",
#' "quadratic", "circular", "triangular", "rational quadratic", "inverse multiquadratic".
#' Default is "exponential".
#' @slot h used for kernel calculation
#' @slot s regularisation paramater for flip flop algorithm
#' @slot lambdaS regularisation for spectra for flip flop algorithm
#' @slot lambdaT regularisation for spectra for flip flop algorithm
#' @slot validation to optimize lambda or not
#' @slot model use in prediction in case of validation is TRUE
#' @slot covMat returning the covariance matrx
#'
#' @name fitSpectra
#' @aliases fitSpectra-class
#' @rdname fitSpectra-class
#'
#' @author Asmita Poddar & Florent Latimier
#'

setClass(
  Class="fitSpectra",
  representation( m                     = "list"
                  , modelname           = "character"
                  , spectra             = "character"
                  , time                = "character"
                  , kerneltypeSpectra   = "character"
                  , kerneltypeTime      = "character"
                  , h                   = "numeric"
                  , s                   = "numeric"
                  , lambdaS             = "numeric"
                  , lambdaT             = "numeric"
                  , validation          = "logical"
                  , model               = "character"
                  , covMat              = "list"
  ),
  prototype( m                   = list(0)
             , modelname         = "full"
             , spectra           = "diag"
             , time              = "diag"
             , kerneltypeSpectra = "exponential"
             , kerneltypeTime    = "exponential"
             , h                 = 10
             , s                 = 0.01
             , lambdaS           = 0.3
             , lambdaT           = 0.3
             , validation        = FALSE
             , model             = "gaussian"
  ),
  # validity function
  validity = function(object)
  {
    if (object@modelname != "full" && object@modelname != "parsimonious")
    { stop("modelname must be either \"full\", \"parsimonious\".")}
    if (object@spectra != "diag" && object@spectra != "unknown" && object@spectra != "kernel")
    { stop("spectra must be either \"diag\", \"unknown\", \"kernel\".")}
    if (object@time != "diag" && object@time != "unknown" && object@time != "kernel")
    { stop("time must be either \"diag\", \"unknown\", \"kernel\".")}
    if (object@kerneltypeSpectra != "epanechnikov" && object@kerneltypeSpectra !="gaussian"
        && object@kerneltypeSpectra != "exponential" && object@kerneltypeSpectra !="uniform"
        &&object@kerneltypeSpectra != "quadratic" && object@kerneltypeSpectra != "circular"
        &&object@kerneltypeSpectra !="triangular" && object@kerneltypeSpectra !="rational quadratic"
        &&object@kerneltypeSpectra !="inverse multiquadratic")
    { stop("wrong kerneltypeSpectra entered. ")}
    if (object@kerneltypeTime != "epanechnikov" && object@kerneltypeTime !="gaussian"
        && object@kerneltypeTime != "exponential" && object@kerneltypeTime !="uniform"
        &&object@kerneltypeTime != "quadratic" && object@kerneltypeTime != "circular"
        &&object@kerneltypeTime !="triangular" && object@kerneltypeTime !="rational quadratic"
        &&object@kerneltypeTime !="inverse multiquadratic")
    { stop("wrong kerneltypeTime entered. ")}
    if (round(object@h) != object@h)
    { stop("h must be an integer.")}
    if (object@lambdaS < 0 | object@lambdaT < 0 )
    { stop("lambdaS and lambdaT must be positiv.")}
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
  signature = "fitSpectra",
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
    if (Object@spectra == "diag")
    {
      covSpectra = parsimoniousSpectra(Object@m)
    }
    if (Object@spectra == "unknown")
    {
      if(Object@time=="unknown")
      {
        if(Object@validation)
        {
          lambda = bestFitLambda(Object,listS = seq(from=0.1,to=10,by=0.1),listT = seq(from=0.1,to=10,by=0.1))
          Object@lambdaS = lambda[[1]]
          Object@lambdaT = lambda[[2]]
        }
        ff = flipflop(data = Object@m, lmean = mean, s=Object@s, lambdaS = Object@lambdaS,lambdaT = Object@lambdaT)
        covSpectra = listof(ff,1)
        covTime = listof(ff,2)
      }
      else
      covSpectra = fullSpectra(Object@m)
    }
    if (Object@spectra == "kernel")
    {
      covSpectra = kernelSpectra(Object@m, Object@kerneltypeSpectra, Object@h)
    }

    if (Object@time == "diag")
    {
      covTime = parsimoniousTime(Object@m)
    }
    if (Object@time == "unknown")
    {
      if(Object@spectra!="unknown")
        covTime = fullTime(Object@m)
    }
    if (Object@time == "kernel")
    {
      covTime = kernelTime(Object@m, Object@kerneltypeTime, Object@h)
    }
  }

  if (Object@modelname=="full")
  {
    Object@covMat = list (covS = covmat, covT = 1, modelname = "full", spectra = Object@spectra
                          , time = Object@time, weight=weight, mean=mean)
  }
  else
  {
    Object@covMat = list (covS = covSpectra, covT = covTime, modelname = Object@modelname
                          , spectra = Object@spectra, time = Object@time, weight=weight, mean=mean)
  }

  Object@covMat
  }
)

#-----------------------------------------------------------------------
#' Initialize an instance of a fitSpectra S4 class.
#'
#' Initialization method of the fitSpectra class.
#'
#' @param .Object object of class fitSpectra
#' @param m spectroscopic data
#' @param modelname name of model to be used for calculating the covariance matrix. Available models are
#' "full", "parsimonious". Default is "full".
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
#' @param s regularisation paramater for flip flop algorithm
#' @param lambdaS regularisation for spectra for flip flop algorithm
#' @param lambdaT regularisation for spectra for flip flop algorithm
#' @param validation to optimize lambda or not
#' @param model use in prediction in case of validation is TRUE
#' @param covMat returning the covariance matrix
#'
#' @name initialize
#' @rdname initialize-method
#' @keywords internal
#'

setMethod(
  "initialize",
  "fitSpectra",
  function(.Object, m = list(0), modelname = "full", spectra = "diag", time = "diag"
           , kerneltypeSpectra = "exponential", kerneltypeTime    = "exponential"
           , h = 10, s = 0.01, lambdaS = 0.3, lambdaT = 0.3,validation = FALSE
           , model = "gaussian")
  { .Object@m = m
  .Object@modelname = modelname
  .Object@spectra = spectra
  .Object@time = time
  .Object@kerneltypeSpectra = kerneltypeSpectra
  .Object@kerneltypeTime = kerneltypeTime
  .Object@h = h
  .Object@s = s
  .Object@lambdaS = lambdaS
  .Object@lambdaT = lambdaT
  .Object@validation = validation
  .Object@model = model
  return(.Object)
  }
)

#' Wrapper function fitSpectra.
#'
#' @param ... any paramaters to be input into the function
#'
#' @name fitSpectra
#' @rdname fitSpectra-class
#'
#' @export
fitSpectra <- function(...)
{
  o = new("fitSpectra", ...)
  fit(o)
}
