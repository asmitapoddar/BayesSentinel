#-----------------------------------------------------------------------
#' Create a list with covariance matrices of the spectra and times
#'
#' Return the covariance matrices
#'
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
#'
#'
#' @examples
#' fittedCov = fit(x, "full")
#'
#' @return A list with the covariance matrices for spectra and time, modelname, spectra,
#' time, weight and mean
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
                  , covMat                 = "list"
  ),
  prototype( m                   = list(0)
             , modelname         = "full"
             , spectra           = "diag"
             , time              = "diag"
             , kerneltypeSpectra = "exponential"
             , kerneltypeTime    = "exponential"
             , h                 = 10
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
    return(TRUE)
  }
)

setGeneric("fit",
           def=function(Object)
           {
             standardGeneric("fit")
           }
)

setMethod(
  f = "fit",
  signature = "fitSpectra",
  definition=function(Object)
  {

  source('~/bayesS4/R/meanData.R')
  source('~/bayesS4/R/full.R')
  source('~/bayesS4/R/kernel.R')
  source('~/bayesS4/R/parsimonious.R')

  #return a list of the list of integer elements
  listof = function(list, int)
  { lapply(list, function(l,int){l[[int]]}, int = int)
  }

  weight = lapply(levels(factor(Object@m[[1]])), function(k,data){length(data[which(data==k)])/length(data)}
                  , data = Object@m[[1]])

  mean = meanData(Object@m)

  if (Object@modelname=="full")
    covmat = full(Object@m)
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
        ff = flipflop(data = Object@m, lmean = mean, s=s)
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
    Object@covMat = list (covS = covmat, covT = 1, modelname = "full", spectra = Object@spectra
                          , time = Object@time, weight=weight, mean=mean)
  else
    Object@covMat = list (covS = covSpectra, covT = covTime, modelname = Object@modelname
                          , spectra = Object@spectra, time = Object@time, weight=weight, mean=mean)
  Object
  }


)
