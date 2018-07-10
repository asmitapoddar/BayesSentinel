#-----------------------------------------------------------------------
#' Create a list with covariance matrices of the spectra and times
#'
#' Return the covariance matrices
#'
#' @param m spectroscopic data
#' @param modelname name of model to be used for calculating the covariance matrix
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
#'

fit <- function(m, modelname = "full", spectra = "diag", time = "diag"
                , kerneltypeSpectra = "exponential",kerneltypeTime = "exponential", h = 10)
{
  source('~/bayes/R/mean.R')
  source('~/bayes/R/full.R')
  source('~/bayes/R/kernelTime.R')
  source('~/bayes/R/parsimonious.R')

  weight = lapply(levels(factor(m[[1]])), function(k,data){length(data[which(data==k)])/length(data)}, data = m[[1]])

  mean = mean(m)

  if (modelname=="full")
    covmat = full(m)
  if (modelname=="parsimonious")
  {
    if (spectra == "diag")
    {
      covSpectra = parsimoniousSpectra(m)
    }
    if (spectra == "unknown")
    {
      covSpectra = fullSpectra(m)
    }
    if (spectra == "kernel")
    {
      covSpectra = kernelSpectra(m, kerneltypeSpectra, h)
    }

    if (time == "diag")
    {
      covTime = parsimoniousTime(m)
    }
    if (time == "unknown")
    {
      covTime = fullTime(m)
    }
    if (time == "kernel")
    {
      covTime = kernelTime(m, kerneltypeTime, h)
    }
  }

  if (modelname=="full")
    covMat = list (covS = covmat, covT = 1, modelname = "full", spectra = spectra, time = time, weight=weight,mean=mean)
  else
    covMat = list (covS = covSpectra, covT = covTime, modelname = modelname, spectra = spectra, time = time, weight=weight,mean=mean)

  covMat
}
