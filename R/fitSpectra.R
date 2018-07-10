#-----------------------------------------------------------------------
#' Create a list with a simulated data set of spectra
#'
#' Simulate one or more Gaussian spectra at regularly sampling time
#'
#' @param x spectroscopic data
#' @param modelname name of model to be used for calculating the covariance matrix
#' @param dist type of distribution
#'
#' @examples
#' fittedCov = fit(x, "full", "gaussian")
#'
#'
#' @return A list with the covaria
#' @author Asmita Poddar & Florent Latimier
#'
#' @example
#' cov = fit(m)
#'

fit <- function(m, modelname = "full", spectra = "diag", time = "diag", kerneltype = "exponential", h = 10)
{

  source('~/bayes/R/full.R')
  source('~/bayes/R/kernelTime.R')
  source('~/bayes/R/parsimonious.R')

  weight = lapply(levels(factor(m[[1]])), function(k,data){length(data[which(data==k)])/length(data)}, data = m[[1]])

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
      #covSpectra = parsimoniousSpectra(m)
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
      covTime = kernelTime(m, kerneltype, h)
    }
  }

  if (modelname=="full")
    covMat = list (covS = covmat, covT = 1, modelname = "full", spectra = spectra, time = time, weight=weight)
  else
    covMat = list (covS = covSpectra, covT = covTime, modelname = modelname, spectra = spectra, time = time, weight=weight)

  covMat
}
