#-----------------------------------------------------------------------
#' Perform regularisation on covariance matrix
#'
#' #-----------------------------------------------------------------------
#' Create a list with a simulated data set of spectra
#'
#' Simulate one or more Gaussian spectra at regularly sampling time
#'
#' @param x spectroscopic data
#' @param modelname name of model to be used for calculating the covariance matrix
#' @param dist type of distribution
#'
#' @examples
#' predModel = predict(m, fittedCov)
#'
#'
#' @return A list with the spectra
#' @author Asmita Poddar & Florent Latimier
#'
#'
#'@example
#'p = predict(m, cov, 50)
#'
#'
#' @param x spectroscopic data
#' @param modelname name of model to be used for calculating the covariance matrix
#' @param dist type of distribution
#'
#' @examples
#' predModel = predict(m, fittedCov)
#'
#'
#' @return A list with the spectra
#' @author Asmita Poddar & Florent Latimier
#'
#'
#'@example
#'p = predict(m, cov, 50)
#'


regularisation = function(fittedCov, lambda)
{
  returnMatrices = function(mat, lambda)
  {
    mat + diag(rep(lambda,nrow(mat)))
  }

  if (fittedCov$modelname == "full")
  {
    reg = lapply(fittedCov$covS, returnMatrices, lambda = lambda )
    #if(model!="mvnorm")
    reg = lapply(reg,inversion)
  }
  else if (fittedCov$modelname == "parsimonious")
  {
    if (fittedCov$spectra == "diag")
    {
      S = lapply(fittedCov$covS,inversion)
    }
    if (fittedCov$spectra == "unknown")
    {
      S = lapply(fittedCov$covS, returnMatrices, lambda = lambda )
      S = lapply(S,inversion)
    }
    if (fittedCov$spectra == "kernel")
    {
      S = lapply(fittedCov$covS,inversion)
    }

    if (fittedCov$time == "diag")
    {
      T = lapply(fittedCov$covT,inversion)
    }
    if (fittedCov$time == "unknown")
    {
      T = lapply(fittedCov$covT, returnMatrices, lambda = lambda )
      T = lapply(T,inversion)
    }
    if (fittedCov$time == "kernel")
    {
      T = lapply(fittedCov$covT,inversion)
    }

    reg = Map('%x%', S , T)
  }

  reg
}


