#-----------------------------------------------------------------------
#' Perform regularisation on covariance matrix
#'
#' Perform regularisation on covariance matrix based on the model selected
#'
#' @param fittedCov covariance matrix
#' @param lambda parameter for regularisation
#'
#' @examples
#' reg = regularisation(fittedCov, 0.3)
#'
#'
#' @return Regularised covariance matrix
#' @author Asmita Poddar & Florent Latimier
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


