#-----------------------------------------------------------------------
#' Perform regularisation on covariance matrix
#'
#' @fittedCov fitted covariance matrix
#' @param lambdaS regularisation parameter for spectra
#' @param lambdaT regularisation parameter for time
#'
#' @examples
#' reg = regularisation(fittedCov)
#'
#' @return A list with regularised covariance matrix
#'
#' @author Asmita Poddar & Florent Latimier
#'

regularisation = function(fittedCov, lambdaS = 0.3 , lambdaT = 0.3)
{
  #print(lambdaT)
  #print(lambdaS)
  returnMatrices = function(mat, lambda)
  {
    regul(mat + diag(rep(lambda,nrow(mat))))
  }

  if (fittedCov$modelname == "full")
  {
    reg = lapply(fittedCov$covS, returnMatrices, lambda = lambdaS )
    reg = lapply(reg,inversion)
  }
  else if (fittedCov$modelname == "parsimonious")
  {
    if (fittedCov$spectra == "diag")
    {
      S = lapply(fittedCov$covS, returnMatrices, lambda = 0)
      S = lapply(S,inversion)
    }
    if (fittedCov$spectra == "unknown")
    {
      S = lapply(fittedCov$covS, returnMatrices, lambda = lambdaS)
      S = lapply(S,inversion)
    }
    if (fittedCov$spectra == "kernel")
    {
      S = lapply(fittedCov$covS, returnMatrices, lambda = 0)
      S = lapply(S,inversion)
    }

    if (fittedCov$time == "diag")
    {
      T = lapply(fittedCov$covT, returnMatrices, lambda = 0 )
      T = lapply(T,inversion)
    }
    if (fittedCov$time == "unknown")
    {
      T = lapply(fittedCov$covT, returnMatrices, lambda =  lambdaT)
      T = lapply(T,inversion)
    }
    if (fittedCov$time == "kernel")
    {
      T = lapply(fittedCov$covT, returnMatrices, lambda = 0 )
      T = lapply(T,inversion)
    }

    reg = Map('%x%', S , T)
  }

  reg
}


# due to the approximation, the matrix can be non-definie positive
regul = function(mat){
  s=0
  while(!is.positive.definite(mat)){
    s = s + 0.000000001
    mat = mat + diag(rep(0.000000001,nrow(mat)))
  }
  mat
}
