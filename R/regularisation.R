#-----------------------------------------------------------------------
#' Regularisation
#'
#' Perform regularisation on covariance matrix
#'
#' @param fittedCov fitted covariance matrix
#' @param lambdaR regularisation parameter for row
#' @param lambdaC regularisation parameter for column
#'
#' @return A list with regularised covariance matrix
#'
#' @author Asmita Poddar & Florent Latimier
#'
#'#' @import matrixcalc is.positive.definite
#'
regularisation = function(fittedCov, lambdaR = 0.3 , lambdaC = 0.3)
{
  # due to the approximation, the matrix can be non-definie positive
  regul = function(mat){
    s=0
    while(!is.positive.definite(mat)){
      s = s + 0.000000001
      mat = mat + diag(rep(0.000000001,nrow(mat)))
    }
    mat
  }

  #do the ridge regularisation
  returnMatrices = function(mat, lambda)
  {
    #to assure the symetry due to approximation
    mat = (mat + t(mat)) /2
    regul((1-lambda)*mat + lambda*mean(diag(mat))*diag(rep(1,nrow(mat))))
  }


  if (fittedCov$modelname == "full")
  {
    reg = lapply(fittedCov$covR, returnMatrices, lambda = lambdaR )
    reg = lapply(reg,inversion)
  }
  else if (fittedCov$modelname == "parsimonious")
  {
    if (fittedCov$row == "diag")
    {
      R = lapply(fittedCov$covR, returnMatrices, lambda = 0)
      R = lapply(R,inversion)
    }
    if (fittedCov$row == "unknown")
    {
      R = lapply(fittedCov$covR, returnMatrices, lambda = lambdaR)
      R = lapply(R,inversion)
    }
    if (fittedCov$row == "kernel")
    {
      R = lapply(fittedCov$covR, returnMatrices, lambda = 0)
      R = lapply(R,inversion)
    }

    if (fittedCov$column == "diag")
    {
      C = lapply(fittedCov$covC, returnMatrices, lambda = 0 )
      C = lapply(C,inversion)
    }
    if (fittedCov$column == "unknown")
    {
      C = lapply(fittedCov$covC, returnMatrices, lambda =  lambdaC)
      C = lapply(C,inversion)
    }
    if (fittedCov$column == "kernel")
    {
      C = lapply(fittedCov$covC, returnMatrices, lambda = 0 )
      C = lapply(C,inversion)
    }

    reg = Map('%x%', R , C)
  }

  reg
}



