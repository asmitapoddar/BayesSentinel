#' kernelCol
#'
#' Calculate the kernel matrix according to the column
#'
#' @param data the 3 dimentional data
#' @param kerneltype a string to represente the type of the kernel.
#' Available kernels are "epanechnikov", "gaussian", "exponential", "uniform",
#' "quadratic", "circular", "triangular", "rational quadratic", "inverse multiquadratic".
#' Default is "exponential".
#' @param h the wide to optimize the kernel
#'
#' @return the kernel square matrix
#'
#' @name kernelCol
#' @export kernelCol
#'
kernelCol = function(data,kerneltype,h)
{
  #do the distancie matrix in Column
  QmatrixCol <- function(data)
  {
    n = length(data[[2]])
    matDiffCol = matrix (0,nrow = n,ncol = n)
    for(i in 1:n)
    {
      matDiffCol[,i] <- abs(data[[2]][i]-data[[2]])
    }
    matDiffCol
  }

  if (kerneltype=="epanechnikov") K = K_E
  if (kerneltype=="gaussian") K = K_G
  if (kerneltype=="exponential") K = K_Exp
  if (kerneltype=="uniform") K = K_U
  if (kerneltype=="quadratic") K = K_Q
  if (kerneltype=="circular") K = K_C
  if (kerneltype=="triangular") K = K_T
  if (kerneltype=="rational quadratic") K = K_RQ
  if (kerneltype=="unverse multiquadratic") K = K_IMQ

  q = K(QmatrixCol(data),h)
  v = diagCol(data)
  for (j in 1:length(unique(data[[1]])) )
  {
    v[[j]] <- sqrt(v[[j]])%*%q%*%sqrt(v[[j]])
  }
  v
}


#' kernelRow
#'
#' Calculate the kernel matrix according to the column
#'
#' @param data the 3 dimentional data
#' @param kerneltype a string to represente the type of the kernel.
#' Available kernels are "epanechnikov", "gaussian", "exponential", "uniform",
#' "quadratic", "circular", "triangular", "rational quadratic", "inverse multiquadratic".
#' Default is "exponential".
#' @param h the wide to optimize the kernel
#'
#' @return the kernel square matrix
#'
#' @name kernelRow
#' @export kernelRow
#'
kernelRow = function(data,kerneltype,h)
{
  #do the distancie matrix in rows
  QmatrixRow <- function(data)
  {
    n = length(data[[6]])
    # change l as the row vector
    l = 1:n
    matDiffRow = matrix (0,nrow = n,ncol = n)
    for(i in 1:n)
    {
      matDiffRow[,i] <- abs(l[i]-l)
    }
    matDiffRow
  }

  if (kerneltype=="epanechnikov") K = K_E
  if (kerneltype=="gaussian") K = K_G
  if (kerneltype=="exponential") K = K_Exp
  if (kerneltype=="uniform") K = K_U
  if (kerneltype=="quadratic") K = K_Q
  if (kerneltype=="circular") K = K_C
  if (kerneltype=="triangular") K = K_T
  if (kerneltype=="rational quadratic") K = K_RQ
  if (kerneltype=="inverse multiquadratic") K = K_IMQ

  q = K(QmatrixRow(data),h)
  v = diagRow(data)
  for (j in 1:length(unique(data[[1]])) )
  {
    v[[j]] <- sqrt(v[[j]])%*%q%*%sqrt(v[[j]])
  }
  v
}

