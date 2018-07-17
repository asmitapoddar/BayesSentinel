source('~/bayesS4/R/kernelTypes.R')

QmatrixTime <- function(data)
{
  n = length(data[[2]])
  matDiffTemps = matrix (0,nrow = n,ncol = n)
  for(i in 1:n)
  {
    matDiffTemps[,i] <- abs(data[[2]][i]-data[[2]])
  }
matDiffTemps
}

kernelTime = function(data,k,h)
{
  if (k=="epanechnikov") K = K_E
  if (k=="gaussian") K = K_G
  if (k=="exponential") K = K_Exp
  if (k=="uniform") K = K_U
  if (k=="quadratic") K = K_Q
  if (k=="circular") K = K_C
  if (k=="triangular") K = K_T
  if (k=="rational quadratic") K = K_RQ
  if (k=="unverse multiquadratic") K = K_IMQ

  q = K(QmatrixTime(data),h)
  v = parsimoniousTime(data)
  for (j in 1:length(unique(data[[1]])) )
  {
   v[[j]] <- sqrt(v[[j]])%*%q%*%sqrt(v[[j]])
  }
  v
}

QmatrixSpectra <- function(data)
{
  n = length(data[[6]])
  # change l as the spectra vector
  l = 1:n
  matDiffTemps = matrix (0,nrow = n,ncol = n)
  for(i in 1:n)
  {
    matDiffTemps[,i] <- abs(l[i]-l)
  }
  matDiffTemps
}

kernelSpectra = function(data,k,h)
{
  if (k=="epanechnikov") K = K_E
  if (k=="gaussian") K = K_G
  if (k=="exponential") K = K_Exp
  if (k=="uniform") K = K_U
  if (k=="quadratic") K = K_Q
  if (k=="circular") K = K_C
  if (k=="triangular") K = K_T
  if (k=="rational quadratic") K = K_RQ
  if (k=="inverse multiquadratic") K = K_IMQ

  q = K(QmatrixSpectra(data),h)
  v = parsimoniousSpectra(data)
  for (j in 1:length(unique(data[[1]])) )
  {
    v[[j]] <- sqrt(v[[j]])%*%q%*%sqrt(v[[j]])
  }()
  v
}



