simulateKernel = function(modelname, kernelSpectra, kernelTime, times, spectra
                          , labels, sigmaL, sigmaS, sigmaT, h)
{
  if(modelname == "full")
  {
    Q = matrix(0, nrow = (length(times)*length(spectra))
                , ncol = (length(times)*length(spectra)) )
    l <- vector("double",length = (length(times)*length(spectra)))
    sigma <- vector("double",length = (length(times)*length(spectra)))
    for(i in 1:length(spectra))
    {
      l[((i-1)*length(times)+1):(i*length(times))] <- spectra[i] * times
      sigma[((i-1)*length(times)+1):(i*length(times))] <- sigmaS[i] * sigmaT
    }
    for (i in 1:length(times))
    {
      Q[,i] <- abs(l[i]-l)
    }
    sigmal = lapply(sigmaL, function(list, int) {int*list}, list = sigma)
    Q = ker(Q, kernelSpectra, h)
  res = lapply(sigmal, function(mat,vect)
    {diag(sqrt(vect)) %*% mat %*% diag(sqrt(vect)) }, mat = Q)
 }


  else
  {
    QS = matrix(0, nrow = length(spectra), ncol = length(spectra))
    QT = matrix(0, nrow = length(times), ncol = length(times))
    if(kernelSpectra != "diag")
    {
      for(i in 1:length(spectra))
      {
        QS[,i] <- abs(spectra[i]-spectra)
      }
      QS = ker(QS,kernelSpectra,h)
    }
    else
    {
      QS = diag(rep(1,length(spectra)))
    }
    if(kernelTime != "diag")
    {
      for(i in 1:length(times))
      {
        QT[,i] <- abs(times[i]-times)
      }
      QT = ker(QT,kernelTime,h)
    }
    else
    {
      QT = diag(rep(1,length(times)))
    }

    res = lapply(sigmaL, function(matS, vectS, matT, vectT, int)
      {(diag(sqrt(int*vectS)) %*% matS %*% diag(sqrt(int*vectS))) %x%
        (diag(sqrt(int*vectT)) %*% matT %*% diag(sqrt(int*vectT)) )
      }, matS = QS , vectS = sigmaS, matT = QT , vectT = sigmaT)
  }

  res
}





ker <- function(mat, kernelType,h)
{
  if (kernelType=="epanechnikov") K = K_E
  if (kernelType=="gaussian") K = K_G
  if (kernelType=="exponential") K = K_Exp
  if (kernelType=="uniform") K = K_U
  if (kernelType=="quadratic") K = K_Q
  if (kernelType=="circular") K = K_C
  if (kernelType=="triangular") K = K_T
  if (kernelType=="rational quadratic") K = K_RQ
  if (kernelType=="inverse multiquadratic") K = K_IMQ

  K(mat,h)
}
