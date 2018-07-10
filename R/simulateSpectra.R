#-----------------------------------------------------------------------
#' Create a list with a simulated data set of spectra
#'
#' Simulate one or more Gaussian spectra at regularly sampling time
#'
#' @param nbPixel number of pixels belonging to class k
#' @param nbCluster number of cluster
#' @param nbSpectrum number of spectra
#' @param nbSampling number of sampling
#' @param sigma a vector of size nbSpectrum giving the variance level of
#' the spectrum
#' @param kernelName [\code{string}] with the kernel to use for the covariance matrix.
#' Available kernels are "gaussian", "tstudent". Default is "gaussian".
#' @param width the width of the kernel to use for Gaussian simulation. Default is 50.
#'              It also signifies the degree of freedom for Student-T simulation.
#'
#' @examples
#' m = simulateSpectra(1000,15,10)
#'
#' @return A list with the labels, times, spectra, clouds, means and sigma
#' @author Serge Iovleff & Asmita Poddar
#'
simulateSpectra<-function(nbPixel, nbCluster, nbSpectrum, kernelName = "gaussian"
                          , nbSampling = 33, sigma = rexp(nbSpectrum) , width = 50
                        )
{
  mean=function(t, nbSpectrum, nbCluster)
  {
    res <- array(0, c(nbCluster, nbSpectrum, length(t)));
    a0 = 100
    b0 = 200

    ak = rexp(length(t))
    lk = rexp(length(t))

    for(i in 1:nbCluster)
    {
      for(j in 1:nbSpectrum)
      {
        s = rep(0, length(t))
        meanLevel = a0*j+b0*i
        #s = meanLevel + colSums(ak*cos((2*pi*lk*t)/365))
        s = meanLevel + ak*cos((2*pi*lk*t)/365)
        res[i,j,]=s
      }
    }
    res
  }
  KernelCov <- function(times, width, sigma, kernelName)
  {
    tLength <- length(times)
    res <- matrix(0, nrow=tLength, ncol=tLength)
    for(i in 1:tLength)
    {
      for(j in 1:tLength)
      {
        s <- times[i]
        t <- times[j]
        if (kernelName == "gaussian")
          res[i,j] <- exp(-abs(s-t)^2/width)
        else if (kernelName == "tstudent")
          res[i,j] <- 1/(1+(abs(s-t))^width)
        else
          stop("Enter a valid kernel name for simulating the data")
      }
    }
    sigma * res
  }
  #times <- seq(from=0, to=365, length.out=nbSampling)
  times = c(0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210
            ,220,230,240,250,260,270,280,290,300,310,321)
  means <- mean(times, nbSpectrum, nbCluster)
  data <- array(0, dim = c(nbPixel, nbSpectrum, nbSampling))
  process <- matrix(0, nrow=nbSpectrum, ncol= nbSampling)

  #creating a vector of size nbPixel containing the labels (number of labels = nbCluster)
  #the probablilty of each cluster being between 0 and 1
  labels <- sample(1:nbCluster, nbPixel , prob = rexp(nbCluster) , replace = T)
  ##prob = rep(1, nbCluster)

  for (i in 1:nbPixel)
    {
      k <- labels[i]
      for ( s in 1:nbSpectrum)
      {
        covariance <- KernelCov(times, width, sigma[s], kernelName)
        if (kernelName == "gaussian")
          process[s,] <- rmvnorm(1, mean = means[k,s,], sigma = covariance )
        if (kernelName == "tstudent")
          process = rt(nbSampling, width, means[k,s,] )
      }
      data[i,,] <- process
    }

  samples <- 1:nbSampling

  if (nbSpectrum == 4)
  {
    spectra = list( spect1  = data[,1,]
                    , spect2  = data[,2,]
                    , spect3  = data[,3,]
                    , spect4  = data[,4,]
                  )
  }

  else if (nbSpectrum == 10)
  {
    spectra = list( spect1  = data[,1,]
                    , spect2  = data[,2,]
                    , spect3  = data[,3,]
                    , spect4  = data[,4,]
                    , spect5  = data[,5,]
                    , spect6  = data[,6,]
                    , spect7  = data[,7,]
                    , spect8  = data[,8,]
                    , spect9  = data[,9,]
                    , spect10 = data[,10,]
                  )
  }

  else
    stop("Number of Spectra must be 4 or 10")


  list(labels=labels , times = times, spectra = spectra
       , clouds = list(years1 = matrix(0, nrow = nbPixel, ncol = length(times) ))
       , means = means, sigma = sigma, process = process
      )
}



