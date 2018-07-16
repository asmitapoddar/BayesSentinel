#-----------------------------------------------------------------------
#' Create a list with a simulated data set of spectra
#'
#' Simulate one or more Gaussian spectra at regularly sampling time
#'
#' @slot nbPixel number of pixels belonging to class k
#' @slot nbCluster number of cluster
#' @slot nbSpectrum number of spectra
#' @slot nbSampling number of sampling
#' @slot sigma a vector of size nbSpectrum giving the variance level of
#' the spectrum
#' @slot kernelName [\code{string}] with the kernel to use for the covariance matrix.
#'
#' @slot width the width of the kernel to use for Gaussian simulation. Default is 50.
#'              It also signifies the degree of freedom for Student-T simulation.
#'
#' @examples
#' m = simulateSpectra(1000,15,10)
#'
#' @return A list with the labels, times, spectra, clouds, means and sigma
#' @author Serge Iovleff & Asmita Poddar
#'

setClass(
  Class="simulateSpectra",
  representation( nbPixel         = "numeric"
                  , nbCluster     = "numeric"
                  , nbSpectrum    = "numeric"
                  , kernelName    = "character"
                  , nbSampling    = "numeric"
                  , sigma         = "numeric"
                  , times         = "numeric"
                  , width         = "numeric"
                  , result        = "list"
  ),
  prototype( nbPixel        = 1000
             , nbCluster    = 15
             , nbSpectrum   = 10
             , kernelName   = "gaussian"
             , nbSampling   = 33
             , sigma        = integer(0)
             , width        = 50
             , times = c(0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210
                               ,220,230,240,250,260,270,280,290,300,310,321)
             , result = list()
  ),
  # validity function
  validity = function(object)
  {
    # check classNumper
    if (round(object@nbPixel) != object@nbPixel)
    { stop("nbPixel must be an integer.")}
    # check classNumper
    if (round(object@nbCluster) != object@nbCluster)
    { stop("nbCluster must be an integer.")}
    if (round(object@nbSpectrum) != object@nbSpectrum)
    { stop("nbSpectrum must be an integer.")}
    if (object@kernelName != "gaussian" && object@kernelName != "tstudent" && object@kernelName != "tskewed")
    { stop("kernelName must be either \"gaussian\", \"tstudent\", \"tskewed\".")}
    if (round(object@nbSampling) != object@nbSampling)
    { stop("nbSampling must be an integer.")}
    if (round(object@width) != object@width)
    { stop("width must be an integer.")}
    return(TRUE)
  }
)

setGeneric("simulate",
           def=function(Object)
           {
             standardGeneric("simulate")
           }
)

setMethod(
  f = "simulate",
  signature = "simulateSpectra",
  definition=function(Object)
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
           s = rep(0, length(Object@times))
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
        else if (kernelName == "tskewed")
          res[i,j] <- 1/(1+(abs(s-t))^width)
        else
          stop("Enter a valid kernel name for simulating the data")
      }
    }
    sigma * res
  }
  #times <- seq(from=0, to=365, length.out=nbSampling)

  Object@sigma = rexp(Object@nbSpectrum)
  Object@times = c(0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210
            ,220,230,240,250,260,270,280,290,300,310,321)
  means <- mean(Object@times, Object@nbSpectrum, Object@nbCluster)
  data <- array(0, dim = c(Object@nbPixel, Object@nbSpectrum, Object@nbSampling))
  process <- matrix(0, nrow=Object@nbSpectrum, ncol= Object@nbSampling)

  #creating a vector of size nbPixel containing the labels (number of labels = nbCluster)
  #the probablilty of each cluster being between 0 and 1
  labels <- sample(1:Object@nbCluster, Object@nbPixel
                   , prob = rexp(Object@nbCluster) , replace = T)
  ##prob = rep(1, nbCluster)

  gamma = 3
  for (i in 1:Object@nbPixel)
  {
    k <- labels[i]
    for ( s in 1:Object@nbSpectrum)
    {
      covariance <- KernelCov(Object@times, Object@width, Object@sigma[s], Object@kernelName)
      if (Object@kernelName == "gaussian")
        process[s,] <- rmvnorm(1, mean = means[k,s,], sigma = covariance )
      if (Object@kernelName == "tstudent")
        process[s,] = rt(Object@nbSampling, Object@width, means[k,s,] )
      if (Object@kernelName == "tskewed")
        process[s,] <- rskt(Object@nbSampling, Object@width, gamma)
    }
    data[i,,] <- process
  }

  samples <- 1:Object@nbSampling

  if (Object@nbSpectrum == 4)
  {
    spectra = list( spect1  = data[,1,]
                    , spect2  = data[,2,]
                    , spect3  = data[,3,]
                    , spect4  = data[,4,]
    )
  }

  else if (Object@nbSpectrum == 10)
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

  Object@result = list(labels=labels)
  Object@result = list(labels=labels , times = Object@times, spectra = spectra
       , clouds = list(years1 = matrix(0, nrow = Object@nbPixel, ncol = length(Object@times) ))
       , means = means, sigma = sigma, process = process
  )
  return(Object)
}

)

setMethod(
  "initialize",
  "simulateSpectra",
  function(.Object, nbPixel = 1000, nbCluster = 15, nbSpectrum = 10, kernelName = "gaussian"
           , nbSampling = 33, sigma = rexp(nbSpectrum)
           , times = c(0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190
                       ,200,210,220,230,240,250,260,270,280,290,300,310,321), width = 50)
  { .Object@nbPixel = nbPixel
    .Object@nbCluster = nbCluster
    .Object@nbSpectrum = nbSpectrum
    .Object@kernelName = kernelName
    .Object@nbSampling = nbSampling
    .Object@sigma = sigma
    .Object@times = times
    .Object@width = width
    return(.Object)
  }
)

#setGeneric("simulateSpectra",
#           function(nbPixel,nbCluster,nbSpectrum,...)
 #            standardGeneric("simulateSpectra")
  #         )
#setMethod("simulateSpectra",
#          signature(a="missing",b="missing"),
#          function(a,b,...) A(as.numeric(1:10),...)

#)



