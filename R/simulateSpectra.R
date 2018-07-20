#-----------------------------------------------------------------------
#' Create a list with a simulated data set of spectra
#'
#' Simulate one or more Gaussian spectra at regularly sampling time
#'
#' @slot nbPixel number of pixels belonging to class k
#' @slot nbCluster number of cluster
#' @slot nbSpectrum number of spectra
#' @slot simulationType type of simulation. Available options are "gaussian" and
#' "tstudent". Default is "gaussian".
#' @slot modelname type of model to be used to build covariance matrix.
#' Available options are "full" and "parsimonious". Default is "full".
#' @slot kernelSpectra type of kernel to be used to simulate  spectra. Available options
#' are "diag", ""epanechnikov", "gaussian", "exponential", "uniform", "quadratic"
#' , "circular", "triangular", "rational quadratic", "inverse multiquadratic".
#' Default is "gaussian".
#' @slot kernelTime type of kernel to be used for simulating time. Available options are
#' "diag", ""epanechnikov", "gaussian", "exponential", "uniform", "quadratic",
#' "circular", "triangular", "rational quadratic", "inverse multiquadratic".
#' Default is "gaussian".
#' @slot nbSampling number of sampling
#' @slot sigma a vector of size nbSpectrum giving the variance level of
#' the spectrum
#' @slot nbSampling number of time intervals of the simulation
#' @slot times time intervals of the simulation
#' @slot width the width of the kernel to use for "gaussian" simulation. Default is 50.
#' @slot gamma degrees of freedom used for simulating "tstudent" distribution of data.
#' Default is 3.
#' @slot labels class labels of the data
#' @slot result return a list of simulated data
#'
#' @examples
#' m = new("simulateSpectra")
#' res = simulate(m)
#'
#' @author Serge Iovleff, Asmita Poddar & Florent Latimier
#'
#' @name simulateSpectra
#' @aliases simulateSpectra-class
#' @rdname simulateSpectra-class
#' @exportClass simulateSpectra
#'

setClass(
  Class="simulateSpectra",
  representation( nbPixel         = "numeric"
                  , nbCluster     = "numeric"
                  , nbSpectrum    = "numeric"
                  , simulationType = "character"
                  , modelname     = "character"
                  , kernelSpectra = "character"
                  , kernelTime    = "character"
                  , nbSampling    = "numeric"
                  , sigma         = "numeric"
                  , times         = "numeric"
                  , width         = "numeric"
                  , gamma         = "numeric"
                  , labels        = "numeric"
                  , result        = "list"
  ),
  prototype( nbPixel        = 1000
             , nbCluster    = 15
             , nbSpectrum   = 10
             , simulationType = "gaussian"
             , modelname     = "full"
             , kernelSpectra = "gaussian"
             , kernelTime    = "gaussian"
             , nbSampling   = 33
             , sigma        = integer(0)
             , width        = 50
             , gamma        = 3
             , times = c(0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170
                         ,180,190,200,210,220,230,240,250,260,270,280,290,300,310,321)
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
    #if (object@kernelName != "gaussian" && object@kernelName != "tstudent"
    # && object@kernelName != "tskewed")
    #{ stop("kernelName must be either \"gaussian\", \"tstudent\", \"tskewed\".")}
    if (round(object@nbSampling) != object@nbSampling)
    { stop("nbSampling must be an integer.")}
    if (round(object@width) != object@width)
    { stop("width must be an integer.")}
    return(TRUE)
  }
)

#' Method num.
#' @name simulate
#' @rdname simulate-methods
#' @exportMethod simulate

setGeneric("simulate",
           def=function(Object)
           {
             standardGeneric("simulate")
           }
)


#' @rdname simulate-methods
#' @aliases simulate

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

  KernelCov <- function(times, spectra, labels, modelname, kernelSpectra, kernelTime
                        , nbCluster, nbSpectrum, nbSampling, h)
  {
    #source('~/bayesS4/R/simulateKernel.R')
    sigmaS = rexp(nbSpectrum)
    sigmaT = rexp(nbSampling)
    sigmaL = rexp(nbCluster)

    simulateKernel(modelname, kernelSpectra, kernelTime, times, spectra, labels
                   , sigmaL, sigmaS, sigmaT, h)
  }

  Object@times = c(0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180
                   ,190,200,210,220,230,240,250,260,270,280,290,300,310,321)
  means <- mean(Object@times, Object@nbSpectrum, Object@nbCluster)

  #creating a vector of size nbPixel containing the labels (number of labels = nbCluster)
  #the probablilty of each cluster being between 0 and 1
  labels <- sample(1:Object@nbCluster, Object@nbPixel
                   , prob = rexp(Object@nbCluster) , replace = T)
  spectra = 1:Object@nbSpectrum
  ##prob = rep(1, nbCluster)

  covariance <- KernelCov(Object@times, spectra, labels, Object@modelname
                          , Object@kernelSpectra, Object@kernelTime, Object@nbCluster
                          , Object@nbSpectrum, Object@nbSampling, Object@width)
  covariance <- lapply(covariance, function(mat){(mat %*% t(mat)) /2}) #to check symmetry

  if (Object@simulationType == "gaussian")
  {
    labels  = sort(labels)
    nb = table(labels)
    process <- lapply(1:Object@nbCluster, function(nb,mean,covariance,label)
    {rmvnorm(nb[label], mean = as.numeric(t(means[label,,]))
             , sigma = covariance[[label]])}
      , nb = nb, mean = means, covariance = covariance)

    data <- do.call("rbind",process)
    process <- lapply(1:Object@nbSpectrum,function(data,spectra,nbSampling)
      {data[,((spectra-1)*nbSampling+1):(spectra*nbSampling)]}
      , data = data, nbSampling = Object@nbSampling)
    names(process) <- paste("spectra", 1:length(process), sep="")
  }
  if (Object@simulationType == "tstudent")
  {
    process <- array(0, dim = c(Object@nbPixel, Object@nbSpectrum, Object@nbSampling))
    d <- matrix(0, nrow = Object@nbSpectrum, ncol = Object@nbSampling)
    for (i in 1:Object@nbPixel)
    {
      k <- labels[i]
      for ( s in 1:Object@nbSpectrum)
      {
        d[s,] = rt(Object@nbSampling, Object@gamma, means[k,s,] )
      }
      process[i,,] <- d
    }
  }


  Object@result = list(labels=labels , times = Object@times, spectra = process
       , clouds = list(years1 = matrix(0, nrow = Object@nbPixel
                                       , ncol = length(Object@times) ))
       , means = means, sigma = sigma
  )
  return(Object@result)
}

)


#' Initialize an instance of a simulateSpectra S4 class.
#'
#' Initialization method of the simulateSpectra class.
#'
#' @slot .Object object of class simulateSpectra
#' @slot nbPixel number of pixels belonging to class k
#' @slot nbCluster number of cluster
#' @slot nbSpectrum number of spectra
#' @slot simulationType type of simulation. Available options are "gaussian" and
#' "tstudent". Default is "gaussian".
#' @slot modelname type of model to be used to build covariance matrix.
#' Available options are "full" and "parsimonious". Default is "full".
#' @slot kernelSpectra type of kernel to be used to simulate  spectra. Available options
#' are "diag", ""epanechnikov", "gaussian", "exponential", "uniform", "quadratic"
#' , "circular", "triangular", "rational quadratic", "inverse multiquadratic".
#' Default is "gaussian".
#' @slot kernelTime type of kernel to be used for simulating time. Available options are
#' "diag", ""epanechnikov", "gaussian", "exponential", "uniform", "quadratic",
#' "circular", "triangular", "rational quadratic", "inverse multiquadratic".
#' Default is "gaussian".
#' @slot nbSampling number of sampling
#' @slot sigma a vector of size nbSpectrum giving the variance level of
#' the spectrum
#' @slot nbSampling number of time intervals of the simulation
#' @slot times time intervals of the simulation
#' @slot width the width of the kernel to use for "gaussian" simulation. Default is 50.
#' @slot gamma degrees of freedom used for simulating "tstudent" distribution of data.
#' Default is 3.
#' @slot labels class labels of the data
#' @slot result return a list of simulated data
#'
#' @examples
#' m = new("simulateSpectra")
#' res = simulate(m)
#'
#' @name initialize
#' @rdname initialize-methods
#' @keywords internal
#'

setMethod(
  "initialize",
  "simulateSpectra",

  function(.Object, nbPixel = 1000, nbCluster = 15, nbSpectrum = 10
           , nbSampling = 33, sigma = rexp(nbSpectrum)
           , times = c(0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170
                       ,180,190,200,210,220,230,240,250,260,270,280,290,300,310,321)
           , width = 50, simulationType = "gaussian", modelname     = "full"
           , kernelSpectra = "gaussian", kernelTime = "gaussian")
  { .Object@nbPixel = nbPixel
  .Object@nbCluster = nbCluster
  .Object@nbSpectrum = nbSpectrum
  .Object@simulationType = simulationType
  .Object@nbSampling = nbSampling
  .Object@sigma = sigma
  .Object@times = times
  .Object@width = width
  .Object@modelname = modelname
  .Object@kernelSpectra = kernelSpectra
  .Object@kernelTime = kernelTime
  return(.Object)
  }

)

#' Wrapper function simulaSpectra.
#'
#' @name simulateSpectra
#' @rdname simulateSpectra-class
#' @export
simulateSpectra <- function(...)
  {
    o = new("simulateSpectra", ...)
    simulate(o)
}

