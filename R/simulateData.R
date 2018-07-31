#-----------------------------------------------------------------------
#' Create a list with a simulated data set of matrix : [nbRow x nbCol]
#'
#'
#' @slot nbSample number of sample belonging to class k
#' @slot nbCluster number of cluster
#' @slot nbRow number of rows
#' @slot simulationType type of simulation. Available options are "gaussian" and
#' "tstudent". Default is "gaussian".
#' @slot modelname type of model to be used to build covariance matrix.
#' Available options are "full" and "parsimonious". Default is "full".
#' @slot kernelRow type of kernel to be used to simulate  rows. Available options
#' are "diag", "epanechnikov", "gaussian", "exponential", "uniform", "quadratic"
#' , "circular", "triangular", "rational quadratic", "inverse multiquadratic".
#' Default is "gaussian".
#' @slot kernelCol type of kernel to be used for simulating columns. Available options are
#' "diag", "epanechnikov", "gaussian", "exponential", "uniform", "quadratic",
#' "circular", "triangular", "rational quadratic", "inverse multiquadratic".
#' Default is "gaussian".
#' @slot sigma a vector of size nbRow giving the variance level of
#' the row
#' @slot nbCol number of columns intervals of the simulation
#' @slot column columns intervals of the simulation
#' @slot width the width of the kernel to use for "gaussian" simulation. Default is 50.
#' @slot gamma degrees of freedom used for simulating "tstudent" distribution of data.
#' Default is 3.
#' @slot labels class labels of the data
#' @slot result return a list of simulated data
#'
#' @examples
#' m = new("simulateData")
#' res = simulate(m)
#'
#' @author Serge Iovleff, Asmita Poddar & Florent Latimier
#'
#' @name simulateData
#' @aliases simulateData-class
#' @rdname simulateData-class
#' @exportClass simulateData
#'

setClass(
  Class="simulateData",
  representation( nbSample        = "numeric"
                  , nbCluster     = "numeric"
                  , nbRow  	  = "numeric"
                  , simulationType = "character"
                  , modelname     = "character"
                  , kernelRow	 = "character"
                  , kernelCol    = "character"
                  , nbCol   	 = "numeric"
                  , sigma         = "numeric"
                  , column         = "numeric"
                  , width         = "numeric"
                  , gamma         = "numeric"
                  , labels        = "numeric"
                  , result        = "list"
  ),
  prototype( nbSample        = 10000
             , nbCluster    = 15
             , nbRow   = 10
             , simulationType = "gaussian"
             , modelname     = "full"
             , kernelRow = "gaussian"
             , kernelCol    = "gaussian"
             , nbCol   = 33
             , sigma        = integer(0)
             , width        = 50
             , gamma        = 3
             , column = c(0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170
                          ,180,190,200,210,220,230,240,250,260,270,280,290,300,310,321)
             , result = list()
  ),
  # validity function
  validity = function(object)
  {
    # check classNumper
    if (round(object@nbSample) != object@nbSample)
    { stop("nbSample must be an integer.")}
    # check classNumper
    if (round(object@nbCluster) != object@nbCluster)
    { stop("nbCluster must be an integer.")}
    if (round(object@nbRow) != object@nbRow)
    { stop("nbRow must be an integer.")}
    #if (object@kernelName != "gaussian" && object@kernelName != "tstudent"
    # && object@kernelName != "tskewed")
    #{ stop("kernelName must be either \"gaussian\", \"tstudent\", \"tskewed\".")}
    if (round(object@nbCol) != object@nbCol)
    { stop("nbCol must be an integer.")}
    if (round(object@width) != object@width)
    { stop("width must be an integer.")}
    return(TRUE)
  }
)

#' Method num.
#'
#' @name simulate
#' @rdname simulate-method
#' @exportMethod simulate

setGeneric("simulate",
           def=function(Object)
           {
             standardGeneric("simulate")
           }
)


#' #' Method num.
#'
#' @param Object object to be input
#'
#' @rdname simulate-method
#' @aliases simulate

setMethod(
  f = "simulate",
  signature = "simulateData",
  definition=function(Object)
  {
    mean=function(t, nbRow, nbCluster)
    {
      res <- array(0, c(nbCluster, nbRow, length(t)));
      a0 = 100
      b0 = 200

      ak = rexp(length(t))
      lk = rexp(length(t))

      for(i in 1:nbCluster)
      {
        for(j in 1:nbRow)
        {
          s = rep(0, length(Object@column))
          meanLevel = a0*j+b0*i
          #s = meanLevel + colSums(ak*cos((2*pi*lk*t)/365))
          s = meanLevel + ak*cos((2*pi*lk*t)/365)
          res[i,j,]=s
        }
      }
      res
    }

    KernelCov <- function(column, spectra, labels, modelname, kernelRow, kernelCol
                          , nbCluster, nbRow, nbCol, h)
    {
      #source('~/bayesS4/R/simulateKernel.R')
      sigmaS = rexp(nbRow)
      sigmaT = rexp(nbCol)
      sigmaL = rexp(nbCluster)

      simulateKernel(modelname, kernelRow, kernelCol, column, spectra, labels
                     , sigmaL, sigmaS, sigmaT, h)
    }

    Object@column = c(0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180
                      ,190,200,210,220,230,240,250,260,270,280,290,300,310,321)
    means <- mean(Object@column, Object@nbRow, Object@nbCluster)

    #creating a vector of size nbSample containing the labels (number of labels = nbCluster)
    #the probablilty of each cluster being between 0 and 1
    labels <- sample(1:Object@nbCluster, Object@nbSample
                     , prob = rexp(Object@nbCluster) , replace = T)
    spectra = 1:Object@nbRow
    ##prob = rep(1, nbCluster)

    covariance <- KernelCov(Object@column, spectra, labels, Object@modelname
                            , Object@kernelRow, Object@kernelCol, Object@nbCluster
                            , Object@nbRow, Object@nbCol, Object@width)
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
      process <- lapply(1:Object@nbRow,function(data,spectra,nbCol)
      {data[,((spectra-1)*nbCol+1):(spectra*nbCol)]}
      , data = data, nbCol = Object@nbCol)
      names(process) <- paste("spectra", 1:length(process), sep="")
    }
    if (Object@simulationType == "tstudent")
    {
      process <- array(0, dim = c(Object@nbSample, Object@nbRow, Object@nbCol))
      d <- matrix(0, nrow = Object@nbRow, ncol = Object@nbCol)
      for (i in 1:Object@nbSample)
      {
        k <- labels[i]
        for ( s in 1:Object@nbRow)
        {
          d[s,] = rt(Object@nbCol, Object@gamma, means[k,s,] )
        }
        process[i,,] <- d
      }
    }


    Object@result = list(labels=labels , column = Object@column, spectra = process
                         , clouds = list(years1 = matrix(0, nrow = Object@nbSample
                                                         , ncol = length(Object@column) ))
                         , means = means, sigma = sigma
    )
    return(Object@result)
  }

)


#' Initialize an instance of a simulateData S4 class.
#'
#' Initialization method of the simulateData class.
#'
#' @param .Object object of class simulateData
#' @param nbSample number of sample belonging to class k
#' @param nbCluster number of cluster
#' @param nbRow number of rows
#' @param simulationType type of simulation. Available options are "gaussian" and
#' "tstudent". Default is "gaussian".
#' @param modelname type of model to be used to build covariance matrix.
#' Available options are "full" and "parsimonious". Default is "full".
#' @param kernelRow type of kernel to be used to simulate  rows. Available options
#' are "diag", "epanechnikov", "gaussian", "exponential", "uniform", "quadratic"
#' , "circular", "triangular", "rational quadratic", "inverse multiquadratic".
#' Default is "gaussian".
#' @param kernelCol type of kernel to be used for simulating columns. Available options are
#' "diag", "epanechnikov", "gaussian", "exponential", "uniform", "quadratic",
#' "circular", "triangular", "rational quadratic", "inverse multiquadratic".
#' Default is "gaussian".
#' @param sigma a vector of size nbRow giving the variance level of
#' the row
#' @param nbCol number of columns intervals of the simulation
#' @param column columns intervals of the simulation
#' @param width the width of the kernel to use for "gaussian" simulation. Default is 50.
#' @param gamma degrees of freedom used for simulating "tstudent" distribution of data.
#' Default is 3.
#' @param labels class labels of the data
#' @param result return a list of simulated data
#'
#' @examples
#' m = new("simulateData")
#' res = simulate(m)
#'
#' @name initialize
#' @rdname initialize-method
#' @keywords internal
#'

setMethod(
  "initialize",
  "simulateData",

  function(.Object, nbSample = 10000, nbCluster = 15, nbRow = 10
           , nbCol = 33, sigma = rexp(nbRow)
           , column = c(0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170
                        ,180,190,200,210,220,230,240,250,260,270,280,290,300,310,321)
           , width = 50, simulationType = "gaussian", modelname     = "full"
           , kernelRow = "gaussian", kernelCol = "gaussian")
  { .Object@nbSample = nbSample
  .Object@nbCluster = nbCluster
  .Object@nbRow = nbRow
  .Object@simulationType = simulationType
  .Object@nbCol = nbCol
  .Object@sigma = sigma
  .Object@column = column
  .Object@width = width
  .Object@modelname = modelname
  .Object@kernelRow = kernelRow
  .Object@kernelCol = kernelCol
  return(.Object)
  }

)

#' Wrapper function simulateData.
#'
#' @param ... any paramaters to be input into the function
#'
#' @name simulateDataSet
#' @rdname simulateData-class
#' @export
simulateDataSet <- function(...)
{
  o = new("simulateData", ...)
  simulate(o)
}

