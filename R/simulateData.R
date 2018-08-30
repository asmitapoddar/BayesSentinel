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
#' @slot a0 the mean distance between two row
#' @slot b0 the mean distance between two cluster
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
#' @importFrom methods new
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
                  , a0            = "numeric"
                  , b0            = "numeric"
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
             , a0           = 7
             , b0           = 225
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

      ak = rexp(length(t))
      lk = rexp(length(t))

      for(i in 1:nbCluster)
      {
        for(j in 1:nbRow)
        {
          s = rep(0, length(Object@column))
          meanLevel = Object@a0*j+Object@b0*i
          #s = meanLevel + colSums(ak*cos((2*pi*lk*t)/365))
          s = meanLevel + ak*cos((2*pi*lk*t)/365)
          res[i,j,]=s
        }
      }
      res
    }

    KernelCov <- function(column, rows, labels, modelname, kernelRow, kernelCol
                          , nbCluster, nbRow, nbCol, h)
    {
      #source('~/bayesS4/R/simulateKernel.R')
      sigmaS = rexp(nbRow)
      sigmaT = rexp(nbCol)
      sigmaL = rexp(nbCluster)

      simulateKernel(modelname, kernelRow, kernelCol, column, rows, labels
                     , sigmaL, sigmaS, sigmaT, h)
    }

    Object@column = c(0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180
                      ,190,200,210,220,230,240,250,260,270,280,290,300,310,321)
    means <- mean(Object@column, Object@nbRow, Object@nbCluster)

    #creating a vector of size nbSample containing the labels (number of labels = nbCluster)
    #the probablilty of each cluster being between 0 and 1
    labels <- sample(1:Object@nbCluster, Object@nbSample
                     , prob = rexp(Object@nbCluster) , replace = T)
    rows = 1:Object@nbRow
    ##prob = rep(1, nbCluster)

    covariance <- KernelCov(Object@column, rows, labels, Object@modelname
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
      process <- lapply(1:Object@nbRow,function(data,rows,nbCol)
      {data[,((rows-1)*nbCol+1):(rows*nbCol)]}
      , data = data, nbCol = Object@nbCol)
      names(process) <- paste("row", 1:length(process), sep="")
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


    Object@result = list(labels=labels , column = Object@column, rows = process
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
#' @param a0 the mean distance between two row
#' @param b0 the mean distance between two cluster
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
           , nbCol = 33, sigma = rexp(nbRow), a0 = 7, b0 = 225
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
  .Object@a0 = a0
  .Object@b0 = b0
  return(.Object)
  }

)

#' Simulate a 3 dimentional matrix
#'
#' Simulate a data that all observation is a matrix : [nbRow : nbCol].
#'
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
#' @param a0 the mean distance between two row
#' @param b0 the mean distance between two cluster
#' @param labels class labels of the data
#'
#' @return simulated data as a list of all observation
#'
#' @name simulateDataMatrix
#' @export simulateDataMatrix
#'
#' @importFrom stats dnorm rnorm rexp rt cov runif sigma
#' @importFrom mvtnorm rmvnorm dmvnorm
#'
simulateDataMatrix <- function(...)
{
  o = new("simulateData", ...)
  simulate(o)
}




#' Simulate a data according to a kernel
#'
#' Simulate the covariance matrix as a diagonal or a kenrel
#'
#' @param modelname type of model to be used to build covariance matrix.
#' Available options are "full" and "parsimonious". Default is "full".
#' @param kernelRow type of kernel to be used to simulate  rows. Available options
#' are "diag", "epanechnikov", "gaussian", "exponential", "uniform", "quadratic"
#' , "circular", "triangular", "rational quadratic", "inverse multiquadratic".
#' Default is "gaussian".
#' @param kernelCol type of kernel to be used to simulate  columns. Available options
#' are "diag", "epanechnikov", "gaussian", "exponential", "uniform", "quadratic"
#' , "circular", "triangular", "rational quadratic", "inverse multiquadratic".
#' Default is "gaussian".
#' @param column columns intervals of the simulation
#' @param rows rows intervals of the simulation
#' @param labels class labels of the data
#' @param sigmaL the variation due to the clusthers
#' @param sigmaR the variation due to the rows
#' @param sigmaC the variation due to the columns
#' @param h the widths parameter to fit the kernel
#'
#' @name simulateKernel
#' @export simulateKernel
#'
simulateKernel = function(modelname, kernelRow, kernelCol, column, rows
                          , labels, sigmaL, sigmaR, sigmaC, h)
{
  if(modelname == "full")
  {
    Q = matrix(0, nrow = (length(column)*length(rows))
               , ncol = (length(column)*length(rows)) )
    l <- vector("double",length = (length(column)*length(rows)))
    sigma <- vector("double",length = (length(column)*length(rows)))
    for(i in 1:length(rows))
    {
      l[((i-1)*length(column)+1):(i*length(column))] <- rows[i] * column
      sigma[((i-1)*length(column)+1):(i*length(column))] <- sigmaR[i] * sigmaC
    }
    for (i in 1:length(column))
    {
      Q[,i] <- abs(l[i]-l)
    }
    sigmal = lapply(sigmaL, function(list, int) {int*list}, list = sigma)
    Q = ker(Q, kernelRow, h)
    res = lapply(sigmal, function(mat,vect)
    {diag(sqrt(vect)) %*% mat %*% diag(sqrt(vect)) }, mat = Q)
  }


  else
  {
    QR = matrix(0, nrow = length(rows), ncol = length(rows))
    QC = matrix(0, nrow = length(column), ncol = length(column))
    if(kernelRow != "diag")
    {
      for(i in 1:length(rows))
      {
        QR[,i] <- abs(rows[i]-rows)
      }
      QR = ker(QR,kernelRow,h)
    }
    else
    {
      QR = diag(rep(1,length(rows)))
    }
    if(kernelCol != "diag")
    {
      for(i in 1:length(column))
      {
        QC[,i] <- abs(column[i]-column)
      }
      QC = ker(QC,kernelCol,h)
    }
    else
    {
      QC = diag(rep(1,length(column)))
    }

    res = lapply(sigmaL, function(matR, vectR, matC, vectC, int)
    {(diag(sqrt(int*vectR)) %*% matR %*% diag(sqrt(int*vectR))) %x%
        (diag(sqrt(int*vectC)) %*% matC %*% diag(sqrt(int*vectC)) )
    }, matR = QR , vectR = sigmaR, matC = QC , vectC = sigmaC)
  }

  res
}




#' ker
#'
#' Apply a kernel type on a matrix of distancies
#'
#' @param mat the matrix of distancies
#' @param kernelType type of kernel to be apply on the matrix
#' @param h the widths parameter to fit the kernel
#'
#' @name ker
#' @export ker
#'
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
