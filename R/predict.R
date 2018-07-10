#-----------------------------------------------------------------------
#' Create a list with a simulated data set of spectra
#'
#' Simulate one or more Gaussian spectra at regularly sampling time
#'
#' @param x spectroscopic data
#' @param modelname name of model to be used for calculating the covariance matrix
#' @param dist type of distribution
#'
#' @examples
#' predModel = predict(m, fittedCov)
#'
#'
#' @return A list with the spectra
#' @author Asmita Poddar & Florent Latimier
#'
#'
#'@example
#'p = predict(m, cov, 50)
#'

predict = function(m, fittedCov, lambda = 0.5, model = "gaussian")

{

  #mvnorm = function(data,reg,mean,weight)
  #{
  #  X = do.call('cbind',m[[6]])
  #  log(weight) + dmvnorm(X,as.vector(mean),reg,log=TRUE)
  #}

  source('~/bayes/R/regularisation.R')

  nbLabel = length(unique(m[[1]]))
  nbPixel = length(m[[1]])

  p = matrix(0, nbLabel, nbPixel)
  weight  = fittedCov$weight
  reg = regularisation(fittedCov, lambda)

  powerLabelG <- function(data,inv,mean,weight){
    X = do.call('cbind',m[[6]]) - t(matrix(rep(as.vector(mean),nbPixel),nrow=length(as.vector(mean))))
    power = rowSums((X %*% inv) * X)
    log(weight * sqrt(abs(det(inv)))) + (-power/2)
  }

  powerLabelF <-function(data,inv,mean,weight){
    X = do.call('cbind',m[[6]]) - t(matrix(rep(as.vector(mean),nbPixel),nrow=length(as.vector(mean))))
    power = rowSums((X %*% inv) * X)
    log(weight) - 1/2*log(abs(det(inv)))-165*log(1+power)

  }


  power <- function(data,invers,mean,weight){
    nbLabel = length(unique(m[[1]]))
    nbPixel = length(m[[1]])
    p = matrix(0, nbLabel, nbPixel)
    if(model == "gaussian"){powerL = powerLabelG}
    if(model == "fisher"){powerL = powerLabelF}
    #if(model == "mvnorm"){powerL = mvnorm}
    for(i in 1:nbLabel){
      p[i,] = powerL(data,invers[[i]],mean[i,,],weight[[i]])
    }
    p
  }

  max.col(t(power(m,reg,m[[3]],weight)))
  # finalp = as.numeric(lapply(1:ncol(p), function(x) {which.max( p[,x] )} ))
}
