#-----------------------------------------------------------------------
#' Predict the label classes of the data
#'
#' Simulate one or more Gaussian spectra at regularly sampling time
#'
#' @param m spectroscopic data
#' @param fittedCov fitted covariance matrix for the data
#' @param lambda parameter for regularisation
#' @param model type of model to be used for prediction of labels
#' Available models are "gaussian", "fisher". Default is "gaussian".
#'
#' @examples
#' p = predict(m, fittedCov)
#'
#' @return A list with the spectra
#'
#' @author Asmita Poddar & Florent Latimier
#'

predict = function(m, fittedCov, lambda = 0.3, model = "gaussian")
{
  mvnorm = function(data,reg,mean,weight)
  {
    X = do.call('cbind',m[[3]])
    weight * dmvnorm(X,mean,reg,log=TRUE)
  }

  source('~/bayes/R/regularisation.R')

  nbLabel = length(unique(m[[1]]))
  nbPixel = length(m[[1]])

  p = matrix(0, nbLabel, nbPixel)
  weight  = fittedCov$weight
  mean = fittedCov$mean
  reg = regularisation(fittedCov, lambda)

  powerLabelG <- function(data,inv,mean,weight){
    X = do.call('cbind',m[[3]]) - mean
    power = rowSums((X %*% inv) * X)
    log(weight) +  log(sqrt(abs(det(inv)))) + (-power/2)
  }

  powerLabelF <-function(data,inv,mean,weight){
    X = do.call('cbind',m[[3]]) - mean
    power = rowSums((X %*% inv) * X)
    log(weight) - 1/2*log(abs(det(inv))) -165*log(1+power)

  }


  power <- function(data,invers,mean,weight){
    nbLabel = length(unique(m[[1]]))
    nbPixel = length(m[[1]])
    p = matrix(0, nbLabel, nbPixel)
    if(model == "gaussian"){powerL = powerLabelG}
    if(model == "fisher"){powerL = powerLabelF}
    if(model == "mvnorm"){powerL = mvnorm}
    for(i in 1:nbLabel){
      p[i,] = powerL(data,invers[[i]],mean[[i]],weight[[i]])
    }
    p
  }

  max.col(t(power(m,reg,mean,weight)))
}
