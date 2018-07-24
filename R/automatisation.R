#-----------------------------------------------------------------------
#' Find percentage of accuracy
#' Find percentage of accuracy when two lists containig the class labels are provided
#'
#' @param list1 first list
#' @param list2 second list
#'
#' @return percentage value - numeric
#'
#' @author Asmita Poddar & Florent Latimier


percent <- function(list1, list2)
{
  diff = list1 == list2
  sum(diff)/length(diff)*100
}


# separera la data en 2 parties -----------------
#'
#' Divide the data into training and test sets
#' Divide the data into training and test sets according to percentage provided by user
#'
#' @param data first list
#' @param pTrain second list
#'
#' @return list containing training and test set
#'

testTrain <- function(data, pTrain)
{
  n = length(data[[1]])
  sampl = runif(n)
  test <- lapply(data[[3]], function(mat, list) {mat[list,]}, list = which(sampl < pTrain))
  train <- lapply(data[[3]], function(mat, list) {mat[list,]}
                  , list = which(sampl >= pTrain))
  test <- list(data[[1]][which(sampl < pTrain)], data[[2]], test
               , data[[4]][which(sampl < pTrain)])
  train <- list(data[[1]][which(sampl >= pTrain)]
                , data[[2]], train,data[[4]][which(sampl >= pTrain)])
  list(test, train)
}

#### inversion matrice ###
#'
#' Inversion of matrix
#'
#' @param mat matrix to be inverted
#' @return inverted matrix

inversion <- function(mat){
  chol2inv(chol(mat))
}

## tous les test possibles :
#'
#' Inversion of matrix
#'
#' @param data matrix to be inverted
#' @param lkernelS list of kernels for spectra
#' @param lkernelT list of kernels for time
#' @param lhS list of widths for spectra
#' @param lhT list of widths for spectra
#' @param listModel list of models
#'
#' @return list containing the accuaracy for different models
#'
#' @export

tablesModel = function( data, lkernelS = list(), lkernelT = list(), lhS = list()
                       , lhT = list(), listModel = list("gaussian"))
{
  #l = testTrain(data,0.8)
  l = list(data,data)
  nr = 2 + length(lkernelS)
  nc = 2 + length(lkernelT)
  res = list()
  name = c()
  for(k in 1:length(listModel))
  {
    #cov <- new("fitSpectra",l[[1]])
    #cov <- fit(cov)
    cov = fitSpectra(l[[1]])
    pred <- new("predictClass", m = l[[2]], fittedCov = cov, model = listModel[[k]]
                , validation = TRUE)
    pred <- predict(pred)
    full <- pred@accuracy
    p = matrix(ncol = nc, nrow = nr)
    for(i in 1:nr){
      if(i==1){
        spectra = "diag"
      }
      if(i==2){
        spectra = "unknown"
      }
      kerneltypeSpectra = ""
      h=0
      if(i>2){
        spectra = "kernel"
        kerneltypeSpectra = lkernelS[[i-2]]
        h = lhS[[i-2]]
      }

      for(j in 1:nc){
        if(j==1){
          time = "diag"
        }
        if(j==2){
          time = "unknown"
        }
        kerneltypeTime = ""
        if(j>2){
          time = "kernel"
          kerneltypeTime = lkernelT[[j-2]]
          h = lhT[[j-2]]
        }
        print( c(spectra, time) )
        cov <- new("fitSpectra", m=l[[1]], modelname = "parsimonious"
                   , spectra = spectra, time = time, kerneltypeSpectra = kerneltypeSpectra
                   , kerneltypeTime = kerneltypeTime, h = h, model = listModel[[k]]
                   , validation = TRUE)
        cov <- fit(cov)
        pred <- new("predictClass", m = l[[2]], fittedCov = cov, model = listModel[[k]]
                    , validation = TRUE)
        pred <- predict(pred)
        p[i,j] <- pred@accuracy
      }
    }

    rownames(p) <- paste("S_",c("diag","unknown", as.character(lkernelS)),sep="")
    colnames(p) <- paste("T_",c("diag", "unknown", as.character(lkernelT)),sep="")
    name = c(name, paste(c("full", "p"), "_", listModel[[k]], sep = ""))
    res = c(res,list(full,p))
  }
  names(res) <- name
  res
}


