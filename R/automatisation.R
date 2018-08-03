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
#' @param data the 3 dimentional data
#' @param pTrain the percent of training sample
#'
#' @return list containing training and test set
#'

testTrain <- function(data, pTrain)
{
  n = length(data[[1]])
  sampl = runif(n)
  train <- lapply(data[[3]], function(mat, list) {mat[list,]}, list = which(sampl < pTrain))
  test <- lapply(data[[3]], function(mat, list) {mat[list,]}
                 , list = which(sampl >= pTrain))
  test <- list(data[[1]][which(sampl < pTrain)], data[[2]], test
               , data[[4]][which(sampl < pTrain)])
  train <- list(data[[1]][which(sampl >= pTrain)]
                , data[[2]], train,data[[4]][which(sampl >= pTrain)])
  list(train = train, test = test)
}


#' inversion
#'
#' Invert a symetrique matrix
#'
#' @param mat matrix to be inverted
#'
#' @return inverted matrix

inversion <- function(mat){
  chol2inv(chol(mat))
}


#' tablesModel
#'
#' Try all the model posible with validation.
#'
#' @param data th 3 dimentinal matrix
#' @param lkernelR list of kernels for row
#' @param lkernelC list of kernels for column
#' @param lhR list of widths for row
#' @param lhC list of widths for row
#' @param listModel list of models
#'
#' @return list containing the accuaracy for different models
#'
#' @export

tablesModel = function( data, lkernelR = list(), lkernelC = list(), lhR = list()
                        , lhC = list(), listModel = list("gaussian"))
{
  #l = testTrain(data,0.8)
  l = list(data,data)
  nr = 2 + length(lkernelR)
  nc = 2 + length(lkernelC)
  res = list()
  name = c()
  for(k in 1:length(listModel))
  {
    #cov <- new("fitData",l[[1]])
    #cov <- fit(cov)
    cov = fitDataMatrix(l[[1]])
    pred <- new("predictClass", m = l[[2]], fittedCov = cov, model = listModel[[k]]
                , validation = TRUE)
    pred <- predictData(pred)
    full <- pred@accuracy
    p = matrix(ncol = nc, nrow = nr)
    for(i in 1:nr){
      if(i==1){
        row = "diag"
      }
      if(i==2){
        row = "unknown"
      }
      kerneltypeRow = ""
      h=0
      if(i>2){
        row = "kernel"
        kerneltypeRow = lkernelR[[i-2]]
        h = lhR[[i-2]]
      }

      for(j in 1:nc){
        if(j==1){
          column = "diag"
        }
        if(j==2){
          column = "unknown"
        }
        kerneltypeCol = ""
        if(j>2){
          column = "kernel"
          kerneltypeCol = lkernelC[[j-2]]
          h = lhC[[j-2]]
        }
        print( c(row, column) )
        cov <- new("fitData", m=l[[1]], modelname = "parsimonious"
                   , row = row, column = column, kerneltypeRow = kerneltypeRow
                   , kerneltypeCol = kerneltypeCol, h = h, model = listModel[[k]]
                   , validation = TRUE)
        cov <- fit(cov)
        pred <- new("predictClass", m = l[[2]], fittedCov = cov, model = listModel[[k]]
                    , validation = TRUE)
        pred <- predict(pred)
        p[i,j] <- pred@accuracy
      }
    }

    rownames(p) <- paste("Row_",c("diag","unknown", as.character(lkernelR)),sep="")
    colnames(p) <- paste("T_",c("diag", "unknown", as.character(lkernelC)),sep="")
    name = c(name, paste(c("full", "p"), "_", listModel[[k]], sep = ""))
    res = c(res,list(full,p))
  }
  names(res) <- name
  res
}


