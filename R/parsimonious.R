### Full Row
vectorRowLabel <- function(data,row,label)
{
  as.numeric(data[[3]][[row]][which(data[[1]]==label),])
}

covRLabel <- function(data,label)
{
  cov(do.call("cbind",lapply(1:length(data[[3]]),vectorRowLabel,label=label,data=data)))
}

#' Unknown row in parsimonious model
#'
#' Return the list of covariance in row for each clusther in case of a parsimonious model with unknown rows and known column.
#'
#' @param data the data with matrix observation
#'
#' @name unknownRow
#' @export unknownRow
#'
unknownRow <- function(data)
{
  lapply(levels(factor(data[[1]])),covRLabel,data=data)
}


### Full column
matrixRowLabel <- function(data,row,label)
{
  data[[3]][[row]][which(data[[1]]==label),]
}

covTLabel <- function(data,label)
{
  cov(do.call("rbind",lapply(1:length(data[[3]]),matrixRowLabel,label=label,data=data)))
}

#' Unknown column in parsimonious model
#'
#' Return the list of covariance in column for each clusther in case of a parsimonious model with unknown column and known rows
#'
#' @param data the data with matrix observation
#'
#' @name unknownCol
#' @export unknownCol
#'
unknownCol <- function(data)
{
  lapply(levels(factor(data[[1]])),covTLabel,data=data)
}



### Parsimonious

#' Diagonal rows in parsimonious model
#'
#' Return the list of covariance in row for each clusther in case of a parsimonious model with diagonal row
#'
#' @param data the data with matrix observation
#'
#' @name diagRow
#' @export diagRow
#'
diagRow <- function(data)
{
  lA = lapply(fullRow(data),diag)
  lapply(lA,diag)
}


#' Diagonal column in parsimonious model
#'
#' Return the list of covariance in column for each clusther in case of a parsimonious model with diagonal column
#'
#' @param data the data with matrix observation
#'
#' @name diagCol
#' @export diagCol
#'
diagCol <- function(data)
{
  lA = lapply(fullCol(data),diag)
  lapply(lA,diag)
}
