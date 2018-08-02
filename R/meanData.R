#' meanData
#'
#' Calculate the mean for each label, in row and in column.
#'
#' @param data object of class fitData
#'
#' @return all the means according to labels
#'
#' @name meanData
#' @export meanData
#'
meanData = function(data)
{
  meanLabel = function(data,label)
  {
    colMeans(do.call("cbind",lapply(1:length(data[[3]]),meanLabelRow,label=label,data=data)))
  }
  meanLabelRow = function(data,label,row)
  {
    data[[3]][[row]][which(data[[1]]==label),]
  }
  lapply(levels(factor(data[[1]])),meanLabel,data=data)
}
