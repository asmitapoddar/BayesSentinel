#' Full model for fitting
#'
#' Return the list of covariance for each clusther in case of a full model.
#'
#' @param data the data with matrix observation
#'
#' @name full
#' @export full
#'
full = function(data)
{
  corLabel = function(data,label)
  {
    cov(do.call("cbind",lapply(1:length(data[[3]]),corLabelRow,label=label,data=data)))
  }
  corLabelRow = function(data,label,row)
  {
    data[[3]][[row]][which(data[[1]]==label),]
  }
  lapply(levels(factor(data[[1]])),corLabel,data=data)
}
