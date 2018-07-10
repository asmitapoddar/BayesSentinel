full = function(data)
{
  corLabel = function(data,label)
  {
    cov(do.call("cbind",lapply(1:length(data[[6]]),corLabelSpectra,label=label,data=data)))
  }
  corLabelSpectra = function(data,label,spect)
  {
    data[[6]][[spect]][which(data[[1]]==label),]
  }
  lapply(levels(factor(data[[1]])),corLabel,data=data)

  #lambda = matrix(0.2, nrow = nrow(covMat[[1]]), ncol = ncol(covMat[[1]]) )   #do it in the regularisation
  #covMat = lapply(covMat, function(x) {x+lambda%*%diag(nrow(x))})
}
