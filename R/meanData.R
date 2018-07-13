meanData = function(data)
{
  meanLabel = function(data,label)
  {
    colMeans(do.call("cbind",lapply(1:length(data[[3]]),meanLabelSpectra,label=label,data=data)))
  }
  meanLabelSpectra = function(data,label,spect)
  {
    data[[3]][[spect]][which(data[[1]]==label),]
  }
  lapply(levels(factor(data[[1]])),meanLabel,data=data)
}
