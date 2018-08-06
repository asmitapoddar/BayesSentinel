#' Plot of the fitted mean
#'
#' Plot all the row means for a label dependind on the columns.
#'
#' @param data 3 dimentional data
#' @param objFitted object give by the fitting
#'
#' @name plot.fitData
#' @export plot.fitData
#'
plot.fitData <- function(data,objFitted){
  means = objFitted$mean

  graphMean = function(data,mean,int)
  {
    color = rainbow(length(data$rows))
    plot(x=data$column,y=rep(0,length(data$column)),ylab="mean",type = "l",col="white",ylim = c(min(mean),max(mean)))
    m = t(matrix(mean,ncol = length(data$row)))
    for(i in 1:nrow(m))
    {
      lines(x=data$column,y=m[i,],col=color[i])
    }
    title(main = paste("Mean for the label : ",data$labels[int]))
    legend(-10,max(mean), legend=paste('row',1:length(data$rows)),col=color, lty=1, cex=0.75)
  }



  for(j in 1:length(means))
  {
    graphMean(data,means[[j]],j)
  }
}




