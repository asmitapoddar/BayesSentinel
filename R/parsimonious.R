### Full Spectra
vectorSpectraLabel <- function(data,spectra,label){
  as.numeric(data[[6]][[spectra]][which(data[[1]]==label),])
}

covSLabel <- function(data,label){
  cov(do.call("cbind",lapply(1:length(data[[6]]),vectorSpectraLabel,label=label,data=data)))
}

fullSpectra <- function(data){
  lapply(levels(factor(data[[1]])),covSLabel,data=data)
}


### Full time
matrixSpectraLabel <- function(data,spectra,label){
  data[[6]][[spectra]][which(data[[1]]==label),]
}

covTLabel <- function(data,label){
  cov(do.call("rbind",lapply(1:length(data[[6]]),matrixSpectraLabel,label=label,data=data)))
}

fullTime <- function(data){
  lapply(levels(factor(data[[1]])),covTLabel,data=data)
}



### Parsionious

parsimoniousSpectra <- function(data){
  lA = lapply(fullSpectra(data),diag)
  lapply(lA,diag)
}

parsimoniousTime <- function(data){
  lA = lapply(fullTime(data),diag)
  lapply(lA,diag)
}
