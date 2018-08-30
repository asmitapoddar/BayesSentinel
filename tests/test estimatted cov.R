estimatedCov = function(fitM){
  res = list()
  for(i in 1:length(fitM[[1]])){
    res[[i]] <- fitM[[1]][[i]] %x% fitM[[2]][[i]]
  }
  res
}

for(i in 1:15){
print(sum((estimatedCov(fitm)[[i]]-m$sigma[[i]]))/330)
}
