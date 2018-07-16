flipflop = function(data,lmean,s,lambdaS,lambdaT){
  X = array(as.numeric(unlist(data[[3]])),dim = c(nrow(data[[3]][[1]]),ncol(data[[3]][[1]]),length(data[[3]])))
  lapply(levels(factor(data[[1]])), function(data,lmean,label,s,lambdaS,lambdaT){flipflopLabel(data = X[which(data[[1]]==label),,],mean = lmean[[which(levels(factor(data[[1]]))==label)]],s = s,lambdaS = lambdaS,lambdaT = lambdaT)},data=data,lmean=lmean,s=s,lambdaS=lambdaS,lambdaT=lambdaT)
}



flipflopLabel = function(data,mean,s, lambdaS=0.1, lambdaT = 0.1)
{
  if(missing(data))  { stop("data is mandatory")}
  ns = dim(data)[3]
  nt = dim(data)[2]
  n = dim(data)[1]

  mean = matrix(mean,ncol=ns)
  X2 = apply(data,1,function(x,mean){x-mean},mean = mean)

  SigmaSold = 0
  SigmaSnew = diag(rep(1,ns))
  SigmaTold = 0
  SigmaTnew = diag(rep(1/sqrt(nt),nt))

  l = 0

  while(l<100 && testFlop(SigmaSold,SigmaSnew) > s && testFlop(SigmaTold,SigmaTnew) > s)
    {
    l <- l+1

    invS = inversion(SigmaSnew+diag(rep(lambdaS,ns) ))
    SigmaTold <- SigmaTnew
    SigmaTnew <- calculFlip(X = data,S = invS,nt = nt)

    invT = inversion(SigmaTnew+diag(rep(lambdaT,nt)))
    SigmaSold <- SigmaSnew
    SigmaSnew <- calculFlop(X=data,S = invT,ns=ns)

  }
  list(SigmaS = SigmaSnew,SigmaT = SigmaTnew)
}




testFlop = function(old, new){
  sqrt(sum((new - old)^2)) / sqrt(sum((old)^2))
}


calculFlip = function(X,S,nt){
  s = matrix(rowMeans(apply(X,1,function(x,S){x %*% S %*% t(x)},S=S)),ncol=nt)
  (s + t(s))/2.
}

calculFlop = function(X,S,ns){
  s = matrix(rowMeans(apply(X,1,function(x,S){t(x) %*% S %*% x},S=S)),ncol=ns)
  (s + t(s))/2.
}
