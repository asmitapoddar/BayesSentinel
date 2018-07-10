percent <- function(list1,list2){
  diff = list1 == list2
  sum(diff)/length(diff)*100
}


# separera la data en 2 parties
testTrain <- function(data,pTrain){
  n = length(data[[1]])
  sampl = runif(n)
  test <- lapply(data[[6]], function(mat,list){mat[list,]},list=which(sampl<pTrain))
  train <- lapply(data[[6]], function(mat,list){mat[list,]},list=which(sampl>=pTrain))
  test <- list(data[[1]][which(sampl<pTrain)],data[[2]],data[[3]],data[[4]],data[[5]],test,data[[7]][which(sampl<pTrain)])
  train <- list(data[[1]][which(sampl>=pTrain)],data[[2]],data[[3]],data[[4]],data[[5]],train,data[[7]][which(sampl>=pTrain)])
  list(test,train)
}


# pourcent predict a partir d'une data
fitPredPrecent <- function(data,modelname = "full", parameter = "time", kernelName = "exponential", h = 10,lambda=0.3,model="gaussian",pTrain=0.1){
  l = testTrain(data,pTrain)
  pourcent(predict(l[[2]],fit(l[[1]],modelname,parameter,kernelName,h),lambda,model),l[[2]][[1]])
}


#### inversion matrice ###
inversion <- function(mat){
  cma <- chol(mat)
  chol2inv(cma)
}

#inversionS <- function(mat){
#  sv <- svd(mat)
#  sv$v %*% diag(1/sv$d) %*% t(sv$u)
#}
