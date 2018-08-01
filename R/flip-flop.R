#-----------------------------------------------------------------------
#' Flip - Flop
#'
#' Implementing the flip flop algorithm
#'
#' @param data sata sent
#' @param list of means
#' @param s parameter for regularisation
#' @param lambdaR list of lambda for row
#' @param lambdaC list of lambda for column
#'
#' @return list of predicted sigma for every label
#'
#' @name flipflop
#' @export flipflop
#'
#' @author Asmita Poddar & Florent Latimier

flipflop = function(data,lmean,s,lambdaR,lambdaC)
{

  flipflopLabel = function(data, mean, s, lambdaR=0.1, lambdaC = 0.1)
  {
    if(missing(data))  { stop("data is mandatory")}
    ns = dim(data)[3]
    nt = dim(data)[2]
    n = dim(data)[1]

    mean = matrix(mean,ncol=ns)
    X2 = apply(data,1,function(x,mean){x-mean},mean = mean)

    sigmaRold = 0
    sigmaRnew = diag(rep(1,ns))
    sigmaCold = 0
    sigmaCnew = diag(rep(1/sqrt(nt),nt))

    l = 0

    while(l<100 && testFlop(sigmaRold,sigmaRnew) > s && testFlop(sigmaCold,sigmaCnew) > s)
    {
      l <- l+1

      invS = inversion(sigmaRnew+diag(rep(lambdaR,ns) ))
      sigmaCold <- sigmaCnew
      sigmaCnew <- calculFlip(X = data,S = invS,nt = nt)

      invT = inversion(sigmaCnew+diag(rep(lambdaC,nt)))
      sigmaRold <- sigmaRnew
      sigmaRnew <- calculFlop(X=data,S = invT,ns=ns)

    }
    list(sigmaR = sigmaRnew, sigmaC = sigmaCnew)
  }

  testFlop = function(old, new)
  {
    sqrt(sum((new - old)^2)) / sqrt(sum((old)^2))
  }

  calculFlip = function(X, S, nt)
  {
    s = matrix(rowMeans( apply(X, 1, function(x,S){x %*% S %*% t(x)}, S = S)), ncol = nt)
    (s + t(s))/2.
  }

  calculFlop = function(X, S, ns)
  {
    s = matrix(rowMeans(apply(X, 1, function(x, S){t(x) %*% S %*% x}, S = S)), ncol = ns)
    (s + t(s))/2.
  }

  X = array(as.numeric(unlist(data[[3]]))
            ,dim = c(nrow(data[[3]][[1]]), ncol(data[[3]][[1]]), length(data[[3]])))
  lapply(levels(factor(data[[1]]))
         , function(data,lmean,label,s,lambdaR,lambdaC)
         {flipflopLabel(data = X[which(data[[1]]==label),,]
                        , mean = lmean[[which(levels(factor(data[[1]]))==label)]]
                        , s = s, lambdaR = lambdaR, lambdaC = lambdaC)}
         , data = data, lmean = lmean, s = s, lambdaR = lambdaR,lambdaC = lambdaC)

}


