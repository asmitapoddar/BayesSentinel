library(BayesSentinel)
sim = simulateSpectra()
cov = fitSpectra(sim)
objPred = predictClass(sim, cov, listLambdaS = (), listLambdaT = () )

predict2 = function(objPred, lambdaS, listLambdaT)
{
  objPred@lambdaS = lambdaS
  lapply(listLambdaT,
         function(objPred,lambdaT)
         { objPred@lambdaT = lambdaT
         predict(objPred)
         }
         , objPred = objPred)
}

p = lapply(objPred@listLambdaS, predict2, objPred = objPred
           , listLambdaT = objPred@listLambdaT)
perc = vapply(p, function(list)
                        { vapply( list, function(pred) {pred@accuracy}
                         , FUN.VALUE = vector('double',length = 1))
                        }
              , FUN.VALUE = vector('double',length = length(p[[1]])))
lambda = c(0)
if(length(objPred@listLambdaS)==1 | length(objPred@listLambdaT)!=1)
{
  lambda = objPred@listLambdaT
  plot(x = lambda, y = perc, type = 'l')
  title(paste(objPred@fittedCov$modelname, objPred@fittedCov$spectra
              ,objPred@fittedCov$time) )
  l = list(lambdaS = lambda[which.max(perc)], lambdaT = lambda[which.max(perc)]
           , predicted = p[[1]][[which.max(perc)]]@predicted_labels, percent = max(perc))
}
if(length(listLambdaS)!=1 | length(listLambdaT)==1)
{
  lambda = objPred@listLambdaS
  plot(x = lambda, y = perc, type = 'l')
  title(paste(objPred@fittedCov$modelname, objPred@fittedCov$spectra
              ,objPred@fittedCov$time) )
  l = list(lambdaS = lambda[which.max(perc)], lambdaT = lambda[which.max(perc)]
           , predicted = p[[which.max(perc)]][[1]]@predicted_labels, percent = max(perc))
}
