library(BayesSentinel)

m = new("simulateData")
res = simulate(m)

plot(res$column,res$means[1,1,],col=1,type='l',ylim=c(100,500))
for(j in 2:10)
{ lines(res$column,res$means[1,j,],col=1)
  lines(res$column,res$means[2,j,],col=2)
  lines(res$column,res$means[3,j,],col=3)
  lines(res$column,res$means[4,j,],col=4)
}

simClass <- res$rows$row1[which(res$labels==1),]
for(j in dim(simClass)[1])
{ lines(res$column,simClass[j,],col="orange")}

simClass <- res$rows$row1[which(res$labels==2),]
for(j in dim(simClass)[1])
{ lines(res$column,simClass[j,],col="orange")}

simClass <- res$rows$row1[which(res$labels==3),]
for(j in dim(simClass)[1])
{ lines(res$column,simClass[j,],col="orange")}

simClass <- res$rows$row1[which(res$labels==4),]
for(j in dim(simClass)[1])
{ lines(res$column,simClass[j,],col="orange")}

