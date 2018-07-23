library(BayesSentinel)
sim = simulateSpectra()
cov = fitSpectra(sim)
pred = predictClass(sim,cov)
pred@accuracy
