library(BayesSentinel)
sim = simulateSpectra()
cov = fitSpectra(sim, modelname = "parsimonious", spectra = "diag", time = "unknown")
objPred = predictClass(sim, cov, validation = TRUE)

