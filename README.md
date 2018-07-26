# BayesSentinel
Package for simulation and classification of multi-dimensional spectroscopic data

## Data
The data is obtained from the [Sentinel-2 satellite] (http://www.cesbio.ups-tlse.fr/us/index_sentinel2.html)

## Environment
R (>= 3.0.2)

## Dependencies
- stats
- graphics
- h5
- methods
- mvtnorm
- matrixcalc

## Usage

##### 1. Cloning the repository.
```
git clone https://github.com/asmitapoddar/BayesSentinel.git
cd BayesSentinel
```

##### 2. Installing the package
```
library(BayesSentinel)  
```
##### 3. Example
```
library(BayesSentinel)  
sim = simulateSpectra()
cov = fitSpectra(sim)
pred = predictClass(sim,cov)
pred@accuracy
```
This will create an S4 object to simulate a dataset with the default values, fit a covarance matrix to the data and use the fitted covarance matrix to classify the data according to the Bayes Classification Rule. An S4 object containing the prediction accuracy is returned.

