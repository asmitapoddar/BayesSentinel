# R Package: BayesSentinel
We have high dimensional temporal satellite spectroscopic data with missing values due to clouds, disturbances etc. taken of the region of France. We make a statistical study of the data to understand the means and covariance between the different spectra as well as the covariance of data between times of sampling. We use the means and covariances to simulate the spectroscopic data.  

We then use different models to estimate the covariance of the spectroscopic data and implement statistical models to classify the data into different classes. We carry out classification using probabilistic models (mixing models) for which we use the conditional independence hypothesis (modeling in the native space). This R Package performs the simulation and classification task as described above.

## Data
The data is obtained from the [Sentinel-2 satellite] (http://www.cesbio.ups-tlse.fr/us/index_sentinel2.html).
The data can be downloaded from: http://www.cesbio.ups-tlse.fr/multitemp/?p=8547

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
Open BayesSentinel.Rproj
Load the package by:
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

