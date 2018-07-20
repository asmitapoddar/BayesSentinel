#-----------------------------------------------------------------------
#     Copyright (C) 2012-2016  Inria
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as
#    published by the Free Software Foundation; either version 2 of the
#    License, or (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public
#    License along with this program; if not, write to the
#    Free Software Foundation, Inc.,
#    59 Temple Place,
#    Suite 330,
#    Boston, MA 02111-1307
#    USA
#
#    Contact : S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
#
#-----------------------------------------------------------------------


#' Read sentinel smoothed data
#'
#' This function read the sentinel smoothed data file
#'
#' @param path the path containing the files
#' @param file the name of the file
#' @param size the number of sample to extract
#'
#' @examples
#' ## the famous sentinel data set
#' \dontrun{
#' multiSentinelSmoothed <- readMultiSentinelSmoothedData()
#' }
#'
#' @return A list with the sentinel data
#' @author Serge Iovleff
#'
readMultiSentinelSmoothedData <- function(path="./data/", file = "gp.hdf5", size=NULL)
{
  f <- h5file(name = paste(path, file, sep = ""))
  y <- f["y"][]
  x <- f["x"][]
  t <- f["dates"][]
  nbSample <- length(as.vector(y))
  if (is.null(size))
  { samples <- 1:nbSample }
  else
  { samples <- sample(1:nbSample, size = size)}

  labels <- list(years1 = as.vector(y)[samples])
  times  <- list(years1 = t)
  spect1 <- list(years1 = as.matrix(x[samples,1:33]))
  spect2 <- list(years1 = as.matrix(x[samples,34:66]))
  spect3 <- list(years1 = as.matrix(x[samples,67:99]))
  spect4 <- list(years1 = as.matrix(x[samples,100:132]))
  spect5 <- list(years1 = as.matrix(x[samples,133:165]))
  spect6 <- list(years1 = as.matrix(x[samples,166:198]))
  spect7 <- list(years1 = as.matrix(x[samples,199:231]))
  spect8 <- list(years1 = as.matrix(x[samples,232:264]))
  spect9 <- list(years1 = as.matrix(x[samples,265:297]))
  spect10 <- list(years1 = as.matrix(x[samples,298:330]))
  clouds <- list(years1 = matrix(0, nrow = length(samples), ncol = length(t) ))

  # return data in a list
  list( labels  = labels
      , times   = times
      , spectra = list( spect1  = spect1
                      , spect2  = spect2
                      , spect3  = spect3
                      , spect4  = spect4
                      , spect5  = spect5
                      , spect6  = spect6
                      , spect7  = spect7
                      , spect8  = spect8
                      , spect9  = spect9
                      , spect10 = spect10
                      )
      , clouds = clouds
      )
}

#' Read sentinel raw data
#'
#' This function read the sentinel raw data files
#'
#' @param path path containing the files
#' @param file name of the file with the raw data
#' @param fileMask name of the file wiht the mask
#' @param size the number of sample to extract
#'
#' @examples
#' ## the famous sentinel data set
#' \dontrun{
#' multiSentinelSmoothed <- readMultiSentinelRawData()
#' }
#'
#' @return A list with the sentinel raw data
#' @author Serge Iovleff
#'
readMultiSentinelRawData <- function(path="./data/", file = "raw.hdf5", fileMask = "mask.hdf5", size=NULL)
{
  f <- h5file(name = paste(path, file, sep = ""))
  times <- f["dates"][]
  labels <- f["y"][]
  spectra <- f["x"][]

  f <- h5file(name = paste(path, fileMask, sep = ""))
  clouds <- f["mask"][]

  nbSample <- length(as.vector(labels))
  if (is.null(size))
  {
    size <- nbSample
    samples <- 1:nbSample
  }
  else
  { samples <- sample(1:nbSample, size = size)}

  times  <- list(years1 = times)
  labels <- list(years1 = as.vector(labels)[samples])
  clouds <- list(years1 = (as.matrix(clouds)[samples,]))
  spectra<- as.matrix(spectra[samples,])

  spect1 <- list(years1 = spectra[,1:30])
  index <- which( spect1$years1 < 0)
  clouds$years1[index] <- 1
  spect2 <- list(years1 = spectra[,31:60])
  index <- which( spect2$years1 < 0)
  clouds$years1[index] <- 1

  spect3 <- list(years1 = spectra[,61:90])
  index <- which( spect3$years1 < 0)
  clouds$years1[index] <- 1

  spect4 <- list(years1 = spectra[,91:120])
  index <- which( spect4$years1 < 0)
  clouds$years1[index] <- 1

  spect5 <- list(years1 = spectra[,121:150])
  index <- which( spect5$years1 < 0)
  clouds$years1[index] <- 1

  spect6 <- list(years1 = spectra[,151:180])
  index <- which( spect6$years1 < 0)
  clouds$years1[index] <- 1

  spect7 <- list(years1 = spectra[,181:210])
  index <- which( spect7$years1 < 0)
  clouds$years1[index] <- 1

  spect8 <- list(years1 = spectra[,211:240])
  index <- which( spect8$years1 < 0)
  clouds$years1[index] <- 1

  spect9 <- list(years1 = spectra[,241:270])
  index <- which( spect9$years1 < 0)
  clouds$years1[index] <- 1

  spect10 <- list(years1 = spectra[,271:300])
  index <- which( spect1$years10 < 0)
  clouds$years1[index] <- 1

  # return data in a list
  list( labels  = labels
      , times   = times
      , spectra = list( spect1  = spect1
                      , spect2  = spect2
                      , spect3  = spect3
                      , spect4  = spect4
                      , spect5  = spect5
                      , spect6  = spect6
                      , spect7  = spect7
                      , spect8  = spect8
                      , spect9  = spect9
                      , spect10 = spect10
                      )
      , clouds = clouds
      )
}
