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
NULL

#' Interface base Class [\code{\linkS4class{ICloHeModel}}] for CloHe models.
#'
#' This class encapsulate the common parameters of all the CloHe models.
#'
#' @slot classNumber the integer with the class attributed to the pixel
#' @slot classLabel the integer attributed to the pixel
#' @slot nk the number of sample in this class
#' @slot pk the proportion of sample in this class
#' @slot nbSpectrum number of spectrum to analyze
#' @slot lnLikelihood log-likelihood of this
#' @slot criterion the value of the criterion
#' @slot nbFreeParameter number of free parameter
#'
#' @examples
#'   getSlots("ICloHeModel")
#'
#' @author Serge Iovleff
#'
#' @name ICloHeModel
#' @aliases ICloHeModel-class
#' @rdname ICloHeModel-class
#'
setClass(
  Class="ICloHeModel",
  representation( classNumber     = "numeric"
                , classLabel      = "numeric"
                , nk              = "numeric"
                , pk              = "numeric"
                , nbSpectrum      = "numeric"
                , lnLikelihood    = "numeric"
                , criterion       = "numeric"
                , nbFreeParameter = "numeric"
                , "VIRTUAL"
                ),
                prototype( classNumber  = integer(0)
                         , classLabel   = integer(0)
                         , nk           = integer(0)
                         , pk           = numeric(0)
                         , lnLikelihood = -.Machine$double.xmax
                         , criterion    =  -.Machine$double.xmax
                         , nbFreeParameter = integer(0)
                 ),
  # validity function
  validity = function(object)
  {
    # check classNumper
    if (round(object@classNumber) != object@classNumber)
    { stop("classNumber must be an integer.")}
    # check classNumper
    if (round(object@classLabel) != object@classLabel)
    { stop("classLabel must be an integer.")}
    if (round(object@nk) != object@nk)
    { stop("nk must be an integer.")}
    if (round(object@nbSpectrum) != object@nbSpectrum)
    { stop("nnSpectrum must be an integer.")}
    # check nbFreeParameter
    if (round(object@nbFreeParameter) != object@nbFreeParameter)
    { stop("nbFreeParameter must be an integer.")}
    return(TRUE)
  }
)

#' Initialize an instance of a CloHe S4 class.
#'
#' Initialization method of the [\code{\linkS4class{ICloHeModel}}] class.
#' Used internally in the 'CloHe' package.
#'
#' @rdname initialize-methods
#' @keywords internal
setMethod(
    f = "initialize",
    signature = c("ICloHeModel"),
    definition=function(.Object)
    {
      .Object@classNumber = 0
      .Object@classLabel  = 0
      .Object@nk          = 0
      .Object@pk          = 0
      .Object@nbSpectrum  = 0
      .Object@nbFreeParameter = 0
      return(.Object)
    }
)


#' @rdname print-methods
#' @aliases print,ICloHeModel-method
#'
setMethod(
  f = "print",
  signature = c("ICloHeModel"),
  function(x,...)
  {
    cat("* Class Number   = ", x@classNumber,"\n")
    cat("* Class Label    = ", x@classLabel,"\n")
    cat("* nk             = ", x@nk,"\n")
    cat("* pk             = ", x@pk,"\n")
    cat("* nbSpectrum     = ", x@nbSpectrum,"\n")
    cat("* lnLikelihood   = ", x@lnLikelihood,"\n")
    cat("* nbFreeParameter= ", x@nbFreeParameter,"\n")
    cat("* criterion      = ", x@criterion, "\n")
  }
)

#' @rdname show-methods
#' @aliases show show,ICloHeModel-method
setMethod(
  f = "show",
  signature = c("ICloHeModel"),
  function(object)
  {
    cat("* Class Number   = ", object@classNumber,"\n")
    cat("* Class Label    = ", object@classLabel,"\n")
    cat("* nk             = ", object@nk,"\n")
    cat("* pk             = ", object@pk,"\n")
    cat("* nbSpectrum     = ", object@nbSpectrum,"\n")
    cat("* lnLikelihood   = ", object@lnLikelihood,"\n")
    cat("* nbFreeParameter= ", object@nbFreeParameter,"\n")
    cat("* criterion      = ", object@criterion, "\n")
  }
)

#' @rdname summary-methods
#' @aliases summary summary,ICloHeModel-method
setMethod(
  f = "summary",
  signature = c("ICloHeModel"),
  function(object,...)
  {
    cat("* Class Number   = ", object@classNumber,"\n")
    cat("* Class Label    = ", object@classLabel,"\n")
    cat("* nk             = ", object@nk,"\n")
    cat("* pk             = ", object@pk,"\n")
    cat("* nbSpectrum     = ", object@nbSpectrum,"\n")
    cat("* lnLikelihood   = ", object@lnLikelihood,"\n")
    cat("* nbFreeParameter= ", object@nbFreeParameter, "\n")
    cat("* criterion      = ", object@criterion, "\n")
  }
)

