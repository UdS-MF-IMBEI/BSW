#' @title Extracting the estimated model parameters of \code{bsw()}
#' @description For objects of class \code{"bsw"}, \code{coef()} extracts the estimated model parameters of \code{bsw()}.
#' @param object An object of class \code{"bsw"}.
#' @return A numeric vector containing the estimated model parameters.
#' @aliases coef,bsw-method
#' @author Adam Bekhit, Jakob Schöpe
#' @export

setMethod(f = "coef",
          signature = "bsw",
          definition = function(object) {
            return(object@coefficients)
          }
)
