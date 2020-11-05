#' Root Mean Square with dB values
#'
#' @name rms.dB
#'
#' @description Function to compute the root mean square (RMS) of values in decibels (dB).
#'
#' @usage rms.dB(x, level="SPL", ref=1, ...)
#'
#' @param x Numerical. A numeric vector or a numeric matrix with dB values.
#' @param level Character. Specify in what scale your data is. \code{SPL} for Sound Pressute Level or \code{IL} for Intensity Level. (By deafault \code{SPL})
#' @param ref Numerical. Reference value for conversion. For Sound in water the ref is 1microPa and on air 20 microPa. (By default 1)
#' @param ... Further arguments passed to \code{\link[base]{mean}}.
#'
#' @details This function converts your dB data to linear values (through \code{\link{dBtoLinear}} function), compute the Root Mean Square (rms), and converts the result back to dB (through \code{\link{LineartodB}} function).
#' @details This function was adapted from \code{\link[seewave]{meandB}} and \code{\link[seewave]{rms}} functions from \code{\link[seewave]{seewave}} package. See their help for more details.
#'
#'
#' @return A numeric value that represents the root mean square of x.
#'
#' @author Cássio Rachid Simões <cassiorachid@@gmail.com>
#'
#' @seealso \code{\link[seewave]{meandB}}, \code{\link[seewave]{rms}}
#'
#' @examples rms.dB(c(80,60,65,62))
#'
rms.dB<-function(x, level="SPL", ref=1, ...){

  x <- dBtoLinear(x, factor=level, ref=ref) #Convert to linear

  x <- sqrt(mean(x^2, ...)) #compute the rms

  x <- LineartodB(x, factor=level, ref=1) #back to decibels

  return(x)
}