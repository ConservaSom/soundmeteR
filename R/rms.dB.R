#' Root Mean Square with dB values
#'
#' @name rms.dB
#'
#' @description Function to compute the root mean square (RMS) of values in decibels (dB).
#'
#' @usage rms.dB(x, level="SPL", ref=1, ...)
#'
#' @param x Numerical. A numeric vector or a numeric matrix with dB values.
#' @param level Character. Specify in what scale your data is. \code{SPL} for Sound Pressute Level or \code{IL} for Intensity Level. (By default \code{SPL})
#' @param ref Numerical. Reference value for conversion. For Sound in water the ref is 1microPa and on air 20 microPa. (By default 1)
#' @param na.rm Logical. Argument passed to \code{\link[base]{mean}}. Should NA be removed? (By default FALSE)
#'
#' @details This function converts your dB data to linear values (through \code{\link{dBtoLinear}} function), compute the Root Mean Square (rms), and converts the result back to dB (through \code{\link{LineartodB}} function).
#' @details This function was adapted from \code{\link[seewave]{meandB}} and \code{\link[seewave]{rms}} functions from \code{\link[seewave]{seewave}} package. See their help for more details.
#'
#'
#' @return A numeric value that represents the root mean square of x.
#'
#'
#' @seealso \code{\link[seewave]{meandB}}, \code{\link[seewave]{rms}}
#'
#' @examples
#' rms.dB(c(80,60,65,62))
#'
#' @export

rms.dB<-function(x, level="SPL", na.rm=FALSE){

  dBtoLinear(x, factor=level) %>%
    .^2 %>%
    mean(na.rm=na.rm) %>%
    sqrt() %>%
    LineartodB(factor=level) %>%
    return()

}
