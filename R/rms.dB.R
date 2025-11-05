#' Root Mean Square with dB values
#'
#' @description Function to compute the root mean square (RMS) of values in decibels (dB).
#'
#' @param x Numeric. A vector or a matrix with decibels values (dB).
#' @param level Character. Specify the factor for conversion. Use \code{SPL} (Sound Pressure Level) for amplitude-base data, with factor 20, or \code{IL} (Intensity Level) for power-base data, with fator 10 (By default \code{SPL}).
#' @param na.rm Logical. Argument passed to \code{\link[base]{mean}}. If \code{TRUE} removes NA (By default \code{FALSE}).
#'
#' @details This function converts the decibels data to linear values using the \code{\link{dBtoLinear}} function, computes the Root Mean Square (rms), and then converts the result back to decibels using the \code{\link{LineartodB}} function.
#' @details This function was adapted from the \code{\link[seewave]{meandB}} and \code{\link[seewave]{rms}} functions from \code{\link[seewave]{seewave}} package. See their documentation for more details.
#'
#' @return A numeric value that representing the root mean square of x.
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
    sqrt(.) %>%
    LineartodB(factor=level) %>%
    return()

}
