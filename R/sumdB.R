#' Sum with dB values
#'
#' @details Function to compute the sum of decibels values (dB).
#'
#' @param x Numerical. A vector or a matrix with decibels values (dB).
#' @param level Character. Specify the factor for conversion. Use \code{SPL} (Sound Pressure Level) for amplitude-base data, with factor 20, or \code{IL} (Intensity Level) for power-base data, with fator 10 (By default \code{IL}).
#' @param na.rm Logical. Argument passed to \code{\link[base]{sum}}. If \code{TRUE} removes NA (by default \code{FALSE}).
#'
#' @details This function converts the decibels data to linear values using the \code{\link{dBtoLinear}} function, computes the sum, and then converts the result back to decibels using the \code{\link{LineartodB}} function.
#'
#' @return A numeric value representing the sum of x.
#'
#' @seealso \code{\link[seewave]{moredB}}
#'
#' @examples
#' sumdB(c(80,60,65,62))
#' sumdB(c(30,30), level="IL")
#' sumdB(c(30,30), level="SPL")
#'
#' @export

sumdB<-function(x, level="IL", na.rm=FALSE){

  dBtoLinear(x, factor=level) %>%
    sum(na.rm=na.rm) %>%
    LineartodB(factor=level) %>%
    return()

}
