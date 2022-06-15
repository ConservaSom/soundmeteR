#' Sum with dB values
#'
#' @details This function computes the sum of dB values.
#'
#' @param x Numerical. A numeric vector or matrix with dB values.
#' @param level Character. Specify in what scale your data is. \code{SPL} for Sound Pressute Level or \code{IL} for Intensity Level. (By default \code{SPL})
#' @param na.rm Logical. Argument passed to \code{\link[base]{sum}}. Should NA be removed? (By default FALSE)
#'
#' @details This function converts your dB data to linear values (through \code{\link{dBtoLinear}} function), compute the sum, and converts the result back to dB (through \code{\link{LineartodB}} function).
#'
#' @return A numeric value that represents the sum of x.
#'
#' @seealso \code{\link[seewave]{moredB}}
#'
#' @examples
#' sumdB(c(80,60,65,62))
#' sumdB(c(30,30), level="IL")
#' sumdB(c(30,30), level="SPL")
#'
#' @export

sumdB<-function(x, level="IL", na.rm=FALSE, ...){

  dBtoLinear(x, factor=level) %>%
    sum(na.rm=na.rm) %>%
    LineartodB(factor=level) %>%
    return()

}
