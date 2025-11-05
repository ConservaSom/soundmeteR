#' Convert deciBels scales to linear
#'
#' @name dBtoLinear
#'
#' @description Function to convert decibels (dB) to linear values (µPa).
#'
#' @usage dBtoLinear(x, factor="IL", ref=1)
#'
#' @param x Numeric. A numeric vector or a numeric matrix with dB values.
#' @param factor Character. Specify the factor to use for conversion. Use \code{SPL} (Sound Pressure Level) for amplitude-base data, with factor \code{20}, or \code{IL} (Intensity Level) for power-base data, with fator \code{10} (by default "IL").
#' @param ref Numeric. The reference value for dB conversion. For sound in water, the common reference is 1 µPa, and for sound in air,it is 20 µPa (by default 1).
#'
#' @details For further details on selecting the appropriate factor, we recommend consulting \href{https://dspillustrations.com/pages/posts/misc/decibel-conversion-factor-10-or-factor-20.html}{this} web page.
#'
#' @return The same input object with the values converted.
#'
#'
#' @seealso \code{\link{rms.dB}}, \code{\link{LineartodB}}, \code{\link[seewave]{convSPL}}
#'
#' @examples dBtoLinear(c(80,60,65,62))
#' LineartodB(dBtoLinear(c(80,60,65,62)))
#'
#' @export

dBtoLinear<-function(x, factor="IL", ref=1){
  if(factor == "IL") {
    fac <- 10
  }else if(factor == "SPL"){
    fac <- 20
  }else{stop("Only 'SPL' or 'IL' acepted for factor argument.")
  }

  return(ref*10^(x/fac))

}
