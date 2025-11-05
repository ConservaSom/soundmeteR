#' Convert linear scales to deciBels
#'
#' @name LineartodB
#'
#' @description Function to convert linear values (µPa) to decibels (dB).
#'
#' @usage LineartodB(x, factor="IL", ref=1)
#'
#' @param x Numeric. A vector or matrix containing linear values (in µPa).
#' @param factor Character. Specify the factor for conversion. Use \code{SPL} (Sound Pressure Level) for amplitude-base data, with factor \code{20}, or \code{IL} (Intensity Level) for power-base data, with factor \code{10} (by default \code{IL}.)
#' @param ref Numeric. The reference value for dB conversion. For sound in water, the common reference is 1 µPa, and for sound in air,it is 20 µPa (by default 1).
#'
#' @details For details about the factor choice, we recommend consulting \href{https://dspillustrations.com/pages/posts/misc/decibel-conversion-factor-10-or-factor-20.html}{this} web page.
#'
#' @return The same input object with the values converted.
#'
#'
#' @seealso \code{\link{rms.dB}}, \code{\link{dBtoLinear}}, \code{\link[seewave]{convSPL}}
#'
#' @examples dBtoLinear(c(80,60,65,62))
#' LineartodB(dBtoLinear(c(80,60,65,62)))
#'
#' @export

LineartodB<-function(x, factor="IL", ref=1){
  if(factor == "IL") {
    fac <- 10
  }else if(factor == "SPL"){
    fac <- 20
  }else{stop("Only 'SPL' or 'IL' acepted for factor argument.")
  }

  return(fac*log10(x/ref))

}
