#' Root Mean Square with dB values
#'
#' @name rms.dB
#'
#' @description Function to compute the root mean square (RMS) of values in decibels (dB).
#'
#' @usage rms.dB(x, level="SPL", ...)
#'
#' @param x Numerical. A numeric vector or a numeric matrix with dB values.
#' @param level Character. Specify in what scale your data is. \code{SPL} for Sound Pressute Level or \code{IL} for Intensity Level.
#' @param ... Further arguments passed to \code{\link[base]{mean}}.
#'
#' @details This function convert your dB data to linear values, compute the Root Mean Square (rms) and converts the result to dB again.
#' @details This function was adapted from \code{\link[seewave]{meandB}} and \code{\link[seewave]{rms}} functions from \code{\link[seewave]{seewave}} package. See their help for more details.
#'
#' @references
#'
#' @return A numeric values that represent the root mean square of x.
#'
#' @author Cássio Rachid Simões <cassiorachid@@gmail.com>
#'
#' @seealso \code{\link[seewave]{meandB}}, \code{\link[seewave]{rms}}
#'
#' @example rms.dB(c(80,60,65,62))
#'
rms.dB<-function(x, level="SPL", ...){
  if(level == "IL") {
    ref <- 10
  }else if(level == "SPL"){
    ref <- 20
  }else{stop("Only 'SPL' or 'IL' acepted for level argument.")
  }

  x <- 10^(x/ref) #Convert to linear

  x <- sqrt(mean(x^2, ...)) #compute the rms

  x <- ref * log10(x) #back to decibels

  return(x)
}
