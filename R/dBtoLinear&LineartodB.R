#' Convert deciBels scales to linear
#'
#' @name dBtoLinear
#'
#' @description Function to convert dB scales. The conversion can be made either from dB to linear (\code{dBtoLinear}) or from linear to dB (\code{\link{LineartodB}}).
#'
#' @usage dBtoLinear(x, factor="IL", ref=1)
#'
#' @param x Numerical. A numeric vector or a numeric matrix with dB values for \code{dBtoLinear} funcion or linear values for \code{\link{LineartodB}}.
#' @param factor Character. Specify in what factor the function should use to convert your data. \code{SPL} (Sound Pressure Level) for amplitude like data (factor \code{20}) or \code{IL} (Intensity Level) for power like (fator \code{10}). (By default \code{IL})
#' @param ref Numerical. Reference value for conversion. For Sound in water the ref is 1 microPa and on air 20 microPa. (By default 1)
#'
#' @details For details about the factor choice, we recommend the reading of \href{https://dspillustrations.com/pages/posts/misc/decibel-conversion-factor-10-or-factor-20.html}{this} web page.
#'
#' @return The same object of the input with the converted values.
#'
#' @author Cássio Rachid Simões <cassiorachid@@gmail.com>
#'
#' @seealso \code{\link{rms.dB}}, \code{\link{LineartodB}}, \code{\link[seewave]{convSPL}}
#'
#' @examples dBtoLinear(c(80,60,65,62))
#' LineartodB(dBtoLinear(c(80,60,65,62)))
#'

dBtoLinear<-function(x, factor="IL", ref=1){
  if(factor == "IL") {
    fac <- 10
  }else if(factor == "SPL"){
    fac <- 20
  }else{stop("Only 'SPL' or 'IL' acepted for factor argument.")
  }

  return(ref*10^(x/fac))

}

#' Convert linear scales to deciBels
#'

#' @name LineartodB
#'
#' @description Function to convert dB scales. The conversion can be made either from dB to linear (\code{\link{dBtoLinear}}) or from linear to dB (\code{LineartodB}).
#'
#' @usage LineartodB(x, factor="IL", ref=1)
#'
#' @param x Numerical. A numeric vector or a numeric matrix with dB values for \code{\link{dBtoLinear}} funcion or linear values for \code{LineartodB}.
#' @param factor Character. Specify in what factor the function should use to convert your data. \code{SPL} (Sound Pressure Level) for amplitude like data (factor \code{20}) or \code{IL} (Intensity Level) for power like (fator \code{10}). (By default \code{IL})
#' @param ref Numerical. Reference value for conversion. For Sound in water the ref is 1 microPa and on air 20 microPa. (By default 1)
#'
#' @details For details about the factor choice, we recommend the reading of \href{https://dspillustrations.com/pages/posts/misc/decibel-conversion-factor-10-or-factor-20.html}{this} web page.
#'
#' @return The same object of the input with the converted values.
#'
#' @author Cássio Rachid Simões <cassiorachid@@gmail.com>
#'
#' @seealso \code{\link{rms.dB}}, \code{\link{dBtoLinear}}, \code{\link[seewave]{convSPL}}
#'
#' @examples dBtoLinear(c(80,60,65,62))
#' LineartodB(dBtoLinear(c(80,60,65,62)))
#'

LineartodB<-function(x, factor="IL", ref=1){
  if(factor == "IL") {
    fac <- 10
  }else if(factor == "SPL"){
    fac <- 20
  }else{stop("Only 'SPL' or 'IL' acepted for factor argument.")
  }

  return(fac*log10(x/ref))

}
