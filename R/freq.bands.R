#' Interval Limits
#'
#' @description Function to compute frequency intervals based on patterns defined by the user. This function is adapted from the \code{\link[seewave]{octaves}} functions from \code{\link[seewave]{seewave}} package.
#'
#' @param x Numerical. Frequency used as the central point from which bands are calculated in kHz (by default 2).
#' @param interval Numeric. Specifies the interval pattern to be applied. See details for further information.
#' @param below Numerical. Number of intervals below x (by default 3).
#' @param above Numerical. Number of intervals above x (by default 3).
#'
#' @details The specified interval is applied as x/interval (for the lower bound) and x×interval (for the upper bound). Some examples of intervals that can be applied are:
#' \itemize{
#'   \item Octaves = 2
#'   \item Third of an octave = 2^(1/3)
#'   \item Perfect fifth (music theory) = 3/2
#'   \item Major third (music theory) = 5/4
#'   \item Minor third (music theory) = 6/5
#'   \item For additional values that can be used, please refer to this \href{https://academics.hamilton.edu/music/spellman/class_notes/music_theory.htm}{link}.
#'}
#'
#' @return A numeric vector containing the frequency limits of each interval.
#'
#'
#' @seealso \code{\link[seewave]{octaves}}
#'
#' @examples
#' freq.bands(1000, interval=2, below = 1, above = 1) #octaves
#' freq.bands(1000, interval=2^(1/3), below = 3, above = 3) #Third of octaves
#'
#' #https://academics.hamilton.edu/music/spellman/class_notes/music_theory.htm
#' freq.bands(440, interval=3/2, below = 0, above = 1) #Perfect fifth (music theory) of A 440
#' freq.bands(440, interval=5/4, below = 0, above = 1) #Major third (music theory) of A 440
#' freq.bands(440, interval=6/5, below = 0, above = 1) #Minor third (music theory) of A 440
#'
#' @export

freq.bands <- function (x, interval=2, below = 3, above = 3){

  y <- numeric(below)
  z <- numeric(above)

  if(below > 0) for(i in 1:below) y[i] <- x/(interval^i)

  if(above > 0) for(i in 1:above) z[i] <- x * interval^i

  res <- c(rev(y), x, z)

  return(res)
}



