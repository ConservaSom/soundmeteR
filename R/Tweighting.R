#' Time Weighting of a Audiofile
#'
#' @description Escrever descrição
#'
#' @param file Wave object
#' @param window Character. Wich time window should be used. 'fast' or 'slow'
#'     are accepted. (by default: "fast")
#' @param ... Further arguments passed to \code{\link{leqbands}}.
#'
#' @details This function split your audiofile in smaller files defined as
#'     \code{fast} (0.125s) and \code{slow} (1s) and analyze each one with
#'     \code{\link{leqbands}} function.
#'
#' @return A numeric vector
#'
#' @seealso \code{\link{leqbands}}, \code{\link{soundmeter}}
#'
#' @examples
#' # creating an example sound file
#' som <- sine(1000, duration = 44500)
#'
#' # default options without calibration (results in dBFS)
#' Tweighting(som)
#'
#' # Simulation of a calib signal with a Leq of 94dB in the field ####
#' # fast
#' Tweighting(som, window = "fast", bands = "octaves", Leq.calib = 94)
#' # slow
#' Tweighting(som, window = "slow", bands = "octaves", Leq.calib = 94)
#'
#' # Using the result of the simulation above to calibrate the sound and output
#' # fast
#' Tweighting(som, window = "fast", bands = "octaves", Calib.value = 309.67)
#' # slow
#' Tweighting(som, window = "slow", bands = "octaves", Calib.value = 309.67)
#'
#' # With tham data
#' data(tham)
#' Tweighting(tham, window = "fast", bands = "octaves", Calib.value = 130.24) # fast
#' Tweighting(tham, window = "slow", bands = "octaves", Calib.value = 130.24) # slow
#'
#' @export


Tweighting <- function(file, window = "fast", ...) {
  if (window == "fast") {
    window <- 0.125
  } else if (window == "slow") {
    window <- 1
  } else {
    stop("Choose a valid window size in seconds ('fast' or 'slow')")
  }

  if (class(file) != "Wave") {
    stop("Only one Wave object accepted on this function")
  }

  res <- sapply(
    1:trunc(duration(file) / window),
    FUN = function(x, file, samp) {
      file %>%
        extractWave(
          from = round((x - 1) * samp),
          to = round(x * samp)
        ) %>%
        leqbands(
          progressbar = F,
          Leq.calib = NULL,
          ...
        ) %>%
        select(-Arquivo) %>%
        return()
    },
    file = file,
    samp = window * file@samp.rate
  ) %>%
    t() %>%
    as.data.frame()

  return(res)
}
