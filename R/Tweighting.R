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


Tweighting <- function(
    file,
    window = "fast",
    channel = "left",
    weighting = "none",
    bands = "thirds",
    ref = 20,
    Calib.value = NULL,
    bandpass = c(0, Inf),
    from = 0,
    to = Inf) {
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

  if (from > 0 | to < Inf) {
    file <- extractWave(file, from = from, to = to, xunit = "time", interact = F)
  }

  if (all(bandpass == c(0, Inf))) {
    res <- sapply(
      1:trunc(duration(file) / window),
      FUN = function(x,
                     file,
                     samp,
                     channel,
                     weighting,
                     bands,
                     Calib.value,
                     ref) {
        file %>%
          extractWave(
            from = round((x - 1) * samp),
            to = round(x * samp)
          ) %>%
          leqbands(
            progressbar = F,
            Leq.calib = NULL,
            channel = channel,
            weighting = weighting,
            bands = bands,
            Calib.value = Calib.value,
            ref = ref
          ) %>%
          select(-Arquivo) %>%
          return()
      },
      file = file,
      samp = window * file@samp.rate,
      channel = channel,
      weighting = weighting,
      bands = bands,
      Calib.value = Calib.value,
      ref = ref
    ) %>%
      t()

    res <- res %>%
      unlist() %>%
      matrix(
        nrow = nrow(res),
        byrow = FALSE,
        dimnames = list(NULL, colnames(res))
      ) %>%
      as.data.frame(check.names = FALSE)
  } else {
    res <- sapply(
      1:trunc(duration(file) / window),
      FUN = function(x, file, channel, samp, bandpass, ref) {
        res <- file %>%
          extractWave(
            from = round((x - 1) * samp),
            to = round(x * samp)
          ) %>%
          pwrspec(
            channel = channel,
            bandpass = bandpass,
            res.scale = "dB",
            ref = ref
          )

        # Implementando curvas de ponderacao ----
        if (any(weighting == c("A", "B", "C", "D", "ITU"))) {
          res$Amp.dB <- dBweight(res$Freq.Hz, dBref = res$Amp.dB)[[weighting]]
        } else if (weighting != "none") {
          stop("Wrong weighting curve. Only 'A', 'B', 'C', 'D', 'ITU', and 'none' accepted. See dBweight()' for details.")
        }

        res <- res %>%
          select(Amp.dB) %>%
          sumdB() %>%
          round(2) %>%
          return()
      },
      file = file,
      channel = channel,
      samp = window * file@samp.rate,
      bandpass = bandpass,
      ref = ref
    ) %>%
      data.frame(Leq = .)

    if (!is.null(Calib.value)) {
      res <- round(res + Calib.value, 2)
    }
  }

  return(res)
}
