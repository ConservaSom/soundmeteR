#' Time Weighting of a Audiofile
#'
#' @description This function split your audiofile into smaller segments, defined as \code{fast} (0.125s) and \code{slow} (1s), and analyze each one using the \code{\link{leqbands}} function.
#'
#' @param file Wave object
#' @param window Character. Specifies the time window to be used. Accepted values are \code{fast} and \code{slow} (default \code{fast}).
#' @param Leq.calib Numeric. Specifies the sound pressure level (in dB SPL) for the signal in the audio file (by default NULL). Cannot be set if Calib.value is specified.
#' @param bandpass Numeric. A vector of length two specifying the lower and upper limits of the bandpass interval in Hz.
#' @param bands Character. Use "octaves" for octave bands or "thirds" for one-third octave bands intervals (by default "thirds").
#' @param Calib.value Numeric. Specifies the calibration value (by default NULL). Cannot be set if Leq.calib is specified.
#' @param channel Character. Choose “left” or “right” channel. Argument passed to mono function from tuneR (by default "left").
#' @param from Numeric. Specifies the start time (in seconds) of the segment to analyze. It can also be relative to the end of the file (in negative values). See examples for details.
#' @param to Numeric. Specifies the end time (in seconds) of the segment to analyze. It can also be relative to the beginning of the file (in negative values). See examples for more details.
#' @param weighting Character. Specifies the end time (in seconds) of the segment to analyze. It can also be relative to the beginning of the file (in negative values). See examples for more details.
#' @param ref Numeric. Defines the weighting curve for the analysis, passed to the dBweight function. Accepted values are "A", "B", "C", "D", "ITU", and "none" (by default "none"). See dBweight for details.
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
