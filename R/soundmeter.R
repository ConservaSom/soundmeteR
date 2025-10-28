#' Sound meter measurements
#'
#' @description Funtion that performs sound meter measurements.
#'
#' @param files Specifies audio file(s) to be analyzed. It can be set to "wd" to select all ".wav" files in the work directory, a single file name, a character vector with multiple file names, a Wave object, or a list of Wave objects (by default "wd"). Only ".wav" files are accepted.
#' @param channel Character. Choose “left” or “right” channel. Argument passed to \link[tuneR]{mono} function from \link[tuneR]{tuneR} (by default "left").
#' @param from Numeric. Specifies the start time (in seconds) of the segment to analyze. It can also be relative to the end of the file (in negative values). See examples for details.
#' @param to Numeric. Specifies the end time (in seconds) of the segment to analyze. It can also be relative to the beginning of the file (in negative values). See examples for more details.
#' @param CalibPosition Numeric. Specifies the calibration position. It can be a negative (relative to the sound file duration) or a positive value, or a data.frame containing these combinations. This parameter is used in conjunction with \code{CalibValue}.
#' @param CalibValue Numeric. Specifies the calibration value (by default NULL). If "CalibPosition" is provided, it serves as the reference value. If "CalibPosition" is absent, "CalibValue" is used as the calibration reference.
#' @param ref Numeric. The reference value for dB conversion. For sound in water, the common reference is 1 µPa, and for sound in air, it is 20 µPa (by default 20).
#' @param fw Character. Defines the weighting curve for the analysis, passed to the dBweight function. Accepted values are "A", "B", "C", "D", "ITU", and "none" (by default "none"). See  \code{\link[seewave]{dBweight}} for details.
#' @param bands Character. Use "octaves" for octave bands or "thirds" for one-third octave bands intervals (by default "thirds").
#' @param tw Character. Specifies the time window to be used. Accepted values are "fast" and "slow" (default "fast").
#' @param saveresults Logical. Set \code{TRUE} to save the results to a .txt file (by default \code{FALSE}).
#' @param outname Character. If \code{saveresults} is \code{TRUE}, specifies a name to append to the file name in the txt file (by default \code{NULL}).
#' @param bandpass Numeric. A vector of length two specifying the lower and upper limits of the bandpass interval in Hz.
#' @param progressbar Logical. Set to \code{TRUE}) to display a progress bar showing elapsed time and the last completed file number, or \code{FALSE} to hide it (default is \code{TRUE}).
#'
#' @details If your reference signal is in a separate file, we recommend obtaining the \code{CalibValue} using the \link{leqbands} function. Examples in the documentation provide further details.
#'
#' @examples
#' data("tham")
#' soundmeter(tham, CalibValue = 130.24, tw = "slow") # slow time window with calib value
#' soundmeter(tham, CalibValue = 130.24, tw = "fast") # fast
#'
#' soundmeter(tham, CalibValue = 130.24, tw = "fast", fw = "A") # fast with frequency weight
#'
#' soundmeter(tham, CalibValue = NULL, tw = "fast", ref = 1) # fast time window in dBFS
#'
#' @export

soundmeter <- function(
    files = "wd",
    channel = "left",
    from = 0,
    to = Inf,
    CalibPosition = NULL,
    CalibValue = NULL,
    ref = 20,
    fw = "none",
    bands = "octaves",
    tw = "fast",
    saveresults = F,
    outname = NULL,
    bandpass = c(0, Inf),
    progressbar = T) {
  if (class(files) == "Wave") {
    files <- list(files)
  } else if (length(files) == 1 && files == "wd") {
    files <- dir(pattern = ".WAV", ignore.case = T)
  } else if (is.data.frame(files)) {
    files <- as.character(files)
  }

  pb <- progress_bar$new(
    format = "[:bar]:percent [:elapsedfull || File :current/:total done]",
    total = length(files),
    complete = "=" # Completion bar character
    , incomplete = "-" # Incomplete bar character
    , current = ">" # Current bar character
    , clear = FALSE # If TRUE, clears the bar when finish
    # , width = 100     # Width of the progress bar
  )

  # organizando a identificação de início e fim do trecho a analizar ####
  from <- c(matrix(from, nrow = length(files)))
  to <- c(matrix(to, nrow = length(files)))

  if (!(channel %in% c("left", "right"))) {
    stop("Only 'left' or 'right' acepted fo channel argument", call. = F)
  }

  if (!is.null(CalibValue) & !is.data.frame(CalibValue)) {
    CalibValue <- matrix(CalibValue, nrow = length(files), ncol = 1, byrow = T)
  } else if (!is.null(CalibValue) & is.data.frame(CalibValue) &&
    nrow(CalibValue) != length(files)) {
    stop("When CalibValue is a data.frame, it must have the number of rows equal to files length.", call. = F)
  }

  if (!is.null(bandpass) & !is.data.frame(bandpass)) {
    bandpass <- matrix(bandpass, nrow = length(files), ncol = 2, byrow = T)
  } else if (!is.null(bandpass) & is.data.frame(bandpass) &&
    nrow(bandpass) != length(files)) {
    stop("When bandpass is a data.frame, it must have the number of rows equal to files length.", call. = F)
  }

  # início do loop maior (por arquivo) ----
  for (i in 1:length(files)) {
    if (!is.null(CalibPosition) && all(CalibPosition < 0)) { # ajustando calibposition

      if (is.character(files[[i]])) { # se for um arquivo para ler
        dur <- readWave(files[[i]], header = T) %>%
          data.frame() %>%
          transmute(dur = samples / sample.rate) %>%
          as.numeric()
      }

      if (class(files[[i]]) == "Wave") { # se for um arquivo já carregado no R
        dur <- seewave::duration(files[[i]])
      }

      calib.ini <- dur + CalibPosition[1]
      calib.fin <- dur + CalibPosition[2]
    } else if (!is.null(CalibPosition)) {
      calib.ini <- CalibPosition[1]
      calib.fin <- CalibPosition[2]
    }

    # calibrando ----
    if (exists("calib.ini") && exists("calib.fin")) {
      if (class(files[[i]]) == "Wave") {
        som <- extractWave(files[[i]],
          from = calib.ini, to = calib.fin,
          xunit = "time", interact = F
        )
      } else {
        som <- readWave(files[[i]],
          from = calib.ini, to = calib.fin,
          units = "seconds"
        )
      }

      CalibValue[i] <- leqbands(som,
        channel = channel, Leq.calib = CalibValue[i],
        weighting = fw, ref = ref, progressbar = F
      )$Calib.value

      rm(som)
    }


    # Reading sound file ####
    if (from[i] < 0 && to[i] < 0) { # ajustando from & to

      if (is.character(files[[i]])) {
        dur <- readWave(files[[i]], header = T) %>% # se for um arquivo para ler
          data.frame() %>%
          transmute(dur = samples / sample.rate) %>%
          as.numeric()
      }

      if (class(files[[i]]) == "Wave"){
        dur <- seewave::duration(files[[i]]) # se for um arquivo já carregado no R
      }

      from[i] <- dur + from[i]
      to[i] <- dur + to[i]
    }

    if (class(files[[i]]) == "Wave") {
      som <- extractWave(files[[i]], from = from[i], to = to[i], xunit = "time", interact = F)
    } else {
      som <- readWave(files[[i]], from = from[i], to = to[i], units = "seconds")
    }

    # Leq & medidas estatísticas
    if (!is.null(CalibValue)) {
      matriz <- Tweighting(
        som,
        channel = channel,
        window = tw,
        bands = bands,
        weighting = fw,
        ref = ref,
        bandpass = c(bandpass[i,1],bandpass[i,2]),
        Calib.value = CalibValue[i]
      )
    } else {
      matriz <- Tweighting(som,
        channel = channel,
        window = tw,
        bands = bands,
        weighting = fw,
        ref = ref,
        bandpass = c(bandpass[i,1],bandpass[i,2])
      )
    }

    # criando e armazenando valores na matriz de resultados ----
    if (i == 1) {
      res <- data.frame(matrix(data = NA, nrow = length(files), ncol = 7 + ncol(matriz[-1])))
      colnames(res) <- c("File", "min", "max", "90", "50", "10", "eq", colnames(matriz[-1]))

      if (fw == "none") {
        colnames(res)[2:7] <- paste0("L", colnames(res)[2:7])
      } else {
        colnames(res)[2:7] <- paste0("L", fw, colnames(res)[2:7])
      }

      if (is.list(files) & !is.null(names(files))) {
        res$File <- names(files)
      } else if (!is.list(files)) {
        res$File <- files
      } else {
        res$File <- 1:length(files)
      }
    }


    # AJUSTAS REFERENCIAS ABAIXO!
    res[i, 2] <- min(matriz$Leq) # Lmin
    res[i, 3] <- max(matriz$Leq) # Lmax

    res[i, 4:6] <- LineartodB(
      quantile(
        dBtoLinear(
          matriz$Leq,
          factor = "SPL",
          ref = ref
        ),
        probs = c(0.1, 0.5, 0.9)
      ),
      factor = "SPL",
      ref = ref
    ) # L90,L50 e L10

    if (all(bandpass[i,] == c(0, Inf))) {
      res[i, 7:ncol(res)] <- leqbands(
        som,
        channel = channel,
        bands = bands,
        weighting = fw,
        ref = ref,
        progressbar = F,
        Calib.value = ifelse(is.null(CalibValue), 0, CalibValue[i])
      )[, -1] # Leq e bandas
    } else {
      espec <- pwrspec(
        som,
        channel = channel,
        bandpass = c(bandpass[i,1],bandpass[i,2]),
        res.scale = "dB",
        ref = ref
      )

      # Implementando curvas de ponderacao ----
      if (any(fw == c("A", "B", "C", "D", "ITU"))) {
        espec$Amp.dB <- dBweight(espec$Freq.Hz, dBref = espec$Amp.dB)[[fw]]
      } else if (fw != "none") {
        stop("Wrong weighting curve. Only 'A', 'B', 'C', 'D', 'ITU', and 'none' accepted. See dBweight()' for details.")
      }

      espec <- espec %>%
        select(Amp.dB) %>%
        sumdB() %>%
        round(2)

      if (!is.null(!is.null(CalibValue))) {
        espec <- round(espec + CalibValue[i], 2)
      }

      res[i, 7] <- espec
    }


    res[i, -1] <- round(res[i, -1], 2) # arredondando valores para duas casas decimais


    if (saveresults) { # Salvando a matriz a cada audio analisado ----
      write.table(res,
        paste("soundmeterResult_", fw, "-weighting_",
          tw,
          ifelse(!is.null(outname), paste("_", outname, sep = ""), paste("")),
          ".txt",
          sep = ""
        ),
        row.names = F, col.names = T, sep = "\t", quote = F
      )
    }

    rm(som)

    if (progressbar) pb$tick()
  } # final do loop maior ####

  return(res)

  pb$terminate()
}
