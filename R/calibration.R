#' Extract calbration factor from reference sinal
#'
#' @name calibration
#'
#' @description ESCREVER
#'
#' @param file Specifies audio file to be analyzed. It can be a single file name on the computer or a Wave object on the environment (by default NULL). Only ".wav" files are accepted.
#' @param channel Character. Choose “left” or “right” channel (by default "left"). Argument passed to \link[tuneR]{mono} function from \link[tuneR]{tuneR} to extract the desired channel.
#' @param from Numeric. The start time (in seconds) of the segment to analyze. It can also be relative to the end of the file (in negative values). See examples for details.
#' @param to Numeric. The end time (in seconds) of the segment to analyze. Could also be relative to the beginning of the file (in negative values). See examples for details.
#' @param weighting Character. Defines the weighting curve for the analysis, passed to the \code{\link[seewave]{dBweight}} function. Accepted values are "A", "B", "C", "D", "ITU", and "none" (by default "none"). See \code{\link[seewave]{dBweight}} for details.
#' #' @param leqsignal Numeric. Specifies the sound pressure level (in dB SPL) for the signal in the audio (by default \code{NA}).
#' @param ref Numeric. The reference value for dB conversion. For sound in water, the common reference is 1 µPa, and for sound in air,it is 20 µPa (by default 20).
#'
#' @export
#' 

calibration <- function(
    file = NULL,
    channel = "left",
    from = 0,
    to = Inf,
    leqsignal = NA,
    weighting = "none",
    ref = 20) {
    if (is.na(leqsignal) | is.null(leqsignal)) {
        stop("leqsignal must be a positive value")
    }

    if(files == "wd"){
      stop("When using calibration function, file must be a single file name or a Wave object.", call. = F)
    }

    calibvalue <- leqbands(
        file,
        channel = channel,
        from = from,
        to = to,
        weighting = weighting,
        ref = ref,
        outname = NULL,
        saveresults = FALSE,
        progressbar = FALSE
    )$Leq %>% # extrai leq em dBFS
        {
            .[] - leqsignal # extrai o valor de calibração para dB absoluto
        } %>%
        abs() # valor como deve ser usado nas demais funções do pacote

    return(calibvalue)
}
