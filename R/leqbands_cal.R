#' leqbands analysis for audiofiles with reference signal
#'
#' @description This function passes the parameters to \code{\link{leqbands}}to automate the calibration process and return spectral analysis results in dB SPL.
#'
#' @param files Specifies audio file(s) to be analyzed. It can be set to "wd" to select all ".wav" files in the work directory, a single file name, a character vector with multiple file names, a Wave object, or a list of Wave objects (by default "wd"). Only ".wav" files are accepted.
#' @param channel Character. Choose “left” or “right” channel. Argument passed to \link[tuneR]{mono} function from \link[tuneR]{tuneR} (by default "left").
#' @param from Numeric. The start time (in seconds) of the segment to analyze. It can also be relative to the end of the file (in negative values). See examples for details.
#' @param to Numeric. The end time (in seconds) of the segment to analyze. It can also be relative to the beginning of the file (in negative values). See examples for details.
#' @param CalibPosition Numeric. Specifies the calibration position. It can be a negative (relative to the sound file duration) or a positive value, or a data.frame containing these combinations. This parameter is used in conjunction with \code{CalibValue}.
#' @param CalibValue Numeric. Specifies the calibration value. If \code{CalibPosition} is provided, it serves as the reference value. If \code{CalibPosition} is absent, \code{CalibValue} is used as the calibration reference.
#' @param ref Numeric. The reference value for dB conversion. For sound in water, the common reference is 1 µPa, and for sound in air,it is 20 µPa (by default 20).
#' @param weighting Character. Defines the weighting curve for the analysis, passed to the \link[seewave]{dBweight} function. Accepted values are "A", "B", "C", "D", "ITU", and "none" (by default "none"). See \link[seewave]{dBweight} for details.
#' @param bands Character. Use "octaves" for octave bands or "thirds" for one-third octave bands intervals (by default "thirds").
#' @param saveresults Logical. Set to \code{TRUE} to save the results to a .txt file (by default \code{FALSE}).
#' @param outname Character. If \code{saveresults} is \code{TRUE}, specifies a name to append to the file name in the txt file (by default \code{NULL}).
#' @param progressbar Logical. Set \code{TRUE} to display a progress bar showing elapsed time and the last completed file number, or \code{FALSE} to hide it (default \code{TRUE}).
#'
#' @details   To use this function, the audio file must follow this structure: 2 seconds of silence, folowed by a reference signal with a known SPL; then another 2 seconds of silence, followed by the sound to be analyzed.
#' @details   The duration of the reference signal must be specified (in seconds) on the \code{SignalDur} argument and his value (in dB SPL) on the \code{refValue} argument.
#'
#' @return ESCREVER
#'
#' @seealso \code{\link{leqbands}}
#'
#' @export

#Coisas para fazer:
#Pensar ao invés de usar um trecho da gravação para calibrar usar um arquivo externo.


leqbands_cal <- function(files="wd", channel="left", from=0, to=Inf,
                     CalibPosition=NULL, CalibValue=NULL, ref=20,
                     weighting="none", bands="thirds", saveresults=F,
                     outname=NULL, progressbar=T){

  if(class(files)=="Wave"){
    files<-list(files)
  }else if(length(files)==1 && files=="wd") {
    files <- dir(pattern=".WAV", ignore.case=T)
  }else if(is.data.frame(files)){
    files <- as.character(files)
  }


  if(progressbar){
    pb <- progress_bar$new(format = "[:bar]:percent [:elapsedfull || File :current/:total done]"
                           , total = length(files)
                           , complete = "="   # Completion bar character
                           , incomplete = "-" # Incomplete bar character
                           , current = ">"    # Current bar character
                           , clear = FALSE    # If TRUE, clears the bar when finish
                           #, width = 100     # Width of the progress bar
    )
  }

  #organizando a identificação de início e fim do trecho a analizar ----
  from=c(matrix(from, nrow=length(files)))
  to=c(matrix(to, nrow=length(files)))

  if(!(channel %in% c("left", "right"))){
    stop("Only 'left' or 'right' acepted fo channel argument", call. = F)
  }

  if(!is.null(CalibValue) & !is.data.frame(CalibValue)){
    CalibValue=matrix(CalibValue, nrow=length(files), ncol=1, byrow=T)
  }else if(!is.null(CalibValue) & is.data.frame(CalibValue) &&
           nrow(CalibValue) != length(files)){
    stop("When CalibValue is a data.frame, it must have the number of rows equal to files length.",call. = F)
  }

  #Loop que analisara os files ----
  for(i in 1:length(files)){

    if(!is.null(CalibPosition) && all(CalibPosition < 0)){ #ajustando calibposition

      if(is.character(files[[i]])){ #se for um arquivo para ler
        dur=readWave(files[[i]], header = T) %>%
          data.frame() %>%
          transmute(dur=samples/sample.rate) %>%
          as.numeric()
      }

      if(class(files[[i]]) == "Wave"){ #se for um arquivo já carregado no R
        dur=duration(files[[i]])
      }

      calib.ini=dur+CalibPosition[1]
      calib.fin=dur+CalibPosition[2]

    }else if(!is.null(CalibPosition)){

      calib.ini=CalibPosition[1]
      calib.fin=CalibPosition[2]

    }

    #calibrando ----
    if(exists("calib.ini") && exists("calib.fin")){

      CalibValue[i]=leqbands(files[[i]], channel=channel, from=calib.ini,
                           to=calib.fin, Leq.calib=CalibValue[i],
                           weighting=weighting, ref=ref,
                           progressbar=F)$Calib.value
    }

    if(i==1){#gerando a matriz de resultados ####
      results<-leqbands(files=files[[i]], from=from[i], to=to[i],
                      channel=channel, Calib.value=CalibValue[i], ref=ref,
                      weighting=weighting, bands=bands, progressbar=F)

    }else {
      results<-rbind(results,
                     leqbands(files=files[[i]], from=from[i], to=to[i],
                            channel=channel, Calib.value=CalibValue[i],
                            ref=ref, weighting=weighting, bands=bands,
                            progressbar=F)
      )
    }

    if(is.list(files) && is.null(names(files))){ # colocando o nome dos files na primeira coluna####
      results[i,1]<-i
    }else if(is.list(files) && !is.null(names(files))){
      results[i,1]<-names(files[i])
    }else {
      results[i,1]<-files[i]
    }

    if(saveresults) { #Salvando a matriz a cada audio analisado ####
      write.table(results,
                  paste("leqbands_calResults_", weighting, "-weighting",
                        ifelse(!is.null(outname),paste("_", outname, sep=""),
                               paste("")), ".txt", sep="")
                  ,row.names = F, col.names = T,sep = "\t", quote=F)
    }

    if(progressbar) pb$tick()

  }

  return(results)

}
