#' Timbre analysis for audiofiles with reference signal
#'
#' @name timbreCal
#'
#' @description This function passes the parameters to \code{\link{timbre}} to automatize the calibration and return spectral analysis with dB SPL results.
#'
#' @param files The audiofile to be analyzed. Can be "wd" to get all ".wav" files on the work directory, a file name (or a character containing a list of filenames) that exist in the work directory (only ".wav" files accepted), or an Wave object (or a list containing more than one Wave object). (By default: "wd")
#' @param channel Argument passed to \link[tuneR]{mono} function from \link[tuneR]{tuneR} to extract the desired channel.
#' @param from Numeric. The start time in seconds of the sample you want to analyze. Could also be relative to the end of the file (in negative values), see examples.
#' @param to Numeric. The end time in seconds of the sample you want to analyze. Could also be relative to the end of the file (in negative values), see examples.
#' @param CalibPosition anda de mãos dadas com calib value. Pode ser negativo, positivo ou data.frame com essas combinações
#' @param CalibValue Anda de mãos dadas com calib position. Quando tem o position, ele é considerado o valor de referência, quando não tem o position, ele é considerado o valor de calibração.
#' @param ref Numerical. The reference value for dB conversion. For sound in water, the common is 1 microPa, and for sound on air 20 microPa. (By default 20)
#' @param weighting Character. Indicate the weighting curve to use on the anlysis. A, B, C and none are supported. (By default: "none")
#' @param bands Character. Choose the type of frequency band of the output. "octaves" to octaves bands intervals or "thirds" to one-third octaves bands intervals. (by deafault: "thirds")
#' @param saveresults Logical. Set \code{TRUE} if you want to save a txt file with the results of the function execution. (By default: \code{FALSE})
#' @param outname Character. If \code{saveresults} is \code{TRUE}, you can specify a name to appear on the txt file name after the default name. (By defaulf: \code{NULL})
#' @param progressbar Logical. Activate or deactivate a progress bar with elapsed time and the last concluded file number. (By default: \code{TRUE})
#'
#' @details   To use this function, the audio file must begin with 2 seconds of silence, followed by a reference signal with known SPL, followed by another 2 seconds of silence, and the following sound to analyze.
#' @details   The duration of the reference signal must be specified (in seconds) on the \code{SignalDur} argument and his value (in dB SPL) on the \code{refValue} argument.
#'
#' @seealso \code{\link{timbre}}
#'
#'
#' @export

#Coisas para fazer:
#Pensar ao invés de usar um trecho da gravação para calibrar usar um arquivo externo.


timbreCal <- function(files="wd", channel="left", from=0, to=Inf,
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
      if(class(files[[i]])=="Wave"){
        som<-extractWave(files[[i]], from=calib.ini, to=calib.fin,
                         xunit = "time", interact = F)
      } else {
        som<-readWave(files[[i]], from = calib.ini, to=calib.fin,
                      units = "seconds")
      }

      CalibValue[i]=timbre(som, channel=channel, Leq.calib=CalibValue[i],
                           weighting=weighting, ref=ref, progressbar=F)$Calib.value

      rm(som)
    }

    if(i==1){#gerando a matriz de resultados ####
      results<-timbre(files=files[[i]], from=from[i], to=to[i],
                      channel=channel, Calib.value=CalibValue[i], ref=ref,
                      weighting=weighting, bands=bands, progressbar=F)

    }else {
      results<-rbind(results,
                     timbre(files=files[[i]], from=from[i], to=to[i],
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
                  paste("TimbreCalResults_", weighting, "-weighting",
                        ifelse(!is.null(outname),paste("_", outname, sep=""),
                               paste("")), ".txt", sep="")
                  ,row.names = F, col.names = T,sep = "\t", quote=F)
    }

    if(progressbar) pb$tick()

  }

  return(results)

}
