#' RMS from a sample of a sound file
#'
#' @description RMS (Root Mean Square) from a sample of a sound file.
#' 
#'
#' @param files Specifies audio file(s) to be analyzed. It can be set to "wd" to select all ".wav" files in the work directory, a single file name, a character vector with multiple file names, a Wave object, or a list of Wave objects (by default "wd"). Only ".wav" files are accepted.
#' @param channel Character. Choose “left” or “right” channel. Argument passed to \link[tuneR]{mono} function from \link[tuneR]{tuneR} (by default "left").
#' @param from Numeric. Specifies the start time (in seconds) of the segment to analyze. It can also be relative to the end of the file (in negative values). See examples for details.
#' @param to Numeric. Specifies the end time (in seconds) of the segment to analyze. It can also be relative to the beginning of the file (in negative values). See examples for more details.
#' @param freq.interval Specifies the frequency interval for computing the RMS. It can be a vector of length two, with the lower and upper frequency bounds (in Hz), or a pattern to calculate a interval (e.g., octaves). For more details, refer to \link{freq.bands}.
#' @param fdom.int A vector of length two specifying the lower and upper frequency bounds (in Hz) to determine the dominant frequency. This frequency will be used as the center of the interval only if a pattern is specified in \code{freq.interval}.
#' @param wl Numeric. Window length for the analysis. It must be a even number of points (by default 512).
#' @param ovlp Numeric. Overlap percentage between two successive windows (by default 50). Argument passed to \link[seewave]{meanspec} function from \link[seewave]{seewave}.
#' @param CalibPosition Numeric. Specifies the calibration position. It can be a negative (relative to the sound file duration) or a positive value, or a data.frame containing these combinations. This parameter is used in conjunction with \code{CalibValue}.
#' @param CalibValue Numeric. Canbe a value to apply of the ref value from a calib signal (specified by CalibPosition).
#' @param freq.weight Character. Defines the weighting curve for the analysis, passed to the \link[seewave]{dBweight} function. Accepted values are "A", "B", "C", "D", "ITU", and "none" (by default "none"). See \link[seewave]{dBweight} for details.
#' @param ref Numeric. The reference value for dB conversion. For sound in water, the common reference is 1 µPa, and for sound in air,it is 20 µPa (by default 20).
#' @param progressbar Logical. Set to \code{TRUE} to display a progress bar showing elapsed time and the last completed file number, or \code{FALSE} to hide it (default \code{TRUE}).
#'
#'
#' @examples
#' song.level(tham, freq.interval=c(22, 20000))
#'
#' song.level(tham, freq.interval=c(22, 20000), CalibValue = 130.24)
#'
#' song.level(tham, freq.interval=c(22, 20000), CalibValue = 130.24, freq.weight="A")
#'
#' song.level(tham, freq.interval=c(22, 20000), CalibValue = 130.24, freq.weight="B")
#'
#' song.level(tham, freq.interval=c(22, 20000), CalibValue = 130.24, freq.weight="C")
#'
#' song.level(tham, from = 3.883035, to=7.044417, freq.interval=c(22.09429, 22627.38), CalibValue = 130.24, freq.weight="A")
#'
#' #Perfect fifth (music theory)
#' song.level(tham, fdom.int = c(800, 2000), from = 3.8, to=7, freq.interval=3/2, CalibValue = 130.24, freq.weight="A")
#'
#' #Perfect fifth (music theory)
#' song.level(tham, fdom.int = c(800, 2000), from = 3.8, to=7, freq.interval=3/2, CalibValue = 130.24, freq.weight="A")
#'
#' @export

song.level<-function(files="wd", channel="left", from=0, to=Inf,
                     freq.interval=c(0, Inf), fdom.int=c(0,Inf), wl=512,
                     ovlp = 50, CalibPosition=NULL, CalibValue=NULL,
                     freq.weight="none", ref=20, progressbar=T){

  if(class(files)=="Wave"){
    files<-list(files)
  }else if(length(files)==1 && files=="wd") {
    files <- dir(pattern=".WAV", ignore.case=T)
  }else if(is.data.frame(files)){
    files <- as.character(files)
  }

  from=c(matrix(from, nrow=length(files)))
  to=c(matrix(to, nrow=length(files)))

  if(!(channel %in% c("left", "right")))stop("Only 'left' or 'right' acepted fo channel argument", call. = F)

  if(!is.null(CalibValue) & !is.data.frame(CalibValue)){
    CalibValue=matrix(CalibValue, nrow=length(files), ncol=1, byrow=T)
  }else if(!is.null(CalibValue) & is.data.frame(CalibValue) && nrow(CalibValue) != length(files)) stop("When CalibValue is a data.frame, it must have the number of rows equal to files length.",call. = F)

  if(!is.data.frame(fdom.int)){
    fdom.int=matrix(fdom.int, nrow=length(files), ncol=2, byrow=T)
  }else if(is.data.frame(fdom.int) & nrow(fdom.int) != length(files)) stop("When fdom.int is a data.frame, it must have the number of rows equal to files length.",call. = F)

  pb <- progress_bar$new(format = "[:bar]:percent [:elapsedfull || File :current/:total done]"
                         , total = length(files)
                         , complete = "="   # Completion bar character
                         , incomplete = "-" # Incomplete bar character
                         , current = ">"    # Current bar character
                         , clear = FALSE    # If TRUE, clears the bar when finish
                         #, width = 100     # Width of the progress bar
  )

  #tabela dos resultados ----
  res=data.frame(File=1:length(files)
                 , Freq.interval=NA
                 , Freq.dom=NA
                 , SongLevel=NA
  )

  if(is.list(files) & !is.null(names(files))){
    res$File = names(files)
  }else if(!is.list(files)){
    res$File=files
  }


  if(length(freq.interval) != 1){
    res=select(res, -Freq.dom)
  }else {
    colnames(res)[2] = paste0("Freq.interval_", freq.interval)
  }


  #Analisando os arquivos ----
  for(i in 1:length(files)){

    if(channel == "right" &&
       !is.list(files) &&
       readWave(files[[i]], header = T)$channels < 2){
      warning(paste0("File ", files[[i]], " doesn't have a right channel."),
              call. = F)
      next
    }

    #Calibração ----
    if(!is.null(CalibPosition) & !is.null(CalibValue)){
        CalibValue[i]=leqbands(files[[i]], channel = channel,
                             from=CalibPosition[1], to=CalibPosition[2],
                             Leq.calib=CalibValue[i], ref=ref,
                             weighting=freq.weight, progressbar=F)$Calib.value
    }


    #localizando Frequencia dominante ----
    if(length(freq.interval) == 1){
      freq.dom=meanspec(readWave(files[[i]], from=from[i], to=to[i],
                                 units = "seconds"),
                        channel = ifelse(channel=="left", 1, 2),
                        wl = wl, ovlp = ovlp, plot = F,
                        dB = ifelse(freq.weight=="none", "max0", freq.weight)) %>% # power spectrum
        .[-1,] %>%
        .[.[,"x"] >= min(fdom.int[i,]/1000) & .[,"x"] <= max(fdom.int[i,]/1000),] %>% #filtro de band pass
        fpeaks(nmax=1, plot = F) #qual o pico?

      if(any(is.na(freq.dom))){
        warning(paste0("File ", files[[i]], "(round ", i, ") doesn't have a peak to find a Dominant Frequency"), call. = F)
        next
      }

      freq.dom=freq.dom %>%
        .[,1]*1000

    }

    #Intervalo para somar ----
    if(length(freq.interval) == 1){
      interval.tosum=freq.bands(freq.dom, interval = freq.interval,
                                below = 1, above = 1) %>%     #intervalo ao redor da dominante ----
      range()
    }else {
      interval.tosum = freq.interval #intervalo fixo, estabelecido pelo usuário ----
    }

    #pwerspec do arquivo ----
    espec=pwrspec(files[[i]], channel = channel, from=from[i], to=to[i],
                  res.scale = "dB", ref=ref)

    #Implementando curvas de ponderacao ----
    if(any(freq.weight == c("A", "B", "C", "D", "ITU"))){
      espec$Amp.dB=dBweight(espec$Freq.Hz, dBref = espec$Amp.dB)[[freq.weight]]
    } else if(freq.weight != "none"){stop("Wrong weighting curve. Only 'A', 'B', 'C', 'D', 'ITU', and 'none' accepted. See dBweight()' for details.")}

    #Energia na banda desejada ----
    espec=espec %>%
      filter(Freq.Hz >= interval.tosum[1] & Freq.Hz < interval.tosum[2])

    song.level=round(sumdB(espec[,"Amp.dB"], level="IL"), 2)

    #Calibrando ----

    if(!is.null(CalibValue)) song.level=song.level+CalibValue[i]


    res[i, "SongLevel"] = song.level
    if(exists("freq.dom")) res[i, "Freq.dom"] = freq.dom
    res[i, 2]=paste0(round(interval.tosum,0), collapse="—") #Freq.interval

    if(progressbar) pb$tick()

  }

  return(res)
  pb$terminate()

}
