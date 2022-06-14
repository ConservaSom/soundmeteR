#' RMS from a sample of a sound file
#'
#'
#' @param files The audiofile to be analyzed. Can be "wd" to get all ".wav" files on the work directory, a file name (or a character containing a list of filenames) that exist in the work directory (only ".wav" files accepted), or an Wave object (or a list containing more than one Wave object). (By default: "wd")
#' @param from Numeric. The start time in seconds of the sample you want to analyze. Could also be relative to the end of the file (in negative values), see examples.
#' @param to Numeric. The end time in seconds of the sample you want to analyze. Could also be relative to the end of the file (in negative values), see examples.
#' @param freq.interval Frequency interval to compute the RMS. Can be a vector with length two with lower and upper interval of frequencies (in Hz), or a pattern to calculate a interval (as octaves). For the last, see \link{freq.bands} function for details.
#' @param fdom.int Vector with length two. Vector with length two with lower and upper interval of frequencies (in Hz) to find the dominant frequency. This frequency will be used as the center the interval only if a pettern is specified.
#' @param CalibPosition
#' @param CalibValue Canbe a value to apply of the ref value from a calib signal (specified by CalibPosition).
#' @param freq.weight Character. Argument passed to dBweight to indicate the weighting curve to use on the anlysis. 'A', 'B', 'C', 'D', 'ITU', and 'none' are supported. See dBweight for details. (By default: "none")
#' @param ref Numerical. The reference value for dB conversion. For sound in water, the common is 1 microPa, and for sound on air 20 microPa. (By default 20)
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

song.level<-function(files="wd", from=0, to=Inf, freq.interval=c(0, Inf), fdom.int=c(0,Inf), CalibPosition=NULL, CalibValue=NULL, freq.weight="none", ref=20){

  if(class(files)=="Wave"){
    files<-list(files)
  }else if(length(files)==1 && files=="wd") {
    files <- dir(pattern=".WAV", ignore.case=T)
  }

  if(!is.null(CalibValue) & !is.data.frame(CalibValue)){
    CalibValue=matrix(CalibValue, nrow=length(files), ncol=1, byrow=T)
  }else if(!is.null(CalibValue) & is.data.frame(CalibValue) && nrow(CalibValue) != length(files)) stop("When CalibValue is a data.frame, it must have the number of rows equal to files length.",call. = F)

  if(!is.data.frame(fdom.int)){
    fdom.int=matrix(fdom.int, nrow=length(files), ncol=2, byrow=T)
  }else if(is.data.frame(fdom.int) & nrow(fdom.int) != length(files)) stop("When fdom.int is a data.frame, it must have the number of rows equal to files length.",call. = F)


  #Analisando os arquivos ----
  for(i in 1:length(files)){

    #Calibração ----
    if(!is.null(CalibPosition) & !is.null(CalibValue)){
        CalibValue[i]=timbre(files[[i]], from=CalibPosition[1], to=CalibPosition[2], Leq.calib=CalibValue[i], ref=ref, weighting=freq.weight)$Calib.value
    }

    espec=pwrspec(files[[i]], from=from, to=to, res.scale = "dB", ref=ref)

    #localizando Frequencia dominante ----
    freq.dom=espec %>%
      filter(Freq.Hz >= fdom.int[i,1] & Freq.Hz <= fdom.int[i,2]) %>%
      slice(which.max(Amp.dB)) %>%
      select(Freq.Hz) %>%
      as.numeric()

    #Intervalo para somar ----
    if(length(freq.interval) == 1){
      interval.tosum=freq.bands(freq.dom, interval = freq.interval, below = 1, above = 1) %>%     #intervalo ao redor da dominante ----
      range()
    }else {
      interval.tosum = freq.interval #intervalo fixo, estabelecido pelo usuário ----
    }

    #Implementando curvas de ponderacao ----
    if(any(freq.weight == c("A", "B", "C", "D", "ITU"))){
      espec$Amp.dB=dBweight(espec$Freq.Hz, dBref = espec$Amp.dB)[[freq.weight]]
    } else if(freq.weight != "none"){stop("Wrong weighting curve. Only 'A', 'B', 'C', 'D', 'ITU', and 'none' accepted. See dBweight()' for details.")}

    #Energia na banda desejada ----
    espec=espec %>%
      filter(Freq.Hz >= interval.tosum[1] & Freq.Hz < interval.tosum[2])

    song.level=round(moredB(espec[,"Amp.dB"], level="IL"), 2)

    #Calibrando ----

    if(!is.null(CalibValue)) song.level=song.level+CalibValue[i]

    if(i == 1){
      res=data.frame(File=1:length(files)
                     , Freq.interval=ifelse(length(freq.interval) == 1, paste0(round(interval.tosum,0), collapse="—"), paste0(freq.interval, collapse="—"))
                     , Freq.dom=freq.dom
                     , SongLevel=NA
                     )

      if(length(freq.interval) != 1) res=select(res, -Freq.dom)
    }

    res[i, "SongLevel"] = song.level

    if(is.list(files) & !is.null(names(files))){
      res[i,"File"] = names(files[i])
    }else if(!is.list(files)){
      res[i,"File"]=files[i]
    }

  }

  return(res)

}
