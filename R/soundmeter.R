#' Function that makes sound meter alike measurements
#'
#' @param files The audiofile to be analyzed. Can be "wd" to get all ".wav" files on the work directory, a file name (or a character containing a list of filenames) that exist in the work directory (only ".wav" files accepted), or an Wave object (or a list containing more than one Wave object). (By default: "wd")
#' @param from Numeric. The start time in seconds of the sample you want to analyze. Could also be relative to the end of the file (in negative values), see examples.
#' @param to Numeric. The end time in seconds of the sample you want to analyze. Could also be relative to the end of the file (in negative values), see examples.
#' @param CalibPosition anda de mãos dadas com calib value. Pode ser negativo, positivo ou data.frame com essas combinações
#' @param CalibValue Anda de mãos dadas com calib position. Quando tem o position, ele é considerado o valor de referência, quando não tem o position, ele é considerado o valor de calibração.
#' @param ref Numerical. The reference value for dB conversion. For sound in water, the common is 1 microPa, and for sound on air 20 microPa. (By default 20)
#' @param fw Character. Argument passed to \code{\link[seewave]{dBweight}} to indicate the frequency weighting curve to use on the anlysis. 'A', 'B', 'C', 'D', 'ITU', and 'none' are supported. See \code{\link[seewave]{dBweight}} for details. (By default: "none")
#' @param bands Character. Choose the type of frequency band of the output. "octaves" to octaves bands intervals or "thirds" to one-third octaves bands intervals. (by deafault: "thirds")
#' @param tw Time weighting
#' @param progressbar Logical. Activate or deactivate a progress bar with elapsed time and the last concluded file number. (By default: \code{TRUE})
#' @param channel Only "left" or "right" acepted. By default "left"
#' @param saveresults Logical. Set \code{TRUE} if you want to save a txt file with the results of the function execution. (By default: \code{FALSE})
#' @param outname Character. If \code{saveresults} is \code{TRUE}, you can specify a name to appear on the txt file name after the default name. (By default: \code{NULL})
#'
#' @details If your reference signal is in a separate file, we recommend you to get the \code{CalibValue} with \link{timbre} function. It's examples provide more details.
#'
#' @examples
#' data("tham")
#' soundmeter(tham, CalibValue = 130.24, tw="slow") #slow time window with calib value
#' soundmeter(tham, CalibValue = 130.24, tw="fast") #fast
#'
#' soundmeter(tham, CalibValue = 130.24, tw="fast", fw="A") #fast with frequency weight
#'
#' soundmeter(tham, CalibValue = NULL, tw="fast", ref=1) #fast time window in dBFS
#'
#' @export

soundmeter <- function(files="wd", from=0, to=Inf, CalibPosition=NULL, CalibValue=NULL, ref=20, fw="none", bands="octaves", tw="fast", progressbar=T, channel="left", saveresults=F, outname=NULL){

  if(class(files)=="Wave"){
    files<-list(files)
  }else if(length(files)==1 && files=="wd") {
    files <- dir(pattern=".WAV", ignore.case=T)
  }else if(is.data.frame(files)){
    files <- as.character(files)
  }

  pb <- progress_bar$new(format = "[:bar]:percent [:elapsedfull || File :current/:total done]"
                           , total = length(files)
                           , complete = "="   # Completion bar character
                           , incomplete = "-" # Incomplete bar character
                           , current = ">"    # Current bar character
                           , clear = FALSE    # If TRUE, clears the bar when finish
                           #, width = 100     # Width of the progress bar
  )

  #organizando a identificação de início e fim do trecho a analizar ####
  from=c(matrix(from, nrow=length(files)))
  to=c(matrix(to, nrow=length(files)))

  if(!(channel %in% c("left", "right")))stop("Only 'left' or 'right' acepted fo channel argument", call. = F)

  if(!is.null(CalibValue) & !is.data.frame(CalibValue)){
    CalibValue=matrix(CalibValue, nrow=length(files), ncol=1, byrow=T)
  }else if(!is.null(CalibValue) & is.data.frame(CalibValue) && nrow(CalibValue) != length(files)) stop("When CalibValue is a data.frame, it must have the number of rows equal to files length.",call. = F)


  #início do loop maior (por arquivo) ----
  for(i in 1:length(files)){

    if(!is.null(CalibPosition) && all(CalibPosition < 0)){ #ajustando calibposition

      if(is.character(files[[i]])) dur=readWave(files[[i]], header = T) %>% #se for um arquivo para ler
          data.frame() %>%
          transmute(dur=samples/sample.rate) %>%
          as.numeric()

      if(class(files[[i]]) == "Wave") dur=duration(files[[i]]) #se for um arquivo já carregado no R

      calib.ini=dur+CalibPosition[1]
      calib.fin=dur+CalibPosition[2]
    }else if(!is.null(CalibPosition)){

      calib.ini=CalibPosition[1]
      calib.fin=CalibPosition[2]

    }

    #calibrando ----
    if(exists("calib.ini") && exists("calib.fin")){
      if(class(files[[i]])=="Wave"){
        som<-extractWave(files[[i]], from=calib.ini, to=calib.fin, xunit = "time", interact = F)
      } else {
        som<-readWave(files[[i]], from = calib.ini, to=calib.fin, units = "seconds")
      }

      CalibValue[i]=timbre(som, channel=channel, Leq.calib=CalibValue[i], weighting=fw, ref=ref, stat.mess = F, time.mess = F)$Calib.value

      rm(som)
    }


    #Reading sound file ####
    if(from[i]<0 && to[i]<0){ #ajustando from & to

      if(is.character(files[[i]])) dur=readWave(files[[i]], header = T) %>% #se for um arquivo para ler
          data.frame() %>%
          transmute(dur=samples/sample.rate) %>%
          as.numeric()

      if(class(files[[i]]) == "Wave") dur=duration(files[[i]]) #se for um arquivo já carregado no R

      from[i]=dur+from[i]
      to[i]=dur+to[i]
    }

    if(class(files[[i]])=="Wave"){
      som<-extractWave(files[[i]], from=from[i], to=to[i], xunit = "time", interact = F)
    } else {
      som<-readWave(files[[i]], from = from[i], to=to[i], units = "seconds")
    }

    #Leq & medidas estatísticas
    if(!is.null(CalibValue)){
      matriz=Tweighting(som, channel=channel, window=tw, bands=bands, weighting=fw, ref=ref, Calib.value=CalibValue[i])
    } else{
      matriz=Tweighting(som, channel=channel, window=tw, bands=bands, weighting=fw, ref=ref)
    }

    #criando e armazenando valores na matriz de resultados ----
    if(i == 1){
      res = data.frame(matrix(data=NA,nrow=length(files),ncol=7+ncol(matriz[-1])))
      colnames(res) = c("File", "min", "max", "90", "50", "10", "eq", colnames(matriz[-1]))

      if(fw == "none"){
        colnames(res)[2:7]=paste0("L", colnames(res)[2:7])
      }else{
        colnames(res)[2:7]=paste0("L", fw, colnames(res)[2:7])
      }

      if(is.list(files) & !is.null(names(files))){
        res$File = names(files)
      }else if(!is.list(files)){
        res$File=files
      }else {
        res$File=1:length(files)
      }
    }


    #AJUSTAS REFERENCIAS ABAIXO!
    res[i,2]=min(matriz$Leq) #Lmin
    res[i,3]=max(matriz$Leq) #Lmax

    res[i,4:6] = LineartodB(quantile(dBtoLinear(matriz$Leq, factor = "SPL", ref = ref), probs=c(0.1, 0.5, 0.9)), factor = "SPL", ref = ref) #L90,L50 e L10

    res[i,7:ncol(res)]=timbre(som, channel=channel, bands=bands, weighting=fw, ref=ref, stat.mess = F, time.mess = F, Calib.value=ifelse(is.null(CalibValue), 0, CalibValue[i]))[,-1] #Leq e bandas

    res[i,-1]=round(res[i,-1],2) #arredondando valores para duas casas decimais


    if(saveresults) { #Salvando a matriz a cada audio analisado ----
      write.table(res,
                  paste("soundmeterResult_", fw, "-weighting_",
                        tw,
                        ifelse(!is.null(outname),paste("_", outname, sep=""),paste("")),
                        ".txt", sep="")
                  ,row.names = F, col.names = T,sep = "\t", quote=F)
    }

    rm(som)

    if(progressbar) pb$tick()

  } #final do loop maior ####

  return(res)

  pb$terminate()

}

