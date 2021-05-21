#' Timbre analysis for audiofiles with reference signal
#'
#' @name timbreCal
#'
#' @description This function passes the parameters to \code{\link{timbre}} to automatize the calibration and return spectral analysis with dB SPL results.
#'
#' @usage timbreCal(files="wd", SignalDur=NULL, RefValue=NULL, ref=20, weighting="none",
#'        bands="thirds", saveresults=F, outname=NULL, time.mess=T, stat.mess=T)
#'
#' @param files The audiofile to be analyzed. Can be "wd" to get all ".wav" files on the work directory, a file name (or a character containing a list of filenames) that exist in the work directory (only ".wav" files accepted), or an Wave object (or a list containing more than one Wave object). (By default: "wd")
#' @param SignalDur Numerical. Specify the reference signal duration (in seconds) on the beggining of the audiofile. (By default: \code{NULL})
#' @param RefValue Numerical. Specify the reference signal sound pressure level (in deciBells SPL) on the beggining of the audiofile. (By default: \code{NULL})
#' @param ref Numerical. The reference value for dB conversion. For sound in water, the common is 1 microPa, and for sound on air 20 microPa. (By default 20)
#' @param weighting Character. Indicate the weighting curve to use on the anlysis. A, B, C and none are supported. (By default: "none")
#' @param bands Character. Choose the type of frequency band of the output. "octaves" to octaves bands intervals or "thirds" to one-third octaves bands intervals. (by deafault: "thirds")
#' @param saveresults Logical. Set \code{TRUE} if you want to save a txt file with the results of the function execution. (By default: \code{FALSE})
#' @param outname Character. If \code{saveresults} is \code{TRUE}, you can specify a name to appear on the txt file name after the default name. (By defaulf: \code{NULL})
#' @param time.mess Logical. Activate or deactivate message of time to complete the function execution. (By default: \code{TRUE})
#' @param stat.mess Logical. Activate or deactivate status message of the function execution. (By default: \code{TRUE})
#'
#' @author Cássio Rachid Simões <cassiorachid@@gmail.com>
#'
#' @details   To use this function, the audio file must begin with 2 seconds of silence, followed by a reference signal with known SPL, followed by another 2 seconds of silence, and the following sound to analyze.
#' @details   The duration of the reference signal must be specified (in seconds) on the \code{SignalDur} argument and his value (in dB SPL) on the \code{refValue} argument.
#'
#' @seealso \code{\link{timbre}}
#'
#'
#'

#Coisas para fazer:
#Pensar ao invés de usar um trecho da gravação para calibrar usar um arquivo externo.


timbreCal<- function(files="wd", SignalDur=NULL, RefValue=NULL, ref=20, weighting="none", bands="thirds", saveresults=F, outname=NULL, time.mess=T, stat.mess=T){

  start.time<-Sys.time()

  require(tuneR)

  if(class(files)=="Wave"){ #vendo o tipo de arquivo usado no imput e armazenando em um objeto ###
    arquivos<-list(files)
  }else if(is.list(files)){
    arquivos<-files
  }else if(length(files)==1 && files=="wd") {
    arquivos <- dir(pattern=".WAV", ignore.case=T)
  }else if(is.character(files)) {
    arquivos <- files
  }else {
    stop("Choose a valid file on the 'files' argument. Could be 'wd', a filename on your work directory, a character object containing filenames, or a list of wave files already loaded in the R environment", call. = F)
  }

  if(!is.numeric(SignalDur)){stop("A numeric value specifying the duration (in seconds) of the referencing signal must be set on the 'SignalDur' argument.")} #rever essa mensagem de erro ----
  if(!is.numeric(RefValue)){stop("A numeric value specifying the intensity (in dB SPL) of the referencing signal must be set on the 'RefValue' argument.")} #rever essa mensagem de erro ----
  if(!saveresults && !is.null(outname)) {stop("You can't set an 'outname' if 'saveresults' is FALSE", call. = F)} #rever essa mensagem de erro ----

  for(i in 1:length(arquivos)){ #Loop que analisara os arquivos####

    if(class(arquivos[[i]])=="Wave"){ #isolando o sinal de refer?ncia do som ####
      som<-extractWave(arquivos[[i]], from=2, to=SignalDur+2, xunit="time", interact=F)
    } else {
      som<-readWave(arquivos[[i]], from=2, to=SignalDur+2, units="seconds"  )
    }

    calib.value<-timbre(files=som, Leq.calib=RefValue, ref=ref, Calib.value=NULL, saveresults=F, outname=NULL, weighting=weighting, time.mess=F, stat.mess=F) #Analisa o sinal de referencia e obtem o valor de calibracao ####
    calib.value<-calib.value[,2]

    if(class(arquivos[[i]])=="Wave"){ #isolando o som a ser analisado ####
      som<-extractWave(arquivos[[i]], from=4+SignalDur, to=Inf, xunit="time", interact=F)
    } else {
      som<-readWave(arquivos[[i]], from=4+SignalDur, units="seconds"  )
    }

    if(i==1){#gerando a matriz de resultados ####
    results<-timbre(files=som, Calib.value=calib.value, ref=ref, Leq.calib=NULL, saveresults=F, outname=NULL, weighting=weighting, bands=bands, time.mess=F, stat.mess=F)
    }else {results<-rbind(results,
                 timbre(files=som, Calib.value=calib.value, ref=ref, Leq.calib=NULL, saveresults=F, outname=NULL, weighting=weighting, bands=bands, time.mess=F, stat.mess=F)
    )}

    if(is.list(arquivos) && is.null(names(arquivos))){ # colocando o nome dos arquivos na primeira coluna####
      results[i,1]<-i
    }else if(is.list(arquivos) && !is.null(names(arquivos))){
      results[i,1]<-names(arquivos[i])
    }else {
      results[i,1]<-arquivos[i]
    }

    if(saveresults) { #Salvando a matriz a cada audio analisado ####
      write.table(results,
                  paste("TimbreCalResults_", weighting, "-weighting",
                        ifelse(!is.null(outname),paste("_", outname, sep=""),paste("")),
                        ".txt", sep="")
                  ,row.names = F, col.names = T,sep = "\t", quote=F)
      }

    if(stat.mess){cat(c("File", i, "of", length(arquivos), "done."), sep=" ", fill = T)}

  }

  if(time.mess){message(cat(c("The code has run in ", format(Sys.time()-start.time), "."), sep = ""))}
  return(results)
}
