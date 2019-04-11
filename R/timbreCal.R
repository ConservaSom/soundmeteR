#Escrevendo funcoes para analise de intensidade por tercas de oitavas com o sinal de refer?ncia imbutido no arquivo de audio
#cassiorachid@gmail.com

#Funcao 'timbreCal' versao 0.2 ####
##2018.05.13
#Novidades da vers?o:
#O??o de escolher intervalo de bandas em oitavas ou ter?as de oitavas implementada. argumento "bands"

#Coisas para fazer:
#Pensar ao invés de usar um trecho da gravação para calibrar usar um arquivo externo.

#In attempt to use this function the audiofile must begin with 1 second of silence, followed by a reference sound with knowed intensity, followed by another 1 sencond of silence and the following sound to analyse.
#The duration of the reference sound must be specified (in seconds) on the 'SignalDur' argument and his intrensity (in dB) on the 'refValue' argument.

##### Arguments ####
#files: The audiofile to be analysed. Can be "wd" to get all wave files on the work directory, a file name (or a character containing a list of filenames) that exist in the work directory, or an Wave object (or a list containing more than one Wave object). (By default: "wd") 
#SignalDur: Numerical. Specify the reference signal duration (in seconds) on the beggining of the audiofile. (By default: NULL)
#RefValue: Numerical. Specify the reference signal intensity (in deciBells) on the beggining of the audiofile. (By default: NULL)
#weighting: Character. Indicate the weighting curve to use on the anlysis. A, B, C and none are supported. (By default: "none") 
#bands: Character. Choose the type of frequency band of the output. "octaves" to octaves bands intervals or "thirds" to one-third octaves bands intervals. (by deafault: "thirds")
#saveresults: Logical. Set TRUE if you want to save a txt file with the results of the function execution. (By default: FALSE)
#outname: Character. If saveresults is TRUE, you can specify a name to appear on the txt file name after the default name. (By defaulf: NULL)
#time.mess: Logical. Activate or deactivate message of time to complete the function execution. (By default: TRUE)
#stat.mess: Logical. Activate or deactivate status message of the function execution. (By default: TRUE)


timbreCal<- function(files="wd", SignalDur=NULL, RefValue=NULL,weighting="none", bands="thirds", saveresults=F, outname=NULL, time.mess=T, stat.mess=T){
  
  start.time<-Sys.time()
  
  require(tuneR)
  
  if(class(files)=="Wave"){ #vendo o tipo de arquivo usado no imput e armazenando em um objeto ###
    arquivos<-list(files)
  }else if(is.list(files)){ 
    arquivos<-files
  }else if(any(files=="wd")) {
    arquivos <- dir(pattern=".WAV", ignore.case=T)
  }else if(is.character(files)) {
    arquivos <- files
  }else {
    stop("Choose a valid file on 'files' argument. Could be 'wd', a valid file name, a character list of files or a list of wavefiles already loaded inde R environment", call. = F)
  }
  
  if(!is.numeric(SignalDur)){stop("A numeric value specifying the duration (in seconds) of the referencing signal must be set on the 'SignalDur' argument.")}
  if(!is.numeric(RefValue)){stop("A numeric value specifying the intensity (in dB) of the referencing signal must be set on the 'RefValue' argument.")}
  if(!exists("timbre") || !is.function(timbre)){stop("'timbre' function is required to run this function.")}
  if(!saveresults && !is.null(outname)) {stop("You can't set an 'outname' if 'saveresults' is FALSE", call. = F)}
  
  for(i in 1:length(arquivos)){ #Loop que analisara os arquivos####
    
    if(class(arquivos[[i]])=="Wave"){ #isolando o som de refer?ncia do som ####
      som<-extractWave(arquivos[[i]], from=2, to=SignalDur+2, xunit="time", interact=F)
    } else {
      som<-readWave(arquivos[[i]], from=2, to=SignalDur+2, units="seconds"  )
    }
    
    calib.value<-timbre(files=som, Leq.calib=RefValue, Calib.value=NULL, saveresults=F, outname=NULL, weighting=weighting, time.mess=F, stat.mess=F) #Analisar o sinal de referencia e obter o valor de calibra??o ####
    calib.value<-calib.value[,2]
    
    if(class(arquivos[[i]])=="Wave"){ #isolando o som a ser analisado ####
      som<-extractWave(arquivos[[i]], from=4+SignalDur, xunit="time", interact=F)
    } else {
      som<-readWave(arquivos[[i]], from=4+SignalDur, units="seconds"  )
    }
    
    if(i==1){#gerando a matriz de resultados ####
    results<-timbre(files=som, Calib.value=calib.value, Leq.calib=NULL, saveresults=F, outname=NULL, weighting=weighting, bands=bands, time.mess=F, stat.mess=F)
    }else {results<-rbind(results,
                 timbre(files=som, Calib.value=calib.value, Leq.calib=NULL, saveresults=F, outname=NULL, weighting=weighting, bands=bands, time.mess=F, stat.mess=F)
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
