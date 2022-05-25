#' Function that makes sound meter alike measurements
#'
#'
#' @param CalibPosition anda de mãos dadas com calib value. Pode ser negativo, positivo ou data.frame com essas combinações
#' @param CalibValue Anda de mãos dadas com calib position. Quando tem o position, ele é considerado o valor de referência, quando não tem o position, ele é considerado o valor de calibração.
#' @param fw Character. Argument passed to \code{\link[seewave]{dBweight}} to indicate the frequency weighting curve to use on the anlysis. 'A', 'B', 'C', 'D', 'ITU', and 'none' are supported. See \code{\link[seewave]{dBweight}} for details. (By default: "none")
#' @param tw Time weighting
#' @param bandpass ainda não implementado
#' @param time.mess Logical. Activate or deactivate message of time to complete the function execution. (By default: \code{TRUE})
#' @param stat.mess Logical. Activate or deactivate status message of the function execution. (By default: \code{TRUE})
#' @param channel Only "left" or "right" acepted. By default "left"
#' @param saveresults Logical. Set \code{TRUE} if you want to save a txt file with the results of the function execution. (By default: \code{FALSE})
#' @param outname Character. If \code{saveresults} is \code{TRUE}, you can specify a name to appear on the txt file name after the default name. (By default: \code{NULL})


soundmeter <- function(files="wd", from=0, to=Inf, CalibPosition=NULL, CalibValue=0, ref=20, fw="none", bands="octaves", banpass=NULL, tw="fast", time.mess=T, stat.mess=T, channel="left", saveresults=F, outname=NULL){

  start.time<-Sys.time()


  require(tuneR)
  require(seewave)

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

  #organizando a identificação de início e fim do trecho a analizar ####
  if(is.list(arquivos)){
    from=data.frame(1:length(arquivos), from)$from
    to=data.frame(1:length(arquivos), to)$to
    channel=data.frame(1:length(arquivos), channel)$channel
  }else {
    from=data.frame(arquivos, from)$from
    to=data.frame(arquivos, to)$to
    channel=data.frame(arquivos, channel)$channel
  }

  #início do loop maior ####
  for(i in 1:length(arquivos)){

    #Reading sound file ####
    if(class(arquivos[[i]])=="Wave"){
      som<-arquivos[[i]]
    } else {
      som<-readWave(arquivos[[i]])
    }

    #Extraindo o canal ####
    if(channel[i] == "right"){
      som=mono(som, "right")
    }

    #Ajustando informações de CalibValue ####
    if(!is.null(CalibValue) & length(CalibValue)>1){

      if(i==1 & length(CalibValue)!=length(arquivos)) stop("When specifying multiples CalibValues, it and files must be the same length")

      cal.val=CalibValue[i]

    }else {cal.val=CalibValue}


    #Análise sem calibração (em dBFS) ####
    if(is.null(CalibValue)){

      matriz=Tweighting(extractWave(som, from=from[i] , to=to[i], xunit = "time", interact=F), window=tw, bands=bands, weighting="none")

    }else if(!is.null(CalibPosition) & !is.null(CalibValue)){#Organizando extraindo valor de calibração ####

      #Verificando de o local da calibração veio num data frame ou é valor único para todos
      if(is.data.frame(CalibPosition)){ #se for data. frame

        if(i==1 & nrow(CalibPosition)!=length(arquivos) & ncol(CalibPosition)!=2) stop("When CalibPosition is a data frame, it must be the same length of files")
        cal.pos = CalibPosition[i,]

      }else if(is.numeric(CalibPosition) & length(CalibPosition) == 2){#se for valores únicos para todos os arquivos

        cal.pos=CalibPosition

      }else {stop("Something went wrong, please check your CalibPosition")}

      #Verificando se o local da calibração é por tempo ou relativo à duração do arquivo
      if(all(cal.pos < 0) & duration(som)>=cal.pos[2]-cal.pos[1]){

        cal.pos = c(duration(som)+cal.pos[1], duration(som)+cal.pos[2])

      }else if(all(!(cal.pos >= 0) | duration(som)<cal.pos[2]-cal.pos[1]) | any(is.na(cal.pos))){stop("Some values of your CalibPosition are wrong")}


      #calibração propriamente dita ####
      cal.val=Tweighting(extractWave(som, from=cal.pos[1] , to=cal.pos[2], xunit = "time", interact=F), window=tw, bands="thirds", Leq.calib=cal.val, weighting="none") #considerando curva Z p calibração

    }

    # if(exists("cal.val")){
      #Usando CalibValue como valor de calibração diretamente ####
      matriz=Tweighting(extractWave(som, from=from[i] , to=to[i], xunit = "time", interact=F), window=tw, bands=bands, weighting=fw, Calib.value = cal.val)
    # }


    #criando e armazenando valores na matriz de resultados ####
    if(i == 1){
      res = data.frame(matrix(data=NA,nrow=length(arquivos),ncol=8+ncol(matriz[-1])))
      colnames(res) = c("File", "Lmin", "Lmax", "L90", "L50", "L10", "Lpeak", "Leq", colnames(matriz[-1]))
    }

    res$File[i]=ifelse(class(arquivos[[i]]) == "Wave", i, arquivos[[i]]) #Nome ou número do arquivo
    res$Lmin[i]=min(matriz$Leq) #Lmin
    res$Lmax[i]=max(matriz$Leq) #Lmax

    res[i,c("L90","L50","L10")] = LineartodB(quantile(dBtoLinear(matriz$Leq, factor = "SPL", ref = ref), probs=c(0.1, 0.5, 0.9)), factor = "SPL", ref = ref) #L90,L50 e L10

    res$Lpeak[i]=LineartodB(max(abs(extractWave(som, from=from[i] , to=to[i], xunit = "time", interact=F)@left/2^(som@bit -1))), factor="SPL", ref=ref) #Lpeak
    if(!is.null(cal.val)) res$Lpeak[i] = res$Lpeak[i]+cal.val

    res[i,8:ncol(res)]=apply(matriz, 2, meandB, level="SPL") #Leq e bandas

    res[i,-1]=round(res[i,-1],2) #arredondando valores para duas casas decimais


    if(saveresults) { #Salvando a matriz a cada audio analisado ----
      write.table(res,
                  paste("soundmeterResult_", fw, "-weighting_",
                        tw,
                        ifelse(!is.null(outname),paste("_", outname, sep=""),paste("")),
                        ".txt", sep="")
                  ,row.names = F, col.names = T,sep = "\t", quote=F)
    }

    if(stat.mess){cat(c("File", i, "of", length(arquivos), "done."), sep=" ", fill = T)} #stat.mess ----

  } #final do loop maior ####

  if(time.mess){message(cat(c("The code has run in ", format(Sys.time()-start.time), "."), sep = ""))}

  return(res)

}

