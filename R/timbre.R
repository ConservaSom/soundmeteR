#' Intensity per octave or octave thirds
#'
#' @name timbre
#'
#' @description
#'
#' @usage timbre(files="wd", weighting="none", bands="thirds", saveresults=F,
#'        outname=NULL, Leq.calib=NULL, Calib.value=NULL, time.mess=T, stat.mess=T)
#'
#' @param files The audiofile to be analyzed. Can be "wd" to get all wave files on the work directory, a file name (or a character containing a list of filenames) that exist in the work directory, or an Wave object (or a list containing more than one Wave object). (By default: "wd")
#' @param weighting Character. Indicate the weighting curve to use on the anlysis. A, B, C and none are supported. (By default: "none")
#' @param bands Character. Choose the type of frequency band of the output. "octaves" to octaves bands intervals or "thirds" to one-third octaves bands intervals. (by deafault: "thirds")
#' @param saveresults Logical. Set \code{TRUE} if you want to save a txt file with the results of the function execution. (By default: \code{FALSE})
#' @param outname Character. If \code{saveresults} is \code{TRUE}, you can specify a name to appear on the txt file name after the default name. (By default: \code{NULL})
#' @param Leq.calib Numerical. The sound intensity level (in dB) that the sound in the audio file must have (by default: NULL). Can not be set if \code{Calib.value} is also set.
#' @param Calib.value Numerical. The calibration value returned from the analysis of a reference sound using \code{Leq.calib} (by default: NULL). Can not be set if \code{Leq.calib} is also set.
#' @param time.mess Logical. Activate or deactivate message of time to complete the function execution. (By default: \code{TRUE})
#' @param stat.mess Logical. Activate or deactivate status message of the function execution. (By default: \code{TRUE})
#'
#' @author Cássio Rachid Simões <cassiorachid@@gmail.com>
#'
#' @details Caution: You need to use an audiofile with entire values of secconds of duration to avoid bugs. Example: 35s, 60s, 19s. By default, the function will trunc your audiofile to the next entire value of seconds.
#' @details These function works only with mono audiofiles.
#' @details The audio files need to have at least 44100Hz of sampling rate.
#'
#' @references Espectro de potencia baseado em http://samcarcagno.altervista.org/blog/basic-sound-processing-r/?doing_wp_cron=1495144982.9675290584564208984375
#' @references Valores de ponderacao das curvas A, B e c baseado em: Bech & Zacharov. 2006. Perceptual Audio Evaluation-Theory, Method and Application

#### Arguments ####


timbre<-function(files="wd", weighting="none", bands="thirds", saveresults=F, outname=NULL, Leq.calib=NULL, Calib.value=NULL, time.mess=T, stat.mess=T){
  start.time<-Sys.time()

  require(tuneR)

  if(class(files)=="Wave"){
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

  if(length(arquivos)==0){stop("There is no wave files on your working directory", call. = F)}

  if(!weighting=="none" && !weighting=="A" && !weighting=="B" && !weighting=="C" ) {
    stop("Wrong weighting curve. Only 'A', 'B', 'C' and 'none' are available.", call. = F)
  }

  if(!saveresults && !is.null(outname)) {stop("You can't set an 'outname' if 'saveresults' is FALSE", call. = F)}

  if((!is.null(Leq.calib)) && (!is.null(Calib.value))) {stop("You can't use 'Leq.calib' and 'Calib.value' in the same function. Please, choose only one.", call. = F)}

  if((!is.null(Leq.calib)) && (!is.numeric(Leq.calib))) {stop("Please, use a numeric value on 'Leq.calib'.", call. = F)}

  if((!is.null(Calib.value)) && (!is.numeric(Calib.value))) {stop("Please, use a numeric value on 'Calib.value'.", call. = F)}

  if(bands!="thirds" && bands!="octaves"){stop("Please, check the 'bands' argument. Only 'octaves' or 'thirds' intervals available.", call. = F)}


  #Defininfo os intervalos de frequência e montando a matriz que amazenara resultados ####
  #
  Freqbands<-c(24.8,31.3,39.4,49.6,62.5,78.7,99.2,125,157.5,198.4,250,315.0,396.9,500,630.0,793.7,1000,1259.9,1587.4,2000,2519.8,3174.8,4000,5039.7,6349.6,8000,10079.4,12699.2,16000,20158.7)
  matriz<-data.frame(matrix(data=NA,nrow=length(arquivos),ncol=2+length(Freqbands)))
  colnames(matriz)<-c("Arquivo","Leq","25","31.5","40","50","63","80","100","125","160","200","250","315","400","500","630","800","1000","1250","1600","2000","2500","3150","4000","5000","6300","8000","10000","12500","16000","20000")

  if(bands=="octaves"){
    matriz.octaves<- data.frame(matrix(data=NA,nrow=length(arquivos),ncol=2+10))
    colnames(matriz.octaves)<-c("Arquivo","Leq","31.5","63","125","250","500","1000","2000","4000","8000","16000")
  }


  if(!is.null(Leq.calib)) {
    calibration<-data.frame(matrix(data=NA,nrow=length(arquivos),ncol=2))
    colnames(calibration)<-c("Arquivo","Calib.value" )}

  for(i in 1:length(arquivos)){

    if(class(arquivos[[i]])=="Wave"){
      som<-arquivos[[i]]
      if(!trunc(length(som)/som@samp.rate)==length(som)/som@samp.rate){
        som<-extractWave(som, to=trunc(length(som)/som@samp.rate)*som@samp.rate, interact=F)
        if(!trunc(length(som)/som@samp.rate)==length(som)/som@samp.rate){stop("Something went wrong with your sound file.")}
      }
    } else {
      som<-readWave(arquivos[[i]])
      if(!trunc(length(som)/som@samp.rate)==length(som)/som@samp.rate){
        som<-readWave(arquivos[[i]], to=trunc(length(som)/som@samp.rate)*som@samp.rate)
        if(!trunc(length(som)/som@samp.rate)==length(som)/som@samp.rate){stop("Something went wrong with your sound file.")}
      }
    }

    if(som@samp.rate<44100){stop("Your audiofiles need to have at least 44100Hz of sampling rate.")}

    s1 <- som@left/2^(som@bit-1) #Mudando a escala da onda para relativo ao maximo (determinado por '2^(bitrate-1)')
    n <- length(s1)
    p <- fft(s1)
    nUniquePts <- ceiling((n+1)/2)
    p <- p[1:nUniquePts] #select just the first half since the second half is a mirror image of the first
    p <- (abs(p)/n)^2 #take the absolute  ('abs')value, or the magnitude, scale by the number of points ('(/n)') so that the magnitude does not depend on the length of the signal or on its sampling frequency, and square it ('^2') to get the power

    if (n %% 2 > 0){  #multiply by two (see technical document for details). odd nfft excludes Nyquist point
      p[2:length(p)] <- p[2:length(p)]*2 # we've got odd number of points fft
    } else {
      p[2: (length(p) -1)] <- p[2: (length(p) -1)]*2 # we've got even number of points fft
    }

    freqArray <- (0:(nUniquePts-1)) * (som@samp.rate / n) #  create the frequency array

    espec<-data.frame(Freq.Hz=freqArray, Int.dB=LineartodB(p, fac="IL", ref=1)) #Tabela contendo frequencias e energia em dBFS. DUVIDA: 20*log10(p) ou 10*log10(p)? usando 10 aqui por que a samcarcacno tmb usa, ver pagina que a função LineartodB() recomenda para entender melhor  ----

    #Calulando a quantidade de energia por banda de frequência ####
    for (j in 1:length(Freqbands)) {


      if(j==1){
        if(is.list(arquivos) && is.null(names(arquivos))){
          matriz[i,1]<-i
        }else if(is.list(arquivos) && !is.null(names(arquivos))){
          matriz[i,1]<-names(arquivos[i])
        }else {
          matriz[i,1]<-arquivos[i]
        }
      }

      sum.int<-LineartodB( sum(#se houver sqrt aqui, eh necessario guargar os resultados sem sqrt para ponderacao ver o link a seguir para compreender a ideia de colocar sqrt aqui #https://www.cirrusresearch.co.uk/blog/2020/03/calculation-of-dba-from-octave-band-sound-pressure-levels/
         dBtoLinear(espec[espec$Freq.Hz>=Freqbands[j]/(2^(1/6)) & espec$Freq.Hz<Freqbands[j]*(2^(1/6)),2], fac="IL", ref=1) #usando fator 10 (IL) por conta da conversao passada ter usado a mesma para chegar em dB
        ), fac="IL", ref=1)

      if(is.finite(sum.int)){
        matriz[i,j+2]<-sum.int
      }else{matriz[i,j+2]<-NA}
    }

    matriz[i,c(-1,-2)]<-round(matriz[i,c(-1,-2)], 1) #arredondando valores para 1 casa decimal

    #Implementando curvas de ponderacao ####
    if(weighting=="A"){ # Implementando Curva A ####
      matriz[i,c(-1,-2)]<-matriz[i,c(-1,-2)]+c(-44.7,-39.4,-34.6,-30.2,-26.2,-22.5,-19.1,-16.1,-13.4,-10.9,-8.6,-6.6,-4.8,-3.2,-1.9,-0.8,0,0.6,1,1.2,1.3,1.2,1,0.5,-0.1,-1.1,-2.5,-4.3,-6.6,-9.3)
      matriz[i,2]<-round(20*log10(sum(10^(matriz[i,c(-1,-2)]/20))),1)

      if(is.numeric(Leq.calib)) {
        calibration[i,1]<-matriz[i,1]
        calibration[i,2]<-Leq.calib-matriz[i,2]

      } else if(is.numeric(Calib.value)) {
        matriz[i,c(-1,-2)]<-matriz[i,c(-1,-2)]+Calib.value
        matriz[i,2]<-round(20*log10(rowSums(10^(matriz[i,c(-1,-2)]/20))),1)
      }

    } else if(weighting=="B") { # Implementando Curva B ####
      matriz[i,c(-1,-2)]<-matriz[i,c(-1,-2)]+c(-20.4,-17.1,-14.2,-11.6,-9.3,-7.4,-5.6,-4.2,-3,-2,-1.3,-0.8,-0.5,-0.3,-0.1,0,0,0,0,-0.1,-0.2,-0.4,-0.7,-1.2,-1.9,-2.9,-4.3,-6.1,-8.4,-11.1)
      matriz[i,2]<-round(20*log10(sum(10^(matriz[i,c(-1,-2)]/20))),1)

      if(is.numeric(Leq.calib)) {
        calibration[i,1]<-matriz[i,1]
        calibration[i,2]<-Leq.calib-matriz[i,2]

      } else if(is.numeric(Calib.value)) {
        matriz[i,c(-1,-2)]<-matriz[i,c(-1,-2)]+Calib.value
        matriz[i,2]<-round(20*log10(rowSums(10^(matriz[i,c(-1,-2)]/20))),1)
      }

    } else if(weighting=="C"){ # Implementando Curva C ####
      matriz[i,c(-1,-2)]<-matriz[i,c(-1,-2)]+c(-4.4,-3,-2,-1.3,-0.8,-0.5,-0.3,-0.2,-0.1,0,0,0,0,0,0,0,0,0,-0.1,-0.2,-0.3,-0.5,-0.8,-1.3,-2,-3,-4.4,-6.2,-8.5,-11.2)
      matriz[i,2]<-round(20*log10(sum(10^(matriz[i,c(-1,-2)]/20))),1)

      if(is.numeric(Leq.calib)) {
        calibration[i,1]<-matriz[i,1]
        calibration[i,2]<-Leq.calib-matriz[i,2]

      } else if(is.numeric(Calib.value)) {
        matriz[i,c(-1,-2)]<-matriz[i,c(-1,-2)]+Calib.value
        matriz[i,2]<-round(20*log10(rowSums(10^(matriz[i,c(-1,-2)]/20))),1)
      }

    } else if(weighting=="none"){ # Implementando Curva Z (nenhuma) ####
      matriz[i,2]<-round(
        LineartodB(  sqrt(sum( #essa equacao resulta igual ao rms do oscilograma
          dBtoLinear(matriz[i,c(-1,-2)], fac="IL", ref=1)
        )) , fac="IL", ref=1)
        ,1)

      if(is.numeric(Leq.calib)) {
        calibration[i,1]<-matriz[i,1]
        calibration[i,2]<-Leq.calib-matriz[i,2]

      } else if(is.numeric(Calib.value)) {
        matriz[i,c(-1,-2)]<-matriz[i,c(-1,-2)]+Calib.value
        matriz[i,2]<-round(20*log10(rowSums(10^(matriz[i,c(-1,-2)]/20))),1)
      }

    }


    #mudando os intervalos para bandas de oitavas ####
    if(bands=="octaves"){
      matriz.octaves[i,1:2]<-matriz[i,1:2] #adicionando o Leq a planilha de oitavas
      matriz[i,c(-1:-2)]<-10^(matriz[i,c(-1:-2)]/20) #linearizando para somar

      #Somando as intensidads das ter?as pertencentes ao intervalo da oitava:
      matriz.octaves[i,"31.5"] = matriz[i,which(colnames(matriz)=="31.5")-1] + matriz[i,which(colnames(matriz)=="31.5")] + matriz[i,which(colnames(matriz)=="31.5")+1]
      matriz.octaves[i,"63"] = matriz[i,which(colnames(matriz)=="63")-1] + matriz[i,which(colnames(matriz)=="63")] + matriz[i,which(colnames(matriz)=="63")+1]
      matriz.octaves[i,"125"] = matriz[i,which(colnames(matriz)=="125")-1] + matriz[i,which(colnames(matriz)=="125")] + matriz[i,which(colnames(matriz)=="125")+1]
      matriz.octaves[i,"250"] = matriz[i,which(colnames(matriz)=="250")-1] + matriz[i,which(colnames(matriz)=="250")] + matriz[i,which(colnames(matriz)=="250")+1]
      matriz.octaves[i,"500"] = matriz[i,which(colnames(matriz)=="500")-1] + matriz[i,which(colnames(matriz)=="500")] + matriz[i,which(colnames(matriz)=="500")+1]
      matriz.octaves[i,"1000"] = matriz[i,which(colnames(matriz)=="1000")-1] + matriz[i,which(colnames(matriz)=="1000")] + matriz[i,which(colnames(matriz)=="1000")+1]
      matriz.octaves[i,"2000"] = matriz[i,which(colnames(matriz)=="2000")-1] + matriz[i,which(colnames(matriz)=="2000")] + matriz[i,which(colnames(matriz)=="2000")+1]
      matriz.octaves[i,"4000"] = matriz[i,which(colnames(matriz)=="4000")-1] + matriz[i,which(colnames(matriz)=="4000")] + matriz[i,which(colnames(matriz)=="4000")+1]
      matriz.octaves[i,"8000"] = matriz[i,which(colnames(matriz)=="8000")-1] + matriz[i,which(colnames(matriz)=="8000")] + matriz[i,which(colnames(matriz)=="8000")+1]
      matriz.octaves[i,"16000"] = matriz[i,which(colnames(matriz)=="16000")-1] + matriz[i,which(colnames(matriz)=="16000")] + matriz[i,which(colnames(matriz)=="16000")+1]

      matriz[i,c(-1:-2)]<-round(20*log10(matriz[i,c(-1:-2)]),1) #voltando a matriz de ter?as para dB para evitar erros nas demais etapas da fun??o

      matriz.octaves[i,c(-1:-2)]<-round(20*log10(matriz.octaves[i,c(-1:-2)]),1) #convertendo a matriz de oitavas pra dB tambem

    }



    if(saveresults) { #Salvando a matriz por som analisado ####

      if(bands=="octaves"){
        matriz.2save<-matriz.octaves
      }else {
        matriz.2save<-matriz
      }
      write.table(matriz.2save,
                  paste("TimbreResults_",weighting,"-weighting",
                        ifelse(!is.null(Leq.calib),paste("_Calib.value"),paste("")),
                        ifelse(!is.null(Calib.value),paste("_Adjusted"),paste("")),
                        ifelse(!is.null(outname),paste("_",outname, sep=""),paste("")),
                        ".txt", sep="")
                  ,row.names = F, col.names = T,sep = "\t", quote=F)}

    if(stat.mess){cat(c("File", i, "of", length(arquivos), "done."), sep=" ", fill = T)}

  }


  if(bands=="octaves"){
    matriz<-matriz.octaves
  }

  if(time.mess){message(cat(c("The code has run in ", format(Sys.time()-start.time), "."), sep = ""))}

  if((is.null(Leq.calib)) & (is.null(Calib.value))) {
    return(matriz)
  }else if(is.numeric(Leq.calib)) {
    return(calibration)
  } else if(is.numeric(Calib.value)) {
    return(matriz)
  } else {warning("Something went wrong. Please, check the arguments and try again.", call. = F)}
}
