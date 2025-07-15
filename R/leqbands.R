#' Equivalent Level per octaves or one-third octaves
#'
#' @name leqbands
#'
#' @description Function to comput the Equivalent Level (Leq) across octaves or one-third octaves.
#'
#' @param files The audio file(s) to be analyzed. Can be set to "wd" to get all ".wav" files in the work directory, or a single file name, or a character containing a list of file names, or an Wave object, or a list containing more than one Wave object. Only ".wav" files are accepted. By default "wd".
#' @param channel Character. Choose “left” or “right” channel. Argument passed to \link[tuneR]{mono} function from \link[tuneR]{tuneR} to extract the desired channel.
#' @param from Numeric. The start time (in seconds) of the segment to analyze. Could also be relative to the end of the file (in negative values). See examples for details.
#' @param to Numeric. The end time (in seconds) of the segment to analyze. Could also be relative to the beginning of the file (in negative values). See examples for details.
#' @param weighting Character. Defines the weighting curve for the analysis, passed to the \code{\link[seewave]{dBweight}} function. Accepted values are 'A', 'B', 'C', 'D', 'ITU', and 'none'. See \code{\link[seewave]{dBweight}} for details. By default "none".
#' @param bands Character. Specifies the type of frequency band for the output. Use "octaves" for octaves bands or "thirds" to one-third octaves bands intervals (default "thirds".
#' @param ref Numerical. The reference value for dB conversion. For sound in water, the common reference is 1 µPa, and for sound in air, 20 µPa. By default 20.
#' @param Leq.calib Numerical. The sound pressure level (in dB SPL) for the signal in the audio file. Cannot be set if \code{Calib.value} is specified. By default \code{NULL}).
#' @param Calib.value Numerical. The calibration value returned from analyzing a reference signal using \code{Leq.calib}. Can not be set if \code{Leq.calib} is specified. By default \code{NULL}.
#' @param saveresults Logical. Set to \code{TRUE} to save the results to a .txt file. By default \code{FALSE}.
#' @param outname Character. If \code{saveresults} is \code{TRUE}, you can specify a name to append to the default file name on the txt file. By default \code{NULL}.
#' @param progressbar Logical. Set \code{TRUE} or \code{FALSE} to display or hide a progress bar showing elapsed time and the last concluded file number. By default \code{TRUE}.
#'
#' @details Caution: Ensure the audio file has a duration that is a whole number of seconds to avoid bugs (e.g., 35s, 60s, 19s). By default, the function will truncate the audio file to the next full second.
#' @details This function works only with mono audio files.
#' @details Audio files must have a minimum sampling rate of 44100Hz.
#' @details If you are working with decibels at full scale (dBFS), we recommend setting \code{ref=1}. This will make the results relative to 0 dBFS.
#'
#' @references Power spectrum adapted from: Carcagno, S. 2013. Basic Sound Processing with R [Blog post]. Retrieved from http://samcarcagno.altervista.org/blog/basic-sound-processing-r/
#' @references Miyara, F. 2017. Software-Based Acoustical Measurements. Springer. 429 pp. DOI: 10.1007/978-3-319-55871-4
#'
#' @return ESCREVER
#'
#' @examples
#' data(tham)
#' leqbands(tham)
#' leqbands(tham, Calib.value=130.24)
#' leqbands(tham, Calib.value=130.24, weighting="A")
#' leqbands(tham, Calib.value=130.24, weighting="A", bands="octaves")
#'
#' leqbands(tham, Leq.calib=48.17, weighting="A")
#'
#' leqbands(tham, from=-3.49, to=-0.9, ref=1) #dB at Full Scale
#' leqbands(tham, from=-3.49, to=-0.9, Calib.value=130.24) #song Sound Pressure Level
#' leqbands(tham, from=0, to=3.8, Calib.value=130.24) #background Sound Pressure Level
#'
#'
#'
#' @export

leqbands<-function(files="wd", channel="left", from=0, to=Inf, weighting="none",
                 bands="thirds", ref=20, Leq.calib=NULL, Calib.value=NULL,
                 saveresults=F, outname=NULL, progressbar=T){

  if(class(files)=="Wave"){
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

  if(length(arquivos)==0){
    stop("There is no wave files on your working directory", call. = F)
  }

  if(!saveresults && !is.null(outname)) {
    stop("You can't set an 'outname' if 'saveresults' is FALSE", call. = F)
  }

  if((!is.null(Leq.calib)) && (!is.null(Calib.value))) {
    stop("You can't use 'Leq.calib' and 'Calib.value' in the same function. Please, choose only one.", call. = F)
  }

  if((!is.null(Leq.calib)) && (!is.numeric(Leq.calib))) {
    stop("Please, use a numeric value on 'Leq.calib'.", call. = F)
  }

  if((!is.null(Calib.value)) && (!is.numeric(Calib.value))) {
    stop("Please, use a numeric value on 'Calib.value'.", call. = F)
  }

  if(bands!="thirds" && bands!="octaves"){
    stop("Please, check the 'bands' argument. Only 'octaves' or 'thirds' intervals available.", call. = F)
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

  #Definindo os intervalos de frequência e montando a matriz que amazenara resultados ####
  Freqbands<-c(24.8,31.3,39.4,49.6,62.5,78.7,99.2,125,157.5,198.4,250,315.0,
               396.9,500,630.0,793.7,1000,1259.9,1587.4,2000,2519.8,3174.8,4000,
               5039.7,6349.6,8000,10079.4,12699.2,16000,20158.7)
  matriz<-data.frame(matrix(data=NA, nrow=length(arquivos),
                            ncol=2+length(Freqbands)
  )
  )
  colnames(matriz)<-c("Arquivo","Leq","25","31.5","40","50","63","80","100",
                      "125","160","200","250","315","400","500","630","800",
                      "1000","1250","1600","2000","2500","3150","4000","5000",
                      "6300","8000","10000","12500","16000","20000")


  for(i in 1:length(arquivos)){

    espec=pwrspec(arquivos[[i]], channel=channel, from=from, to=to,
                  res.scale = "dB", ref=ref)

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

      sum.int<-
        sumdB(
          espec[espec$Freq.Hz>=Freqbands[j]/(2^(1/6)) &
                  espec$Freq.Hz<Freqbands[j]*(2^(1/6)),2]
          , level="IL")

      if(is.finite(sum.int)){
        matriz[i,j+2]<-sum.int
      }else{matriz[i,j+2]<-NA}
    }

    matriz[i,c(-1,-2)]<-round(matriz[i,c(-1,-2)], 2) #arredondando valores para 2 casas decimais

    #Implementando curvas de ponderacao ####
    if(any(weighting == c("A", "B", "C", "D", "ITU"))){
      matriz[i,c(-1,-2)]<-matriz[i,c(-1,-2)]+round(seewave::dBweight(as.numeric(colnames(
        matriz[,-1:-2])))[[weighting]],2)

    } else if(weighting != "none"){
      stop("Wrong weighting curve. Only 'A', 'B', 'C', 'D', 'ITU', and 'none' accepted. See dBweight()' for details.")
    }

    #Calculando Leq ####
    #Equation 1.83 from Miraya (2017) to sum one-third octave bands
    #It needs to be factor 10 (same as 'IL') and without reference (same as 1)
    matriz[i,2]<-round(sumdB(matriz[i,c(-1,-2)], level="IL", na.rm=T),2)

    #Gerando valor de calibração ####
    if(is.numeric(Leq.calib)) {
      if(i==1){
        calibration<-data.frame(matrix(data=NA,nrow=length(arquivos),ncol=2))
        colnames(calibration)<-c("Arquivo","Calib.value" )
      }

      calibration[i,1]<-matriz[i,1]
      calibration[i,2]<-Leq.calib-matriz[i,2]


    } else if(is.numeric(Calib.value)) { #calibrando ####
      matriz[i,c(-1,-2)]<-matriz[i,c(-1,-2)]+Calib.value
      matriz[i,2]<-round(sumdB(matriz[i,c(-1,-2)], level="IL", na.rm = T),2)
    }

    #mudando os intervalos para bandas de oitavas ####
    if(bands=="octaves"){
      if(i==1){
        matriz.octaves<- data.frame(matrix(data=NA,nrow=length(arquivos),ncol=2+10))
        colnames(matriz.octaves) <- c("Arquivo","Leq","31.5","63","125","250",
                                      "500","1000","2000","4000","8000","16000")
      }

      matriz.octaves[i,1:2]<-matriz[i,1:2] #adicionando o Leq a planilha de oitavas
      matriz[i,c(-1:-2)]<-dBtoLinear(matriz[i,c(-1:-2)], factor="IL", ref=1) #Assumming that eq. 1.83 from Miyara (2017) can be applyed here too. It needs to be factor 10 (same as 'IL') and without reference (same as 1)

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

      matriz[i,c(-1:-2)]<-round(LineartodB(matriz[i,c(-1:-2)], factor = "IL",
                                           ref=1),2) #here IL and ref is relative to eq 1.83 from Miyara (2017). It needs to be factor 10 (same as 'IL') and without reference (same as 1)

      matriz.octaves[i,c(-1:-2)]<-round(LineartodB(matriz.octaves[i,c(-1:-2)],
                                                   factor = "IL", ref=1),2) #here IL and ref is relative to eq 1.83 from Miyara (2017). It needs to be factor 10 (same as 'IL') and without reference (same as 1)

    }



    if(saveresults) { #Salvando a matriz por som analisado ####

      if(bands=="octaves"){
        write.table(matriz.octaves,
                    paste("leqbandsResults_",weighting,"-weighting",
                          ifelse(!is.null(Leq.calib),paste("_Calib.value"),paste("")),
                          ifelse(!is.null(Calib.value),paste("_Adjusted"),paste("")),
                          ifelse(!is.null(outname),paste("_",outname, sep=""),paste("")),
                          ".txt", sep="")
                    ,row.names = F, col.names = T,sep = "\t", quote=F)
      }else {
        write.table(matriz,
                    paste("leqbandsResults_",weighting,"-weighting",
                          ifelse(!is.null(Leq.calib),paste("_Calib.value"),paste("")),
                          ifelse(!is.null(Calib.value),paste("_Adjusted"),paste("")),
                          ifelse(!is.null(outname),paste("_",outname, sep=""),paste("")),
                          ".txt", sep="")
                    ,row.names = F, col.names = T,sep = "\t", quote=F)
      }}

    if(progressbar) pb$tick()

  }


  if(bands=="octaves"){
    matriz<-matriz.octaves
  }

  if((is.null(Leq.calib)) & (is.null(Calib.value))) {
    return(matriz)
  }else if(is.numeric(Leq.calib)) {
    return(calibration)
  } else if(is.numeric(Calib.value)) {
    return(matriz)
  } else {warning("Something went wrong. Please, check the arguments and try again.",
                  call. = F)
  }
}
