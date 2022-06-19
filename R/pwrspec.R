#' Power Spectrum of a sound file
#'
#' @description  This functions computes a power spectrum from a sound file.
#'
#' @param file A wave file path in your computer or a class Wave object already loaded into R environment.
#' @param channel Argument passed to \link[tuneR]{mono} function from \link[tuneR]{tuneR} to extract the desired channel.
#' @param from Numeric. The start time in seconds of the sample you want to analyze. Could also be relative to the end of the file (in negative values), see examples.
#' @param to Numeric. The end time in seconds of the sample you want to analyze. Could also be relative to the end of the file (in negative values), see examples.
#' @param bandpass A vector with length two with lower and upper limits of the band pass interval in Hz.
#' @param res.scale Character. Specify the kind of scale the power spectrum amplitude should be adjusted. \code{microPa} for linear values in µPa, and \code{dB} fo decibells values in dB-SPL. (By default "microPa")
#' @param ref Numerical. Reference value for dB conversion. For Sound in water the ref is 1 µPa and on air 20 µPa. (By default 1)
#'
#' @return This function returns a data.frame with Frequency (Hz) and Intensity (microPa) of the file.
#'
#' @references Power spectrum adapted from: Carcagno, S. 2013. Basic Sound Processing with R [Blog post]. Retrieved from http://samcarcagno.altervista.org/blog/basic-sound-processing-r/
#' @references Miyara, F. 2017. Software-Based Acoustical Measurements. Springer. 429 pp. DOI: 10.1007/978-3-319-55871-4
#'
#' @examples
#' data(tham)
#' pwrspec(tham)
#'
#' #with from, to and bandpass
#' pwrspec(tham, from=3.8, to=7.1, bandpass=c(900,2000))
#'
#' #from and to relative to duration (negatives values from the end of the soundfiles)
#' pwrspec(tham, from=-3.49, to=-0.9, bandpass=c(900,2000))
#'
#'
#' @export

pwrspec <- function(file, channel="left", from=0, to=Inf, bandpass=c(0,Inf), res.scale="microPa", ref=1){

  if(from < 0 & to < 0){ #ajustando para from e to relativo ao tamanho do arquivo

    if(is.character(file)) dur=readWave(file, header = T) %>% #se for um arquivo para ler
        data.frame() %>%
        transmute(dur=samples/sample.rate) %>%
        as.numeric()

    if(class(file) == "Wave") dur=duration(file) #se for um arquivo já carregado no R

    from=dur+from
    to=dur+to
  }


  #Loading the sound file ####
  if(class(file) == "Wave") som=extractWave(mono(file, channel), from = from, to=to, xunit="time", interact=F)
  if(is.character(file)) som=mono(readWave(file, from = from, to=to, units = "seconds"), channel)

  #freee unused memore ####
  rm(file)

  #Trunc samples to duration with 3 decimal places if file bigger than 1s####
  if(length(som)/som@samp.rate > 1){ #
    if(!trunc((length(som)/som@samp.rate)*10^3)/10^3==length(som)/som@samp.rate){
      som<-extractWave(som, to=(trunc((length(som)/som@samp.rate)*10^3)/10^3)*som@samp.rate, interact=F)
      if(!trunc((length(som)/som@samp.rate)*10^3)/10^3==length(som)/som@samp.rate){message("Warning: It may take a little longer than usual to analyze this file (file duration with more than three decimal places)")}
    }
  }


  s1 <- som@left/2^(som@bit-1) #scaled to the maximum possible (as result of '/2^(bitrate-1)')
  n <- length(s1)
  p <- fft(s1)
  nUniquePts <- ceiling((n+1)/2)
  p <- p[1:nUniquePts] #select just the first half since the second half is a mirror image of the first
  p <- 2*(abs(p/n)) #changed here and next if/else in 2020.11.03 to match Miyara (2017) code routine in topic 8.6.6

  if (n %% 2 > 0){  #Routine to remove Nyquist point. Odd nfft excludes Nyquist point
    p[2:length(p)] <- p[2:length(p)]
  } else {
    p[2: (length(p) -1)] <- p[2: (length(p) -1)]
  }

  freqArray <- (0:(nUniquePts-1)) * (som@samp.rate / n) #create the frequency array

  rm(som) #to free memory usage

  espec<-data.frame(Freq.Hz=freqArray, Amp.microPa=p^2) #p^2 is part of Miyara (2017) code routine in topic 8.6.6

  espec=espec %>%
    filter(Freq.Hz >= bandpass[1] & Freq.Hz <= bandpass[2])

  if(res.scale == "dB"){
    espec$Amp.dB=LineartodB(
      sqrt( #equation based on Miyara 2017 8.6.6 topic
        espec$Amp.microPa
      )/sqrt(2) #/sqrt(2) to be able to apply calibration (Miyara 2017, topic 8.6.6 codes)
      , factor="SPL", ref=ref)

    espec=select(espec, -Amp.microPa)
  }

  return(espec)

}
