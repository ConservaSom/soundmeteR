#' Time Weighting of a Audiofile
#'
#' @name Tweighting
#'
#' @description Escrever descrição
#'
#' @param file Wave object
#' @param window Character. Wich time
#' @param Leq.calib Numeric. The sound pressure level (in dB SPL) that the signal in the audio file must have (by default: NULL). This parameter is passed to \code{\link{timbre}} function.
#' @param ... Further arguments passed to \code{\link{timbre}}.
#'
#' @details This function split your audiofile in smaller files defined as \code{fast} (0.125s) and \code{slow} (1s) and analyze each one with \code{\link{timbre}} function.
#'
#' @return A numeric vector
#'
#' @seealso \code{\link{timbre}}, \code{\link{soundmeter}}
#'
#' @examples
#' #creating an example sound file
#' som=sine(1000, duration = 44500)
#'
#' #default options without calibration (results in dBFS)
#' Tweighting(som)
#'
#' #Simulation of a calib signal with a Leq of 94dB in the field ####
#' #fast
#' Tweighting(som, window = "fast", bands="octaves", Leq.calib=94)
#' #slow
#' Tweighting(som, window = "slow", bands="octaves", Leq.calib=94)
#'
#' #Using the result of the simulation above to calibrate the sound and output
#' #fast
#' Tweighting(som, window = "fast", bands="octaves", Calib.value=309.67)
#' #slow
#' Tweighting(som, window = "slow", bands="octaves", Calib.value=309.67)
#'
#'
#'
#'
#'
#'
#'


Tweighting <- function(file, window="fast", Leq.calib=NULL,...){

  require(tuneR)

  if(window == "fast"){
    window = 0.125
  }else if(window == "slow"){
    window = 1
  }else stop("Choose a valid window size in seconds ('fast' or 'slow')")

  if(class(file) != 'Wave') stop("Only one Wave object accepted on this function")

  res=sapply(1:trunc(duration(file)/window)
             , FUN=function(x, file, samp){
               return(timbre(extractWave(file, from = round((x-1)*samp), to=round(x*samp)), stat.mess = F, time.mess = F, Leq.calib=NULL, ...)[,-1])
             }
             , file=file, samp=window*file@samp.rate)

  if(!is.numeric(res)){
    res=t(res)
    res=as.data.frame(matrix(unlist(res), nrow = nrow(res), byrow = F, dimnames = list(NULL, colnames(res))), check.names=F) #convertendo em data.frame p facilitar manuseio
  }

  if(!is.null(Leq.calib)) res=Leq.calib - seewave::meandB(res$Leq, level = "SPL") #Extraindo valor de calibração a partir do da média dos valores de RMS (exatamente o Leq)

  return(res)

}
