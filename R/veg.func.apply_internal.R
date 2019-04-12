#' Internal noiseR function
#'
#' @description Internal \code{noiseR} function.
#' @description Function to improve the time computation of \code{\link{propag.ruido}}. It's only to be used inside an \code{\link[base]{apply}} family function.
#'
#' @note These function is not to be called by the user.
#'
#' @author Cássio Rachid Simões <cassiorachid@@gmail.com>
#'
#Funcao 'veg.func.apply' versao 0.4 ####
#2018.05.31


veg.func.apply<-function(x, pi.x, li.x=NULL , contour.veg){

  if(!is.null(li.x)){
    dists.line.cell<-gDistance(SpatialPoints(matrix(x, ncol=2)), li.x, byid=T)
    pi.x<-data.frame(li.x[which.min(dists.line.cell)])
  }


  linhax<-SpatialLines(list(Lines(Line(rbind(x, pi.x)), ID=1))) #Criando linha fonte-celula
  crs(linhax)<-crs(contour.veg)
  compr.linhax<-SpatialLinesLengths(linhax, longlat = F) #comprimento da linha fonte-celula

  if(compr.linhax>0){
    intersec.x<-gIntersection(linhax, contour.veg, byid = T) #intersection linha-borda da mina

    if(!is.null(intersec.x)){
      if(class(intersec.x)=="SpatialLines"){
        compr.intersec.x<-SpatialLinesLengths(intersec.x, longlat = F)  #distancia fonte-borda da mina
      } else {compr.intersec.x<-0}
    }
  }

  if(ifelse( exists("compr.intersec.x"), !compr.linhax<=compr.intersec.x, FALSE)){
    SpatialLinesLengths(linhax, longlat = F)-compr.intersec.x
  } else {NA}

}
