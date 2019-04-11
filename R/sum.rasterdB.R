#' Sum raster with dB values
#'
#' @name sum.rasterdB
#'
#' @description Function to sum rasters with decibel values (Intensity Level).
#' @description Last upate: 2018.06.01
#'
#' @usage sum.rasterdB(x, na.rm=FALSE)
#'
#' @author Cássio Rachid Simões <cassiorachid@@gmail.com>
#'
#' @param x A list of raster objects with values in decibels (Intensity Level) to sum.
#' @param na.rm Logical. `NA` action, if `TRUE`, ignore NA values during the sum.
#'
#' @details Include details.
#'
#' @references Include citations.
#'
#' @return One raster with the sum of all other raster of the input.
#'
#Funcao 'sum.rasterdB' versao 0.1 ####
sum.rasterdB<-function(x, na.rm=FALSE){

  require(raster)

  if(!is.list(x)){stop("The input need to be a list with multiples rasters.")}

  compareRaster(x, showwarning=F, stopiffalse=T
                , extent=TRUE, rowcol=TRUE, crs=TRUE, res=TRUE, orig=FALSE, rotation=TRUE, values=FALSE
                )

  r<-x[[1]]
  r[]<-NA
  names(r)<-NA

  results<-data.frame(matrix(NA, ncol=length(x), nrow=length(values(x[[1]]))))
  for(i in 1:length(x)){
    results[,i]<-values(x[[i]])
  }

  r[]<-10*log10( #logaritmizando de novo
    rowSums( #somando as linas
      10^(results/10) #linearizando tudo
      , na.rm=na.rm) #fim: somando as linhas
  )#fim: logaritmizando de novo

  return(r)
}
