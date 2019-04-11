
#Função para somar rasters de valores em decibeis (Nivel de Intensidade)
#cassiorachid@gmail.com
#2018.06.01

####Arguments ####
#x: A list of raster objects with values in dB (Intensity Level) to sum.
#na.rm: Logical. NA action, if TRUE, ignore NA values during the sum.

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
