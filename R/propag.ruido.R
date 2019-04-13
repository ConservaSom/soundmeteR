#' Sound decay models for points on GIS
#'
#' @name propag.ruido
#'
#' @description
#'
#' @usage propag.ruido(points=NA, ID=NULL, NIS0=NA, Alfa=NA, Beta=NA, elev.r=NA,
#'        contour.veg=NA, save.meta=T, name.meta=NULL, multicore=F,
#'        cl=detectCores()-1)
#'
#' @param points Matrix with longitude (first column) and latitude (second column) of each point to be modeled.
#' @param ID Matrix with IDs for each point the user wants to model. Do not use repeated names.
#' @param NIS0 Sound intensity level at 1m. Parameter of the sound decay model. See details for explanations of the equation.
#' @param Alfa Matrix with values of  \emph{alfa} constant of the decay equation (see details for more information). Need to be one row for each point, even if they are the same.
#' @param Beta Same as the \emph{alfa} parameter, but representing values of \emph{beta} constant.
#' @param elev.r Raster with elevation values. OBS.: This raster will be used to extract all relevant GIS information (\emph{i. e.} reslution, size, crs).
#' @param contour.veg Object type 'SpatialPolygonDataFrame' with the contour of the open areas (areas without vegetagion cover). Still in improvement, need to work with multiple open areas.
#' @param save.meta Logical. If \code{TRUE} all steps of the process of modeling will be saved in your work directory. Objects will be saved with 'Results.propag_' at the beggining of the file name. (By default = \code{TRUE})
#' @param name.meta Optional. Define a custom name to be added after de default name if \code{save.meta} is \code{TRUE}.
#' @param multicore Logical. If \code{TRUE} multicore processing will be active. This parameter improve the processing time. (by deafult = \code{TRUE})
#' @param cl Numerical. Cluster object. Number of processes to be used in the \code{propag.ruido} computation. Only considered if \code{multicore} is \code{TRUE}. See documentation of package \code{\link{parallel}} and \code{\link{snow}} for details. (By default = number of CPU cores of the current host less one)
#'
#' @details Include here the equations used and their explanation.
#'
#' @references
#'
#' @return
#'
#' @author Cássio Rachid Simões <cassiorachid@@gmail.com>
#' @author Carlos Barros de Araújo <cabarau@@gmail.com>
#'
#' @seealso \code{\link{propag.ruido.linha}}
#'
#' @example
#'
#'

#####Script de automatizacao de analise de propagacao de ruido com decaimento de vegetacao e relevo por ponto #####
#2018.05.28
#versao 0.3

propag.ruido <-function(points=NA, ID=NULL, NIS0=NA, Alfa=NA, Beta=NA, elev.r=NA, contour.veg=NA, save.meta=T, name.meta=NULL, multicore=F, cl=detectCores()-1){

  start.time<-Sys.time()

  #Carregando pacotes necessarios ####
  require(sp)
  require(raster)
  require(rgeos)
  if(multicore){require(parallel)}

  #Criando camada em branco para referencia ####
  r <- raster(extent(elev.r@extent[1:4]), ncol=elev.r@ncols, nrow=elev.r@nrows, crs=crs(elev.r))
  r.points<-rasterToPoints(r) #convertendo camada para pontos


  results.propag<-list()
  results.propag[["decai"]]<-data.frame(matrix(data = NA,nrow = length(values(r)), ncol=nrow(points)))
  colnames(results.propag[["decai"]])<-paste("ID",ID, sep = "")

  results.propag[["veg"]]<-data.frame(matrix(data = NA,nrow = length(values(r)), ncol=nrow(points)))
  colnames(results.propag[["veg"]])<-paste("ID",ID, sep = "")

  results.propag[["veg_spread"]]<-data.frame(matrix(data = NA,nrow = length(values(r)), ncol=nrow(points)))
  colnames(results.propag[["veg_spread"]])<-paste("ID",ID, sep = "")

  results.propag[["rel.h.fonte"]]<-data.frame(matrix(data = NA,nrow = length(values(r)), ncol=nrow(points)))
  colnames(results.propag[["rel.h.fonte"]])<-paste("ID",ID, sep = "")

  results.propag[["rel.h.recep"]]<-data.frame(matrix(data = NA,nrow = length(values(r)), ncol=nrow(points)))
  colnames(results.propag[["rel.h.recep"]])<-paste("ID",ID, sep = "")

  results.propag[["final.h.fonte"]]<-data.frame(matrix(data = NA,nrow = length(values(r)), ncol=nrow(points)))
  colnames(results.propag[["final.h.fonte"]])<-paste("ID",ID, sep = "")

  results.propag[["final.h.recep"]]<-data.frame(matrix(data = NA,nrow = length(values(r)), ncol=nrow(points)))
  colnames(results.propag[["final.h.recep"]])<-paste("ID",ID, sep = "")
  ########################################################################################################################################
  #Fazendo os calculos para o ponto i ####
  for (i in 1:nrow(points)) {

    #Definindo localizacao do emissor
    X=points[i,1]
    Y=points[i,2]
    pi.x<-cbind(X,Y)
    S.pi.x<-SpatialPoints(pi.x) #Criando spatial point e pi.x
    crs(S.pi.x)<-crs(r) #atribuindo referencia geografica

    #calulando relevo relativo a altitude do emissor
    elev.rel<-elev.r
    values(elev.rel)<-values(elev.rel)-extract(elev.rel, S.pi.x)

    #####################################################################################################################################
    #Decaimento####
    di<-distanceFromPoints(r, S.pi.x) #criando camada de distancia

    NIS1 <- NIS0[i] - Beta[i]*log10(di) - Alfa[i]*(di) #calculando o decaimento com a distancia

    results.propag[["decai"]][,i]<-values(NIS1)



    #####################################################################################################################################
    #Correcao da vegetacao####

    #correc.veg<-data.frame(matrix(data=NA,nrow=nrow(r.points),ncol=1))
    #colnames(correc.veg)<-c("d2")


    if(multicore){ #implementando multicore ####
      cl.int <- makeCluster(cl, type = "PSOCK")
      clusterExport(cl.int, list("SpatialLines", "Lines", "Line", "crs","crs<-","SpatialLinesLengths", "gIntersection")) #enviando as fun??es necess?rias para os processos

      correc.veg<-data.frame(d2=parRapply(cl.int, r.points, veg.func.apply, pi.x=pi.x, contour.veg=contour.veg))

      stopCluster(cl.int)

    }else {
      correc.veg<-data.frame(d2=apply(r.points, 1, veg.func.apply, pi.x=pi.x, contour.veg=contour.veg))
    }

    correc.veg$perda.dB<-5.2504*log(correc.veg$d2)-9.8094
    correc.veg[correc.veg$perda.dB<=0 & !is.na(correc.veg$perda.dB), 2]<-NA
    results.propag[["veg_spread"]][,i]<-correc.veg$perda.dB

    correc.veg$perda.dB_carol<-2.3902*log10(correc.veg$d2)
    correc.veg[correc.veg$perda.dB_carol<=0 & !is.na(correc.veg$perda.dB_carol), 2]<-NA
    results.propag[["veg"]][,i]<-correc.veg$perda.dB_carol

    #####################################################################################################################################
    #Correcao do relevo####

    #correc.rel<-data.frame(matrix(data=NA,nrow=nrow(r.points),ncol=4))
    #colnames(correc.rel)<-c("X","R","h","h.recep")


    if(multicore){
      cl.int <- makeCluster(cl, type = "PSOCK")
      clusterExport(cl.int, list("SpatialLines", "Lines", "Line", "crs", "crs<-", "SpatialLinesLengths", "extract", "SpatialPoints", "rasterToPoints", "gDistance"))

      correc.rel<-parRapply(cl.int, r.points, rel.func.apply, pi.x=pi.x, r=r, elev.rel=elev.rel, S.pi.x=S.pi.x)
      correc.rel<-data.frame(matrix(correc.rel,ncol = 4, byrow = T))
      colnames(correc.rel)<-c("X","R","h","h.recep")

      stopCluster(cl.int)

    }else {
      correc.rel<-data.frame(t(apply(r.points, 1, rel.func.apply, pi.x=pi.x, r=r, elev.rel=elev.rel, S.pi.x=S.pi.x)))
      colnames(correc.rel)<-c("X","R","h","h.recep")
    }


    correc.rel.pes<-correc.rel*0.3048 #convertendo de metros para pes

    #correcao altura pela fonte
    BPD1<-with(correc.rel.pes, (sqrt((h^2)+(R^2))+sqrt((h^2)+((X-R)^2))-X))
    results.propag[["rel.h.fonte"]][,i]<-13.547*BPD1^0.2293

    #correcaoaltura pelo receptor
    BPD2<-with(correc.rel.pes, (sqrt((h.recep^2)+(R^2))+sqrt((h.recep^2)+((X-R)^2))-X))
    results.propag[["rel.h.recep"]][,i]<-13.547*BPD2^0.2293

    #####################################################################################################################################
    #Criando valores finais ####
    results.propag[["final.h.fonte"]][,i]<-rowSums(cbind(results.propag[["decai"]][,i],-results.propag[["veg"]][,i],-results.propag[["rel.h.fonte"]][,i]), na.rm = T)

    results.propag[["final.h.recep"]][,i]<-rowSums(cbind(results.propag[["decai"]][,i],-results.propag[["veg"]][,i],-results.propag[["rel.h.recep"]][,i]), na.rm = T)

    #####################################################################################################################################
    #salvando a lista com os resultados de pi.x ####
    if(save.meta){
      save(results.propag, file=paste("Results.propag_", ifelse(!is.null(name.meta), name.meta, ""),".RData", sep = ""))
    }
  } #encerando loop ponto i

  #######################################################################################################################################
  #Criando os dois rasters finais ####
  f.h.fonte<-10^(results.propag[["final.h.fonte"]]/10)
  f.h.fonte<-10*log10(rowSums(f.h.fonte, na.rm=T))

  f.h.recep<-10^(results.propag[["final.h.recep"]]/10)
  f.h.recep<-10*log10(rowSums(f.h.recep, na.rm=T))

  final.rasters.x<-list()
  final.rasters.x[["final.h.fonte"]]<-r
  values(final.rasters.x[["final.h.fonte"]])<-f.h.fonte
  final.rasters.x[["final.h.recep"]]<-r
  values(final.rasters.x[["final.h.recep"]])<-f.h.recep

  if(save.meta){
    writeRaster(final.rasters.x[["final.h.fonte"]], filename=paste("Results.propag_raster.h.emis_", ifelse(!is.null(name.meta), name.meta, ""),".tif", sep = ""),format="GTiff", overwrite=T)
    writeRaster(final.rasters.x[["final.h.recep"]], filename=paste("Results.propag_raster.h.recep_", ifelse(!is.null(name.meta), name.meta, ""),".tif", sep = ""),format="GTiff", overwrite=T)
  }

  print(Sys.time()-start.time)
  return(final.rasters.x)

}



