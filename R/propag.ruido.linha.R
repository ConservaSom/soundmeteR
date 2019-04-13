#' Sound decay models for lines on GIS
#'
#' @name propag.ruido.linha
#'
#' @description
#'
#' @usage propag.ruido.linha(line.origin=NA, line.end=NA, ID=NULL, NIS0=NA, Alfa=NA,
#'        Beta=NA, elev.r=NA, contour.veg=NA, save.meta=T, name.meta=NULL, multicore=F,
#'        cl=detectCores()-1)
#'
#' @param line.origin Matrix with longitude (first column) and latitude (second column) of the beggining  of each line to be modeled (only straight lines are acepted).
#' @param line.end Matrix with longitude (first column) and latitude (second column) of the end  of each line to be modeled (only straight lines are acepted)
#' @param ID Matrix with IDs for each line the user wants to model. Do not use repeated names.
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
#' @note Function modified from \code{\link{propag.ruido}}.
#'
#' @seealso \code{\link{propag.ruido}}
#'
#' @example
#'
#'
#2018.05.31
#versao 0.1

propag.ruido.linha <-function(line.origin=NA, line.end=NA, ID=NULL, NIS0=NA, Alfa=NA, Beta=NA, elev.r=NA, contour.veg=NA, save.meta=T, name.meta=NULL, multicore=F, cl=detectCores()-1){

  start.time<-Sys.time()

  if(nrow(line.origin)!=nrow(line.end)){stop("'line.origin' and 'line.end' have different lengths, please check your input data.")}

  #Carregando pacotes necessarios ####
  require(sp)
  require(raster)
  require(rgeos)
  if(multicore){require(parallel)}

  #Criando camada em branco para referencia ####
  r <- raster(extent(elev.r@extent[1:4]), ncol=elev.r@ncols, nrow=elev.r@nrows, crs=crs(elev.r))
  r.points<-rasterToPoints(r) #convertendo camada para pontos
  S.r.points<-as(r, "SpatialPoints") #convertendo camada para pontos espaciais


  results.propag<-list()
  results.propag[["decai"]]<-data.frame(matrix(data = NA,nrow = length(values(r)), ncol=nrow(line.origin)))
  colnames(results.propag[["decai"]])<-paste("ID",ID, sep = "")

  results.propag[["veg"]]<-data.frame(matrix(data = NA,nrow = length(values(r)), ncol=nrow(line.origin)))
  colnames(results.propag[["veg"]])<-paste("ID",ID, sep = "")

  results.propag[["veg_spread"]]<-data.frame(matrix(data = NA,nrow = length(values(r)), ncol=nrow(line.origin)))
  colnames(results.propag[["veg_spread"]])<-paste("ID",ID, sep = "")

  results.propag[["rel.h.fonte"]]<-data.frame(matrix(data = NA,nrow = length(values(r)), ncol=nrow(line.origin)))
  colnames(results.propag[["rel.h.fonte"]])<-paste("ID",ID, sep = "")

  results.propag[["rel.h.recep"]]<-data.frame(matrix(data = NA,nrow = length(values(r)), ncol=nrow(line.origin)))
  colnames(results.propag[["rel.h.recep"]])<-paste("ID",ID, sep = "")

  results.propag[["final.h.fonte"]]<-data.frame(matrix(data = NA,nrow = length(values(r)), ncol=nrow(line.origin)))
  colnames(results.propag[["final.h.fonte"]])<-paste("ID",ID, sep = "")

  results.propag[["final.h.recep"]]<-data.frame(matrix(data = NA,nrow = length(values(r)), ncol=nrow(line.origin)))
  colnames(results.propag[["final.h.recep"]])<-paste("ID",ID, sep = "")
  ########################################################################################################################################
  #Fazendo os calculos para o ponto i ####
  for (i in 1:nrow(line.origin)) {

    #Definindo localizacao do emissor
    S.li.x<-SpatialLines(list(Lines(Line(rbind(as.numeric(line.origin[i,]), as.numeric(line.end[i,]))),ID=1))) #criando a linha
    crs(S.li.x)<-crs(r) #atribuindo referencia geografica

    #Transformando a linha em pontos ####
    li.x<-as(rasterize(S.li.x, r), "SpatialPoints")


    #calulando relevo relativo a altitude do emissor
    elev.rel<-elev.r #Resolvido por hora, pensar numa mudan?a junto com o apply depois ####
    #values(elev.rel)<-values(elev.rel)-extract(elev.rel, S.pi.x)

    #####################################################################################################################################
    #Decaimento cilindrico####
    di<-gDistance(S.li.x, S.r.points, byid=T) #criando camada de distancia

    NIS1 <- NIS0[i] - Beta[i]*log10(di) - Alfa[i]*(di) #calculando o decaimento com a distancia

    results.propag[["decai"]][,i]<-NIS1



    #####################################################################################################################################
    #Correcao da vegetacao####

    #correc.veg<-data.frame(matrix(data=NA,nrow=nrow(r.points),ncol=1))
    #colnames(correc.veg)<-c("d2")

    if(multicore){ #implementando multicore ####
      cl.int <- makeCluster(cl, type = "PSOCK")
      clusterExport(cl.int, list("SpatialLines", "Lines", "Line", "crs","crs<-","SpatialLinesLengths", "gIntersection", "gDistance", "SpatialPoints")) #enviando as fun??es necess?rias para os processos

      correc.veg<-data.frame(d2=parRapply(cl.int, r.points, veg.func.apply, pi.x=NULL, li.x=li.x, contour.veg=contour.veg))

      stopCluster(cl.int)

    }else {
      correc.veg<-data.frame(d2=apply(r.points, 1, veg.func.apply, pi.x=NULL, li.x=li.x, contour.veg=contour.veg))
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
      clusterExport(cl.int, list("SpatialLines", "Lines", "Line", "crs", "crs<-", "SpatialLinesLengths", "extract", "SpatialPoints", "rasterToPoints", "gDistance", "SpatialPoints"))

      correc.rel<-parRapply(cl.int, r.points, rel.func.apply, pi.x=NULL, li.x=li.x, r=r, elev.rel=elev.rel, S.pi.x=NULL)
      correc.rel<-data.frame(matrix(correc.rel,ncol = 4, byrow = T))
      colnames(correc.rel)<-c("X","R","h","h.recep")

      stopCluster(cl.int)

    }else {
      correc.rel<-data.frame(t(apply(r.points, 1, rel.func.apply, pi.x=NULL, li.x=li.x, r=r, elev.rel=elev.rel, S.pi.x=NULL)))
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
      save(results.propag, file=paste("Results.propag.linha_", ifelse(!is.null(name.meta), name.meta, ""),".RData", sep = ""))
    }
  } #encerando loop linha i

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

