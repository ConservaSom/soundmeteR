###################################################################################################################
#
#
#####Script de automatizacao de analise de propagacao de ruido com decaimento de vegetacao e relevo por linha #####
#Modificado a partir da propag.ruido para poder ser implementada em fondes de ruido em forma de linha (ex. estradas e esteiras)
#2018.05.31
#versao 0.1
###################################################################################################################

#OBSERVA??ES ####
#implementado apenas para linhas retas at? o momento. Com apenas duas coordenadas: Inicio e fim da linha.
#No futuro, implementar para linhas a partir de um shape


####Argumentos ####
#line.origin: Coordenadas da origem da linha (obrigat?riamente uma reta).
#
##line.end: Coordenadas do final da linha (obrigat?riamente uma reta).
#
#ID: matriz contendo a identificacao dos pontos que serao modelado. nao podem existir nome iguais.
#
#NIS0: Valor de NIS a 1 metro para a equacao de decaimento do ponto.
#
#Alfa: Matriz contendo os valores da constante alfa da equacao de decaimente. Deve ser uma matriz com um valor para cada ponto, mesmo que seja repetido.
#
#Beta: Semelhante ao argumento alfa, porem apresentando os valores da contante beta.
#
#elev.r: Raster contendo informacoes de elevacao. NOTA: Esse raster vai ser a referencia de tamanho,  resolucao e coordenadas de toda a analise.
#
#contour.veg: objeto do tipo SpatialPolygonDataFrame que contem o contor da regiao que nao apresenta vegetacao.
#
#save.meta: Argumento que define se deseja salvar ou nao as matrizes armazenadas em lista com os resultados de todas as etapas do calculo no seu diretorio de trabalho. O padrao e TRUE. Objetos serao salvos como "Results.propag_*" no seu diretorio de trabalho.
#
#name.meta: Opcional. o nome definido aqui serah adicionado ao fim do nome padrao de seva.meta.
#
#multicore: Logico. Se TRUE, processamento dee multiplos n?cleos ativada. By default: FALSE
#
#cl: Numerical. Cluster object. See package 'parallel' and 'snow' for details. By default: number of CPU cores less one on the current host.
######cl for dummies: o numero de Cluster significa o numero de processos que ser?o iniciados para executar o dado processamento. Fazer mais testes para entender.
#


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

