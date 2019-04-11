#Escrevendo funcoes para analise de propagação de rúido
#cassiorachid@gmail.com

#Funcao 'rel.func.apply' versao 0.2 ####
#2018.05.30

#Novidades da versão:
#Implementado opção para trablahar com linhas.

#Internal. This function is not to be called by the user.

#Function to improve the time computation of 'propag.ruido' function.
#It's only to be used inside an apply family function.

rel.func.apply<-function(x, pi.x, li.x=NULL, r, elev.rel, S.pi.x){
  
  if(!is.null(li.x)){
    dists.line.cell<-gDistance(SpatialPoints(matrix(x, ncol=2)), li.x, byid=T)
    S.pi.x<-li.x[which.min(dists.line.cell)]
    pi.x<-as.numeric(data.frame(S.pi.x))
  }
  
  linhax<-SpatialLines(list(Lines(Line(rbind(x, pi.x)),ID=1))) #Criando linha fonte-celula
  crs(linhax)<-crs(r)
  compr.linhax<-SpatialLinesLengths(linhax, longlat = F) #distancia linha fonte-celula
  
  if(compr.linhax>0){
    
    rel.linha.x<-as.data.frame(extract(elev.rel,linhax, cellnumbers=T)) #Valores de relevo ao longo da linha fonte-celula
    
    #deixando o relevo relativo à altura do ponto mais próximo da linha.
    if(!is.null(li.x)){ #teste lógico paliativo até ver uma estratégia melhor 
      rel.linha.x[,2]<-rel.linha.x[,2]-extract(elev.rel, S.pi.x)
      h.recep<-extract(elev.rel,SpatialPoints(rbind(x)))-extract(elev.rel, S.pi.x)
    } else {
      h.recep<-extract(elev.rel,SpatialPoints(rbind(x)))
    }
    
    rel.linha.x<-rel.linha.x[which.max(rel.linha.x[,2]),] #valor maximo de relevo e qual celula ele esta
    
    h.rel.recep.x<-rel.linha.x[,2]-h.recep #extraindo altura do ponto maximo relativo ao receptor
    
    coord.rel.max<-SpatialPoints(t(matrix(rasterToPoints(r)[rel.linha.x[,1],]))) #Localizacao no mapa da coordenada maxima
    crs(coord.rel.max)<-crs(r)
    dist.bar.x<-gDistance(S.pi.x, coord.rel.max) #distancia fonte-barreira
    
  } else {
    
    dist.bar.x = 0
    rel.linha.x = data.frame(0,0)
    h.rel.recep.x = 0
    
  }
  
  
  cbind(
    compr.linhax
    , dist.bar.x
    , rel.linha.x[,2]
    , h.rel.recep.x
  )
}
