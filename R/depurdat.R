#depurdat.R.- Depuración de series climatológicas y relleno de lagunas.

#depudm.- Depuración de datos mensuales, con opción de visualización.
depudm <- function(varcli, anyi, anyf, nm=12, wa=100, dz.max=2, difumb=0.05,
  leer=TRUE, a=0, b=1, wz=0.001, sqrtrans=FALSE, ttip=3, refglob=FALSE, ndec=1,
  pval=0.05, graf=FALSE, auto=FALSE, verb=TRUE) {

  mes3 <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
# Lectura (opcional) de los datos desde ficheros varcli_anyi-anyf.dat y
# varcli_anyi-anyf.est, a una matriz dat (datos) y una tabla est.c
# (coordenadas, indicativos y nombres de las estaciones). Las constantes
# a y b son para transformar los datos que se leen mediante y=a+bx :
  if(leer) leerdm(varcli,anyi,anyf,nm=nm,b=b,a=a)
# Obtención de la matriz de pesos (est.w).  wa: 0(no ponderar) o >>5 para
# depuración; <5 para relleno de lagunas (si los datos son de alta calidad).
  matpesos(wa=wa, wz=wz, verb=verb) #matriz de pesos (est.w)
  if(auto==TRUE || graf==FALSE) {
    sink("climatol.log",append=TRUE)
    cat("\n=========== HOMOGENIZATION OF ",varcli,", ",anyi,"-",anyf,". (",
      date(),")\n",sep="")
    cat("depudm: nm=",nm,", wa=",wa,", dz.max=",dz.max,", difumb=",difumb,
      ", leer=",leer,", a=",a,", b=",b,",\n",sep="")
    cat("  wz=",wz,", sqrtrans=",sqrtrans,", ttip=",ttip,", refglob=",refglob,
      ", ndec=",ndec,", pval=",pval,",\n",sep="")
    cat("  graf=",graf,", auto=",auto,", verb=",verb,"\n\n",sep="")
    sink()
  }
# Transformación raíz cuadrada? Aplicar a difumb si >1:
  if(sqrtrans && difumb>1) difumb <- sqrt(difumb)
# proceso de depuración iterativo, por meses (u otras series, si nm!=12):
  me <- 1
  while(me <= nm) {
    if(graf && !refglob) { #posibilidad de saltar a otro mes
      cat("Month ",me," ? (1-",nm,"; RETURN=confirm; f=end)\n",sep="")
      zz <- scan(what=character(),n=1,quiet=TRUE)
      if(length(zz)==0) zz <- "0"
      if(zz=="f") break
      if(substr(zz,1,1)!="1") zz <- "0"
      me2 <- as.integer(zz)
      if(me2 > 0 && me2 <= nm) me <- me2
    }
    if(refglob) { #todos los meses a la vez
      dat.d <<- dat
      attr(dat.d,"dim") <<- c(na*nm,ne)
    }
    else {
      cat("\n================== MONTH",me,"\n")
      dat.d <<- dat[me,,] #datos de 1 mes (dat.d)
    }
    #transformación raíz cuadrada opcional (sólo valores > 1):
    if(sqrtrans) {
      for(i in 1:na) {
        for(j in 1:ne) {
          if(!is.na(dat.d[i,j]) && dat.d[i,j]>1) dat.d[i,j] <<- sqrt(dat.d[i,j])
        }
      }
    }
    
    dat.na <- is.na(dat.d) #índice de datos ausentes
    tipif(ttip=ttip) #tipificación (dat.m, [dat.s,] dat.z)
    dat.m0 <- dat.m #copia de las medias
    #proceso iterativo:  
    ite <- 0
    repeat { #iterar hasta estabilizar las medias:
      ite <- ite+1
      cat("Iteration",ite,"...  ")
      datest(verb=verb) #obtención de series estimadas (dat.e)
      if(auto) { #eliminación automática de los datos anómalos:
        elim <- abs(dat.z-dat.e) > dz.max #datos a eliminar
        elim[is.na(elim)] <- FALSE #eliminar los molestos NA
        nelim <- length(elim[elim==TRUE]) #no. de datos a eliminar
        destipif(ttip=ttip) #destipificación de dat.e
        if(nelim>0) {
          sink("climatol.log",append=TRUE)
          for(i in 1:ne) { #listado de los datos que se van a eliminar
            for(j in 1:na) {
              if(elim[j,i]) cat(as.character(est.c[i,4]),
                anyi+j-1,me,":",round(dat.d[j,i],ndec),
                "->",round(dat.e[j,i],ndec),"\n")
            }
          }
          sink()
          dat.na[elim] <- TRUE #actualización índice de datos ausentes
          cat(nelim,"data rejected.   ")
        }
      }
      else destipif(ttip=ttip) #destipificación de dat.e
      dat.d[dat.na] <<- dat.e[dat.na] #relleno de los datos ausentes
      tipif(ttip=ttip) #actualizar dat.m y dat.z (y dat.s, si ttip=3)
      maxdif <- max(abs(dat.m-dat.m0))
      cat("Max. average difference: ",maxdif,"\n")
      if(maxdif<=difumb) break
      dat.m0 <- dat.m #copia de las medias
    }
    #grabar los datos depurados de errores puntuales y con lagunas rellenadas:
    arch <- paste(varcli,"_",substr(as.character(anyi),3,4),"-",substr(as.character(anyf),3,4),".dep",sep="") # nombre del archivo
    #deshacer transformación raíz cuadrada?:
    if(sqrtrans) dat.d[dat.d>1] <<- dat.d[dat.d>1] * dat.d[dat.d>1]
    if(refglob) {
      dat <<- round(dat.d,ndec)
      attr(dat,"dim") <<- c(nm,na,ne)
      for(j in 1:nm) write(dat[j,,],arch,ncolumns=10,append=(j!=1))
    }
    else write(round(dat.d,ndec),arch,ncolumns=10,append=(me!=1))
    #no usar datos rellenados para gráficas ni pruebas de homogeneidad:
    dat.d[dat.na] <<- NA #no representar los datos rellenados!
    dat.z[dat.na] <<- NA # "      "       "    "       "
    #gráficas interactivas de anomalías y datos (incopatible con refglob!):
    if(graf && !refglob) {
      #etiqueta del mes (si nm!=12, poner sólo el número):
      if(nm==12) labm <- mes3[me] else labm <- me
      grafanom(me=me,labm=labm,ttip=ttip,ndec=ndec,sqrtrans=sqrtrans,pval=pval)
    }
    else if(pval>0) { #aplicar los tests de homogeneidad de grafanom()
      sink("climatol.log",append=TRUE)
      i <- 1
      while(i>0 && i<=ne) { #para cada estación
        switch(ttip,   #tipificar los datos estimados
          y <- dat.e[,i] - dat.m[i],
          if(dat.m[i]>1) y <- dat.e[,i] / dat.m[i]
            else y <- dat.e[,i],
          y <- (dat.e[,i]-dat.m[i]) / dat.s[i],
          y <- dat.e[,i]
        )
        y <- dat.z[,i] - y #ahora como anomalía (tipificada)
        indme <- paste(est.c[i,4],sprintf("%02d",as.integer(me))) #indicativo y mes
        #prueba móvil de la t (ventana de 10 términos):
        movttest(y,indme,pval=pval,verb=verb)
        #prueba móvil de la t (ventana de 10 términos):
        movttest(y,indme,nterm=20,pval=pval,verb=verb)
        tpv <- cor.test(y,anyi:anyf)$p.value
        if(tpv<=pval) cat(indme,"Trend p-value:",round(tpv,3),"\n")
        i <- i+1
      }
      sink()
    }
    me <- me+1
  }
}

# Lectura de los datos mensuales, y transformación a+bx:
leerdm <- function(varcli,anyi,anyf,nm=12,b=1,a=0) {
  varcli <<- varcli
  anyi <<- anyi
  anyf <<- anyf
  na <<- anyf-anyi+1 #no. de años
  arch <- paste(varcli,"_",substr(as.character(anyi),3,4),"-",substr(as.character(anyf),3,4),".",sep="") #raíz de los nombres de archivo 
  arche <- paste(arch,"est",sep="") #nombre del archivo de estaciones
  est.c <<- read.table(arche) #leer coordenadas y nombres de las estaciones
  ne <<- nrow(est.c) #no. de estaciones
  nd <<- na*nm #no. de datos por estación
  archd <- paste(arch,"dat",sep="") #nombre del archivo de datos
  dat <<- scan(archd) #lectura de los datos
  dat <<- a+dat*b #transformación lineal (ej.: si estaban grabados en décimas)
 #attr(dat,"dim") <<- c(nd,ne) #conversión de vector a matriz
  attr(dat,"dim") <<- c(nm,na,ne) #conversión de vector a matriz
}

#cálculo de la matriz de pesos:
matpesos <- function(wa=5,wz=.001,verb=TRUE) {
#wa: parámetro de suavización de la función de peso (> 0). Asignarle
# un valor cero (o negativo) para que no aplicar ponderaciones (todos los
# pesos valdrán la unidad).
#wz: factor de escala de z. El valor por defecto es apropiado si z se da en m
# y x,y en km. También sirve para sobreponderar z, o para hallar las distancias
# únicamente en el plano horizontal (wz=0).
  ne <- nrow(est.c) #no. de estaciones
  est.w <<- mat.or.vec(ne,ne)
  if(verb) cat("Computing inter-station weight matrix:\n")
  for(i in 1:(ne-1)) {
    if(verb) cat(" ",i)
    for(j in (i+1):ne) {
      if(wa<=0) est.w[i,j] <- 1
      else {
        dx <- est.c[i,1]-est.c[j,1]
        dy <- est.c[i,2]-est.c[j,2]
        dz <- (est.c[i,3]-est.c[j,3])*wz
        d2 <- dx*dx+dy*dy+dz*dz #distancia cuadrática
        est.w[i,j] <<- wa/(wa+d2)
      }
      est.w[j,i] <<- est.w[i,j]
    }
  }
  if(verb) cat("\n")
}

#tipificación de los valores (genera dat.m, dat.s y dat.z):
tipif <- function(ttip=3) {
  dat.m <<- apply(dat.d,2,mean,na.rm=TRUE) #medias
  if(ttip==3) dat.s <<- apply(dat.d,2,sd,na.rm=TRUE) #desviaciones típicas
  nd <- attr(dat.d,"dim")[1] #no. de datos
  ne <- attr(dat.d,"dim")[2] #no. de estaciones
  dat.z <<- mat.or.vec(nd,ne)
#tipos de tipificación: ttip=1(diferencias), 2(proporciones), 3(estandari-
#  zación); 0(ninguna).
  switch(ttip,
    for(i in 1:ne) dat.z[,i] <<- dat.d[,i] - dat.m[i],
    #en proporciones, no dividir por menos de 1:
    for(i in 1:ne) { if(dat.m[i]>1) dat.z[,i] <<- dat.d[,i] / dat.m[i]
      else dat.z[,i] <<- dat.d[,i] },
    for(i in 1:ne) dat.z[,i] <<- (dat.d[,i]-dat.m[i]) / dat.s[i],
    dat.z <<- dat.d
  )
}

#destipificación de los valores estimados (destipifica dat.e):
destipif <- function(ttip=3) {
  if(!ttip) return()
  ne <- attr(dat.d,"dim")[2] #no. de estaciones
#tipos de tipificación: ttip=1(diferencias), 2(proporciones), 3(estandari-
#  zación); 0(ninguna).
  switch(ttip,
    for(i in 1:ne) dat.e[,i] <<- dat.e[,i] + dat.m[i],
    for(i in 1:ne) if(dat.m[i]>1) dat.e[,i] <<- dat.e[,i] * dat.m[i],
    for(i in 1:ne) dat.e[,i] <<- dat.e[,i] * dat.s[i] + dat.m[i],
    return()
  )
}

#obtención de series estimadas en función de las vecinas:
datest <- function(verb=FALSE) {
# dat.m <- apply(dat.d,2,mean,na.rm=TRUE) #medias
# dat.s <- apply(dat.d,2,sd,na.rm=TRUE) #medias
  nd <- attr(dat.d,"dim")[1] #no. de datos
  ne <- attr(dat.d,"dim")[2] #no. de estaciones
  dat.e <<- mat.or.vec(nd,ne) #datos estimados <<matrix(nd,ne)?
  if(verb) cat("Computing estimated series for each station:\n")
  for(i in 1:ne) { #para cada estación
    if(verb) cat(" ",i)
    for(j in 1:nd) { #para cada dato
      se <- 0
      sw <- 0
      for(ir in 1:ne) { #para cada estación de referencia
        if(ir==i) next #es la misma estación
        if(is.na(dat.z[j,ir])) next #sin dato de referencia
        se <- se + dat.z[j,ir]*est.w[i,ir]
        sw <- sw + est.w[i,ir]
      }
      dat.e[j,i] <<- se / sw #dato estimado (tipificado)
    }
  }
  if(verb) cat("\n")
  #si hay NaN, convertirlos en NA. (Esto sucede si faltan datos de referencia):
  dat.e[is.nan(dat.e)] <<- NA
}

#Gráficas de anomalías y datos originales y estimados (series mensuales)
grafanom <- function(me=me, labm="", ttip=3, ndec=1, sqrtrans=FALSE, pval=.05){
  old.par <- par(no.readonly = TRUE)
  on.exit(par(old.par))
# ne <- attr(dat.d,"dim")[2] #no. de estaciones
  par(pty="s",las=1,tck=1,fg=grey(.7)) #etiquetas, grid de color gris
  i <- 1
  while(i>0 && i<=ne) {
    switch(ttip,   #tipificar los datos estimados antes de representarlos
      y <- dat.e[,i] - dat.m[i],
      if(dat.m[i]>1) y <- dat.e[,i] / dat.m[i]
        else y <- dat.e[,i],
      y <- (dat.e[,i]-dat.m[i]) / dat.s[i],
      y <- dat.e[,i]
    )
    y <- dat.z[,i] - y #ahora como anomalía (tipificada)
    ym <- mean(y, na.rm=TRUE)
    y <- y-ym #centramos las anomalías respecto a su propia media
    tit <- paste(i,": ",est.c[i,4],"-",est.c[i,5],sep="")
    if(ttip==3) tit <- paste(tit,"  (s=",round(dat.s[i],ndec+1),")",sep="")
    tit <- paste(tit,"       ",varcli," (",labm,")",sep="")
    plot(anyi:anyf,y,type="h",lty=1,lwd=3,main=tit,ylim=c(-5,5),xlab="Years",ylab="Centered normalized anomalies (original - estimated)",col=hsv(.7,1,.9))
    grid()
    abline(1,0,col=hsv(.3,1,.8),lty=2); abline(-1,0,col=hsv(.3,1,.8),lty=2)
    abline(2,0,col="red"); abline(-2,0,col="red")
    #pruebas de homogeneidad (saltos en la media y tendencias):
    if(pval>0) {
      indme <- paste(est.c[i,4],sprintf("%02d",as.integer(me))) #indicativo y mes
      #prueba móvil de la t (ventana de 10 términos):
      movttest(y,indme,pval=pval)
      lines(anyi:anyf+.5,pv*2/pval,col="red") #p-valores de la prueba anterior
      #prueba móvil de la t (ventana de 10 términos):
      movttest(y,indme,nterm=20,pval=pval)
      lines(anyi:anyf+.5,pv*2/pval,col=hsv(.3,1,.8)) #p-valores prueba anterior
      mtext("p-val",4,las=1,line=1,adj=0,at=5,col="red")
      mtext(pval,4,las=1,line=1,adj=0,at=2,col="red")
      mtext(0,4,las=1,line=1,adj=0,at=0,col="red")
      tpv <- cor.test(y,anyi:anyf)$p.value
      cat(indme,"Trend p-value:",round(tpv,3),"\n")
      if(tpv<pval) {
        areg <- lm(y~c(anyi:anyf))
        abline(areg,col="blue")
      }
    }
    cat("Station (1-",ne,"; RETURN=next, v=view data, g=graphic save, f=end)\n",sep="")
    zz <- scan(what=character(),n=1,quiet=TRUE)
    if(length(zz)==0) zz <- "0"
    if(zz=="f") break
    if(zz=="g") grabeps()  #grabar gráfico en EPS
    else if(zz=="v") {
      yrotu <- paste(varcli,"data")
      y <- dat.e[,i]; yd <- dat.d[,i]
      if(sqrtrans) y[y>1] <- y[y>1] * y[y>1]
      yx <- max(max(y,na.rm=TRUE),max(yd,na.rm=TRUE))
      yn <- min(min(y,na.rm=TRUE),min(yd,na.rm=TRUE))
      yx <- yx + 0.1*(yx-yn)
      plot(anyi:anyf,y,ylim=c(yn,yx),type="l",lty=1,main=tit,xlab="Years",ylab=yrotu,col=hsv(.6))
      points(anyi:anyf,yd,type="b",col="red")
#     points(anyi:anyf,dat.d[,i],col="red")
#     segments(anyi:anyf,dat.d[,i],anyi:anyf,dat.e[,i],col="red")
      par(fg="black")
      legend(anyi,yx,legend=c("Originals","Estimated"),lty=c(1,1),col=c("red",hsv(.6)),bg="white",horiz=TRUE)
      par(fg=grey(.7))
      cat("RETURN=next station; g=graph save; f=end)\n",sep="")
      zz <- scan(what=character(),n=1,quiet=TRUE)
      if(length(zz)>0) {
        if(zz=="f") break
        if(zz=="g") grabeps()  #grabar gráfico en EPS
        z0 <- substr(zz,1,1)
        if(z0>="0" && z0<="9") {
          qi <- as.numeric(zz)
          if(qi>0 && qi<=ne) i <- qi-1
        }
      }
    }
    else {
      z0 <- substr(zz,1,1)
      if(z0>="0" && z0<="9") {
        qi <- as.numeric(zz)
        if(qi>0 && qi<=ne) i <- qi-1
      }
    }
    i <- i+1
  }
  invisible()   #reset old.par (restablecemos parámetros gráficos anteriores)
}

#depstat.- Cálculo de promedios depurados del periodo deseado (-> *.med)
depstat <- function(varcli, anyi,anyf, anyip=anyi,anyfp=anyf, nm=12, ndec=1,
  vala=2) {
  if(anyip<anyi) {
    cat("First year former than initial year of data!\n")
    return()
  }
  if(anyfp>anyf) {
    cat("Last year later than final year of data!\n")
    return()
  }
  na <- anyf-anyi+1 #no. de años
  arch <- paste(varcli,"_",substr(as.character(anyi),3,4),"-",substr(as.character(anyf),3,4),".",sep="") #raíz de los nombres de archivo 
  arche <- paste(arch,"est",sep="") #nombre del archivo de estaciones
  est.c <<- read.table(arche) #leer coordenadas y nombres de las estaciones
  ne <<- nrow(est.c) #no. de estaciones
 #nd <<- na*nm #no. de datos por estación
  archd <- paste(arch,"dep",sep="") #nombre del archivo de datos depurados
  dep <<- scan(archd) #lectura de los datos depurados
  attr(dep,"dim") <<- c(na,ne,nm) #conversión a matriz tridimensional
  ars <- paste(varcli,"_",substr(as.character(anyip),3,4),"-",substr(as.character(anyfp),3,4),".med",sep="") #nombre del archivo de salida
  sink(ars) #abrir archivo de salida
  cat("Average values of ",varcli," (",anyip,"-",anyfp,")\n\n",sep="")
  cat("Station    Jan    Feb    Mar    Apr    May    Jun    Jul    Aug    Sep    Oct    Nov    Dec Annual\n")
  cat("--------------------------------------------------------------------------------------------------\n")
  for(i in 1:ne) {
    #datos mensuales:
    dmens <- round(apply(dep[(anyip-anyi+1):(anyfp-anyi+1),i,],2,mean),ndec)
    #valor anual:
    switch(vala,
      vanu <- sum(dmens),
      vanu <- round(mean(dmens),ndec),
      vanu <- max(dmens),
      vanu <- min(dmens),
      vanu <- round(mean(dmens),ndec)
    )
    #imprimir los datos:
    cat(sprintf("%-8s",as.character(est.c[i,4]))) #indic. de la estación
    cat(formatC(dmens,ndec,6,"f"),formatC(vanu,ndec,6,"f"),"\n")
  }
  sink()  #cerrar archivo de salida
}

#Prueba de la t para una ventana móvil
movttest <- function(x,indme,nterm=10,pval=.05,verb=TRUE) {
  n2 <- nterm/2; n1 <- n2-1; n3 <- 2*n2-1
  pmin <- 1; imin <- 0
  pv <<- rep(NA,na) #vector para guardar los p-valores
#if(verb) cat("Prueba de la t para una ventana móvil de",nterm,"términos:\n")
  for(i in 1:(na-n3)) {
    n=0
    for(j in 0:n1) if(!is.na(x[i+j])) n <- n+1
    if(n<5) next #faltan datos
    n=0
    for(j in n2:n3) if(!is.na(x[i+j])) n <- n+1
    if(n<5) next #faltan datos
    pv[i+n1] <<- t.test(x[i:(i+n1)],x[(i+n2):(i+n3)])$p.value
#if(verb && pv[i+n1]<pval) cat(anyi+i+n1-1,"-",anyi+i+n1,": ",round(pv[i+n1],3),"\n",sep="")
    if(pmin>pv[i+n1]) { pmin=pv[i+n1]; imin=i }
  }
  if(verb && pmin<=pval) {
    cat(indme," Minimum p-value (",nterm," t.): ",sep="")
    if(pmin==1) cat("Not enough correlative data\n")
    else cat(round(pmin,3)," (",anyi+imin+n1-1,"-",anyi+imin+n1,")\n",sep="")
  }
}

#Grabar el gráfico en eps
grabeps <- function() {
  if(dev.next()==dev.cur()) {
    ar <- system("date +R%y%m%d%H%M-%%02d.eps",TRUE)
    postscript(ar,onefile=FALSE,horizontal=FALSE,height=8,width=8,pointsize=12)
    dev.set(dev.prev())
  }
  dev.copy(which=dev.next())
  dev.set(dev.prev())
}

