#depurdat.R.- Depuración y homogeneización de series climatológicas.
#(Comments in Spanish, sorry!)

#depudm.- Depuración de datos mensuales (o de otra periodicidad).
depudm <- function(varcli, anyi, anyf, nm=12, a=0, b=1,  wz=.001, deg=FALSE,
  std=3, nref=10, wa=10000, dz.max=4, mxdif=.005, xf=0, nmd=5, rtrans=0, ndec=1,
  graf=2, aref=FALSE, maxite=50, vmin=NA, vmax=NA) {
  #deg: indica si las coordenadas están en grados
  #nmd: no. mínimo de datos para poder seguir. (Se borrará la estación que
  #tenga menos datos). El mínimo absoluto es 5 términos:
  nmd <- max(5,nmd)
  #nref: limitar el no. de datos de referencia
  #xf: vector de fechas (opción para datos diarios)
  if(class(xf)=="Date") fechas <- TRUE else fechas <- FALSE
  #aref: autorreferencia (pasado por homogen, para la última fase).
  #maxite: limitar el no. máximo de iteraciones

  #Obtener matrices de distancias (est.d), pesos (est.w) y proximidades (est.p):
  matpesos(wa=wa, wz=wz, deg=deg)

  #Generar las matrices de datos estimados y calculados:
  dat.z <<- matrix(NA,nd,ne) #datos observados (estandarizados)
  dat.e <<- matrix(NA,nd,ne) #datos estimados (estandarizados)
  dat.c <<- matrix(NA,nd,ne) #datos calculados (estimados, sin estandarizar)
  oneref <<- matrix(FALSE,nd,ne) #sólo 1 referencia?
  anom <<- matrix(NA,nd,ne) #anomalías
  sanom <<- matrix(NA,nd,ne) #anomalías estandarizadas
  mindist <<- matrix(NA,nd,ne) #distancias mínimas
  used <<- matrix(FALSE,ne,ne) #flags de estaciones usadas

  dat.d <<- dat #copia de trabajo de los datos a procesar
  #comprobar si hay algún término con pocos o ningún dato:
  numdat <- apply(!is.na(dat.d),1,sum)
  nmin=min(numdat)
  if(!nmin) {
    cat("Some terms have become void of data!:\n")
    stop("Cannot continue! Shorten the study period, add series with data in the empty terms, or loose the outlier rejection.")
  }
  #comprobar y eliminar las estaciones con menos de nmd datos:
  numdat <- apply(!is.na(dat.d),2,sum)
  nmin=min(numdat)
  if(nmin<nmd) {
    cat("There are stations with less than ",nmd," data:\n",sep="")
    cat("Stations...:",formatC(which(numdat<nmd),0,4),"\n")
    cat("Nr. of data:",formatC(numdat[numdat<nmd],0,4),"\n")
    cat("These stations will be deleted in order to proceed:\n")
    for(idel in which(numdat<nmd)) {
      z <- paste(est.c[idel,4],"(",idel,") ",est.c[idel,5],sep="")
      cat(z,"  DELETED\n")
    }
    #eliminar las estaciones con pocos datos:
    estelim(which(numdat<nmd),std=std,wa=wa,wz=wz,deg=deg)
    dat.d <<- dat #actualizar los datos a procesar
  }
  dat.na <<- is.na(dat.d) #índice de datos ausentes
  fac <- nd/(1-apply(dat.na,2,sum)) #factor acelerador de la convergencia
  fac[fac>5] <- 5 #valor máximo del factor (para evitar divergencias)
  #primera estima de medias y desviaciones típicas, o conservar anteriores:
  if(exists("dat.m0")) { #conservar las anteriores si ya existían
    dat.m <<- dat.m0
    if(std==3) dat.s <<- dat.s0
  }
  else { #estimar medias por diferencias o proporciones
    datmed <<- apply(dat.d,1,mean,na.rm=TRUE) #serie media global
    refmed <<- mean(datmed) #media global de referencia
    dat.m <<- apply(dat.d,2,mean,na.rm=TRUE) #medias de partida
    if(std==3) {
      refstd <<- sd(datmed) #desv. típica global de referencia
      dat.s <<- apply(dat.d,2,sd,na.rm=TRUE) #desv. típ. de partida
    }
    switch(std,
      for(i in 1:ne) dat.m[i] <<- dat.m[i] + refmed - mean(datmed[!is.na(dat.d[,i])]),
      for(i in 1:ne) dat.m[i] <<- dat.m[i] * refmed / mean(datmed[!is.na(dat.d[,i])]),
      {
        for(i in 1:ne) dat.m[i] <<- dat.m[i] + refmed - mean(datmed[!is.na(dat.d[,i])])
        for(i in 1:ne) dat.s[i] <<- dat.s[i] + refstd - sd(datmed[!is.na(dat.d[,i])])
      },
      for(i in 1:ne) dat.m[i] <<- dat.m[i] + refmed - mean(datmed[!is.na(dat.d[,i])])
    )
    dat.m0 <<- dat.m #copia de las medias
  }
  #estandarizar los datos (obtener dat.z únicamente):
  stdrz(std=std,mscalc=FALSE)
  if(std==3) dat.s0 <<- dat.s #copia de las desv. típicas
  #proceso iterativo:  
  ite <- 0
  cat("Iterative computation of missing data, with (optional) outlier removal:\n")
  cat("Max. data change (station)\n")
  repeat { #iterar hasta estabilizar los datos estimados y las medias:
    ite <- ite+1
    #cálculo series estimadas estandarizadas (dat.e, dat.c, mindist):
    #(aquí forzamos aref=FALSE para que no haga autocorrección durante el
    #cómputo iterativo de la media y la desviación típica)
    datest(std=std,nr=nref,aref=FALSE)
    #eliminación automática de los datos anómalos:
    anom <<- dat.z-dat.e #anomalías
    anom[dat.na] <<- NA  #no arrastrar anomalías de datos rellenados!
    #estandarizar las anomalías:
    anomm <- apply(anom,2,mean,na.rm=TRUE) #anomalías medias
    anoms <- apply(anom,2,sd,na.rm=TRUE) #desv. típicas de las anomalías
    for(j in 1:ne) sanom[,j] <<- (anom[,j]-anomm[j])/anoms[j] #an. estand.
    elim <- abs(sanom)>dz.max #datos a eliminar
    elim[is.na(elim)] <- FALSE #eliminar los molestos NA
    nelim <- sum(elim) #no. de datos a eliminar
    if(nelim>0) { #eliminar los datos originales anómalos
      #listado de los datos a eliminar:
      for(i in 1:ne) {
        for(j in 1:nd) if(elim[j,i] & !is.na(oneref[j,i])) {
          outan[j,iest[i]] <<- sanom[j,i] #guardar la anomalía del outlier
          do <- dat.d[j,i] #dato original
          dc <- dat.c[j,i] #dato calculado
          if(rtrans>0) {
            if(do>1) do <- do^rtrans
            if(dc>1) dc <- dc^rtrans
          }
          cat(as.character(est.c[i,4]),"(",i,") ",sep="")
          if(nm>1) {
            me <- j%%nm; if(me==0) me <- 12
            cat(anyi+j%/%nm,me)
          }
          else if(fechas) cat(format(xf[j]))
          else cat(j)
          cat(": ",do," -> ",round(dc,ndec)," (stan=",round(sanom[j,i],2),")",sep="")
          if(oneref[j,i]) { #no eliminar si sólo tenían una referencia!
            cat(" Only 1 reference! (Unchanged)")
            elim[j,i] <- FALSE
            oneref[j,i] <<- NA #evitar repetir este mismo mensaje
          }
          cat("\n")
        }
      }
      dat[elim] <<- NA #eliminación de los datos anómalos
      dat.na[elim] <<- TRUE #actualización índice de datos ausentes
    }
    dat.d[dat.na] <<- dat.c[dat.na] #relleno de los datos ausentes
    if(ite>1) {
      maxddif <- max(abs(dat.d-dat.d0),na.rm=TRUE) #máx. dat. dif.
      kmaxdif <- ceiling(which.max(abs(dat.d-dat.d0))/nd) #estación máx. dif.
    }
    dat.d0 <- dat.d #copia de los datos
    stdrz(std=std) #actualizar dat.m[, dat.s] y dat.z
    if(ite>1) {
      cat(round(maxddif,ndec+2)," (",est.c[kmaxdif,4],")\n",sep="")
      if(maxddif<=mxdif || ite==maxite) {
        if(ite==maxite) cat("\nAverage computation skipped after",ite,"iterations\n")
        else cat("\n")
        break
      }
    }
  }
  dat.m0 <<- dat.m #copia de las medias
  if(std==3) dat.s0 <<- dat.s #copia de las desv. típicas
  #restablecer valores de oneref:
  oneref[is.na(oneref)] <<- TRUE
  #si aref=TRUE, repetir el relleno de lagunas con autocorrección de las
  #series fragmentadas:
  if(aref==TRUE) {
    datest(std=std,nr=nref,aref=TRUE)
    dat.d[dat.na] <<- dat.c[dat.na] #relleno de los datos ausentes
    stdrz(std=std) #actualizar dat.m[, dat.s] y dat.z
  }
  #grabar los datos depurados de errores puntuales y con lagunas rellenadas:
  arch <- paste(varcli,"_",anyi,"-",anyf,".dah",sep="") #nombre archivo
  z <- dat.d
  if(rtrans>0) z[z>1] <- z[z>1]^rtrans
  #corregir valores fuera de límites?:
  if(!is.na(vmin)) z[z<vmin] <- vmin
  if(!is.na(vmax)) z[z<vmax] <- vmax
  if(nm>0) ncols <- nm else ncols <- 10
  write(round(z,ndec),arch,ncolumns=ncols)
  #calcular los valores finales de las anomalías y sus estandarizaciones:
  anom <<- dat.z-dat.e #anomalías
  anom[dat.na] <<- NA  #no arrastrar anomalías de datos rellenados!
  anomm <- apply(anom,2,mean,na.rm=TRUE) #anomalías medias
  anoms <- apply(anom,2,sd,na.rm=TRUE) #desv. típicas de las anomalías
  for(j in 1:ne) sanom[,j] <<- (anom[,j]-anomm[j])/anoms[j] #an. estand.
}

# Lectura de los datos mensuales (con transformación a+bx opcional):
leerdm <- function(varcli,anyi,anyf,b=1,a=0,na.strings="NA") {
  arch <- paste(varcli,"_",anyi,"-",anyf,sep="") #raíz nombres de archivo 
  arche <- paste(arch,".est",sep="") #nombre del archivo de estaciones
  #leer coordenadas y nombres de las estaciones:
  est.c <<- read.table(arche,colClasses=c("numeric","numeric","numeric","character","character"))
  ne <<- nrow(est.c) #no. de estaciones
  archd <- paste(arch,".dat",sep="") #nombre del archivo de datos
  dat <<- scan(archd,na.strings=na.strings) #lectura de los datos
  dat <<- a+dat*b #transformación lineal (ej.: si estaban grabados en décimas)
  numdat <- length(dat) #no. de datos leídos
  nd <<- numdat/ne #no. de datos por estación
  if(nd-floor(nd)>1e-16) {
    cat(ne,"stations read from",arche,"\n")
    cat(numdat,"data read from",archd,"\n")
    stop("The number of data is not multiple of the number of stations!")
  }
  dim(dat) <<- c(nd,ne) #conversión de vector a matriz
}

#cálculo de las matrices de distancias, pesos, y rangos de proximidad:
matpesos <- function(wa=100, wz=.001, deg=FALSE) {
#wz: factor de escala de z. El valor por defecto es apropiado si z se da en m
# y x,y en km. También sirve para sobreponderar z, o para hallar las distancias
# únicamente en el plano horizontal (wz=0).
  est.d <<- matrix(NA,ne,ne) #matriz de distancias
  est.w <<- matrix(NA,ne,ne) #matriz de pesos
  est.p <<- matrix(NA,ne,ne) #matriz de rangos de proximidad
  cat("Computing inter-station weights")
  if(ne>100) cat(" (",date(),")",sep="")
  else cat(":")
  cat("\n")
  for(i in 1:(ne-1)) {
    cat(" ",i)
    for(j in (i+1):ne) {
      dx <- est.c[i,1]-est.c[j,1]
      dy <- est.c[i,2]-est.c[j,2]
      if(deg) {  #convertir grados a km
        dx <- dx*111*cos((est.c[i,2]+est.c[j,2])*pi/360)
        dy <- dy*111
      }
      dz <- (est.c[i,3]-est.c[j,3])*wz
      d2 <- dx*dx+dy*dy+dz*dz #distancia cuadrática
      #evitar autocorrección en series con las mismas coordenadas:
      if(d2==0 && iest[i]!=iest[j]) est.d[i,j] <<- 0.01 #evitar autocorrección 
      else est.d[i,j] <<- sqrt(d2) #distancia
      est.d[j,i] <<- est.d[i,j]  #matriz simétrica
      if(wa<=0) est.w[i,j] <<- 1 #todas las estaciones pesan lo mismo
      else est.w[i,j] <<- wa/(wa+d2) #peso
      est.w[j,i] <<- est.w[i,j]  #matriz simétrica
    }
  }
  cat("\nComputing proximity ranks...\n")
  #matriz de rangos de proximidad:
  for(i in 1:ne) est.p[i,] <<- order(est.d[i,])
  if(ne>100) cat("Done! (",date(),")\n",sep="")
}

#estandarización de los valores (genera dat.m, dat.s y dat.z):
stdrz <- function(std=3,mscalc=TRUE) {
  if(mscalc) { #calcular medias y, opcionalmente, desviaciones típicas
    dat.m <<- apply(dat.d,2,mean,na.rm=TRUE) #medias
    if(std==3) dat.s <<- apply(dat.d,2,sd,na.rm=TRUE) #desviaciones típicas
  }
  #tipos de estandarización: std=1(diferencias), 2(proporciones),
  #  3(estandarización); 0(ninguna).
  switch(std,
    for(i in 1:ne) dat.z[,i] <<- dat.d[,i] - dat.m[i],
    #en proporciones, no dividir por menos de 1:
    for(i in 1:ne) { if(dat.m[i]>1) dat.z[,i] <<- dat.d[,i] / dat.m[i]
      else dat.z[,i] <<- dat.d[,i] },
    for(i in 1:ne) dat.z[,i] <<- (dat.d[,i]-dat.m[i]) / dat.s[i],
    dat.z <<- dat.d
  )
}

#datest.- Obtención de series estimadas (dat.e, dat.c) en función de las
#vecinas. También actualiza used[ne,ne] y mindist[ne,ne].
datest <- function(std=3,nr=10,aref=FALSE) {
  for(i in 1:ne) { #para cada estación
    for(j in 1:nd) { #para cada dato
      se <- 0
      sw <- 0
      nref <- 0
      for(ir in 1:ne) { #para cada estación de referencia
        kr <- est.p[i,ir]
        if(kr==i) next #es la misma estación
        if(dat.na[j,kr]) next #sin dato original de referencia
        nref <- nref+1
        used[i,kr] <<- TRUE #marca de estación usada
        #distancia mínima (distancia al dato más próximo, de 1 a 1000 km):
        if(nref==1) mindist[j,i] <<- max(est.d[i,kr],1)
        se <- se + est.w[i,kr] * dat.z[j,kr]
        sw <- sw + est.w[i,kr]
        #si no. máx. de referencias, o autoreferencia, terminar:
        if(nref==nr | (aref & !est.d[i,kr])) break
      }
      if(!nref) dat.e[j,i] <<- dat.z[j,i] #sin referencia! conservar original
      else {
        #si sólo hay una referencia, marcar para no corregir la serie:
        if(nref==1 & !is.na(oneref[j,i])) oneref[j,i] <<- TRUE
        #no permitir datos negativos si std=2 (precipitación, etc):
        if(std==2 & se<0) se <- 0
        dat.e[j,i] <<- se / sw #dato estimado (estandarizado)
      }
    }
  }
  #si hay NaN, convertirlos en NA. (Sucede a veces con std=2):
  n <- sum(is.nan(dat.e))
  if(n>0) {
    cat(n,"NaN's in dat.e ! (changing them to NA's...)\n")
    dat.e[is.nan(dat.e)] <<- NA
  }
  #valores calculados por desestandarización de dat.e:
  if(!std) dat.c <<- dat.e
  else {
    switch(std,
      for(i in 1:ne) dat.c[,i] <<- dat.e[,i] + dat.m[i],
      for(i in 1:ne) if(dat.m[i]>1) dat.c[,i] <<- dat.e[,i] * dat.m[i]
                     else dat.c[,i] <<- dat.e[,i],
      for(i in 1:ne) dat.c[,i] <<- dat.e[,i] * dat.s[i] + dat.m[i],
      return()
    )
  }
  #no permitir datos negativos si std=2 (precipitación, etc):
  if(std==2) dat.c[dat.c<0] <- 0
}

#dahstat.- Estadísticas de datos homogeneizados.
dahstat <- function(varcli, anyi,anyf, anyip=anyi, anyfp=anyf, nm=12, ndec=1,
  vala=2, mnpd=0, mnsh=0, out="med", prob=.5, func=FALSE, pernum=100,
  eshcol=4, sep=" ", eol="\n") {
#out="med"(valores medios), "mdn"(medianas), "max"(máximos), "min"(mínimos),
# "std"(desv.típ.), "q"(cuantiles), "tnd"(tendencias)
#vala= 0(ninguno), 1(suma), 2(media), 3(máximo), 4(mínimo)
#mnpd=Mínimo porcentaje de datos originales
#mnsh=Máximo SNHT
#func=Listar sólo las estaciones funcionando al final del periodo.
#pernum=No. de años sobre los que se expresan las tendencias (100 por defecto)
#eshcol=Columnas de est.c seleccionadas para el listado
  if(anyip<anyi) stop("Asked initial year before first year of data!")
  if(anyfp>anyf) stop("Asked final year beyond last year of data!")
  estvar <- c("X","Y","Z","Code","Name","PD","io","op","SNHT")
  mes3 <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
  na <- anyf-anyi+1 #no. de años
  if(nm==1) vala <- 0 #no es necesario el valor anual
  else if(vala<1 || vala>4) vala <- 2 #valor medio en caso de vala erróneo
  funa <- c("sum","mean","max","min")[vala] #función para el valor anual
  #función elegida para el cálculo de los valores mensuales:
  fun <- c("mean","median","max","min","sd","quantile")[which(c("med","mdn","max","min","std","q","tnd")==out)]
  arch <- paste(varcli,"_",anyi,"-",anyf,sep="") #raíz nombres de archivo 
  arche <- paste(arch,".esh",sep="") #archivo de estaciones
  est.c <- read.table(arche,as.is=TRUE) #leer coordenadas y nombres estaciones
  ne <- nrow(est.c) #no. de estaciones
  archd <- paste(arch,".dah",sep="") #archivo de datos homogeneizados
  dah <- scan(archd) #lectura de los datos homogeneizados
  dim(dah) <- c(nm,na,ne) #conversión a tres dimensiones
  #si no se reconoce la opción out, terminar aquí:
  if(!length(fun)) {
    cat("Data available in dah[",nm,",",na,",",ne,"] for the ",ne," stations in est.c\n",sep="")
    cat("(No further output, since option out='",out,"' is not recognized)\n",sep="")
    dah <<- dah; est.c <<- est.c
    return(invisible(ne))
  }
  #seleccionar las estaciones con un mínimo de mnpd % de datos originales:
  esel <- est.c[,6]>=mnpd #vector de la selección
  #seleccionar las estaciones con un SNHT menor o igual a mnsh:
  if(mnsh) esel <- est.c[,9]<=mnsh #vector de la selección
  if(func) esel <- esel & as.logical(est.c[,8]) #sólo últimos fragmentos?
  if(sum(esel)==0) stop("No station selected: No output")
  dah <- dah[,,esel]
  if(nm==1) dim(dah) <- c(1,na,length(esel))
  est.c <- est.c[esel,]
  ne <- sum(esel) #no. de estaciones seleccionadas
  if(vala) { #calcular los valores anuales
    aval <- as.vector(apply(dah,2:3,funa))
    dim(dah) <- c(nm,na*ne)
    dah <- rbind(dah,aval)
    nm <- nm+1
    dim(dah) <- c(nm,na,ne)
  }
  #dimensionar valores a calcular:
  val <- matrix(NA,ne,nm)
  #calcular los valores deseados:
  if(out=="tnd") { #tendencias del periodo escogido
    x <- anyip:anyfp #variable independiente
    for(i in 1:ne) {
      if(nm==1) { #una sola subserie
        aj <- lm(dah[1,(anyip-anyi+1):(anyfp-anyi+1),i]~x) #regresión lineal
        val[i,] <- round(aj$coefficients[2]*pernum,ndec)
      }
      else {
        for(j in 1:nm) {
          aj <- lm(dah[j,(anyip-anyi+1):(anyfp-anyi+1),i]~x) #regresión lineal
          val[i,j] <-  round(aj$coefficients[2]*pernum,ndec)
        }
      }
    }
  }
  else { #aplicar la función deseada al periodo escogido
    for(i in 1:ne) {
      if(nm==1) {
        if(out=="q") val[i,]<-round(eval(call(fun,dah[,(anyip-anyi+1):(anyfp-anyi+1),i]),prob),ndec)
        else val[i,]<-round(eval(call(fun,dah[,(anyip-anyi+1):(anyfp-anyi+1),i])),ndec)
      }
      else { #datos mensuales:
        if(out=="q") val[i,]<-round(apply(dah[,(anyip-anyi+1):(anyfp-anyi+1),i],1,fun,prob),ndec)
        else val[i,]<-round(apply(dah[,(anyip-anyi+1):(anyfp-anyi+1),i],1,fun),ndec)
      }
    }
  }
  #imprimir mensaje:
  if(out=="med") cat("Mean")
  else if(out=="mdn") cat("Median")
  else if(out=="max") cat("Maximum")
  else if(out=="min") cat("Minimum")
  else if(out=="std") cat("Standard deviation")
  else if(out=="q") cat(prob,"prob. quantile")
  else if(out=="tnd") cat("Trend")
  cat(" values of ",varcli," (",anyip,"-",anyfp,")",sep="")
  if(out=="tnd") cat(", expressed in units per ",pernum," years,",sep="")
  dahs <<- data.frame(cbind(est.c[eshcol],val))
  if(nm==12) ndf <- c(estvar[eshcol],mes3)
  else if(nm==13) ndf <- c(estvar[eshcol],mes3,"Annual")
  else ndf <- c(estvar[eshcol],1:nm)
  names(dahs) <<- ndf
  #grabar los valores calculados:
  ars <- paste(varcli,"_",anyip,"-",anyfp,".",out,sep="") #nombre del archivo de salida
  if(out=="q") ars <- paste(ars,formatC(100*prob,width=2,flag="0"),sep="") #añadir percentil a la q
  write.table(dahs,ars,row.names=FALSE,sep=sep,eol=eol)
  cat("\n  written to",ars,"(and remain in memory as 'dahs')\n")
}

#snhtw.- SNHT para ventanas solapadas de 2*nt términos válidos.
snhtw <- function(x,nt=48) {
  ntt <- length(x) #no. total de términos de la serie
  ntv <- sum(!is.na(x)) #no. de términos válidos de la serie
  if(2*nt>ntv) return(c(0,0)) #no hay suficientes datos válidos para la prueba
  tV <- 0 #inicialización del tV máximo a devolver
  pk <- 0 #inicialización de la posición a devolver
  #inicialización de los límites muestrales (a1-b1, a2-b2):
  k <- 1; while(k<ntt && is.na(x[k])) k <- k+1; a1 <- k
  n<-1; while(n<nt && k<ntt) { k <- k+1; if(!is.na(x[k])) n <- n+1; }
  b1 <- k
  k <- k+1; while(k<ntt && is.na(x[k])) k <- k+1; a2 <- k
  n<-1; while(n<nt && k<ntt) { k <- k+1; if(!is.na(x[k])) n <- n+1; }
  b2 <- k
  #aplicación de SNHT a las ventanas solapadas:
  repeat {
    st <- snht(x[a1:b2])
    stx <- max(st,na.rm=TRUE)
    if(stx>tV) { tV <- stx; pk <- which.max(st)+a1-1 }
    if(b2==ntt) return(c(tV,pk))
    #desfasar las ventanas hacia adelante:
    a1 <- a2; b1 <- b2
    k <- b2+1; while(k<ntt && is.na(x[k])) k <- k+1
    if(is.na(x[k])) return(c(tV,pk)) else a2 <- k
    n<-1; while(n<nt && k<ntt) { k <- k+1; if(!is.na(x[k])) n <- n+1; }
    b2 <- k
  }
}

#Standard Normal Homogeneity Test (Alexandersson)
snht <- function(x,nmt=3) {
#nmt: no. mínimo de términos de cada muestra
  n <- length(x)
  Tsnht <- rep(NA,n)
  if(n<nmt*2) return(Tsnht) #insuficientes datos
  z <- (x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)
  yc <- FALSE #ya calculado?
  for(i in (nmt+1):(n-nmt)) { #(despreciar los primeros y últimos nmt términos)
    if(is.na(x[i]) && yc) next #test ya calculado
    n1 <- sum(!is.na(x[1:(i-1)])) #no. de términos de la muestra 1
    n2 <- sum(!is.na(x[i:n])) #no. de términos de la muestra 2
    if(n1<nmt | n2<nmt) next #al menos una muestra es demasiado pequeña
    z1 <- mean(z[1:(i-1)],na.rm=TRUE)
    z2 <- mean(z[i:n],na.rm=TRUE)
    Tsnht[i] <- n1*z1*z1 + n2*z2*z2
    if(is.na(x[i])) yc <- TRUE #marca de test ya calculado
    else if(yc) yc <- FALSE    #quitar la marca
  }
  return(Tsnht)
}

#estelim.- Eliminar estaciones (por tener muy pocos datos)
estelim <- function(x,std,wa=100,wz=.001,deg=FALSE) {
  est <- setdiff(1:ne,x) #estaciones a mantener
  est.c <<- est.c[est,]
  dat <<- dat[,est]
  if(exists("dat.d")) dat.d <<- dat.d[,est]
  if(exists("dat.z")) dat.z <<- dat.z[,est]
  if(exists("dat.e")) dat.e <<- dat.e[,est]
  if(exists("dat.c")) dat.c <<- dat.c[,est]
  if(exists("anom")) anom <<- anom[,est]
  if(exists("sanom")) sanom <<- sanom[,est]
  if(exists("mindist")) mindist <<- mindist[,est]
  if(exists("dat.m")) dat.m <<- dat.m[est]
  if(std==3 && exists("dat.s")) dat.s <<- dat.s[est]
  if(exists("iest")) iest <<- iest[est]
  if(exists("oneref")) oneref <<- oneref[,est]
  #actualizar el no. de estaciones
  ne <<- ne-length(x)
  #recalcular la matriz de pesos
  matpesos(wa=wa, wz=wz, deg=deg)
}

#Gráfico de anomalías estandarizadas
plotstan <- function(i,varcli,x,xlab,swa,lw) {
  y <- sanom[,i] #vector de anomalías de la estación
  ylab="Standardized anomalies (observed - computed)"
  tit <- paste(varcli," at ",est.c[i,4],"(",i,"), ",est.c[i,5],sep="")
  plot(x,y,type="h",lwd=lw,ylim=c(-5,5),main=tit,xlab=xlab,ylab=ylab,col=hsv(.7,1,.9))
  grid(col=grey(.4))
  abline(-3,0,lty=3,col=grey(.4)); abline(-5,0,lty=3,col=grey(.4))
  lines(x,log10(mindist[,i])-5,col=verde)
  mtext(" 1",4,las=1,adj=0,at=-5,col=verde)
  mtext(" 10",4,las=1,adj=0,at=-4,col=verde)
  mtext(" 100",4,las=1,adj=0,at=-3,col=verde)
  mtext("min.d.",4,las=1,adj=0,at=-5.4,col=verde)
  mtext(" (km)",4,las=1,adj=0,at=-2,col=verde)
}

#outrename.- Append a suffix to the output files, to avoid overwrites.
outrename <- function(varcli, anyi, anyf, suffix) {
  fbn <- paste(varcli,"_",anyi,"-",anyf,sep="") #file base name
  for(ext in c(".txt",".dah",".esh",".pdf"))
    file.rename(paste(fbn,ext,sep=""),paste(fbn,"-",suffix,ext,sep=""))
}

#cerrar.- Cerrar los archivos de salida.
cerrar <- function() {
  sink()
  graphics.off()
}

#dd2m.- Cálculo de valores mensuales a partir de datos diarios.
dd2m <- function(varcli, anyi, anyf, ini, anyip=anyi, anyfp=anyf, ndec=1,
  valm=2, nmin=15, na.strings="NA") {
  arch <- paste(varcli,"-d_",anyi,"-",anyf,sep="") #raíz nombres de archivo
  do <- scan(paste(arch,".dat",sep=""),na.strings=na.strings) #datos originales
  dd <- scan(paste(arch,".dah",sep="")) #datos homogeneizados
  est <- read.table(paste(arch,".est",sep="")) #coordenadas de las estaciones
  ne <- nrow(est) #no. de estaciones
  nd <- length(dd)/ne #no. de datos por estación
  dim(do) <- c(nd,ne)
  dim(dd) <- c(nd,ne)
  f <- as.Date(0:(nd-1),origin=ini) #fechas
  fun <- c("sum","mean","max","min")[valm] #función para el valor mensual
  na <- anyfp-anyip+1 #no. de años
  dm <- array(NA,c(12,na,ne)) #datos mensuales
  for(ie in 1:ne) { #para cada estación
    for(aa in anyip:anyfp) { #para cada año solicitado
      ka <- aa-anyip+1 #índice del año
      for(me in 1:12) { #para cada mes
        if(me<10) aamm <- paste(aa,"-0",me,sep="")
        else aamm <- paste(aa,"-",me,sep="")
        sel <- substring(f,1,7)==aamm #índice de datos a seleccionar
        d <- dd[sel,ie]
        if(length(d)<nmin) next #pocos datos!
        if(sum(!is.na(do[sel,ie]))<nmin) next #pocos datos originales!
        dm[me,ka,ie] <- eval(call(fun,d))
        dm[is.nan(dm)] <- NA #si no hay datos, poner NA
      }
    }
  }
  #grabar los datos mensuales:
  fichsal <- paste(varcli,"-m_",anyi,"-",anyf,".dah",sep="")
  write(round(dm,ndec),fichsal,ncolumns=12)
  cat("Monthly",fun,"values output to file",fichsal,"\n")
}

#homogeneización automática iterativa de un conjunto de datos
homogen <- function(varcli, anyi, anyf, nm=12, nref=10, dz.max=5,
  wd=c(0,0,100), tVt=25, tVf=.02, swa=60, snhtt=50, mxdif=.05, force=FALSE,
  a=0, b=1, wz=.001, deg=FALSE, rtrans=0, std=3, ndec=1, mndat=0, leer=TRUE,
  gp=3, na.strings="NA", nclust=100, maxite=50, ini="", vmin=NA, vmax=NA,
  verb=TRUE) {
  #nref: no. (máximo) de estaciones de referencia
  #wd: Weight distance (km; distancia a la que el peso se reduce a la mitad;
  #    wd=0: todas las estaciones de referencia pesan lo mismo).
  #tVt: tV tolerado en las series (se cortarán con un tV mayor; si tVt=0, no
  #    cortar las series, sólo rellenar las lagunas).
  #tVf: factor de tolerancia (por referencia disponible) para cortes en cadena.
  #swa: Semi-Window Amplitude (no. de términos)
  #mndat: Mínimo no. de datos para fragentar las series
  #gp: Parámetro de gráficos (0=ninguno; 1=anomalías globales e histogramas;
  #    2=id+gráficos mensuales de diagnóstico; 3=id+gráficos de medias anuales
  #    móviles y correcciones; 4=id., pero con sumas anuales móviles.
  #nclust: no. máximo de estaciones a usar en el análisis de agrupamiento
  #ini: fecha inicial (para datos diarios, en formato 'AAAA-MM-DD').
  #vmin, vmax: rango de valores permitidos en la variable climática.
  #force: forzar cortes aun con una sóla referencia.
  #------------------------------------------------------------------  
  #en caso de error, cerrar archivos de salida:
  options(error=cerrar)
  verde <<- hsv(.33,1,.6) #color muy usado
  #si son datos diarios añadir '-d' al nombre de la variable:
  if(nm<1) varcli <- paste(varcli,"-d",sep="")
  #abrir archivo de bitácora y escribir cabecera:
  archlog <- paste(varcli,"_",anyi,"-",anyf,".txt",sep="")
  sink(archlog,split=verb)
  cat("\nHOMOGEN() APPLICATION OUTPUT  (From R's contributed package 'climatol')\n\n")
  cat("\n=========== Homogenization of ",varcli,", ",anyi,"-",anyf,". (",
        date(),")\n",sep="")
  cat("homogen:")
  arg <- names(formals()) #lista de los argumentos de la función
  for(i in 1:length(arg)) {
    cat(" ",arg[i],"=",sep="")
    cat(eval(as.symbol(arg[i])),sep=",")
  }
  cat("\n\n")
  #establecer valor de mndat, si no se especificó como parámetro:
  if(mndat<=0) { 
    if(nm>5) mndat <- nm
    else if(nm<=0) mndat <- swa/2
    else mndat <- 5
  }
  #parámetros que han de tener 3 valores:
  k <- length(wd); if(k<3) wd <- c(wd,rep(wd[k],3-k))
  k <- length(nref); if(k<3) nref <- c(nref,rep(nref[k],3-k))
  k <- length(dz.max); if(k<3) dz.max <- c(dz.max,rep(dz.max[k],3-k))
  #etiquetas mensuales (de tres letras):
  mes3 <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
  #lectura inicial de datos:
  if(leer) leerdm(varcli,anyi,anyf,b=b,a=a,na.strings=na.strings)
  else {
    ne <<- nrow(est.c) #no. de estaciones
    nd <<- length(dat[,1]) #no. de datos por estación
  }
  if(!is.na(vmin)) { #comprobar si hay valores inferiores al mínimo permitido
    n <- sum(dat<vmin,na.rm=TRUE)
    if(n) {
      dat[dat<vmin] <<- vmin #corrección de los valores erróneos
      cat("Number of data forced to",vmin,":",n,"\n")
    }
  }
  if(!is.na(vmax)) { #comprobar si hay valores superiores al máximo permitido
    n <- sum(dat>vmax,na.rm=TRUE)
    if(n) {
      dat[dat>vmax] <<- vmax #corrección de los valores erróneos
      cat("Number of data forced to",vmax,":",n,"\n")
    }
  }
  if(gp>2) dat.o <- dat #copia de los datos originales
  if(nd<100) lw=3 #anchura de las barras de anomalías
  else if(nd<300) lw=2
  else lw=1
  #otras inicializaciones:
  if(swa>=nd) swa <- floor(swa/2) #evitar semiventana demasiado grande
  if(nm>0) {
    x <- (anyi*nm):(anyf*nm+nm-1)/nm #abscisas en años
    xlab <- "Years"
    na <<- anyf-anyi+1 #no. de años
    nsy <- rep(0,na)   #no. de cortes por año
    if(nd != nm*na) {
      cat(nd,"data *",ne,"stations =",nd*ne,"data, but\n")
      cat(nd,"data /",na,"years = nm =",nd/na,"data per year!\n")
      stop("Inconsistent number of data!")
    }
  }
  else {
    if(ini>0) {
      x <- as.Date(0:(nd-1),origin=ini) #vector de fechas
      x2 <- 1:nd
      xlab <- "Dates"
    }
    else {
      x <- 1:nd
      xlab <- "Index"
    }
  }
  nei <<- ne  #no. inicial de estaciones
  est.i <<- est.c #datos iniciales de las estaciones
  nsp <- rep(0,nei)  #no. de cortes de cada estación original
  iest <<- 1:ne   #índice señalando la serie original de cada subserie
  outan <<- matrix(NA,nd,ne) #anomalías de los outliers
  #gráficos iniciales:
  if(gp>0) {
    #activar salida gráfica a documento pdf, con rótulo inicial:
    pdf(paste(varcli,"_",anyi,"-",anyf,".pdf",sep=""))
    plot(-1:1,-1:1,type="n",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
    text(0,0.4,"CLIMATOL",cex=4)
    text(0,-0.45,paste("Homogenization\ngraphic output of\n",varcli,"\n",anyi,"-",anyf,sep=""),cex=3)
    #no. de datos de cada término
    numdat <- apply(!is.na(dat),1,sum)
    plot(x,numdat,type="l",ylim=c(0,ne),col="blue",xlab=xlab,ylab="Nr. of data",main=paste("Nr. of",varcli,"data in all stations"))
    abline(h=5,lty=2,col="green")
    abline(h=3,lty=2,col="red")
    grid()
    #si hay algún término sin datos no podremos continuar:
    if(!min(numdat)) {
      stop("At least one term has missing data in all stations! (See the PDF graph)\nCannot continue. (Shorten the study period, or add series with data in those void terms)")
    }
    #boxplots de los datos de cada estación:
    if(nm>1) { #boxplots para cada nm (mensuales, etc)
      dim(dat) <<- c(nm,na,ne) #dimensiones provisionales
      for(me in 1:nm) { #para cada mes
        z <- data.frame(dat[me,,])
        names(z) <- 1:ne
        #etiqueta del mes (si nm!=12, poner sólo el número):
        if(nm==12) labm <- mes3[me] else labm <- me
        labm <- paste(" (",labm,")",sep="")
        if(ne>nclust) hist(as.matrix(z),xlab=varcli,main=paste("Data values of ",varcli,labm,sep=""),col="wheat")
        else {
          boxplot(z,xlab="Stations",ylab="Values",main=paste("Data values of ",varcli,labm,sep=""),col="wheat",border=hsv(.7,1,.9))
          grid(col=grey(.4))
          abline(h=0)
        }
      }
      dim(dat) <<- c(nd,ne) #restablecer dimensiones de trabajo
    }
    else { #un sólo gráfico con los boxplots de cada estación
      z <- data.frame(dat)
      names(z) <- 1:ne
      if(ne>nclust) hist(as.matrix(z),xlab=varcli,main=paste("Data values of",varcli),col="wheat")
      else {
        boxplot(z,xlab="Stations",ylab="Values",main=paste("Data values of",varcli),col="wheat",border=hsv(.7,1,.9))
        grid(col=grey(.4))
        abline(h=0)
      }
    }
  }
  if(rtrans>0) { #transformar los datos?
    dat[!is.na(dat)&dat>1] <<- dat[!is.na(dat)&dat>1]^(1/rtrans)
    #disminuir mxdif si se transforman los datos. Mantener la precisión de
    #mxdif para una precipitación de 100 mm:
    z <- (100+mxdif)^(1/rtrans)-100^(1/rtrans)
    if(mxdif > z) mxdif <- z
  }
  if(gp>0) { #continuamos con los gráficos iniciales
    #histograma de todos los datos (distribución quasi-normal?)
    if(rtrans>0) main="Histogram of all (transformed) data"
    else main="Histogram of all data"
    hist(dat,xlab=varcli,main=main,col=hsv(.4,1,.8))
    #correlograma de series diferenciadas a nivel mensual (r <-> distancia)
    #(si hay más de nclust estaciones, sólo de una muestra aleatoria de nclust)
    if(ne>nclust) { splc <- sample(1:ne,nclust); nec <- nclust }
    else { splc <- 1:ne; nec <- ne }
    est.d <<- matrix(NA,nec,nec) #matriz de distancias
    cat("Computing inter-station distances:")
    for(i in 1:(nec-1)) {
      cat(" ",i)
      for(j in (i+1):nec) {
        dx <- est.c[splc[i],1]-est.c[splc[j],1]
        dy <- est.c[splc[i],2]-est.c[splc[j],2]
        if(deg) {  #convertir grados a km
          dx <- dx*111*cos((est.c[splc[i],2]+est.c[splc[j],2])*pi/360)
          dy <- dy*111
        }
        dz <- (est.c[splc[i],3]-est.c[splc[j],3])*wz
        d2 <- dx*dx+dy*dy+dz*dz #distancia cuadrática
        est.d[i,j] <<- sqrt(d2) #distancia
        est.d[j,i] <<- est.d[i,j]  #matriz simétrica
      }
    }
    cat("\n")
    data <- dat[,splc] #copia de los datos
    if(nm>1) { #calcular las series diferenciadas por meses
      dim(data) <- c(nm,na,nec) #dimensionamos por meses
      difd <- array(NA,c(nm,na-1,nec)) #matriz de datos diferenciados
      for(i in 1:nm) { #diferenciar las series, mes a mes
        difd[i,,] <- diff(data[i,,])
      }
      dim(difd) <- c(nd-nm,nec) #redimensionar
    }
    else difd <- diff(data) #series diferenciadas globalmente
    corm <- cor(difd,use="p") #matriz de correlaciones
    #eliminar |r|==1 (debidos a estaciones con sólo 2 datos en común):
    corm[corm==1] <- NA; corm[corm==-1] <- NA
    if(ne>nclust) main <- paste("Correlogram of first difference",nclust,"sampled series")
    else main <- "Correlogram of first difference series"
    if(rtrans>0) main <- paste(main," (root ",rtrans," transformed)",sep="")
    xd <- as.vector(est.d); y <- as.vector(corm)
    plot(xd,y,xlim=c(0,max(est.d,na.rm=TRUE)),xlab="Distance (km)",ylab="Correlation coefficient",main=main,col="blue")
    grid(col=gray(.4))
    if(ne>2) {  #dendrograma de las estaciones
      dism <- dist(corm) #matriz de disimilaridad
      #si hay NA's en la matriz de disimilaridad, no intentar clustering
      if(!sum(is.na(dism))) {
        hc <- hclust(dism,method="ward")
        if(ne>nclust) main <- paste("Dendrogram of",nclust,"sampled stations")
        else main <- "Dendrogram of station clusters"
        plot(hc,xlab="Stations",sub="",ylab="Dissimilarity",main=main)
        #clasificación de las estaciones cortando por la mitad del rango de las
        #disimilaridades (salvo que el número de clases sea superior a 9, en cuyo
        #caso se incrementará el corte en .1):
        if(ne<=nclust) { #(no hacer grupos de una muestra)
          cutlev <- mean(range(hc$height))
          repeat {
            ct <- cutree(hc,h=cutlev)
            nc <- length(levels(factor(ct)))
            if(nc<10) break
            cutlev <- cutlev + .1
          }
          if(nc>1) abline(h=cutlev,col="red",lty=2)
        }
      } else { nc <- 1; ct <- 1 } #un sólo grupo si dism tiene NA's
      if(ne>nclust) { nc <- 1; ct <- 1 } #un sólo grupo si >nclust estaciones
      #mapa de las estaciones:
      if(nc==1) { col="blue"; main=paste(varcli,"station locations") }
      else {
        col=rainbow(nc,1,.6)[ct]
        main=paste(varcli," station locations (",nc," clusters)",sep="")
      }
      if(ne>nclust) { #dibujar símbolos si hay más de nclust estaciones
        if(deg) plot(est.c[,1:2],asp=1/(cos(mean(est.c[,2])*pi/180)),xlab="Longitude (deg)",ylab="Latitude (deg)",col=col,pch=ct,main=main)
        else plot(est.c[,1:2],asp=1,xlab="X (km)",ylab="Y (km)",col=col,pch=ct,main=main)
      }
      else { #hasta nclust estaciones, poner el número
        if(deg) plot(est.c[,1:2],type="n",asp=1/(cos(mean(est.c[,2])*pi/180)),xlab="Longitude (deg)",ylab="Latitude (deg)",main=main)
        else plot(est.c[,1:2],type="n",asp=1,xlab="X (km)",ylab="Y (km)",main=main)
        text(est.c[,1:2],labels=1:ne,col=col)
      }
      grid(col=gray(.4))
    }
    rm(data,difd,corm) #borrar objetos temporales
  }
  #----------- Continuar con el proceso:
  if(gp==1) { #si sólo se deseaban los gráficos iniciales, terminar
    graphics.off() #volcar el buffer del último gráfico
    cat("\nOnly the initial exploratory graphics were demanded.\nSee them in ",varcli,"_",anyi,"-",anyf,".pdf\n",sep="")
  }
  else { #aplicación iterativa de depudm():
    if(exists("dat.m0")) rm(dat.m0,pos=1) #borrar posible dat.m0 anterior
    for (k in 1:3) { #proceso en 3 etapas: 1) cortes por ventanas móviles;
      #           2) cortes por SNHT en toda la serie; 3) relleno de lagunas
      if(tVt==0 & k<3) next #no realizar cortes, sólo rellenar lagunas
      if(k==2) {
        if(snhtt>0) tVt <- snhtt #umbral para SNHT en toda la series
        else next #no realizar cortes por SNHT
      }
      cat("\n\n============ STAGE",k,"=======================\n\n")
      #gráfico separador de niveles:
      if(gp>0) {
        plot(-1:1,-1:1,type="n",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
        text(0,0.4,paste("Stage",k),cex=4)
        if(k==1) text(0,-0.3,paste("Binary splits on",swa,"term\nstepped windows\nwith tV >",tVt,"\nand wd =",wd[k],"km"),cex=3)
        else if(k==2) text(0,-0.3,paste("Binary splits on\nwhole series\nwith SNHT >",snhtt,"\nand wd =",wd[k],"km"),cex=3)
        else text(0,-0.3,paste("Anomalies after\nmissing data\nrecalculation\nwith wd =",wd[k],"km\n( swa =",swa,")"),cex=2.5)
      }
      repeat { #repetir hasta que no se corte ninguna serie
        #empezamos por invocar depudm(), que calculará las anomalías y
        #eliminará los datos anómalos:
        depudm(varcli, anyi, anyf, nm=nm, a=a, b=b, wz=wz, deg=deg, std=std,
          nref=nref[k], wa=wd[k]*wd[k], dz.max=dz.max[k], mxdif=mxdif, xf=x,
          nmd=mndat, rtrans=rtrans, ndec=ndec, aref=(k==3 & tVt>0),
          maxite=maxite, vmin=vmin, vmax=vmax)
        #en la última etapa, limitarse al último relleno de lagunas:
        if(k>2) break
        #analizar los saltos en la media de las series, cortándolas cuando
        #el máximo tV del test supere el umbral (tVt):
        nn <- 0 #inic. no. de nuevas estaciones
        tVx <<- rep(0,ne) #máximos valores del shift test (por estación)
        kpx <<- rep(NA,ne) #posiciones de los máximos tV (por estación)
        splt <- rep(0,ne)  #tV con que se cortaron las estaciones
        modif <- FALSE #inicialización modificación series
        cat("Performing shift analysis for the",ne,"stations.")
        if(k==1) cat("  tV values:\n") else cat("  SNHT values:\n")
        for(i in 1:ne) { #análisis de saltos en la media para cada estación
          cat(" ",i,":",sep="")
          y <- sanom[,i] #anomalías estandarizadas de la estación
          if(k==1) { #análisis de saltos en ventanas móviles
            st <- snhtw(y,swa) #prueba SNHT en ventanas solapadas
            tVx[i] <<- st[1]
            kpx[i] <<- st[2]
            cat(round(tVx[i],1)) #tV's de corte
          }
          else { #análisis de saltos en toda la serie
            st <- snht(y)
            if(sum(!is.na(st))>0) {
              tVx[i] <<- max(st,na.rm=TRUE)
              kpx[i] <<- which.max(st)
            }
            cat(round(tVx[i],1)) #SNHT
          }
          if(!i%%5) cat("\n") #evitar líneas muy largas
        }
        if(i%%5) cat("\n")
        #cortar las series cuyo tV máximo supere el umbral, de mayor a menor,
        #siempre que no se hayan usado series recién cortadas con un tV similar
        tVxx <- max(tVx,na.rm=TRUE) #máximo tVx de todas las estaciones
        while(tVxx > tVt) {
          i <- which.max(tVx) #estación con el máximo tVx
          #si i usó referencias cortadas con un tV demasiado grande, iniciar
          #una nueva iteración:
          if(max(splt[used[i,]])>tVxx*(1+tVf*min(nref,sum(used[i,])))) break
          kp <- kpx[i] #posición del tVx en la estación i
          cat("\n",as.character(est.c[i,4]),"(",i,")",sep="")
          if(oneref[kp,i] & !force) cat(" could break at ") #sólo 1 referencia?
          else cat(" breaks at ")
          if(nm>1) cat(anyi+floor((kp-1)/nm)," ",(kp-1)%%nm+1,sep="")
          else if(nm==1) cat(anyi+kp-1)
          else cat(x[kp])
          cat(" (",round(tVx[i],1),")",sep="")
          if(oneref[kp,i] & !force) { #no cortar con una sóla referencia!:
            cat(", but it has only one reference\n")
            tVx[i] <<- -1 #pasar el tVx de esta estación a -1
            tVxx <- max(tVx,na.rm=TRUE) #máximo tVx de las estaciones restantes
            next
          }
          #gráfico de anomalías con la posición del corte:
          if(gp>1) {
            plotstan(i,varcli,x,xlab,swa,lw)
            lines(rep(x[kp],2),c(-5,4.8),col="red",lty=2) #marca del corte
            text(x[kp],5,floor(tVxx))
          }
          if(sum(!is.na(dat[1:(kp-1),i])) < mndat) {
            dat[1:(kp-1),i] <<- NA
            cat(" Fragment with less than",mndat,"data DELETED\n")
          }
          else if(sum(!is.na(dat[kp:nd,i])) < mndat) {
            dat[kp:nd,i] <<- NA
            cat(" Fragment with less than",mndat,"data DELETED\n")
          }
          else {
            nn <- nn+1 #incrementamos el no. de nuevas series
            iest <<- c(iest,iest[i]) #añadir índice a la serie original
            nsp[iest[i]] <- nsp[iest[i]]+1 #y también su no. de saltos
            if(nm>0) { #contar no. de saltos por año
              z <- floor((kp-1)/nm) #término anual del salto
              nsy[z] <- nsy[z]+1 #no. de saltos por año (base 0)
            }
            dat <<- cbind(dat,rep(NA,nd)) #nueva columna de datos
            #pasar los datos a la nueva serie:
            dat[kp:nd,ne+nn] <<- dat[kp:nd,i]
            dat[kp:nd,i] <<- NA #borrar los datos pasados a la nueva serie
            #copiar las coordenadas y poner sufijo a indicativo y nombre:
            #(Usamos la lista original de estaciones, por si se borra alguna)
            z <- data.frame(est.i[iest[i],1:3],V4=paste(est.i[iest[i],4],"-",1+nsp[iest[i]],sep=""),V5=paste(est.i[iest[i],5],"-",1+nsp[iest[i]],sep=""))
            est.c <<- rbind(est.c,z)
            attr(est.c,"row.names")[ne+nn] <<- as.integer(ne+nn)
            #estimar las nuevas medias (y desv. típicas?):
            switch(std,
              { dat.m0[i] <<- mean(dat[,i],na.rm=TRUE) + refmed - mean(datmed[!is.na(dat[,i])])
                dat.m0 <<- c(dat.m0, mean(dat[,ne+nn],na.rm=TRUE)+refmed-mean(datmed[!is.na(dat[,ne+nn])])) },
              { dat.m0[i] <<- mean(dat[,i],na.rm=TRUE) * refmed / mean(datmed[!is.na(dat[,i])])
                dat.m0 <<- c(dat.m0, mean(dat[,ne+nn],na.rm=TRUE)*refmed/mean(datmed[!is.na(dat[,ne+nn])])) },
              { dat.m0[i] <<- mean(dat[,i],na.rm=TRUE) + refmed - mean(datmed[!is.na(dat[,i])])
                dat.m0 <<- c(dat.m0, mean(dat[,ne+nn],na.rm=TRUE)+refmed-mean(datmed[!is.na(dat[,ne+nn])]))
                dat.s0[i] <<- sd(dat[,i],na.rm=TRUE) + refstd - sd(datmed[!is.na(dat[,i])])
                dat.s0 <<- c(dat.s0, sd(dat[,ne+nn],na.rm=TRUE)+refstd-sd(datmed[!is.na(dat[,ne+nn])])) },
              { dat.m0[i] <<- mean(dat[,i],na.rm=TRUE) + refmed - mean(datmed[!is.na(dat[,i])])
                dat.m0 <<- c(dat.m0, mean(dat[,ne+nn],na.rm=TRUE)+refmed-mean(datmed[!is.na(dat[,ne+nn])])) }
            )
          }
          #actualizar tVx y banderas para continuar el bucle:
          modif <- TRUE #marcar si se han modificado series
          splt[i] <- tVx[i] #tV de corte de la estación i
          tVx[i] <<- 0 #anular el tVx de esta estación
          tVxx <- max(tVx,na.rm=TRUE) #máximo tVx de las estaciones restantes
        }
        if(nn) {
          cat("\n\nUpdate number of stations: ",ne,"+",nn,"new = ")
          ne <<- ne+nn  #actualizar el no. de estaciones
          cat(ne,"\n\n")
        }
        if(!nn && !modif) {
          #máximos tV, por estaciones:
          w <- ceiling(200/ne) #anchura de las barras
          if(w>10) w <- 10
          ylim <- c(0,80)
          fcol <- 4 #factor para los colores de las barras
          if(k==1) {
            xlab2 <- "tV"
            ylab <- "Maximum tV" 
            main <- "Station's maximum tV"
          }
          else {
            xlab2 <- "SNHT"
            ylab <- "Maximum SNHT"
            main <- "Station's maximum SNHT"
          }
          plot(tVx,type="n",ylim=ylim,ylab=ylab,xlab="Stations",main=main)
          z <- which(tVx>=15*fcol)
          lines(z,tVx[z],type="h",lwd=w,col="red")
          z <- which(tVx>=10*fcol & tVx<15*fcol)
          lines(z,tVx[z],type="h",lwd=w,col=hsv(.08,.9,.9))
          z <- which(tVx>=5*fcol & tVx<10*fcol)
          lines(z,tVx[z],type="h",lwd=w,col=hsv(.16,1,.8))
          z <- which(tVx<5*fcol & tVx>0)
          lines(z,tVx[z],type="h",lwd=w,col=hsv(.33,1,.8))
          grid(col=grey(.4))
          #histograma de máximos tV, globales (sin 0's, que no son reales):
          z <- tVx[!is.na(tVx) & tVx!=0]
          if(k==1) main <- "Histogram of maximum tV"
          else main <- "Histogram of maximum SNHT"
          if(sum(!is.na(z))) hist(z,breaks=min(0,floor(min(z))):max(20,ceiling(max(z))),xlab=xlab2,col="purple",main=main)
          #si hay posibles cortes con 1 sóla referencia, colorear de rojo:
          if(min(z,na.rm=TRUE)<0) hist(z[z<0],breaks=1,col=2,add=TRUE)
          if(k==2 | snhtt<1) {
            #histograma de no. de cortes por estación:
            hist(nsp,breaks=0:max(9,max(nsp)+1)-.5,col="orange",xlab="Number of splits",ylab="Number of stations",main="Number of splits per station")
            if(nm>0) { #frecuencias de fragmentación por años:
              w <- min(5,ceiling(400/na)) #anchura de las barras
              plot(anyi:anyf,nsy,type="h",lwd=w,col=2,ylim=c(0,max(10,max(nsy))),xlab="Years",ylab="Number of splits",main="Number of splits per year")
              grid(col=grey(.4))
            }
          }
          break #salir del bucle para ir al siguiente nivel
        }
      }
    }
    #gráficos de anomalías de las series homogeneizadas (con tVx máximos,
    #ordenados por series originales):
    tVx <<- rep(NA,ne) #(guardaremos los máximos tV finales)
    snhx <- rep(NA,ne) #(guardaremos los máximos SNHT finales)
    if(gp>1) for(io in 1:nei) { #para cada serie original
      wi <- which(iest==io) #estaciones derivadas de la estación io
      lwi <- length(wi)
      if(!lwi) next #(estación totalmente borrada!)
      for(i in wi) { #para cada serie derivada de la original
        y <- sanom[,i] #anomalías estandarizadas de la estación
        plotstan(i,varcli,x,xlab,swa,lw)
        #aplicar SNHTw y marcar su tV máximo (si >=1):
        st <- snhtw(y,swa); zz <- floor(st[1]); tVx[i] <<- st[1]
        if(zz) {
          kp <- st[2]
          lines(rep(x[kp],2),c(-5,4.8),col=verde,lty=2) #marca máximo SNHTw
          text(x[kp],5,zz,col=verde) #valor
        }
        #aplicar SNHT y marcar su máximo:
        st <- snht(y)
        if(sum(!is.na(st))>0) {
          kp <- which.max(st)
          snhx[i] <- floor(max(st,na.rm=TRUE))
          zz <- floor(snhx[i])
          lines(rep(x[kp],2),c(-5,4.8),lty=4) #marca máximo SNHT
          text(x[kp],-5.2,zz) #valor
        }
        #dibujar la tendencia si es significativa:
        if(nm<1) zz <- cor.test(y,x2)$p.value
        else zz <- cor.test(y,x)$p.value
        if(zz <= 0.05) {
          if(nm<1) aj <- lm(y~x2)
          else aj <- lm(y~x)
          abline(aj,col=4)
        }
      }
    }
    #anular los tVx que no se han podido calcular:
    tVx[tVx==0] <- NA
    #gráficas de las series homogeneizadas y sus correcciones:
    if(gp>2 && tVt>0) {
      plot(-1:1,-1:1,type="n",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
      text(0,0.4,"Final graphics",cex=3.5)
      text(0,-0.3,"Recalculated series\nand applied corrections",cex=2.5)
      if(nm>0) xlab <- "Years" else xlab <- "Dates"
      #lectura datos homogeneizados:
      dh <- scan(paste(varcli,"_",anyi,"-",anyf,".dah",sep="")) 
      dim(dh) <- c(nd,ne) #asignación de dimensiones
      old.par <- par(no.readonly=TRUE)
      layout(matrix(1:2,2,1,byrow=TRUE))
      par(las=1,cex=.8)
      for(i in 1:nei) { #para cada estación original
        wi <- which(iest==i) #estaciones derivadas de la estación i
        lwi <- length(wi)
        if(!lwi) next #(estación totalmente borrada!)
        if(lwi>1) vi <- TRUE else vi <- FALSE
        if(nm>1) { #representar valores anuales
          fltr <- rep(1,nm) #filtro para sumas anuales móviles
          if(gp>3) ylab <- "Running annual totals"
          else { ylab <- "Running annual means"; fltr <- fltr/nm }#medias móv.
        }
        else { fltr <- 1; ylab <- "Values" }
        tit <- paste(varcli," at ",est.c[i,4],"(",i,"), ",est.c[i,5],sep="")
        yo <- as.vector(dat.o[,i]) #datos originales
        y <- dh[,wi] #datos homogeneizados
        par(mar=c(0,4,4,2),xaxt="n")
        matplot(x,filter(y,fltr),type="l",lty=1,col=2:20,ylab=ylab,main=tit)
        lines(x,filter(yo,fltr))
        grid(col=grey(.4))
        par(mar=c(5,4,0,2),xaxt="s")
        #correcciones:
        #(no se usa matplot porque no maneja bien las fechas en el eje X)
        if(std==2) {
          yo[yo==0] <- NA; yd <- y/yo; ylab <- "Correction factors"
          if(!vi) ylim <- c(0,2)
        }
        else {
          yd <- y-yo; ylab <- "Correction terms"
          if(!vi) ylim <- c(-1,1)
        }
        if(vi) {
          ylim <- c(floor(min(yd,na.rm=TRUE)),ceiling(max(yd,na.rm=TRUE)))
          plot(x,yd[,1],type="n",ylim=ylim,ylab=ylab,xlab=xlab)
        }
        else plot(x,yd,type="n",ylim=ylim,ylab=ylab,xlab=xlab)
        matlines(x,yd,type="l",lty=1,col=2:20)
        grid(col=grey(.4))
      }
      par(old.par)
    }
    cat("\n======== End of the homogenization process (",date(),")\n",sep="")
    cat("\n----------- Final computations:\n")
    #autocorrelaciones de las anomalías de cada estación
    cat("\nACmx: Station maximum absolute autocorrelations of anomalies\n")
    sac <- rep(NA,ne) #vector de máximas autocorrelaciones
    for(i in 1:ne) {
      zz <- acf(anom[,i],plot=FALSE,na.action=na.pass)$acf
      zz[1] <- 0 #anulamos la autocorrelación trivial
      sac[i] <- max(abs(zz)) #máxima autocorrelación con diferentes desfases
    }
    print(summary(round(sac,2)))
    #prueba SNHT-w de cada estación
    cat("\ntVx: Stepped window SNHT (on anomaly series)\n")
    print(summary(round(tVx,1)))
    #prueba SNHT de cada estación
    cat("\nSNHT: Standard normal homogeneity test (on anomaly series)\n")
    print(summary(round(snhx,1)))
    #errores típicos de las estimas (sin estandarizar):
    cat("\nRMSE: Root mean squared error of the estimated data\n")
    #desestandarizamos las anomalías si std>1 :
    if(std==2 && dat.m[i]>1) { #no multiplicar por menos de 1
      for(i in 1:ne) anom[,i] <- anom[,i] * dat.m[i]
    }
    else if(std==3) {
      for(i in 1:ne) anom[,i] <- anom[,i] * dat.s[i]
    }
    if(rtrans) { #deshacer transformación, mediante límites de confianza:
      rmse <- rep(NA,ne) #rmse medios de cada estación
      sdz <- sd(anom,na.rm=TRUE)
      for(i in 1:ne) {
        l1 <- (dat.m[i]-sdz[i])^rtrans
        l2 <- (dat.m[i]+sdz[i])^rtrans
        rmse[i] <- (l2-l1)/2 #rmse de la estación i
      }
    }
    else rmse <- sd(anom,na.rm=TRUE)
    zz <- summary(rmse)
    print(zz)
    sedec <- max(1,2-ceiling(log10(zz[4]))) #no. de decimales de RMSE
    pod <- floor(100*(nd-apply(dat.na,2,sum))/nd) #porcentaje de datos originales
    cat("\nPD: Percentage of original data\n")
    print(summary(pod))
    #imprimir resumen de resultados
    cat("\n")
    print(data.frame(ACmx=round(sac,2),tVx=round(tVx,2),SNHT=snhx,RMSE=round(rmse,sedec),PD=pod,Code=est.c[,4],Name=est.c[,5]),right=FALSE)
    #averiguar qué estaciones derivadas funcionan al final del periodo:
    cur <- apply(!is.na(dat[(nd-mndat+1):nd,]),2,sum) #últimos mndat términos
    cur[cur>0] <- 1
    #escribir la nueva tabla de estaciones (con extensión 'esh', añadiendo
    #el porcentaje de datos originales, la estación original, si funciona
    #actualmente, y el SNHT):
    arch <- paste(varcli,"_",anyi,"-",anyf,".esh",sep="") #nombre del archivo
    write.table(cbind(est.c,pod,iest,cur,snhx),arch,row.names=FALSE,col.names=FALSE)
    cat("\n----------- Generated outputs: --------------------------------\n\n")
    cat(varcli,"_",anyi,"-",anyf,".txt : This text output\n",sep="")
    cat(varcli,"_",anyi,"-",anyf,".dah : Homogenized data (postprocess with 'dahstat()')\n",sep="")
    cat(varcli,"_",anyi,"-",anyf,".esh : List of homogenized stations (original and split)\n",sep="")
    if(gp>1) { #últimos gráficos
      #histograma de las anomalías de los outliers:
      main <- "Histogram of normalized anomalies"
      if(sum(!is.na(outan))) {
        xlim <- range(sanom,na.rm=TRUE)
        xlim[1] <- min(xlim[1], min(outan,na.rm=TRUE))
        xlim[2] <- max(xlim[2], max(outan,na.rm=TRUE))
        hist(sanom,breaks=30,xlim=xlim,xlab="Anomaly",col="green",main=main)
        hist(outan,col="red",add=TRUE)
      }
      else hist(sanom,breaks=30,xlab="Anomaly",col="green",main=main)
      abline(v=-dz.max,col=2); abline(v=dz.max,col=2)
      #histograma de tVx:
      z <- tVx; main <- "Histogram of maximum tV"
#     if(sum(!is.na(z))) hist(z,breaks=min(0,floor(min(z))):max(20,ceiling(max(z))),xlab="tVx",col=verde,main=main)
      if(sum(!is.na(z))) hist(z,breaks=20,xlab="tVx",col=verde,main=main)
      #histograma de SNHT:
      z <- snhx; main <- "Histogram of maximum SNHT"
#     if(sum(!is.na(z))) hist(z,breaks=min(0,floor(min(z))):max(20,ceiling(max(z))),xlab="SNHT",col="purple",main=main)
      if(sum(!is.na(z))) hist(z,breaks=20,xlab="SNHT",col="purple",main=main)
      #gráfico de calidad/singularidad:
      plot(rmse,snhx,type="n",xlim=c(0,max(1,max(rmse,na.rm=TRUE))),ylim=c(0,max(50,max(snhx,na.rm=TRUE))),xlab="RMSE",ylab="SNHT",main="Station's quality/singularity")
      grid(col=grey(.4))
      text(rmse,snhx,col=hsv(.7,1,.9))
    }
    if(gp>0) { #cerrar la salida gráfica
      graphics.off() #volcar el buffer del último gráfico
      cat(varcli,"_",anyi,"-",anyf,".pdf : Diagnostic graphics\n",sep="")
    }
    cat("\n")
  }
  sink() #cerrar el archivo de bitácora
}

