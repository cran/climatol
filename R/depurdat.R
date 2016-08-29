#depurdat.R.- Depuración y homogeneización de series climatológicas.
#(Most comments in Spanish; sorry!)

climatol.version <- '3.0'

#- cerrar.- Cerrar los archivos de salida.
cerrar <- function() {
  sink()
  graphics.off()
}

#- dahgrid.- Generación de grids de datos homogeneizados.
dahgrid <- function(varcli, anyi, anyf, anyip=anyi, anyfp=anyf, grid,
  mh=FALSE, std=NA, ini=NA, obsonly=TRUE) {
  #anyip: año inicial de referencia (para el cálculo de las anomalías)
  #anyfp: año final de referencia
  #grid: grid a interpolar, de clase SpatialPixel
  #mh: Si TRUE, leer datos mensuales de la homogeneización diaria (*-mh_*.dat)
  #std: Incluido en el fichero *.rda, pero si mh=TRUE ese fichero no se lee,
  #     y entonces conviene especificarlo. (=3 por defecto)
  #ini: Incluido en el fichero *.rda, pero si mh=TRUE ese fichero no se lee,
  #     y entonces puede especificarse aquí. (Por defecto, 1 de enero de anyi))
  #obsonly: do not interpolate missing data estimated by homogen()
  if(!requireNamespace("sp", quietly=TRUE)
   | !requireNamespace("gstat", quietly=TRUE)
   | !requireNamespace("raster", quietly=TRUE)
   | !requireNamespace("ncdf4", quietly=TRUE)
  ) stop('This function requires packages sp, gstat, raster and ncdf4.\nPlease, install the lacking packages an re-run the function')
  if(anyip<anyi) stop("Asked initial reference year before first year of data!")
  if(anyfp>anyf) stop("Asked final reference year beyond last year of data!")
  #- lectura de los datos originales y homogeneizados
  if(!mh) {
    fbas <- sprintf('%s_%d-%d',varcli,anyi,anyf) #raíz nombres de fichero
    load(sprintf('%s.rda',fbas))
  } else { #leer los datos mensuales *-mh_*.dat (generados por dd2m)
    fbas <- sprintf('%s-mh_%d-%d',varcli,anyi,anyf) #raíz nombres de fichero
    est.c <- read.table(sprintf('%s.est',fbas),colClasses=c('numeric','numeric','numeric','character','character','numeric','numeric','numeric','numeric'))
    ne <- nrow(est.c)
    dah <- scan(sprintf('%s.dat',fbas))
    nd <- length(dah)/ne
    est.b <- read.table(sprintf('%s-m_%d-%d.est',varcli,anyi,anyf),colClasses=c("numeric","numeric","numeric","character","character"))
    nei <- nrow(est.b)
    if(obsonly) { #leer datos originales y asegurar el orden correcto:
      dat <- scan(sprintf('%s-m_%d-%d.dat',varcli,anyi,anyf))
      if(!identical(est.b[1:nei,4],est.c[1:nei,4])) {
        dim(dat) <- c(nd,nei)
        dat <- dat[,match(est.c[1:nei,4],est.b[1:nei,4])]
      }
    }
    if(is.na(std)) std <- 3 #valor por defecto si mh=TRUE
    nm <- 12 #mh se usa solo con valores mensuales
  }
  dim(dah) <- c(nd,ne) #conversión a 2 dimensiones
  dim(dat) <- c(nd,nei) #conversión a 2 dimensiones
  #- seleccionar las series de los fragmentos más largos
  sel <- tapply(est.c[,6],est.c[,7],which.max)
  ksel <- rep(NA,nei)
  for(k in 1:nei) ksel[k] <- which(est.c[,7]==k)[sel[k]]
  #- retener solo los datos homogeneizados (dah) de las subseries más largas
  dah <- dah[,ksel]
  #- calcular sus medias y desviaciones típicas en el periodo escogido
  if(anyip==anyi & anyfp==anyf) { ki <- 1; kf <- nd } else {
  ki <- (anyip-anyi)*nm+1; kf=ki+(anyfp-anyip+1)*nm-1 } #pos. inicial y final
  m <- apply(dah[ki:kf,],2,mean)
  if(std>2) s <- apply(dah[ki:kf,],2,sd)
  #- grabarlas, con sus coordenadas, para su uso con GIS
  if(std<3) {
    df <- data.frame(est.c[ksel,1:4],m)
    names(df) <- c('X','Y','Z','Code','Means')
  } else {
    df <- data.frame(est.c[ksel,1:4],m,s)
    names(df) <- c('X','Y','Z','Code','Means','Std.Dev.')
  }
  fmeans <- sprintf('%s_%d-%d_means.csv',varcli,anyip,anyfp)
  write.csv(df,fmeans,row.names=FALSE)
  #- normalizar las series
  switch(std,
    daz <- scale(dah,center=m,scale=FALSE), #std=1
    if(min(m)<1) { z <- which(m > 1)
      daz <- dah
      daz[,z] <- scale(dah[,z],center=FALSE,scale=m[z]) }
    else daz <- scale(dah,center=FALSE,scale=m),
    daz <- scale(dah,center=m,scale=s) #std=3 (default)
  )
  #- if(obsonly), blanquear los datos ausentes en las series originales
  if(obsonly) daz[is.na(dat)] <- NA
  rg <- range(daz,na.rm=TRUE) #rango de valores
  #- interpolar las medias (y desv. típicas), y grabarlas en NetCDF
  df <- data.frame(est.c[ksel,1:2],m)
  names(df) <- c('x','y','z')
  sp::coordinates(df) <- ~x+y
  m <- gstat::idw(z~1, df, grid, debug.level=0) #medias interpoladas
  dimLon <- ncdf4::ncdim_def(name='lon', units='deg.E', vals=unique(grid@coords[,1]))
  dimLat <- ncdf4::ncdim_def(name='lat', units='deg.N', vals=rev(unique(grid@coords[,2])))
  varCli.m <- ncdf4::ncvar_def(name=sprintf('%s.m',varcli), units='', dim=list(dimLon,
    dimLat), missval=NA, longname=sprintf('%s %d-%d means',varcli,anyip,anyfp))
  listvar <- list(varCli.m)
  nc <- ncdf4::nc_create(sprintf('%s_m.nc',fbas), listvar) #abrir el fichero netcdf
  zz <- raster::rasterFromXYZ(m)
  ncdf4::ncvar_put(nc,varCli.m,zz@data@values)
  ncdf4::nc_close(nc)
  if(std>2) {
    df <- data.frame(est.c[ksel,1:2],s)
    names(df) <- c('x','y','z')
    sp::coordinates(df) <- ~x+y
    s <- gstat::idw(z~1, df, grid, debug.level=0) #desv. típicas interpoladas
    varCli.s <- ncdf4::ncvar_def(name=sprintf('%s.s',varcli), units='',
      dim=list(dimLon, dimLat), missval=NA,
      longname=sprintf('%s %d-%d std. deviations',varcli,anyip,anyfp))
    listvar <- list(varCli.s)
    nc <- ncdf4::nc_create(sprintf('%s_s.nc',fbas), listvar) #abrir el fichero netcdf
    zz=raster::rasterFromXYZ(s)
    ncdf4::ncvar_put(nc,varCli.s,zz@data@values)
    ncdf4::nc_close(nc)
  }
  #- === crear un netcdf con los grids interpolados en cada paso de tiempo
  if(is.na(ini)) ini <- sprintf('%d-01-01',anyi) #fecha inicial por defecto
  if(nm>0) x <- seq(as.Date(ini),length.out=nd,by=sprintf('%d months',12/nm))
  else x <- seq(as.Date(ini),length.out=nd,by='1 day')
  dimTime <- ncdf4::ncdim_def(name='Date', units='days since 1970-01-01',
    vals=as.numeric(x), calendar='standard')
  varCli <- ncdf4::ncvar_def(name=varcli, units='', dim=list(dimLon, dimLat,
    dimTime), missval=NA)
  listvar <- list(varCli)
  nc <- ncdf4::nc_create(sprintf('%s.nc',fbas), listvar) #abrir el fichero netcdf
  #- para cada paso de tiempo:
  cat(sprintf('Interpolating %d grids...:      ',nd))
  kz <- max(10,round(nd/100))
  for(k in 1:nd) {
    if(!k%%kz) cat('\b\b\b\b\b',sprintf('%2s %%',round(k*100/nd)))
    #- interpolar (IDW) la variable estandarizada a los puntos del grid
    df <- data.frame(est.c[ksel,1:2],daz[k,])
    if(obsonly) df <- df[!is.na(df[,3]),]
    names(df) <- c('x','y','z')
    sp::coordinates(df) <- ~x+y
    z <- gstat::idw(z~1, df, grid, debug.level=0) #desv. típicas interpoladas
    #pasar de SpatialPointsDataFrame a RasterLayer:
    zz=raster::rasterFromXYZ(z)
    #- grabar los valores en el netcdf
    ncdf4::ncvar_put(nc,varCli,zz@data@values,start=c(1,1,k),count=c(-1,-1,1))
  }
  cat(' (done)\n\n')
  #- cerrar el netcdf y terminar
  ncdf4::nc_close(nc)
  cat(sprintf('Normalized grids (%f to %f) saved to file %s.nc',rg[1],rg[2],fbas),'\n')
  cat('Means')
  if(std>2) cat(' and standard deviations')
  cat(' (of the whole series) saved to files\n')
  cat(sprintf('%s_m.nc',fbas))
  if(std>2) cat(',',sprintf('%s_s.nc',fbas))
  cat(' and',fmeans,'\n\n')
}

#- dahstat.- Estadísticas de datos homogeneizados.
dahstat <- function(varcli, anyi,anyf, anyip=anyi, anyfp=anyf, stat="me",
  ndc=1, vala=2, cod=NULL, mnpd=0, mxsh=0, prob=.5, last=FALSE, long=FALSE,
  mh=FALSE, pernys=100, ini=NA, estcol=4, sep=',', dec='.', eol="\n") {
#stat='tnd' genera también los p-valores, en *.pval
#stat="me"(valores medios), "mdn"(medianas), "max"(máximos), "min"(mínimos),
# "std"(desv.típ.), "q"(cuantiles), "tnd"(tendencias), "series"(todos los datos)
#ndc=no. de decimales (con prioridad sobre el ndec guardado en el *.rda)
#vala= 0(ninguno), 1(suma), 2(media), 3(máximo), 4(mínimo)
#cod=lista de códigos de las estaciones a listar
#mnpd=Mínimo porcentaje de datos originales
#mxsh=Máximo SNHT
#last=Listar sólo las estaciones que operaban al final del periodo.
#long=Listar sólo las estaciones con mayor fragmento original
#mh=Si TRUE, leer datos mensuales de la homogeneización diaria (*-mh_*.dat)
#pernys=No. de años sobre los que se expresan las tendencias (100 por defecto)
#estcol=Columnas de est.c seleccionadas para el listado
#ini=Fecha inicial, sólo si mh=TRUE: si no se fija se supone que es el 1 de
#    enero de anyi. (Si mh=FALSE, ini se toma del fichero *.rda)
#Los parámetros sep, dec y eol permite personalizar los formatos de salida
  #- inicializaciones
  if(anyi>anyf) stop ('First year of data greater than the last year!')
  if(anyip<anyi) stop("Asked initial year before first year of data!")
  if(anyfp>anyf) stop("Asked final year beyond last year of data!")
  #función elegida para el cálculo de los valores mensuales:
  fun <- c("mean","median","max","min","sd","quantile")[which(c("me","mdn","max","min","std","q","tnd")==stat)]
  #- si no se reconoce la opción stat, terminar aquí
  if(!length(fun) & stat!='series')
    stop(sprintf("Option stat='%s' not recognized!",stat))
  estvar <- c('X','Y','Z','Code','Name','pod','ios','ope','snht')
  mes3 <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
##mes3 <- c("Ene","Feb","Mar","Abr","May","Jun","Jul","Ago","Sep","Oct","Nov","Dic")
  #- leer los datos de entrada
  if(!mh) load(sprintf('%s_%d-%d.rda',varcli,anyi,anyf))
  else { #leer los datos mensuales *-mh_*.dat (generados por dd2m)
    dah <- scan(sprintf('%s-mh_%d-%d.dat',varcli,anyi,anyf))
    est.c <- read.table(sprintf('%s-mh_%d-%d.est',varcli,anyi,anyf),colClasses=c("numeric","numeric","numeric","character","character","numeric","numeric","numeric","numeric"))
    ne <- nrow(est.c)
    nd <- length(dah)/ne
    est.b <- read.table(sprintf('%s-m_%d-%d.est',varcli,anyi,anyf),colClasses=c("numeric","numeric","numeric","character","character"))
    nei <- nrow(est.b)
    nm <- 12 #sólo datos mensuales si mh=TRUE
    na <- anyf-anyi+1
    if(nd != na*nm) stop(sprintf('Number of data (%d) is not number of months (%d) times number of years (%d)',nd,nm,na))
    dim(dah) <- c(nm,na,ne)
    if(stat=='series') {
      dat <- scan(sprintf('%s-m_%d-%d.dat',varcli,anyi,anyf))
      dim(dat) <- c(nm,na,nei)
    } else rm(est.b)
  }
  ndec <- ndc #actualizar el no. de decimales (ndec estaba en el fichero *.rda)
  if(nm>0) na <- anyf-anyi+1 else na <- nd #no. de años o de datos diarios
  if(nm<2) {
    dim(dah) <- c(1,na,ne)
    if(stat=='series') dim(dat) <- c(1,na,nei)
  }
  if(nm<2 | stat=='series') vala <- 0 #valor anual innecesario
  else {
    if(vala<0 | vala>4) vala <- 2 #valor medio en caso de vala erróneo
    funa <- c("sum","mean","max","min")[vala] #función para el valor anual
  }
  #- seleccionar estaciones solicitadas:
  if(!is.null(cod)) {
    ksel <- which(est.c[,4] %in% cod) #estaciones originales solicitadas
    ksel <- which(est.c[,7] %in% ksel) #id. madres e hijas
    dah <- dah[,,ksel] #datos de esas series
    if(nm<2) dim(dah) <- c(1,na,length(ksel))
    est.c <- est.c[ksel,] #coordenadas de esas estaciones
  }
  #seleccionar las estaciones con un mínimo de mnpd % de datos originales:
  esel <- est.c[,6]>=mnpd #vector de la selección
  #seleccionar las estaciones con un SNHT menor o igual a mxsh:
  if(mxsh>0) esel <- esel & est.c[,9]<=mxsh #vector de la selección
  if(last) esel <- esel & as.logical(est.c[,8]) #sólo últimos fragmentos?
  else if(long) {
    lsel <- rep(TRUE,length(esel)) #inicializar vector
    for(ko in 1:nei) { #para cada estación original
      kest <- which(est.c[,7]==ko) #series de la misma estación ko
      if(length(kest)>1) { #si hay más de un fragmento...
        kmax <- which.max(est.c[kest,6]) #cuál tiene mayor % de datos orig.?
        lsel[kest[-kmax]] <- FALSE #desseleccionamos los fragmentos menores
      }
    }
    esel <- esel & lsel
  }
  ne <- sum(esel) #no. de estaciones seleccionadas
  if(ne==0) stop("No station selected: No output")
  dah <- dah[,,esel]
  dim(dah) <- c(max(c(1,nm)),na,sum(esel))
  est.c <- est.c[esel,] #lista de estaciones seleccionadas
  iest <- est.c[,7] #índice de las correspondientes estaciones originales
  #- if(vala), calcular los valores anuales
  if(vala) { #calcular los valores anuales
    aval <- as.vector(apply(dah,2:3,funa))
    dim(dah) <- c(nm,na*ne)
    dah <- rbind(dah,aval)
    nm <- nm+1
    dim(dah) <- c(nm,na,ne)
  }
  #dimensionar valores a calcular:
  if(stat!="series") val <- matrix(NA,ne,max(1,nm))
  if(nm>0) { #valores mensuales, estacionales, anuales, ...
    x <- anyip:anyfp #vector de años solicitados
    xk <- x-anyi+1 #posiciones de x dentro de anyi:anyf
    pernum <- pernys
  } else { #valores diarios
    xk <- 1:nd
    if(is.na(ini)) ini <- sprintf('%d-01-01',anyi) #fecha inicial por defecto
    x <- seq(as.Date(ini),length.out=nd,by='1 day')
    pernum <- pernys*365.25
  }
  #- if(stat=="tnd"), calcular las tendencias del periodo escogido
  if(stat=="tnd") { #tendencias del periodo escogido
    pval=val #matriz para almacenar los p-valores
    for(i in 1:ne) {
      if(nm<2) { #una sola subserie
        aj <- lm(dah[1,xk,i]~x) #regresión lineal
        val[i,] <- round(aj$coefficients[2]*pernum,ndec)
        pval[i,] <- round(summary(aj)$coefficients[2,4],3)
      }
      else {
        for(j in 1:nm) {
          aj <- lm(dah[j,xk,i]~x) #regresión lineal
          val[i,j] <-  round(aj$coefficients[2]*pernum,ndec)
          pval[i,j] <- round(summary(aj)$coefficients[2,4],3)
        }
      }
    }
  }
  #- else if(stat=="series"), listar los valores (arch. individuales)
  else if(stat=="series") { #listar en ficheros sueltos con formato CSV
    # los valores y sus flags (0:=original, 1:rellenado, 2:corregido)
    for(kest in 1:ne) { #para cada estación homogeneizada seleccionada
      dh <- dah[,xk,kest] #datos homogeneizados (del periodo solicitado)
      ik <- iest[kest]
      if(mh) ik <- which(est.b[,4]==est.c[ik,4]) #(podría cambiar el orden...)
      do <- dat[,xk,ik] #datos originales (del periodo solicitado)
      #comparación datos homogeneizados con originales (no usar '=='!):
      df <- abs(dh-do) < 1e-9
      df <- as.numeric(df) #TRUE=1, FALSE=0
      df[df==0] <- 2 #datos distintos a los originales
      df[df==1] <- 0 #datos iguales a los originales
      df[is.na(df)] <- 1 #datos rellenados (originales ausentes)
      dim(df) <- dim(dh)
      #nombre de los archivos de salida:
      ard <- sprintf('%s_%d-%d_%s.csv',varcli,anyip,anyfp,est.c[kest,4])
      arf <- sprintf('%s_%d-%d_%s-flg.csv',varcli,anyip,anyfp,est.c[kest,4])
      if(nm>1) { dh <- cbind(x,t(dh)); df <- cbind(x,t(df)) }
      else { dh <- cbind(format(x),dh); df <- cbind(format(x),df) }
      if(nm==12) colnames(dh) <- c('Year',mes3)
      else if(nm>1) colnames(dh) <- c('Year',as.character(1:nm))
      else colnames(dh) <- c('Date','Value')
      colnames(df) <- colnames(dh)
      if(stat=="series") {
        write.csv(dh,ard,row.names=FALSE,quote=FALSE)
        write.csv(df,arf,row.names=FALSE,quote=FALSE)
      }
    }
    cat(sprintf('Homogenized values written to %s_%d-%d_*.csv,\nwith flags in %s_%d-%d_*-flg.csv:\n',varcli,anyip,anyfp,varcli,anyip,anyfp))
    cat('  0: Observed data\n')
    cat('  1: Missing data (filled)\n')
    cat('  2: Corrected data\n')
    ars <- paste(varcli,"_",anyip,"-",anyfp,".pval",sep="") #fichero de salida
    return(invisible())
  }
  #- else, aplicar la función deseada al periodo escogido
  else { #aplicar la función deseada al periodo escogido
    for(i in 1:ne) {
      if(nm<2) {
        if(stat=="q") val[i,] <- round(eval(call(fun,dah[,xk,i],prob)),ndec)
        else val[i,] <- round(eval(call(fun,dah[,xk,i])),ndec)
      }
      else { #datos mensuales:
        if(stat=="q") val[i,] <- round(apply(dah[,xk,i],1,fun,prob),ndec)
        else val[i,] <- round(apply(dah[,xk,i],1,fun),ndec)
      }
    }
  }
  #- imprimir mensaje con los ficheros generados
  if(stat=="me") cat("Mean")
  else if(stat=="mdn") cat("Median")
  else if(stat=="max") cat("Maximum")
  else if(stat=="min") cat("Minimum")
  else if(stat=="std") cat("Standard deviation")
  else if(stat=="q") cat(prob,"prob. quantile")
  else if(stat=="tnd") cat("Trend")
  cat(" values of ",varcli," (",anyip,"-",anyfp,")",sep="")
  if(stat=="tnd") cat(", expressed in units per ",pernys," years,",sep="")
  dahs <- data.frame(cbind(est.c[estcol],val))
  if(nm==12) ndf <- c(estvar[estcol],mes3)
  else if(nm==13) ndf <- c(estvar[estcol],mes3,"Annual")
  else if(nm<2) ndf <- c(estvar[estcol],"Value")
  else ndf <- c(estvar[estcol],1:nm)
  names(dahs) <- ndf
  #- grabar los valores en los ficheros
  #fichero de salida:
  if(stat=="q") ars <- sprintf('%s_%d-%d_%s%d.csv',varcli,anyip,anyfp,stat,100*prob)
  else ars <- sprintf('%s_%d-%d_%s.csv',varcli,anyip,anyfp,stat)
  write.csv(dahs,ars,row.names=FALSE,quote=FALSE)
  cat("\n  written to",ars,"\n")
  if(stat=="tnd") { #grabar los p-valores
    dahs2 <- data.frame(cbind(est.c[estcol],pval))
    names(dahs2) <- ndf
##  ars <- paste(varcli,"_",anyip,"-",anyfp,".pval",sep="") #fichero de salida
    ars <- sprintf('%s_%d-%d_pval.csv',varcli,anyip,anyfp)
    write.csv(dahs2,ars,row.names=FALSE,quote=FALSE)
    cat("P-values written to",ars,"\n")
  }
}

#- db2dat.- Get data from a database and build input files *.dat and *.est for
#the homogen() function. (ODBC must be intalled and properly configured.)
#----------------------------------------------------------------------
#Example for a database called "climate", with user "USER" and password "PASS":
# R  #start R (version 3 or higher)
# library(RODBC)
# ch <- odbcConnect("climate",uid="USER",pwd="PASS") #connect to database
# db2dat('HRel',1961,2015,10,FALSE,ch,'%Y-%m-%d','monthly_relhum','Station',
# 'Date','Value','stations','Station','Name','Longitude','Latitude','Elevation')
# odbcClose(ch) #close connection to mcheng
#----------------------------------------------------------------------
# This example will compile monthly average relative humidity for the period
# 1961-2015 excluding series with less than 10 years of data (120 monthly data)
# in files HRel_1961-2015.dat and HRel_1961-2015.est, which you can
# homogenize later with the Climatol R package with, e.g.:
# library(climatol)
# homogen('HRel',1961,2015,vmin=0,vmax=100)
# -------------------------------------------------------------------
db2dat <- function(varcli,anyi,anyf,minny=5,daily=TRUE,ch,
  dformat='%Y-%m-%d',vtable,vcode,vdate,vval,stable,scode,sname,sx,sy,sz) {
  #varcli: Achronym of the climatic variable under study
  #anyi:   Fist year of the study period
  #anyf:   Last year of the study period
  #minny:  Minimum number of years with data in the series to study
  #ch:     Name of the ODBC conexion to the database
  #dformat:Format of dates in the database
  #vtable: Name of the table containing our climatic variable
  #vcode:  Name of the variable containing station codes in the database
  #vdate:  Name of the variable containing dates in the database
  #vval:   Name of the climatic variable in the database
  #stable: Name of the table containing station information (metadata)
  #scode:  Name of the variable containing station codes
  #sname:  Name of the variable containing station names
  #sx:     Name of the variable containing longitudes (degrees with decimals!)
  #sy:     Name of the variable containing latitudes (degrees with decimals!)
  #sz:     Name of the variable containing elevations (meters)
  #- initializations
  na <- anyf-anyi+1 #nr. of years
  if(na<=0) stop('Last year must be greater than the first year')
  if(daily) {
    x <- seq(as.Date(sprintf('%d-01-01',anyi)),as.Date(sprintf('%d-12-31',anyf)),by='1 day')
    ndmin <- round(minny*365.25) #min. nr. of daily data
  } else {
    x <- seq(as.Date(sprintf('%d-01-01',anyi)),as.Date(sprintf('%d-12-01',anyf)),by='1 month')
    ndmin <- minny*12 #min. nr. of monthly data
  }
  nd <- length(x) #nr. of data per station
  #- read station names and coordinates
  cat('Getting station names and coordinates...\n')
  ds <- RODBC::sqlQuery(ch,sprintf("SELECT %s, %s, %s, %s, %s FROM %s", sx,sy,sz,scode,sname,stable,scode))
  ds[,1:2] <- round(ds[,1:2],5) #round coordinates to 5 decimals
  ds[,3] <- round(ds[,3],1) #round elevations to 1 decimal
  ds[,4] <- as.character(ds[,4]) #force codes as character strings
  ds[,5] <- as.character(ds[,5]) #force names as character strings
  ds <- ds[order(ds[,4]),] #order stations by code
  ns <- nrow(ds); ndat <- rep(0,nd)
  #- open data and stations files
  dfile <- sprintf('%s_%d-%d.dat',varcli,anyi,anyf)
  efile <- sprintf('%s_%d-%d.est',varcli,anyi,anyf)
  Fd <- file(dfile,'w')
  Fe <- file(efile,'w')
  #- get data from the ODBC connection, station by station
  cat('Getting data for every station...\n')
  ne <- 0
  for(i in 1:ns) { #for every station
    cat(unlist(ds[i,]),'\n')
    if(sum(is.na(ds[i,]))>0) {
      cat('Warning: Incomplete metadata (station skipped)\n')
      next
    }
    dd <- RODBC::sqlQuery(ch,sprintf("SELECT %s,%s FROM %s WHERE %s >= '%s' AND %s <= '%s' AND %s = '%s'",vdate,vval,vtable,vdate,strftime(x[1],dformat),vdate,strftime(x[nd],dformat),vcode,ds[i,4]))
    if(sum(!is.na(dd[,2])) < ndmin) next #not enough data
    dd[,1] <- as.Date(dd[,1],format=dformat,tz='') #force vdate to class Date
    k <- match(dd[,1],x) #match data time steps
    if(sum(is.na(k))>0) {
      cat('Warning: Station skipped because some or all of its dates do not match the expected values\n')
      next
    }
    dat <- rep(NA,nd) #initialize data vector
    dat[k] <- dd[,2] #assign data
    write(dat,Fd,ncolumns=ifelse(daily,10,12)) #write into data file
    write.table(ds[i,],Fe,row.names=FALSE,col.names=FALSE) #write metadata
    ne <- ne + 1 #count nr. of saved series
    ndat <- ndat + !is.na(dat) #count nr. of data at every time step
  }
  #close files:
  close(Fe); close(Fd)
  cat(sprintf('\nFiles %s and %s successfully generated.',dfile,efile))
  #check data availability along time:
  if(min(ndat)==0) {
    ks <- which(ndat==0)
    cat(sprintf(' BUT:\nNo data available in any station for %s',ifelse(daily,'day','month')))
    if(length(ks)>1) cat('s:\n') else cat(' ')
    print(x[ks])
    cat(sprintf('Add stations or shorten the study period to avoid %s without data\n',ifelse(daily,'days','months')))
  } else cat('\n')

}

#- dd2m.- Cálculo de valores mensuales a partir de datos diarios.
dd2m <- function(varcli, anyi, anyf, anyip=anyi, anyfp=anyf, ndec=1, suf=NA,
  valm=2, namax=10, na.strings="NA", homog=FALSE, ini=NA) {
  #suf: sufijo opcional a añadir al nombre de la variable para leer los datos.
  #valm: Valor mensual (1=suma, 2=media, 3=máximo, 4=mínimo)
  #namax: Máximo no. permitido de datos diarios originalmente ausentes
  #homog: Usar datos ya homogeneizados? (poner homog=TRUE)
  #ini: Fecha inicial. Si es NA se supone que es el 1 de enero de anyi
  fbas <- sprintf('%s_%d-%d',varcli,anyi,anyf) #raíz nombres de fichero
  if(is.na(suf)) fntr <- fbas #raíz nombres ficheros de entrada
    else fntr <- sprintf('%s-%s_%d-%d',varcli,suf,anyi,anyf)
  if(homog) {
    load(sprintf('%s.rda',fntr)) #leer datos homogeneizados y originales
  } else {
    dah <- scan(sprintf('%s.dat',fntr),na.strings=na.strings) #datos originales
    est.c <- read.table(sprintf('%s.est',fntr),colClasses=c("numeric","numeric","numeric","character","character")) #coord. estaciones
    ne <- nrow(est.c) #no. de estaciones originales
    nd <- length(dah)/ne #no. de datos por estación
    dim(dah) <- c(nd,ne)
    na <- anyf-anyi+1 #no. de años
    nm <- nd/na #no. de "meses" (no. de datos por año y estación)
    if(nm < 180) stop(sprintf("These data does not seem daily (%d items per year)",round(nm)))
  }
  if(!is.na(ini)) fech <- as.Date(0:(nd-1),origin=ini) #fechas
  else fech <- as.Date(0:(nd-1),origin=sprintf('%4d-01-01',anyi))
  me <- strftime(fech,"%m")
  anyo <- strftime(fech,"%Y")
  fun <- c("sum","mean","max","min")[valm] #función para el valor mensual
  na <- anyfp-anyip+1 #no. de años
  dm <- matrix(NA,na*12,ne) #datos mensuales
  for(ie in 1:ne) { #para cada estación
    cat(' ',ie)
    z <- aggregate(dah[,ie],list(me,anyo),fun,na.rm=TRUE) #valores mensuales
    z[,3] <- round(z[,3],ndec) #redondear
    z2 <- aggregate(is.na(dah[,ie]),list(me,anyo),sum) #no. de datos ausentes
    #conservar solo el periodo deseado:
    zp <- z[,2]>=anyip & z[,2]<=anyfp
    z <- z[zp,]; z2 <- z2[zp,] 
    zz <- z2[,3] <= namax #meses con suficientes datos
    dm[zz,ie] <- z[zz,3] #asignar los datos mensuales a la matriz general
  }
  dm[is.nan(dm)] <- NA #si no hay datos, poner NA
  #grabar los datos mensuales:
  if(homog) {
    fichsal <- sprintf("%s-mh_%d-%d.dat",varcli,anyip,anyfp)
    fichest <- sprintf("%s-mh_%d-%d.est",varcli,anyip,anyfp)
  } else {
    fichsal <- sprintf("%s-m_%d-%d.dat",varcli,anyip,anyfp)
    fichest <- sprintf("%s-m_%d-%d.est",varcli,anyip,anyfp)
  }
  write(round(dm,ndec),fichsal,ncolumns=12)
  write.table(est.c,fichest,row.names=FALSE,col.names=FALSE)
  cat("\n\nMonthly",fun,"values output to file",fichsal,"\n")
  if(namax>0 & !homog) cat('  (Months with more than',namax,'missing original daily data\n  have also been set to missing)\n')
}

#- homogen.- homogeneización automática de un conjunto de series de datos.
homogen <- function(varcli, anyi, anyf, suf=NA, nm=NA, nref=c(10,10,4), std=3,
swa=NA, ndec=1, dz.max=5, dz.min=-dz.max, wd=c(0,0,100), snht1=25, snht2=snht1,
tol=.02, mxdif=NA, force=FALSE, wz=.001, trf=0, mndat=NA, gp=3, ini=NA,
na.strings="NA", maxite=50, vmin=NA, vmax=NA, nclust=100,
clustmethod='ward.D2', grdcol=grey(.5), mapcol=grey(.65), hires=TRUE,
expl=FALSE, metad=FALSE, sufbrk='m', verb=TRUE) {
  #varcli: variable climática (acrónimo usado)
  #anyi: año inicial
  #anyf: año final
  #nm: número de meses. (Si no se fija, se calcula por el no. de datos)
  #suf: sufijo opcional a añadir al nombre de la variable para leer los datos.
  #nref: no. (máximo) de estaciones de referencia
  #dz.max: límite superior de tolerancia de anomalías 
  #dz.min: límite inferior de tolerancia de anomalías
  #wd: Weight distance (km; distancia a la que el peso se reduce a la mitad;
  #    wd=0: todas las estaciones de referencia pesan lo mismo).
  #snht1: Umbral del SNHT en la primera pasada. (snht1=0 sólo rellena lagunas)
  #snht2: Umbral del SNHT en la segunda pasada. (=snht1 por defecto;
  #                                              =0 para saltarse la pasada)
  #tol: factor de tolerancia (por referencia disponible) para cortes en cadena
  #swa: Semi-Window Amplitude (no. de términos; 60 por defecto, o 365 si nm=0). 
  #trf: Transformar los datos? (0:no transformar; 1:log(x+1); >1:raíz trf)
  #mndat: Mínimo no. de datos para fragentar las series
  #gp: Parámetro de gráficos. 0=ninguno; 1=anomalías globales e histogramas;
  #    2=id+gráficos mensuales de diagnóstico; 3=id+gráficos de medias anuales
  #    móviles y correcciones; 4=id., pero con sumas anuales móviles
  #nclust: no. máximo de estaciones a usar en el análisis de agrupamiento
  #clustmethod: hierarchical clustering method
  #ini: fecha inicial (para datos diarios, en formato 'AAAA-MM-DD').
  #vmin, vmax: rango de valores permitidos en la variable climática.
  #force: forzar cortes aun con una sola referencia.
  #wz: factor de escala de z. El valor por defecto es apropiado si z se da en m
  #   y x,y en km. También sirve para sobreponderar z, o para hallar las
  #   distancias únicamente en el plano horizontal (wz=0).
  #grdcol: color de las retículas de los gráficos.
  #mapcol: color del mapa de fondo.
  #hires: mapa de fondo en alta resolución. (Poner a FALSE si no hace falta).
  #expl: Pasada exploratoria? (FALSE por defecto).
  #metad: Usar metadatos? En ese caso se fragmentarán las series en los lugares
  #   indicados en *brk.csv y solo se rellenarán lagunas. (FALSE por defecto).
  #sufbrk: sufijo a añadir al nombre de la variable para leer los metadatos.
  #   ('m' por defecto, para leer breaks detectados a escala mensual).
  #verb: Ver mensajes del proceso por pantalla (además de en el fichero *.txt).
  # -----------------------------------------------------------------  
  #- inicializaciones
  warnlt1=FALSE #flag to warn if there were means lower than 1
  #funciones auxiliares:
  datmed.mean <- function(x) mean(datmed[x])
  datmed.sd <- function(x) sd(datmed[x])
  #en caso de error, cerrar archivos de salida:
  options(error=cerrar)
  #establecer mxdif en función de la precisión elegida:
  if(is.na(mxdif)) mxdif=10^(-ndec)/2 #0.05 para un decimal
  verde <- hsv(.33,1,.6) #color muy usado
  #establecer altos umbrales de snht si es una pasada exploratoria:
  if(expl) snht1=snht2=9999
  #dz.min ha de ser negativo!:
  z=dz.min>0; if(sum(z)>0) dz.min[z] <- -dz.min[z]
  #- abrir fichero de bitácora y escribir cabecera
  archlog <- paste(varcli,"_",anyi,"-",anyf,".txt",sep="")
  sink(archlog,split=verb)
  cat("\nHOMOGEN() APPLICATION OUTPUT  (From R's contributed package 'climatol' ",climatol.version,")\n",sep='')
  cat("\n=========== Homogenization of ",varcli,", ",anyi,"-",anyf,". (",
    date(),")\n",sep="")
  time1 <- Sys.time() #tiempo al inicio del proceso
  cat("\nParameters:")
  arg <- names(formals()) #lista de los argumentos de la función
  for(i in 1:length(arg)) {
    cat(" ",arg[i],"=",sep="")
    cat(eval(as.symbol(arg[i])),sep=",")
  }
  cat("\n\n")
  #parámetros que han de tener 3 valores:
  k <- length(wd); if(k<3) wd <- c(rep(0,3-k),wd)
  k <- length(nref); if(k<3) nref <- c(nref,rep(nref[k],3-k))
  k <- length(dz.max); if(k<3) dz.max <- c(dz.max,rep(dz.max[k],3-k))
  k <- length(dz.min); if(k<3) dz.min <- c(dz.min,rep(dz.min[k],3-k))
  #etiquetas mensuales (de tres letras):
  mes3 <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
  #- lectura inicial de datos
  fbas <- sprintf('%s_%d-%d',varcli,anyi,anyf) #raíz nombres de fichero 
  if(is.na(suf)) fntr <- fbas #raíz nombres ficheros de entrada
    else fntr <- sprintf('%s-%s_%d-%d',varcli,suf,anyi,anyf)
  fiche <- sprintf('%s.est',fntr) #nombre del fichero de estaciones
  #leer coordenadas y nombres de las estaciones:
  est.c <- read.table(fiche,colClasses=c("numeric","numeric","numeric","character","character"))
  ne <- nrow(est.c) #no. de estaciones
  #comprobar si las coordenadas no están en grados:
  if(max(abs(est.c[,1]))>180 | max(abs(est.c[,2]))>90) deg <- FALSE
  else {
    deg <- TRUE
    if(gp>0) {
      if(requireNamespace("maps",quietly=TRUE)) mapok <- TRUE else mapok <- FALSE
      if(hires & requireNamespace("mapdata",quietly=TRUE)) maphr <- TRUE else maphr <- FALSE
    }
  }
  fichd <- sprintf('%s.dat',fntr) #nombre del fichero de datos
  dat <- scan(fichd,na.strings=na.strings) #lectura de los datos
  numdat <- length(dat) #no. de datos leídos
  nd <- numdat/ne #no. de datos por estación
  if(nd-floor(nd)>1e-16) {
    cat(ne,"stations read from",fiche,"\n")
    cat(numdat,"data read from",fichd,"\n")
    stop("The number of data is not multiple of the number of stations!")
  }
  dim(dat) <- c(nd,ne) #conversión de vector a matriz
  #- si hay algún término sin datos, emitir aviso y terminar
  numdat <- apply(!is.na(dat),1,sum) #no. de datos de cada término
  if(!min(numdat)) {
    z <- range(which(numdat==0))
    cat('\n',sum(numdat==0),' time steps (between terms ',z[1],' and ',z[2],') have missing data in all stations!\n',sep='')
    cat("Cannot continue. (Shorten the study period, or add series with data in those void terms)\n\n")
    stop()
  }
  na <- anyf-anyi+1 #no. de años
  nsy <- rep(0,na)   #no. de saltos por año
  if(is.na(nm)) { #calcular no. de datos por año y estación:
    z <- nd/na
    if(z>=1) nm <- ceiling(z)
    if(nm > 12) nm <- 0 #datos diarios
    else if(!nm%in%c(1,2,3,4,6,12)) {
      cat(sprintf('Computed nr. of data per year/station: %d.\n',nm))
      stop('Please set manually the right value of nm (one of 1,2,3,4,6,12)')
    }
  }
  if(is.na(swa)) { if(nm==0) swa=365 else swa=60 } #valores de swa por defecto
  else if(swa>=nd) swa <- ceiling(nd/4) #evitar semiventana demasiado grande
  #establecer valor de mndat, si no se especificó:
  if(is.na(mndat)) { 
    if(nm<=0) mndat <- swa/2 else mndat <- max(5,nm)
  }
  #- comprobar si hay valores fuera del rango permitido
  if(!is.na(vmin)) { #hay valores inferiores al mínimo permitido?
    n <- sum(dat<vmin,na.rm=TRUE)
    if(n) {
      dat[dat<vmin] <- vmin #corrección de los valores erróneos
      cat(n,"data forced to",vmin,":\n")
    }
  }
  if(!is.na(vmax)) { #hay valores superiores al máximo permitido?
    n <- sum(dat>vmax,na.rm=TRUE)
    if(n) {
      dat[dat>vmax] <- vmax #corrección de los valores erróneos
      cat(n,"data forced to",vmax,":\n")
    }
  }
  dat.o <- dat #copia de los datos originales
  if(nd<100) lw=3 #anchura de las barras de anomalías
  else if(nd<300) lw=2
  else lw=1
  #- generar vector temporal x
  if(is.na(ini)) ini <- sprintf('%d-01-01',anyi) #fecha inicial por defecto
  if(nm>0) x <- seq(as.Date(ini),length.out=nd,by=sprintf('%d months',12/nm))
  else x <- seq(as.Date(ini),length.out=nd,by='1 day')
  nei <- ne  #no. inicial de estaciones
  est.i <- est.c #datos iniciales de las estaciones
  nsp <- rep(0,nei)  #no. de cortes de cada estación original
  iest <- 1:ne   #índice señalando la serie original de cada subserie
  outan <- matrix(NA,nd,ne) #anomalías de los outliers
  #- if(gp>0), generar gráficos iniciales
  if(gp>0) {
    #activar salida gráfica a documento pdf, con rótulo inicial:
    pdfname <- paste(varcli,"_",anyi,"-",anyf,".pdf",sep="")
    pdf(pdfname,title=pdfname)
    plot(-1:1,-1:1,type="n",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
    text(0,0.4,sprintf("CLIMATOL %s",climatol.version),cex=4)
    text(0,-0.45,paste("Homogenization\ngraphic output of\n",varcli,"\n",anyi,"-",anyf,sep=""),cex=3)
    #datos disponibles en cada serie. (Solo si ne<=5*nclust, pues cuando el no.
    #de estaciones es muy grande, el gráfico ocupa mucho espacio):
    if(ne<=5*nclust) {
      image(x,1:ne,dat,xlab='Time',ylab='Stations',main=paste(varcli,'data availability'),col=4)
      grid(col=grdcol)
    } else cat(sprintf('\nMore than 5*%d stations: Per station data availability graph skipped\n\n',nclust))
    #no. de datos de cada término
    numdat <- apply(!is.na(dat),1,sum)
    plot(x,numdat,type="l",ylim=c(0,ne),col="blue",xlab='Time',ylab="Nr. of data",main=paste("Nr. of",varcli,"data in all stations"))
    grid(col=grdcol)
    abline(h=5,lty=2,col="green")
    abline(h=3,lty=2,col="red")
    #si hay algún término sin datos no podremos continuar:
    if(!min(numdat)) {
      stop("At least one term has missing data in all stations! (See the PDF graph)\nCannot continue. (Shorten the study period, or add series with data in those void terms)")
    }
    #boxplots de los datos de cada estación:
    if(nm>1) { #boxplots para cada nm (mensuales, etc)
      dim(dat) <- c(nm,na,ne) #dimensiones provisionales
      for(me in 1:nm) { #para cada mes
        z <- data.frame(dat[me,,])
        names(z) <- 1:ne
        #etiqueta del mes (si nm!=12, poner sólo el número):
        if(nm==12) labm <- mes3[me] else labm <- me
        labm <- paste(" (",labm,")",sep="")
        if(ne>nclust) hist(as.matrix(z),xlab=varcli,main=paste("Data values of ",varcli,labm,sep=""),col="wheat")
        else {
          boxplot(z,xlab="Stations",ylab="Values",main=paste("Data values of ",varcli,labm,sep=""),col="wheat",border=hsv(.7,1,.9))
          grid(col=grdcol)
          abline(h=0)
        }
      }
      dim(dat) <- c(nd,ne) #restablecer dimensiones de trabajo
    }
    else if(ne<=nclust) { #un sólo gráfico con los boxplots de cada estación
      z <- data.frame(dat)
      names(z) <- 1:ne
      boxplot(z,xlab="Stations",ylab="Values",main=paste("Data values of",varcli),col="wheat",border=hsv(.7,1,.9))
      grid(col=grdcol)
      abline(h=0)
    }
  }
  #- Transformar los datos? :
  if(trf>=1) {
    if(min(dat,na.rm=TRUE)<0) stop('Your data has negative values: cannot apply transformations!')
    if(trf==1) dat <- log1p(dat) else dat <- dat^(1/trf)
    if(is.na(vmin) | vmin<0) vmin <- 0 #evitar valores negativos
    if(mxdif>0.01) mxdif <- 0.01 #rebajar mxdif
  }
  if(gp>0) { #continuamos con los gráficos iniciales
    #histograma de todos los datos (distribución quasi-normal?)
    if(trf) main="Histogram of all (transformed) data"
    else main="Histogram of all data"
    hist(dat,xlab=varcli,main=main,col=hsv(.4,1,.8))
    #correlograma de series diferenciadas a nivel mensual (r <-> distancia)
    #(si hay más de nclust estaciones, sólo de una muestra aleatoria de nclust)
    if(ne>nclust) { splc <- sample(1:ne,nclust); nec <- nclust }
    else { splc <- 1:ne; nec <- ne }
    est.d <- matrix(NA,nec,nec) #matriz de distancias
    for(i in 1:(nec-1)) {
      for(j in (i+1):nec) {
        dx <- est.c[splc[i],1]-est.c[splc[j],1]
        dy <- est.c[splc[i],2]-est.c[splc[j],2]
        if(deg) {  #convertir grados a km
          dx <- dx*111*cos((est.c[splc[i],2]+est.c[splc[j],2])*pi/360)
          dy <- dy*111
        } else { #convertir m a km
          dx <- dx/1000.; dy <- dy/1000.
        }
        dz <- (est.c[splc[i],3]-est.c[splc[j],3])*wz
        d2 <- dx*dx+dy*dy+dz*dz #distancia cuadrática
        est.d[i,j] <- sqrt(d2) #distancia
        est.d[j,i] <- est.d[i,j]  #matriz simétrica
      }
    }
    data <- dat[,splc] #copia de los datos
    if(nm>1) { #calcular las series diferenciadas por meses
      dim(data) <- c(nm,na,nec) #dimensionamos por meses
      difd <- apply(data,c(1,3),diff)
      dim(difd) <- c(nd-nm,nec) #redimensionar
    }
    else difd <- diff(data) #series diferenciadas globalmente
    corm <- cor(difd,use="p") #matriz de correlaciones
    #eliminar |r|==1 (debidos a estaciones con sólo 2 datos en común):
    corm[corm==1] <- NA; corm[corm==-1] <- NA
    if(ne>nclust) main <- sprintf('Correlogram of %d sampled series (first differences)',nclust)
    else if(nm>0) main <- "Correlogram of first difference series"
    else main <- "Correlogram of the daily series"
    if(trf) main <- paste(main," (transformed)",sep="")
    xd <- as.vector(est.d); y <- as.vector(corm)
    plot(xd,y,xlim=c(0,max(est.d,na.rm=TRUE)),xlab="Distance (km)",ylab="Correlation coefficient",main=main,col="blue")
    grid(col=gray(.4))
    if(ne>2) {  #dendrograma de las estaciones
      dism <- dist(corm) #matriz de disimilaridad
      #si hay NA's en la matriz de disimilaridad, no intentar clustering
      if(!sum(is.na(dism))) {
        hc <- hclust(dism,method=clustmethod)
        if(ne>nclust) main <- paste("Dendrogram of",nclust,"sampled stations")
        else main <- "Dendrogram of station clusters"
        plot(hc,xlab="Stations",sub="",ylab="Dissimilarity",main=main)
        #clasificación de las estaciones cortando por la media más 1 desv.típ.
        #de las disimilaridades (salvo que el número de clases sea superior
        #a 9, en cuyo caso se incrementará el corte en .1):
        cutlev <- mean(hc$height)+sd(hc$height)
        repeat {
          ct <- cutree(hc,h=cutlev)
          nc <- length(levels(factor(ct)))
          if(nc<10) break
          cutlev <- cutlev + .1
        }
        if(nc>1) abline(h=cutlev,col="red",lty=2)
      } else { #un sólo grupo si dism tiene NA's
        cat("\nNA's in similarity matrix: NO CLUSTER ANALYSIS\n\n")
        nc <- 1; ct <- 1
      }
      #mapa de las estaciones:
      if(nc==1) { col="blue"; main=paste(varcli,"station locations") }
      else {
        col=rainbow(nc,1,.55)[ct]
        main=paste(varcli," station locations (",nc," clusters)",sep="")
      }
      if(deg) asp=1/(cos(mean(est.c[,2])*pi/180)) #relación de aspecto
      if(ne>nclust) { #dibujar símbolos si hay más de nclust estaciones
        #dibujar primero en negro las estaciones que no están en la muestra:
        if(deg) plot(est.c[-splc,1:2],asp=asp,xlab="Longitude (deg)",ylab="Latitude (deg)",pch='+',cex=.5,main=main)
        else plot(est.c[-splc,1:2],asp=1,xlab="X (km)",ylab="Y (km)",pch='+',cex=.5,main=main)
        #y ahora las estaciones de la muestra, con sus colores:
        points(est.c[splc,1:2],col=col,pch=ct)
      }
      else { #hasta nclust estaciones, poner el número
        if(deg) plot(est.c[,1:2],type="n",asp=1/(cos(mean(est.c[,2])*pi/180)),xlab="Longitude (deg)",ylab="Latitude (deg)",main=main)
        else plot(est.c[,1:2],type="n",asp=1,xlab="X (km)",ylab="Y (km)",main=main)
        text(est.c[,1:2],labels=1:ne,col=col)
      }
      grid(col=gray(.4))
      if(deg==TRUE) {
        if(maphr) try(maps::map('worldHires',col=mapcol,add=TRUE))
        else if(mapok) try(maps::map('world',col=mapcol,add=TRUE))
      }
    }
    rm(data,difd,corm) #borrar objetos temporales
  }
  #- if(gp==1), terminar (solo se deseaban los gráficos iniciales)
  if(gp==1) {
    graphics.off() #volcar el buffer del último gráfico
    cat("\nOnly the initial exploratory graphics were demanded.\nSee them in ",varcli,"_",anyi,"-",anyf,".pdf\n",sep="")
    sink() #cerrar el fichero de bitácora
    return(invisible())
  }
  #  Proceso de homogeneización en tres etapas:
  #  1) cortes por SNHT en ventanas móviles
  #  2) cortes por SNHT en toda la serie
  #  3) relleno de lagunas
  #- abrir archivos de outliers y breaks
  Fout <- file(sprintf('%s_out.csv',fbas),'w')
  write('"Code","Date","Observed","Suggested","Stand. deviations"',Fout)
  if(!metad) {
    Fbrk <- file(sprintf('%s_brk.csv',fbas),'w')
    write('"Code","Date","SNHT"',Fbrk)
  }
  #- cálculo de las matrices de distancias y rangos de proximidad
  est.d <- matrix(0,ne,ne) #matriz de distancias
  cat("Computing inter-station distances:")
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
      if(d2==0 & iest[i]!=iest[j]) est.d[i,j] <- 0.01 #evitar autocorrección 
      else est.d[i,j] <- sqrt(d2) #distancia
      est.d[j,i] <- est.d[i,j]  #matriz simétrica
    }
  }
  cat("\n")
  est.p <- t(apply(est.d,1,order)) #matriz de rangos de proximidad
  refhom <- substr(est.c[,4],1,1)=='*' #referencias homogéneas
  #- Estima inicial de medias y desv. típ. por diferencias o proporciones
  datmed <- apply(dat,1,mean,na.rm=TRUE) #serie media global
  refmed <- mean(datmed) #media global de referencia
  dat.m <- apply(dat,2,mean,na.rm=TRUE) #medias de partida
  if(std==3) {
    refstd <- sd(datmed) #desv. típica global de referencia
    dat.s <- apply(dat,2,sd,na.rm=TRUE) #desv. típ. de partida
  }
  switch(std,
    dat.m <- dat.m + refmed - apply(!is.na(dat),2,datmed.mean),
    dat.m <- dat.m * refmed / apply(!is.na(dat),2,datmed.mean),
    {dat.m <- dat.m + refmed - apply(!is.na(dat),2,datmed.mean)
     dat.s <- dat.s + refstd - apply(!is.na(dat),2,datmed.sd)},
    dat.m <- dat.m + refmed - apply(!is.na(dat),2,datmed.mean)
  )
  #- metad==TRUE? Leer *_brk.csv y cortar las series por donde indica
  if(metad) {
    cat('\nSplitting the series following the metadata file...:\n')
    if(sufbrk=='') fichbrk <- sprintf('%s_%d-%d_brk.csv',varcli,anyi,anyf)
    else fichbrk <- sprintf('%s-%s_%d-%d_brk.csv',varcli,sufbrk,anyi,anyf)
    brk <- read.csv(fichbrk,colClasses=c("character","character","numeric"))
    brk[,2] <- as.Date(brk[,2])
    nbrk <- nrow(brk); nn <- 0
    if(nbrk<1) break #sin breaks en el fichero!
    for(kb in nbrk:1) { #para cada break (en orden inverso):
      i <- match(brk[kb,1],est.c[,4]) #estación a cortar
      if(is.na(i)) {
        cat(sprintf('\nCode %s not found in station list; break skipped',brk[kb,1]))
        next
      }
      kp <- match(brk[kb,2],x) #posición de corte
      cat(sprintf('\n%s(%d) breaks at %s',est.c[i,4],i,format(x[kp])))
      if(sum(!is.na(dat[1:(kp-1),i])) < mndat) {
        dat[1:(kp-1),i] <- NA
        cat(" Fragment with less than",mndat,"data DELETED")
      }
      else if(sum(!is.na(dat[kp:nd,i])) < mndat) {
        dat[kp:nd,i] <- NA
        cat(" Fragment with less than",mndat,"data DELETED")
      }
      else {
        nn <- nn+1 #incrementamos el no. de nuevas series
        iest <- c(iest,iest[i]) #añadir índice a la serie original
        nsp[iest[i]] <- nsp[iest[i]]+1 #y también su no. de saltos
        if(nm>0) { #contar no. de saltos por año
          z <- 1 + floor((kp-1)/nm) #término anual del salto
          nsy[z] <- nsy[z] + 1 #no. de saltos por año
        }
        dat <- cbind(dat,rep(NA,nd)) #nueva columna de datos
        #pasar los datos a la nueva serie:
        dat[kp:nd,ne+nn] <- dat[kp:nd,i]
        dat[kp:nd,i] <- NA #borrar los datos pasados a la nueva serie
        #copiar las coordenadas y poner sufijo a indicativo y nombre:
        #(Usamos la lista original de estaciones, por si se borra alguna)
        z <- data.frame(est.i[iest[i],1:3],paste(est.i[iest[i],4],"-",1+nsp[iest[i]],sep=""),paste(est.i[iest[i],5],"-",1+nsp[iest[i]],sep=""))
        names(z) <- names(est.i)
        est.c <- rbind(est.c,z)
        switch(std,
          { dat.m[i] <- mean(dat[,i],na.rm=TRUE) + refmed - mean(datmed[!is.na(dat[,i])])
            dat.m <- c(dat.m, mean(dat[,ne+nn],na.rm=TRUE)+refmed-mean(datmed[!is.na(dat[,ne+nn])])) },
          { dat.m[i] <- mean(dat[,i],na.rm=TRUE) * refmed / mean(datmed[!is.na(dat[,i])])
            dat.m <- c(dat.m, mean(dat[,ne+nn],na.rm=TRUE)*refmed/mean(datmed[!is.na(dat[,ne+nn])])) },
          { dat.m[i] <- mean(dat[,i],na.rm=TRUE) + refmed - mean(datmed[!is.na(dat[,i])])
            dat.m <- c(dat.m, mean(dat[,ne+nn],na.rm=TRUE)+refmed-mean(datmed[!is.na(dat[,ne+nn])]))
            dat.s[i] <- sd(dat[,i],na.rm=TRUE) + refstd - sd(datmed[!is.na(dat[,i])])
            dat.s <- c(dat.s, sd(dat[,ne+nn],na.rm=TRUE)+refstd-sd(datmed[!is.na(dat[,ne+nn])])) },
          { dat.m[i] <- mean(dat[,i],na.rm=TRUE) + refmed - mean(datmed[!is.na(dat[,i])])
            dat.m <- c(dat.m, mean(dat[,ne+nn],na.rm=TRUE)+refmed-mean(datmed[!is.na(dat[,ne+nn])])) }
        )
      }
    }
    cat("\n\nUpdate number of series: ",ne,"+",nn,"= ")
    ne <- ne+nn  #actualizar el no. de estaciones
    cat(ne,"\n\n")
    snht1=0 #pasar directamente a relleno de lagunas
  }
  #- for (ks in 1:3) #(snht en ventanas, snht total, y relleno de lagunas)
  for (ks in 1:3) { #para cada etapa:
    #- if(snht1==0 & ks<3) next #no realizar cortes, sólo rellenar lagunas
    if(snht1==0 & ks<3) next
    if(ks==2) snht1 <- snht2 #umbral para SNHT en toda la serie
    cat("\n\n========== STAGE",ks)
    switch(ks,
      cat(" (SNHT on overlapping temporal windows) ===========\n\n"),
      cat(" (SNHT on the whole series) =======================\n\n"),
      cat(" (Final computation of all missing data) ==========\n\n")
    )
    #- calcular la matriz de pesos
    est.w <- matrix(1,nei,nei) #matriz de pesos
    if(wd[ks]>0) { #pesos diferentes de 1
      cat("Computing inter-station weights...")
      wd2 <- wd[ks]*wd[ks]
      for(i in 1:(nei-1)) {
        for(j in (i+1):nei) {
          est.w[i,j] <- wd2/(wd2+est.d[i,j]*est.d[i,j])
          est.w[j,i] <- est.w[i,j]  #matriz simétrica
        }
      }
      cat(' (done)\n\n')
    }
    #- if(gp>0), pintar rótulo separador de niveles
    if(gp>0) {
      plot(-1:1,-1:1,type="n",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
      text(0,0.4,paste("Stage",ks),cex=4)
      if(ks==1) text(0,-0.3,paste("Binary splits on",swa,"term\nstepped windows\nwith SNHT >",snht1,"\nand wd =",wd[ks],"km"),cex=3)
      else if(ks==2) text(0,-0.3,paste("Binary splits on\nwhole series\nwith SNHT >",snht1,"\nand wd =",wd[ks],"km"),cex=3)
      else text(0,-0.3,paste("Anomalies after\nmissing data\nrecalculation\nwith wd =",wd[ks],"km\n( swa =",swa,")"),cex=2.5)
    }
    #valores dependientes de la pasada:
    if(ks==3) aref <- TRUE else aref <- FALSE
    nrefk <- nref[ks]
    dz.maxk <- dz.max[ks]; dz.mink <- dz.min[ks]
    #- repetir hasta que no se corte ninguna serie
    repeat {
      #---------- Cálculo de anomalías y eliminación de datos anómalos:
      #- inicializar matrices dat.z|e|c oneref anom sanom mindist nrefs used
      dat.z <- matrix(NA,nd,ne) #datos observados (estandarizados)
      dat.e <- matrix(NA,nd,ne) #datos estimados (estandarizados)
      dat.c <- matrix(NA,nd,ne) #datos calculados (estimados, sin estand.)
      oneref <- matrix(FALSE,nd,ne) #sólo 1 referencia?
      anom <- matrix(NA,nd,ne) #anomalías
      sanom <- matrix(NA,nd,ne) #anomalías estandarizadas
      mindist <- matrix(NA,nd,ne) #distancias mínimas
      nrefs <- matrix(NA,nd,ne) #no. de referencias
      used <- matrix(FALSE,ne,ne) #flags de estaciones usadas
      #copia de trabajo de los datos:
      dat.d <- dat
      #- si hay algún término sin ningún dato, avisar y terminar
      numdat <- apply(!is.na(dat.d),1,sum)
      nmin=min(numdat)
      if(!nmin) {
        cat("\nThere are terms with NO DATA!:\n")
        for(j in which(numdat==0)) cat(format(x[j]),"\n")
        stop("Cannot continue! Shorten the study period, add series with data in the empty terms, or be more tolerant to outliers.")
      }
      #- eliminar las estaciones con menos de mndat datos
      numdat <- apply(!is.na(dat.d),2,sum)
      numdat[iest==0] <- NA #evitar repetir mensajes anteriores
      nmin=min(numdat,na.rm=TRUE)
      if(nmin<mndat) {
        cat("There are stations with less than ",mndat," data:\n",sep="")
        cat("Stations...:",formatC(which(numdat<mndat),0,4),"\n")
        cat("Nr. of data:",formatC(numdat[numdat<mndat],0,4),"\n")
        cat("These stations will be deleted in order to proceed:\n")
        z <- which(numdat<mndat); 
        for(idel in z) {
          cat(paste(est.c[idel,4],"(",idel,") ",est.c[idel,5],"  DELETED\n",
            sep=""))
        }
        cat('\n')
        #eliminar las estaciones con pocos datos:
        iest[z] <- 0 #anular los índices de estación original
      }
      dat.na <- is.na(dat.d) #índice de datos ausentes
      #- primera estima de medias y desv. típicas, o usar las anteriores
      #usar las medias y desviaciones típicas anteriores si existen:
      if(exists('dat.m0')) {
        dat.m <- dat.m0
        if(std==3) dat.s <- dat.s0
      }
      #- estandarizar los datos dat.d (obtener dat.z)
      switch(std,
        dat.z <- scale(dat.d,center=dat.m,scale=FALSE),     #std=1
        if(min(dat.m)<1) { warnlt1=TRUE
          z <- which(dat.m > 1)
          dat.z <- dat.d
          dat.z[,z] <- scale(dat.d[,z],center=FALSE,scale=dat.m[z]) }
        else dat.z <- scale(dat.d,center=FALSE,scale=dat.m), #std=2
        dat.z <- scale(dat.d,center=dat.m,scale=dat.s),     #std=3
        dat.z <- dat.d
      )
      #- ite=0 y repetir hasta estabilizar los datos estimados
      #proceso iterativo de estima de las medias de cada serie:  
      ite <- 0
      cat("Computation of missing data with outlier removal\n")
      cat('(Suggested data replacements are provisional)\n')
      if(length(dat)>10000000) cat('This process may take a very long time (many days)\n')
      else if(length(dat)>1000000) cat('This process may take a long time (many hours)\n')
      if(ks==3) cat("\nThe following lines will have one of these formats:\n")
      cat("Station(rank) Date: Observed -> Suggested (Standard dev.)\n")
      if(ks==3) cat("Iteration Max.mean.difference (Station_code)")
      repeat {
        ite <- ite+1
        #- ite+=1 y obtener las series estimadas (dat.e|c) con las vecinas
        #  actualizando used, nrefs y mindist:
        for(i in 1:ne) { #para cada estación
          if(!iest[i]) next #estación borrada
          ik <- iest[i] #índice de referencia estación inicial
          for(j in 1:nd) { #para cada dato
            se <- 0
            sw <- 0
            nr <- 0
            for(ir in 1:nei) { #para cada estación (posible referencia)
              kr <- est.p[ik,ir]
              krf <- which(iest==kr) #fragmentos de la referencia
              k <- which(!dat.na[j,krf]) #cuál tiene dato observado?
              if(length(k)!=1) next #ningún fragmento con dato
              k <- krf[k] #índice del fragmento con dato
              if(i==k) next #es la misma estación
              nr <- nr+1 #no. de referencias
              used[i,k] <- TRUE #marca de estación usada
              #distancia mínima (distancia al dato más próximo):
              if(nr==1) mindist[j,i] <- max(est.d[ik,kr],1)
              w <- est.w[ik,kr]
              se <- se + w * dat.z[j,k]
              sw <- sw + w
              if(nr>=nrefk) break #si no. máx. de referencias, terminar
            }
            if(!nr) { #sin referencia!
              dat.e[j,i] <- dat.z[j,i] #conservar el dato original
              nrefs[j,i] <- NA
            } else {
              nrefs[j,i] <- nr
              #si sólo hay una referencia, marcar para no corregir la serie:
              if(nr==1 & !is.na(oneref[j,i])) oneref[j,i] <- TRUE
              #no permitir datos negativos si std=2 (precipitación, etc):
              if(std==2 & se<0) se <- 0
              dat.e[j,i] <- se / sw #dato estimado (estandarizado)
            }
          }
        }
        #si hay NaN, convertirlos en NA. (Sucede a veces con std=2):
        n <- sum(is.nan(dat.e))
        if(n>0) {
          cat(n,"NaN's in dat.e ! (changing them to NA's...)\n")
          dat.e[is.nan(dat.e)] <- NA
        }
        #valores calculados por desestandarización de dat.e:
        switch(std,
          dat.c <- scale(dat.e,center=-dat.m,scale=FALSE),
          if(warnlt1) { z <- which(dat.m > 1)
            dat.c <- dat.e
            dat.c[,z] <- scale(dat.e[,z],center=FALSE,scale=1/dat.m[z]) }
          else dat.c <- scale(dat.e,center=FALSE,scale=1/dat.m), #std=2
          { dat.c <- scale(dat.e,center=FALSE,scale=1/dat.s)
            dat.c <- scale(dat.c,center=-dat.m,scale=FALSE) },
          dat.c <- dat.e
        )
        #- cálculo de anomalías (anom, sanom) y eliminación de outliers
        anom <- dat.z-dat.e #anomalías
        anom[dat.na] <- NA  #no arrastrar anomalías de datos rellenados
        #estandarizar las anomalías:
        anomm <- apply(anom,2,mean,na.rm=TRUE) #anomalías medias
        anoms <- apply(anom,2,sd,na.rm=TRUE) #desv. típicas de las anomalías
#       for(j in 1:ne) sanom[,j] <- (anom[,j]-anomm[j])/anoms[j] #an. estand.
        sanom <- scale(anom,center=anomm,scale=anoms)
        if(!expl & dz.maxk>.1) { #eliminar outliers
          elim <- sanom<dz.mink | sanom>dz.maxk #datos a eliminar
          elim[is.na(elim)] <- FALSE #eliminar los molestos NA
          nelim <- sum(elim) #no. de datos a eliminar
          if(nelim>0) { #eliminar los datos originales anómalos
            #listado de los datos a eliminar:
            for(i in 1:ne) {
              for(j in 1:nd) if(elim[j,i] & !is.na(oneref[j,i])) {
                outan[j,iest[i]] <- sanom[j,i] #guardar la anomalía del outlier
                do <- dat.d[j,i] #dato original
                dc <- dat.c[j,i] #dato calculado
                if(trf==1) { do <- expm1(do); dc <- expm1(dc) }
                else if(trf>1) { do <- do^trf; dc <- dc^trf }
                cat(sprintf('%s(%d) %s',est.c[i,4],i,format(x[j])))
                cat(": ",do," -> ",round(dc,ndec)," (",round(sanom[j,i],2),")",sep="")
                if(oneref[j,i]) { #no eliminar si sólo tenían una referencia!
                  cat(" Only 1 reference! (Unchanged)")
                  elim[j,i] <- FALSE
                  oneref[j,i] <- NA #evitar repetir este mismo mensaje
                }
                else { #escribir en Fout
                  write(c(est.c[iest[i],4],format(x[j]),round(do,ndec),
                  round(dc,ndec),round(sanom[j,i],2)),Fout,ncolumns=5,sep=',')
                }
                cat("\n")
              }
            }
            dat[elim] <- NA #eliminación de los datos anómalos
            dat.na[elim] <- TRUE #actualización índice de datos ausentes
          }
          else if(!aref) cat('(No detected outliers)\n')
        }
        #- relleno de las lagunas de datos
        dat.d[dat.na] <- dat.c[dat.na] 
        if(ite>1) {
          maxddif <- max(abs(dat.d-dat.d0),na.rm=TRUE) #máx. dat. dif.
          kmaxdif <- ceiling(which.max(abs(dat.d-dat.d0))/nd) #estación
        }
        dat.d0 <- dat.d #copia de los datos
        #- actualizar dat.m|s|z
        dat.m <- apply(dat.d,2,mean,na.rm=TRUE)
        if(std==3) dat.s <- apply(dat.d,2,sd,na.rm=TRUE)
        switch(std,
          dat.z <- scale(dat.d,center=dat.m,scale=FALSE), #std=1
          if(min(dat.m)<1) {
            warnlt1=TRUE
            z <- which(dat.m > 1)
            dat.z <- dat.d
            dat.z[,z] <- scale(dat.d[,z],center=FALSE,scale=dat.m[z])
          } else dat.z <- scale(dat.d,center=FALSE,scale=dat.m), #std=2
          dat.z <- scale(dat.d,center=dat.m,scale=dat.s), #std=3
          dat.z <- dat.d
        )
        #- if(!aref) break (no afinar los datos ausentes hasta el final)
        #  porque no parece necesario y alarga el tiempo de proceso 
        if(!aref) break
        #- if(ite>1), si los datos ya no varían, break
        if(ite>1) {
          cat(ite,' ',round(maxddif,ndec+2)," (",est.c[kmaxdif,4],")\n",sep="")
          if(maxddif<=mxdif | ite==maxite) {
            if(ite==maxite) cat("\nAverage computation skipped after",ite,"iterations\n")
            else cat("\n")
            break
          }
        }
      }
      #- guardar dat.m|s en dat0.m|s
      dat.m0 <- dat.m #copia de las medias
      if(std==3) dat.s0 <- dat.s #copia de las desv. típicas
      #restablecer valores de oneref:
      oneref[is.na(oneref)] <- TRUE
      #- if(aref==TRUE), repetir relleno de lagunas con autocorrección
      if(aref==TRUE) { # de las series fragmentadas:
        cat('Last series readjustment (please, be patient...)\n')
        #- obtener las series estimadas (dat.e, dat.c) con las vecinas
        #  y actualizar used[ne,ne], nrefs[nd,ne] y mindist[ne,ne]:
        for(i in 1:ne) { #para cada estación
          if(!iest[i]) next #estación borrada
          ik <- iest[i] #índice de referencia estación inicial
          for(j in 1:nd) { #para cada dato
            se <- 0
            sw <- 0
            nr <- 0
            for(ir in 1:nei) { #para cada estación (posible referencia)
              kr <- est.p[ik,ir]
              krf <- which(iest==kr) #fragmentos de la referencia
              k <- which(!dat.na[j,krf]) #cuál tiene dato observado?
              if(length(k)!=1) next #ningún fragmento con dato
              k <- krf[k] #índice del fragmento con dato
              if(i==k) next #es la misma estación
              nr <- nr+1 #no. de referencias
              used[i,k] <- TRUE #marca de estación usada
              #distancia mínima (distancia al dato más próximo):
              if(nr==1) mindist[j,i] <- max(est.d[ik,kr],1)
              w <- est.w[ik,kr]
              se <- se + w * dat.z[j,k]
              sw <- sw + w
              #si no. máx. de referencias, o autoreferencia, terminar:
              if(nr>=nrefk | (aref & ir==1)) break
            }
            if(!nr) { #sin referencia!
              dat.e[j,i] <- dat.z[j,i] #conservar el dato original
              nrefs[j,i] <- NA
            } else {
              nrefs[j,i] <- nr
              #si sólo hay una referencia, marcar para no corregir la serie:
              if(nr==1 & !is.na(oneref[j,i])) oneref[j,i] <- TRUE
              #no permitir datos negativos si std=2 (precipitación, etc):
              if(std==2 & se<0) se <- 0
              dat.e[j,i] <- se / sw #dato estimado (estandarizado)
            }
          }
        }
        #si hay NaN, convertirlos en NA. (Sucede a veces con std=2):
        n <- sum(is.nan(dat.e))
        if(n>0) {
          cat(n,"NaN's in dat.e ! (changing them to NA's...)\n")
          dat.e[is.nan(dat.e)] <- NA
        }
        #valores calculados por desestandarización de dat.e:
        switch(std,
          dat.c <- scale(dat.e,center=-dat.m,scale=FALSE),
          if(warnlt1) { 
            z <- which(dat.m > 1)
            dat.c <- dat.e
            dat.c[,z] <- scale(dat.e[,z],center=FALSE,scale=1/dat.m[z])
          } else dat.c <- scale(dat.e,center=FALSE,scale=1/dat.m), #std=2
          { dat.c <- scale(dat.e,center=FALSE,scale=1/dat.s)
            dat.c <- scale(dat.c,center=-dat.m,scale=FALSE) },
          dat.c <- dat.e
        )
        #- relleno de los datos ausentes
        dat.d[dat.na] <- dat.c[dat.na]
        if(!is.na(vmax)) dat.d[dat.d > vmax] <- vmax
        if(!is.na(vmin)) dat.d[dat.d < vmin] <- vmin
        #- cálculo final de dat.m|s|z
        dat.m <- apply(dat.d,2,mean,na.rm=TRUE) #medias
        if(std==3) dat.s <- apply(dat.d,2,sd,na.rm=TRUE) #desv. típicas
        switch(std, #datos estandarizados
          dat.z <- scale(dat.d,center=dat.m,scale=FALSE), #std=1
          if(min(dat.m)<1) { 
            warnlt1=TRUE
            z <- which(dat.m > 1)
            dat.z <- dat.d
            dat.z[,z] <- scale(dat.d[,z],center=FALSE,scale=dat.m[z])
          } else dat.z <- scale(dat.d,center=FALSE,scale=dat.m), #std=2
          dat.z <- scale(dat.d,center=dat.m,scale=dat.s), #std=3
          dat.z <- dat.d
        )
      }
      #- calcular los valores finales de las anomalías (anom, sanom)
      anom <- dat.z-dat.e #anomalías
      anom[dat.na] <- NA  #no arrastrar anomalías de datos rellenados!
      anomm <- apply(anom,2,mean,na.rm=TRUE) #anomalías medias
      anoms <- apply(anom,2,sd,na.rm=TRUE) #desv. típicas de las anomalías
      sanom <- scale(anom,center=anomm,scale=anoms) #anomalías estandarizadas
      #- ----------- Análisis de saltos en la media (binary split):
      #- if(ks>2) break (en la última etapa, solo relleno final de lagunas)
      if(ks>2) break
      #analizar los saltos en la media de las series, cortándolas cuando
      #el máximo snht del test supere el umbral (snht1):
      nn <- 0 #inic. no. de nuevas estaciones
      tVx <- rep(NA,ne) #máximos valores del shift test (por estación)
      kpx <- rep(NA,ne) #posiciones de los máximos tV (por estación)
      splt <- rep(0,ne)  #tV con que se cortaron las estaciones
      modif <- FALSE #inicialización modificación series
      cat("\nPerforming shift analysis on the",ne,"stations...\n")
      for(i in 1:ne) { #análisis de saltos en la media para cada estación
        if(refhom[i]) next
        y <- sanom[,i] #anomalías estandarizadas de la estación
        if(ks==1) { #análisis de saltos en ventanas móviles
          st <- snhtw(y,swa) #prueba SNHT en ventanas solapadas
          if(st[1]>0) tVx[i] <- st[1] else tVx[i] <- NA
          kpx[i] <- st[2]
        }
        else { #análisis de saltos en toda la serie
          st <- snht(y)
          if(sum(!is.na(st))>0) {
            tVx[i] <- max(st,na.rm=TRUE)
            kpx[i] <- which.max(st)
          }
        }
      }
      #- cortar las series cuyo snht máximo supere el umbral, de mayor a menor
      #  siempre que no se hayan usado series recién cortadas con snht similar
      tVxx <- max(tVx,na.rm=TRUE) #máximo tVx de todas las estaciones
      while(tVxx > snht1) {
        i <- which.max(tVx) #estación con el máximo snht
        #si i usó referencias cortadas con un snht demasiado grande, iniciar
        #una nueva iteración:
        if(max(splt[used[i,]])>tVxx*(1+tol*min(nr,sum(used[i,])))) break
        kp <- kpx[i] #posición del tVx en la estación i
        if(oneref[kp,i] & !force) { #no cortar con una sola referencia
          tVx[i] <- -1 #pasar el tVx de esta estación a -1
          tVxx <- max(tVx,na.rm=TRUE) #máximo tVx de las estaciones restantes
          next
        }
        cat(sprintf('\n%s(%d) breaks at %s (%.1f)',est.c[i,4],i,
          format(x[kp]),tVx[i]))
        write(sprintf('%s,%s,%.1f',est.c[iest[i],4],format(x[kp]),tVx[i]),
          Fbrk,ncolumns=3)
        #gráfico de anomalías con la posición del corte:
        if(gp>1) {
          y <- sanom[,i] #vector de anomalías de la estación
          ylab="Standardized anomalies (observed - computed)"
          tit <- paste(varcli," at ",est.c[i,4],"(",i,"), ",est.c[i,5],sep="")
          plot(x,y,type="h",lwd=lw,ylim=c(-5,5),main=tit,xlab='Time',ylab=ylab,col=hsv(.7,1,.9))
          grid(col=grdcol)
          abline(-3,0,lty=3,col=grdcol); abline(-5,0,lty=3,col=grdcol)
          lines(x,log10(nrefs[,i])-5,col='orange2')
          lines(x,log10(mindist[,i])-5,col=verde)
          mtext(" 1",4,las=1,adj=0,at=-5,col=verde)
          mtext(" 10",4,las=1,adj=0,at=-4,col=verde)
          mtext(" 100",4,las=1,adj=0,at=-3,col=verde)
          mtext("min.d.",4,las=1,adj=0,at=-5.4,col=verde)
          mtext(" (km)",4,las=1,adj=0,at=-2,col=verde)
          mtext("n.ref.",4,las=1,adj=0,at=-5.8,col='orange2')
          lines(rep(x[kp],2),c(-5,4.8),col="red",lty=2) #marca del corte
          text(x[kp],5,floor(tVxx))
        }
        if(nm>0) { #contar no. de saltos por año
          z <- 1 + floor((kp-1)/nm) #término anual del salto
          nsy[z] <- nsy[z] + 1 #no. de saltos por año
        }
        if(sum(!is.na(dat[1:(kp-1),i])) < mndat) {
          dat[1:(kp-1),i] <- NA
          cat(" Fragment with less than",mndat,"data DELETED\n")
        }
        else if(sum(!is.na(dat[kp:nd,i])) < mndat) {
          dat[kp:nd,i] <- NA
          cat(" Fragment with less than",mndat,"data DELETED\n")
        }
        else {
          nn <- nn+1 #incrementamos el no. de nuevas series
          iest <- c(iest,iest[i]) #añadir índice a la serie original
          nsp[iest[i]] <- nsp[iest[i]]+1 #y también su no. de saltos
          dat <- cbind(dat,rep(NA,nd)) #nueva columna de datos
          #pasar los datos a la nueva serie:
          dat[kp:nd,ne+nn] <- dat[kp:nd,i]
          dat[kp:nd,i] <- NA #borrar los datos pasados a la nueva serie
          #copiar las coordenadas y poner sufijo a indicativo y nombre:
          #(Usamos la lista original de estaciones, por si se borra alguna)
          z <- data.frame(est.i[iest[i],1:3],paste(est.i[iest[i],4],"-",1+nsp[iest[i]],sep=""),paste(est.i[iest[i],5],"-",1+nsp[iest[i]],sep=""))
          names(z) <- names(est.i)
          est.c <- rbind(est.c,z)
          switch(std,
            { dat.m0[i] <- mean(dat[,i],na.rm=TRUE) + refmed - mean(datmed[!is.na(dat[,i])])
              dat.m0 <- c(dat.m0, mean(dat[,ne+nn],na.rm=TRUE)+refmed-mean(datmed[!is.na(dat[,ne+nn])])) },
            { dat.m0[i] <- mean(dat[,i],na.rm=TRUE) * refmed / mean(datmed[!is.na(dat[,i])])
              dat.m0 <- c(dat.m0, mean(dat[,ne+nn],na.rm=TRUE)*refmed/mean(datmed[!is.na(dat[,ne+nn])])) },
            { dat.m0[i] <- mean(dat[,i],na.rm=TRUE) + refmed - mean(datmed[!is.na(dat[,i])])
              dat.m0 <- c(dat.m0, mean(dat[,ne+nn],na.rm=TRUE)+refmed-mean(datmed[!is.na(dat[,ne+nn])]))
              dat.s0[i] <- sd(dat[,i],na.rm=TRUE) + refstd - sd(datmed[!is.na(dat[,i])])
              dat.s0 <- c(dat.s0, sd(dat[,ne+nn],na.rm=TRUE)+refstd-sd(datmed[!is.na(dat[,ne+nn])])) },
            { dat.m0[i] <- mean(dat[,i],na.rm=TRUE) + refmed - mean(datmed[!is.na(dat[,i])])
              dat.m0 <- c(dat.m0, mean(dat[,ne+nn],na.rm=TRUE)+refmed-mean(datmed[!is.na(dat[,ne+nn])])) }
          )
        }
        #actualizar tVx y banderas para continuar el bucle:
        modif <- TRUE #marcar si se han modificado series
        splt[i] <- tVx[i] #tV de corte de la estación i
        tVx[i] <- 0 #anular el tVx de esta estación
        tVxx <- max(tVx,na.rm=TRUE) #máximo tVx de las estaciones restantes
      }
      if(nn) {
        cat("\n\nUpdate number of series: ",ne,"+",nn,"= ")
        ne <- ne+nn  #actualizar el no. de estaciones
        cat(ne,"\n\n")
        refhom <- c(refhom,rep(FALSE,nn)) #actualizar referencias homogéneas
      }
      #- sin nuevos cortes? histogramas de snht y break
      if(!nn & !modif) {
        if(gp>1) {
          #histograma de máximos tV, globales (sin 0's, que no son reales):
          z <- tVx[!is.na(tVx) & tVx!=0]
          main <- paste("Histogram of maximum SNHT (Stage ",ks,")",sep='')
##        if(sum(!is.na(z))) hist(z,breaks=min(0,floor(min(z))):max(20,ceiling(max(z))),xlab='SNHT',col="purple",main=main) 
          if(sum(!is.na(z))) hist(z,breaks=20,xlab='SNHT',col="purple",main=main)
          #si hay posibles cortes con 1 sóla referencia, colorear de rojo:
          if(min(z,na.rm=TRUE)<0) hist(z[z<0],breaks=1,col=2,add=TRUE)
          if(ks==2 | snht2<1) {
            #histograma de no. de cortes por estación:
            hist(nsp,breaks=0:max(9,max(nsp)+1)-.5,col="orange2",xlab="Number of splits",ylab="Number of stations",main="Number of splits per station")
            if(nm>0) { #frecuencias de fragmentación por años:
              w <- min(5,ceiling(400/na)) #anchura de las barras
              plot(anyi:anyf,nsy,type="h",lwd=w,col=2,ylim=c(0,max(10,max(nsy))),xlab="Years",ylab="Number of splits",main="Number of splits per year")
              grid(col=grdcol)
            }
          }
        }
        #lista de posibles cortes que solo tienen una referencia:
        z <- which(tVx<0)
        if(length(z)>0) {
          cat('Stations that could break but had only one reference:\n')
          print(est.c[z,4])
        }
        break #salir del bucle para ir al siguiente nivel
      }
    }
  }
  #------------ Fin de las tres fases de la homogeneización ---------
  #RMSE de los datos calculados:
  if(trf==1) z <- expm1(dat.c)
  else if(trf>1) z <- dat.c^trf
  else z <- dat.c
  zo <- dat.o[,iest]; zo[dat.na] <- NA
  rmse <- apply(z-zo,2,sd,na.rm=TRUE)
  #- gráficos de anomalías de las series homogeneizadas
  #  (con tVx máximos, ordenados por series originales):
  tVx <- rep(NA,ne) #(guardaremos los máximos tV finales)
  snhx <- rep(NA,ne) #(guardaremos los máximos SNHT finales)
  for(io in 1:nei) { #para cada serie original
    wi <- which(iest==io) #estaciones derivadas de la estación io
    lwi <- length(wi)
    if(!lwi) next #(estación totalmente borrada!)
    for(i in wi) { #para cada serie derivada de la original
      y <- sanom[,i] #anomalías estandarizadas de la estación
      if(gp>1) {
        ylab="Standardized anomalies (observed - computed)"
        tit <- paste(varcli," at ",est.c[i,4],"(",i,"), ",est.c[i,5],sep="")
        plot(x,y,type="h",lwd=lw,ylim=c(-5,5),main=tit,xlab='Time',ylab=ylab,col=hsv(.7,1,.9))
        grid(col=grdcol)
        abline(-3,0,lty=3,col=grdcol); abline(-5,0,lty=3,col=grdcol)
        lines(x,log10(nrefs[,i])-5,col='orange2')
        lines(x,log10(mindist[,i])-5,col=verde)
        mtext(" 1",4,las=1,adj=0,at=-5,col=verde)
        mtext(" 10",4,las=1,adj=0,at=-4,col=verde)
        mtext(" 100",4,las=1,adj=0,at=-3,col=verde)
        mtext("min.d.",4,las=1,adj=0,at=-5.4,col=verde)
        mtext(" (km)",4,las=1,adj=0,at=-2,col=verde)
        mtext("n.ref.",4,las=1,adj=0,at=-5.8,col='orange2')
      }
      #aplicar SNHTw y marcar su tV máximo (si >=1):
      st <- snhtw(y,swa); tVx[i] <- st[1]; zz <- floor(st[1])
      if(zz) {
        kp <- st[2]
        if(gp>1) {
          lines(rep(x[kp],2),c(-5,4.8),col=verde,lty=2) #marca máximo SNHTw
          text(x[kp],5,zz,col=verde) #valor
        }
      }
      #aplicar SNHT y marcar su máximo:
      st <- snht(y)
      if(sum(!is.na(st))>0) {
        kp <- which.max(st)
        snhx[i] <- round(max(st,na.rm=TRUE),1)
        zz <- floor(snhx[i])
        if(gp>1) {
          lines(rep(x[kp],2),c(-5,4.8),lty=4) #marca máximo SNHT
          text(x[kp],-5.2,zz) #valor
        }
      }
    }
  }
  #datos homogeneizados:
  if(trf==1) dah <- expm1(dat.d) #deshacer transformación logarítmica
  else if(trf>1) dah <- dat.d^trf #deshacer transformación raíz
  else dah <- dat.d
  dah <- round(dah,ndec) #redondear con el no. de decimales deseado
  #- gráficos de las series homogeneizadas y sus correcciones
  if(gp>2) {
    plot(-1:1,-1:1,type="n",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
    text(0,0.4,"Final graphics",cex=3.5)
    text(0,-0.3,"Adjusted series and\napplied corrections",cex=2.5)
    if(nm>0) xlab <- "Years" else xlab <- "Dates"
    old.par <- par(no.readonly=TRUE)
    layout(matrix(1:2,2,1,byrow=TRUE))
    par(las=1,cex=.8)
    for(i in 1:nei) { #para cada estación original
      wi <- which(iest==i) #estaciones derivadas de la estación i
      lwi <- length(wi)
      if(!lwi) next #(estación totalmente borrada!)
      if(lwi>1) vi <- TRUE else vi <- FALSE
      #filtros para valores anuales:
      if(nm>0) fltr <- rep(1,nm) else fltr <- rep(1,365)
      if(gp>3) ylab <- "Running annual totals"
      else {
        ylab <- "Running annual means"
        if(nm>0) fltr <- fltr/nm else fltr <- fltr/365
      }
      tit <- paste(varcli," at ",est.i[i,4],"(",i,"), ",est.i[i,5],sep="")
      yo <- as.vector(dat.o[,i]) #datos originales
      y <- dah[,wi] #datos homogeneizados
      par(mar=c(0,4,4,2),xaxt="n")
      matplot(x,filter(y,fltr),type="l",lty=1,col=2:20,ylab=ylab,main=tit)
      lines(x,filter(yo,fltr))
      grid(col=grdcol)
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
        plot(x,yd[,1],type="n",ylim=ylim,ylab=ylab,xlab='Time')
      }
      else {
        if(trf) ylim <- c(floor(min(yd,na.rm=TRUE)),ceiling(max(yd,na.rm=TRUE)))
        plot(x,yd,type="n",ylim=ylim,ylab=ylab,xlab='Time')
      }
      matlines(x,yd,type="l",lty=1,col=2:20)
      grid(col=grdcol)
    }
    par(old.par)
  }
  if(snht1>0) cat("\n======== End of the homogenization process, after ")
  else cat("\n======== End of the missing data filling process, after ")
  cat(format(round(Sys.time()-time1,2)),'\n')
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
  #prueba SNHT de cada estación
  cat("\nSNHT: Standard normal homogeneity test (on anomaly series)\n")
  print(summary(round(snhx,1)))
  #errores típicos de las estimas (sin estandarizar):
  cat("\nRMSE: Root mean squared error of the estimated data\n")
  zz <- summary(rmse)
  print(zz)
  sedec <- max(1,2-ceiling(log10(zz[4]))) #no. de decimales de RMSE
  pod <- floor(100*(nd-apply(dat.na,2,sum))/nd) #porcentaje de datos originales
  cat("\nPOD: Percentage of original data\n")
  print(summary(pod))
  #- imprimir resumen de resultados
  cat("\n")
  print(data.frame(ACmx=round(sac,2),SNHT=snhx,RMSE=round(rmse,sedec),POD=pod,Code=est.c[,4],Name=est.c[,5]),right=FALSE)
  #averiguar qué estaciones derivadas funcionan al final del periodo:
  cur <- apply(!is.na(dat[(nd-mndat+1):nd,]),2,sum) #últimos mndat términos
  cur[cur>0] <- 1
  #añadir cuatro nuevas columnas a la tabla de estaciones (porcentaje de datos
  #originales, estación original, si funciona actualmente, y SNHT):
  est.c <- cbind(est.c,pod,iest,cur,snhx)
  #- if(gp>1), últimos gráficos (hist. de anomalías y snht; calidad/singular.)
  if(gp>1) {
    #histograma de las anomalías (las de los outliers, en rojo):
    main <- "Histogram of normalized anomalies"
    z <- hist(c(sanom,outan),plot=FALSE)
    zx <- z$breaks
    zy <- z$counts; zy[zy==0] <- NA
#Esto da Error en if (logy && !is.null(ylim) && min(ylim) <= 0) stop("log scale error: 'ylim' <= 0") : valor ausente donde TRUE/FALSE es necesario
#   ymax <- max(zy)
#   barplot(zy,log='y',space=0,ylab='Frequency',col='green',main=main,xlab='Anomalies (standard deviations)',ylim=c(.1,ymax))
    barplot(zy,log='y',space=0,ylab='Frequency',col='green',main=main,xlab='Anomalies (standard deviations)')
    axis(1,1:length(zx)-1,labels=as.character(zx),las=2)
    if(sum(!is.na(outan))) { #repintar las frec. de outan en rojo
      zy <- hist(outan,breaks=zx,plot=FALSE)$counts; zy[zy==0] <- NA
#     barplot(zy,log='y',space=0,col=hsv(0,.75),add=TRUE,ylim=c(.1,ymax))
      barplot(zy,log='y',space=0,col=hsv(0,.75),add=TRUE)
    }
    #histograma de tVx:
    z <- tVx; main <- "Histogram of maximum windowed SNHT"
##    if(sum(!is.na(z))) hist(z,breaks=min(0,floor(min(z))):max(20,ceiling(max(z))),xlab="tVx",col=verde,main=main)
    if(sum(!is.na(z))) hist(z,breaks=20,xlab="SNHT",col=verde,main=main)
    #histograma de SNHT:
    z <- snhx; main <- "Histogram of maximum global SNHT"
##    if(sum(!is.na(z))) hist(z,breaks=min(0,floor(min(z))):max(20,ceiling(max(z))),xlab="SNHT",col="purple",main=main)
    if(sum(!is.na(z))) hist(z,breaks=20,xlab="SNHT",col="purple",main=main)
    #gráfico de calidad/singularidad:
    plot(rmse,snhx,type="n",xlim=c(0,max(1,max(rmse,na.rm=TRUE))),ylim=c(0,max(50,max(snhx,na.rm=TRUE))),xlab="RMSE",ylab="SNHT",main="Station's quality/singularity")
    grid(col=grdcol)
    text(rmse,snhx,col=hsv(.7,1,.9))
  }
  if(gp>0) graphics.off() #cerrar la salida gráfica
  #- grabar los resultados en un fichero rda
  kelim <- rev(which(iest==0)) #estaciones eliminadas (en orden inverso)
  nelim <- length(kelim) #no. de estaciones eliminadas
  if(nelim>0) { #si se eliminaron estaciones:
    #ajustar los índices de estación original:
    for(ke in 1:nelim) iest[iest>kelim[ke]] <- iest[iest>kelim[ke]] - 1
    dat <- dat.o[,iest[1:nei]>0] #datos originales sin series eliminadas
    nei <- sum(iest[1:nei]>0) #no. de estaciones originales no eliminadas
    dah <- dah[,iest>0] #datos homogeneizados sin series eliminadas
    ne <- sum(iest>0) #no. de estaciones homogeneizadas no eliminadas
    est.c[,7] <- iest #actualizar los índices de estación original
    est.c <- est.c[iest>0,] #lista de estaciones homogeneizadas no eliminadas
  }
  else dat <- dat.o
  if(nm>0) {
    dim(dat) <- c(nm,na,nei)
    dim(dah) <- c(nm,na,ne)
  }
  names(est.c) <- c('X','Y','Z','Code','Name','pod','ios','ope','snht')
  rownames(est.c) <- 1:ne
  save(dat,dah,est.c,nd,ne,nei,nm,ndec,std,ini, file=sprintf('%s.rda',fbas))
  #ordenar archivos de outliers y breaks:
  if(!metad) { 
    close(Fbrk)
    brk <- read.csv(sprintf('%s_brk.csv',fbas),colClasses=c("character","character","numeric"))
    brk <- brk[order(brk[,1],brk[,2]),]
    write.csv(brk,sprintf('%s_brk.csv',fbas),row.names=FALSE)
  }
  close(Fout)
  out <- read.csv(sprintf('%s_out.csv',fbas),colClasses=c("character","character","numeric","numeric","numeric"))
  out <- out[order(out[,1],out[,2]),]
  write.csv(out,sprintf('%s_out.csv',fbas),row.names=FALSE)
  cat("\n----------- Generated output files: -------------------------\n\n")
  cat(sprintf('%s.txt :  This text output',fbas),'\n')
  cat(sprintf('%s_out.csv :  List of corrected outliers',fbas),'\n')
  cat(sprintf('%s_brk.csv :  List of corrected breaks',fbas),'\n')
  if(gp>0) cat(sprintf('%s.pdf :  Diagnostic graphics',fbas),'\n')
  cat(sprintf('%s.rda :  Homogenization results.',fbas))
  cat(' Postprocess with (examples):\n')
  cat(sprintf('   dahstat(\'%s\',%d,%d) #get averages in file %s-me.csv',varcli,anyi,anyf,fbas),'\n')
  cat(sprintf('   dahstat(\'%s\',%d,%d,stat=\'tnd\') #get OLS trends and their p-values',varcli,anyi,anyf),'\n')
  cat(sprintf('   dahgrid(\'%s\',%d,%d,grid=YOURGRID) #get homogenized grids',varcli,anyi,anyf),'\n')
  cat('   ... (See other options in the package documentation)\n\n')
  if(warnlt1) cat(' \n=====>> Warning!!! There were series with means lower than 1.\nThis may jeopardize the homogenization results when std=2 has been set.\nPlease, revise these results and, if unsatisfactory, try to run homogen()\nagain after having multiplied all your input data by a suitable factor\nin order to prevent any mean from being lower than 1.\n\n')
  sink() #cerrar bitácora
}

#- homogsplit.- Apply homogen() on overlapping split areas.
homogsplit <- function(varcli, anyi, anyf, xc=NULL, yc=NULL, xo=.5, yo=.38,
  maponly=FALSE, suf=NA, nm=NA, nref=c(10,10,4), swa=NA, std=3, ndec=1,
  dz.max=5, dz.min=-dz.max, wd=c(0,0,100), snht1=25, snht2=snht1, tol=.02,
  mxdif=NA, force=FALSE, wz=.001, trf=0, mndat=NA, gp=3, ini=NA,
  na.strings="NA", maxite=50, vmin=NA, vmax=NA, nclust=100,
  clustmethod='ward.D2', grdcol=grey(.5), mapcol=grey(.65), hires=TRUE,
  expl=FALSE, metad=FALSE, sufbrk='m', verb=TRUE) {
  
  #output files:
  f.bas <- sprintf('%s_%d-%d',varcli,anyi,anyf) #base name
  f.txt <- sprintf('%s.txt',f.bas) #console output
  f.rda <- sprintf('%s.rda',f.bas) #homogenization objects
  if(file.exists(f.txt)) file.rename(f.txt,sprintf('%s.bak',f.txt))
  if(file.exists(f.rda)) file.rename(f.rda,sprintf('%s.bak',f.rda))
  #................ process:
  est.c <- read.table(sprintf('%s.est',f.bas),colClasses=c("numeric","numeric","numeric","character","character"))
  nei <- nrow(est.c); na <- anyf-anyi+1
  dat <- scan(sprintf('%s.dat',f.bas))
  nd <- length(dat)/nei
  dim(dat) <- c(nd,nei)
  #check whether coordinates are in degrees:
  if(max(abs(est.c[,1]))>180 | max(abs(est.c[,2]))>90) deg <- FALSE
  else {
    deg <- TRUE
    if(gp>0) {
      if(requireNamespace("maps",quietly=TRUE)) mapok <- TRUE else mapok <- FALSE
      if(hires & requireNamespace("mapdata",quietly=TRUE)) maphr <- TRUE else maphr <- FALSE
    }
  }
  if(is.null(xc) | is.null(yc)) { #plot sites to help choosing split borders:
    plot(est.c[,1:2],xlab='X',ylab='Y'); grid(col=grdcol)
    cat('mean(x)=',mean(est.c[,1]),';   mean(y)=',mean(est.c[,2]),'\n')
    if(deg) cat('Aspect ratio=',1/cos(mean(est.c[,2])*pi/180),'\n')
    cat('Choose the cut x and y values and call homogsplit() again specifying\nthem in vectors xc and yc.\n')
  }
  else { #apply homogen() on the selected areas:
    est0 <- est.c; dat0 <- dat #keep whole original data for selections
    pne <- 0 #previous number of stations
    nxc <- length(xc) #nr. of x cut borders
    nyc <- length(yc) #nr. of y cut borders
    #save a map of available stations and split areas:
    f.map <- sprintf('%s-map.pdf',f.bas)
    main=paste('Split areas of the',nei,'available',varcli,'stations')
    pdf(f.map)
    if(deg) {
      asp=1/(cos(mean(range(est.c[,2]))*pi/180)) #aspect ratio
      plot(est.c[,1:2],pch='+',col=hsv(.6,.7,1),asp=asp,xlab="Longitude (deg)",ylab="Latitude (deg)",main=main)
      if(maphr) try(maps::map('worldHires',add=TRUE))
      else if(mapok) try(maps::map('world',add=TRUE))
    } else plot(est.c[,1:2],pch='+',asp=1,xlab="X (km)",ylab="Y (km)",main=main)
    grid(col=gray(.4))
    abline(h=yc,col=2); abline(v=xc,col=2)
    abline(h=yc+yo,col=3); abline(h=yc-yo,col=3)
    abline(v=xc+xo,col=3); abline(v=xc-xo,col=3)
    mtext(paste('xc=',paste(xc,collapse=','),'   xo=',xo,sep=''),3)
    mtext(paste('yc=',paste(yc,collapse=','),'   yo=',yo,sep=''),4)
    graphics.off()
    if(maponly) {
      cat('Split areas and station map saved as',f.map,'\n')
      cat('No further action required\n')
      return()
    }
    #-------------------- homogenize overlapping areas: -------------
    #open console output file:
    sink(f.txt,split=TRUE)
    cat("\nHOMOGSPLIT() APPLICATION OUTPUT  (From R's contributed package 'climatol' ",climatol.version,")\n",sep='')
    cat("\n=========== Homogenization of ",varcli,", ",anyi,"-",anyf,". (",
      date(),")\n",sep="")
    time1 <- Sys.time() #time when starting homogenization
    cat("\nParameters:")
    arg <- names(formals()) #lista de los argumentos de la función
    for(i in 1:length(arg)) {
      cat(" ",arg[i],"=",sep="")
      cat(eval(as.symbol(arg[i])),sep=",")
    }
    cat("\n\n")
    noa <- 0 #nr. of overlapping areas
    ioa <- rep(0,nei) #index of overlapping areas assigned to each station
    listest <- list()
    #select target and overlapping stations:
    for(i in 1:(nyc+1)) { #for every y interval
      if(i==1) { #first y interval
        itg <- est0[,2] < yc[i] # y target stations
        iov <- est0[,2] < yc[i]+yo # y overlapping+target stations
      } else if(i>nyc) { #last y interval
        itg <- est0[,2] > yc[nyc] # y target stations
        iov <- est0[,2] > yc[nyc]-yo # y overlapping+target stations
      } else { #intermediate y intervals
        itg <- est0[,2]>yc[i-1] & est0[,2]<yc[i] # y target stations
        iov <- est0[,2]>yc[i-1]-yo & est0[,2]<yc[i]+yo # y ov.+tg. stations
      }
      for(j in 1:(nxc+1)) { #for every x interval
        cat('\n==================================================\n\n')
        cat('              AREA ',i,j,'\n')
        cat('\n==================================================\n\n')
        if(j==1) { #first x interval
          jtg <- est0[,1] < xc[j] # x target stations
          jov <- est0[,1] < xc[j]+xo # x overlapping+target stations
        } else if(j>nxc) { #last x interval
          jtg <- est0[,1] > xc[nxc] # x target stations
          jov <- est0[,1] > xc[nxc]-xo # x overlapping+target stations
        } else { #intermediate x intervals
          jtg <- est0[,1]>xc[j-1] & est0[,1]<xc[j] # x target stations
          jov <- est0[,1]>xc[j-1]-xo & est0[,1]<xc[j]+xo # x ov.+tg. sta.
        }
        #define station selections based on current x and y intervals:
        est.tg <- itg & jtg
        if(sum(est.tg)==0) { #no target stations in current selection
          cat('No target stations in this area\n')
          next
        }
        est.ov <- iov & jov
        if(exists('est.ov0')) { #previous selection with too few stations?
          est.tg <- est.tg | est.tg0 #add previous few stations
          est.ov <- est.ov | est.ov0 #add previous few stations
          rm(est.ov0,est.tg0)
        }
        if(sum(est.ov)<10) { #too few stations in current selection?
          est.ov0 <- est.ov #keep selection to be added to the next
          est.tg0 <- est.tg #keep selection to be added to the next
          cat('Only',sum(est.ov),'stations in this area:\n')
          if(i==nyc+1 & j==nxc+1) cat('As this is the last area, they will not be homogenized.\nPlease, choose new cutting lines to include them in a broader area,\nsince this homogenization results will have inconsistent number of stations.\n')
          else cat('They will be included in the next selection.\n')
          next #go to next selection
        }
        #generate temporary files for current overlapping area:
        noa <- noa + 1 #increase nr. of overlapping areas
        ioa[est.tg] <- noa #index of assigned area
        var <- sprintf('%s-%d',varcli,noa)
        basef <- sprintf('%s_%d-%d',var,anyi,anyf) #base file name
        write.table(est0[est.ov,],sprintf('%s.est',basef),row.names=FALSE,
          col.names=FALSE)
        write(dat0[,est.ov],sprintf('%s.dat',basef),ncolumns=max(c(10,nm),na.rm=TRUE))
        #homogenize current overlapping area:
        homogen(var, anyi, anyf, suf=suf, nm=nm, nref=nref, std=std, swa=swa,
        ndec=ndec, dz.max=dz.max, dz.min=dz.min, wd=wd, snht1=snht1,
        snht2=snht2, tol=tol, mxdif=mxdif, force=force, wz=wz, trf=trf,
        mndat=mndat, gp=gp, ini=ini, na.strings=na.strings, maxite=maxite,
        vmin=vmin, vmax=vmax, nclust=nclust, clustmethod=clustmethod,
        grdcol=grdcol, mapcol=mapcol, hires=hires, expl=expl, metad=metad,
        sufbrk=sufbrk, verb=verb)
      }
    }
    cat("\n======== End of homogenization of overlapping areas, after ")
    cat(format(round(Sys.time()-time1,2)),'\n')
    #joint areal homogenization files:
    for(i in 1:noa) {
      f.baux <- sprintf('%s-%d_%d-%d',varcli,i,anyi,anyf)
      load(sprintf('%s.rda',f.baux))
      dim(dah) <- c(nd,ne)
      #keep only data from inner area (not from overlapping margins):
      qcod <- est0[ioa==i,4] #original station codes
      qest <- match(qcod,est.c[,4]) #indexes of original stations
      sel <- est.c[,7]%in%qest #selected stations and their derivatives
      dah <- dah[,sel]; est.c <- est.c[sel,]; ne <- sum(sel)
      if(i==1) {
        zdah <- dah; zest.c <- est.c; zne <- ne
        zout <- read.csv(sprintf('%s_out.csv',f.baux),colClasses=c("character","character","numeric","numeric","numeric"))
        zbrk <- read.csv(sprintf('%s_brk.csv',f.baux),colClasses=c("character","character","numeric"))
      } else {
        zdah <- cbind(zdah,dah)
        zest.c <- rbind(zest.c,est.c)
        zne <- zne+ne
        zout <- rbind(zout,read.csv(sprintf('%s_out.csv',f.baux),colClasses=c("character","character","numeric","numeric","numeric")))
        zbrk <- rbind(zbrk,read.csv(sprintf('%s_brk.csv',f.baux),colClasses=c("character","character","numeric")))
      }
    }
    dat <- dat0; dah <- zdah; est.c <- zest.c; ne <- zne 
    if(nm>0) {
      dim(dat) <- c(nm,nd/nm,nei)
      dim(dah) <- c(nm,nd/nm,ne)
    }
    save(dat,dah,est.c,nd,ne,nei,nm,ndec,std, file=f.rda)
    out <- zout[order(zout[,1],zout[,2],zout[,3]),]
    brk <- zbrk[order(zbrk[,1],zbrk[,2],zbrk[,3]),]
    out <- out[row.names(unique(out[,1:3])),]
    brk <- brk[row.names(unique(brk[,1:3])),]
    write.csv(out,sprintf('%s_out.csv',f.bas),row.names=FALSE)
    write.csv(brk,sprintf('%s_brk.csv',f.bas),row.names=FALSE)
    cat("\n----------- Generated output files: -------------------------\n\n")
    cat(sprintf('%s.txt :  This text output\n',f.bas))
    cat(sprintf('%s_out.csv :  List of corrected outliers\n',f.bas))
    cat(sprintf('%s_brk.csv :  List of corrected breaks\n',f.bas))
    if(gp>0) cat(sprintf('%s-*_%d-%d.pdf :  Diagnostic graphics (one file per area)\n',varcli,i,anyi,anyf))
    cat(sprintf('%s :  Map of specified areas\n',f.map))
    cat(sprintf('%s.rda :  Homogenization results.',f.bas))
    cat(' Postprocess with (examples):\n')
    cat(sprintf('   dahstat(\'%s\',%d,%d) #get averages in file %s-me.csv\n',varcli,anyi,anyf,f.bas))
    cat(sprintf('   dahstat(\'%s\',%d,%d,stat=\'tnd\') #get OLS trends and their p-values\n',varcli,anyi,anyf))
    cat(sprintf('   dahgrid(\'%s\',%d,%d,grid=YOURGRID) #get homogenized grids\n',varcli,anyi,anyf))
    cat('   ... (See other options in the package documentation)\n\n')
    sink() #close console output file
  }
}

#- outrename.- Append a suffix to the output files, to avoid overwrites.
outrename <- function(varcli, anyi, anyf, suffix, restore=FALSE) {
  #if restore=TRUE, the suffix will be removed! 
  fbn <- paste(varcli,"_",anyi,"-",anyf,sep="") #original file base name
  #destination file base name:
  fbn2 <- paste(varcli,"-",suffix,"_",anyi,"-",anyf,sep="")
  for(ext in c(".txt",".rda",".pdf")) {
    if(restore) file.rename(paste(fbn2,ext,sep=""),paste(fbn,ext,sep=""))
    else file.rename(paste(fbn,ext,sep=""),paste(fbn2,ext,sep=""))
  }
  if(restore) {
    file.rename(sprintf('%s_out.csv',fbn2),sprintf('%s_out.csv',fbn))
    file.rename(sprintf('%s_brk.csv',fbn2),sprintf('%s_brk.csv',fbn))
  } else {
    file.rename(sprintf('%s_out.csv',fbn),sprintf('%s_out.csv',fbn2))
    file.rename(sprintf('%s_brk.csv',fbn),sprintf('%s_brk.csv',fbn2))
  }
  return(invisible())
}

#- snht.- Standard Normal Homogeneity Test (Alexandersson)
snht <- function(x,nmt=3) {
#nmt: no. mínimo de términos de cada muestra
  n <- length(x)
  Tsnht <- rep(NA,n)
  if(n<nmt*2) return(Tsnht) #insuficientes datos
  z <- (x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)
  yc <- FALSE #ya calculado?
  for(i in (nmt+1):(n-nmt)) { #(despreciar los primeros y últimos nmt términos)
    if(is.na(x[i]) & yc) next #test ya calculado
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

#- snhtw.- SNHT para ventanas solapadas de 2*nt términos válidos.
snhtw <- function(x,nt=48) {
  ntt <- length(x) #no. total de términos de la serie
  ntv <- sum(!is.na(x)) #no. de términos válidos de la serie
  if(2*nt>ntv) return(c(0,0)) #no hay suficientes datos válidos para la prueba
  tV <- 0 #inicialización del tV máximo a devolver
  pk <- 0 #inicialización de la posición a devolver
  #inicialización de los límites muestrales (a1-b1, a2-b2):
  k <- 1; while(k<ntt & is.na(x[k])) k <- k+1; a1 <- k
  n<-1; while(n<nt & k<ntt) { k <- k+1; if(!is.na(x[k])) n <- n+1; }
  b1 <- k
  k <- k+1; while(k<ntt & is.na(x[k])) k <- k+1; a2 <- k
  n<-1; while(n<nt & k<ntt) { k <- k+1; if(!is.na(x[k])) n <- n+1; }
  b2 <- k
  #aplicación de SNHT a las ventanas solapadas:
  repeat {
    st <- snht(x[a1:b2])
    stx <- max(st,na.rm=TRUE)
    if(stx>tV) { tV <- stx; pk <- which.max(st)+a1-1 }
    if(b2==ntt) return(c(tV,pk))
    #desfasar las ventanas hacia adelante:
    a1 <- a2; b1 <- b2
    k <- b2+1; while(k<ntt & is.na(x[k])) k <- k+1
    if(is.na(x[k])) return(c(tV,pk)) else a2 <- k
    n<-1; while(n<nt & k<ntt) { k <- k+1; if(!is.na(x[k])) n <- n+1; }
    b2 <- k
  }
}

