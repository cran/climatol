#  quizá se pueda arreglar tratando esas series como si fueran diarias...
#  nref[3]=nref[2] maxdif=mxdif

#depurdat.R.- Depuración y homogeneización de series climatológicas.
#(Most comments in Spanish; sorry!)

climatol.version <- '3.1'

#- cerrar.- Cerrar los archivos de salida.
cerrar <- function() {
  sink()
  graphics.off()
}

#- daily2climatol.- Convert DAILY data files to CLIMATOL input format.
daily2climatol <- function(stfile, stcol=1:6, datcol=1:4, varcli,
  anyi=NA, anyf=NA, mindat=365, sep='', dec='.', na.strings='NA',
  header=FALSE) {
  #stfile: file with the file names and station coordinates, codes and names
  #stcol: columns in stfile holding file names, longitudes, latitudes, elevations and station codes and names (0 when missing)
  #datcol: columns in data files holding year,month,day,value (default to 1:4) 
  #varcli: acronym of the climatic variable under study
  #anyi: first year to study (defaults to the first year available in data)
  #anyf: last year to study (defaults to the last year available in data)
  #mindat: minimum required number of data per station
  #sep: data separator ('' by default, meaning any white space)
  #dec: decimal point ('.' by default)
  #na.strings: strings coding missing data ('NA' by default)
  #header: does files have a header? (FALSE by default)
  #read stations file:
  st <- read.table(stfile,as.is=TRUE,sep=sep,dec=dec,header=header)
  ne <- nrow(st) #nr. of stations
  if(is.na(anyi) | is.na(anyf)) { #check the time period of the data:
    cat('\nChecking the period covered by the data...\n')
    inidate <- as.Date('3000-12-31'); enddate <- as.Date('0001-01-01')
    for(i in 1:ne) { #for every station
      cat(' ',i)
      d <- read.table(st[i,stcol[1]],sep=sep,dec=dec,header=header,na.strings=na.strings)
      dates <- as.Date(sprintf('%d-%02d-%02d',d[,datcol[1]],d[,datcol[2]],d[,datcol[3]]))
      rdates <- range(dates,na.rm=TRUE) #range of dates in the file
      nadates <- is.na(dates)
      if(sum(nadates)>0) {
        cat('Abnormal dates found in file',st[i,stcol[1]],':\n')
        print(d[nadates,])
      }
      dates[is.na(d[,datcol[4]])] <- NA #remove dates without data
      rdates <- range(dates,na.rm=TRUE) #range of dates with data
      if(rdates[1]<inidate) inidate <- rdates[1]
      if(rdates[2]>enddate) enddate <- rdates[2]
    }
    cat('\n')
  } else {
    if(anyf<anyi) stop('Set initial year (anyi) lower or equal than final year (anyf)')
    inidate <- as.Date(sprintf('%d-01-01',anyi))
    enddate <- as.Date(sprintf('%d-12-31',anyf))
  }
  dates <- seq(inidate,enddate,by='1 day') #vector of dates
  nd <- length(dates) #number of dates (=data per station)
  cat(sprintf('%d days between %s and %s\n',nd,inidate,enddate))
  dat <- matrix(NA,nd,ne)
  #populate data matrix:
  cat('\nGenerating',varcli,'input files for Climatol from daily files...:\n\n')
  for(i in 1:ne) { #for every station
    cat(sprintf('%3d %s\n',i,st[i,stcol[1]]))
    d <- read.table(st[i,stcol[1]],sep=sep,dec=dec,na.strings=na.strings) #data
    ddates <- as.Date(as.Date(sprintf('%d-%02d-%02d',d[,datcol[1]],d[,datcol[2]],d[,datcol[3]])))
    kd <- match(ddates,dates) #match data dates with the dates vector
    #avoid "NAs are not allowed in subscripted assignments" error:
    if(sum(is.na(kd))>0) { d <- d[!is.na(kd),]; kd <- kd[!is.na(kd)] }
    ddat <- d[,datcol[4]]
    dat[kd,i] <- ddat
  }
# dat[dat==mis] <- NA #use R missing data code
  #remove stations without mindat data:
  ndat <- apply(!is.na(dat),2,sum)
  sel <- ndat < mindat
  if(sum(sel)>0) { dat <- dat[,!sel]; st <- st[!sel,] }
  #write data file:
  anyi <- format(inidate,'%Y'); anyf <- format(enddate,'%Y')
  fich <- sprintf('%s_%s-%s.dat',varcli,anyi,anyf)
  write(dat,fich,ncolumns=10)
  cat('\nData saved to file',fich,':\n')
  print(summary(as.vector(dat)))
  #write stations file:
  nc <- ncol(st)
  st <- st[,2:nc]; nc <- nc-1 #remove first column (file names)
  if(stcol[5]==0) {
    cod <- as.character(1:ne) #assign station codes
    if(nc>3) st <- cbind(st[,1:3],cod,st[,4:ncol(st)])
    else st <- cbind(st,cod) 
  }
  if(stcol[6]==0) st <- cbind(st,sprintf('st%04d',1:ne)) #assign names
  fich <- sprintf('%s_%s-%s.est',varcli,anyi,anyf)
  write.table(st,fich,row.names=FALSE,col.names=FALSE)
  cat('\nStation coordinates and names saved to file',fich,':\n')
  names(st) <- c('X (lon)','Y (lat)','Z (elev)','Code','Name')
  print(summary(st))
}

#- rclimdex2climatol.- Convert DAILY data from RClimDex to CLIMATOL.
rclimdex2climatol <- function(stfile, kvar, varcli='', chrcod=c(6,10),
  anyi=NA, anyf=NA, mis=-99.9, mindat=365, names=FALSE) {
  #stfile: file with the station codes and coordinates (in HOMER format)
  # (data files will be named as in stfile)
  #kvar.- RClimDex variable to extract: 1(RR), 2(TX), 3(TN)
  #chrcod: initial and final characters of data file names to use as codes
  #anyi: initial year to study (defaults to the first year available in data)
  #anyf: final year to study (defaults to the last year available in data)
  #mindat: minimum number of data per station
  #names: station names in the 9th column? (FALSE by default)
  if(varcli=='') varcli=c('RR','TX','TN')[kvar] #acronym of the variable
  cat('\nGenerating',varcli,'Climatol input files from RClimDex files...:\n\n')
  st <- read.table(stfile,as.is=TRUE) #stations
  ne <- nrow(st)
  if(is.na(anyi) | is.na(anyf)) { #check the time period of the data:
    inidate <- as.Date('3000-12-31'); enddate <- as.Date('0001-01-01')
    for(i in 1:ne) { #for every station
      d <- read.table(st[i,1])
      dates <- as.Date(sprintf('%d-%02d-%02d',d[,1],d[,2],d[,3]))
      rdates <- range(dates,na.rm=TRUE) #range of dates with data
      nadates <- is.na(dates)
      if(sum(nadates)>0) {
        cat('Abnormal dates found in file',st[i,1],':\n')
        print(d[nadates,])
      }
      if(rdates[1]<inidate) inidate <- rdates[1]
      if(rdates[2]>enddate) enddate <- rdates[2]
    }
  }else {
    if(anyf<anyi) stop('Set initial year (anyi) lower or equal than final year (anyf)')
    inidate <- as.Date(sprintf('%d-01-01',anyi))
    enddate <- as.Date(sprintf('%d-12-31',anyf))
  }
  dates <- seq(inidate,enddate,by='1 day') #vector of dates
  nd <- length(dates) #number of dates (=data per station)
  dat <- matrix(NA,nd,ne)
  #populate data matrix:
  for(i in 1:ne) { #for every station
    cat(st[i,1],'\n')
    d <- read.table(st[i,1]) #data
    ddates <- as.Date(sprintf('%d-%02d-%02d',d[,1],d[,2],d[,3]))
    kd <- match(ddates,dates) #match data dates with the dates vector
    #avoid "NAs are not allowed in subscripted assignments" error:
    if(sum(is.na(kd))>0) { d <- d[!is.na(kd),]; kd <- kd[!is.na(kd)] }
    ddat <- d[,kvar+3]
    dat[kd,i] <- ddat
  }
  dat[dat==mis] <- NA #use R missing data code
  #remove stations without mindat data:
  ndat <- apply(!is.na(dat),2,sum)
  sel <- ndat < mindat
  if(sum(sel)>0) { dat <- dat[,!sel]; st <- st[!sel,] }
  #write data file:
  anyi <- format(inidate,'%Y'); anyf <- format(enddate,'%Y')
  fich <- sprintf('%s_%s-%s.dat',varcli,anyi,anyf)
  write(dat,fich,ncolumns=10)
  cat('\nData from',format(inidate),'to',format(enddate),'saved to file',fich,'\n')
  #find longest period without concurrent missing data in all stations:
  avd=apply(!is.na(dat),1,sum)>0
  if(sum(!avd)>0) {
    rle=rle(avd)
    maxrle=which.max(rle$lengths)
    ki=diffinv(rle$lengths)[maxrle]+1
    kf=diffinv(rle$lengths)[maxrle+1]
    cat('The longest period without concurrent missing data in all stations\n')
    cat('  goes from',format(dates[ki]),'to',format(dates[kf]),'\n')
  }
  #write stations file:
  neg <- st[,5]<0; st[neg,5] <- -st[neg,5]
  X <- round(st[,5]+st[,6]/60.+st[,7]/3600.,6)
  X[neg] <- -X[neg]
  neg <- st[,2]<0; st[neg,2] <- -st[neg,2]
  Y <- round(st[,2]+st[,3]/60.+st[,4]/3600.,6)
  Y[neg] <- -Y[neg]
  cod <- substr(st[,1],chrcod[1],chrcod[2])
  if(names) df <- data.frame(X,Y,st[,8],cod,st[,9])
  else df <- data.frame(X,Y,st[,8],cod,cod)
  fich <- sprintf('%s_%s-%s.est',varcli,anyi,anyf)
  write.table(df,fich,row.names=FALSE,col.names=FALSE)
  cat('Station coordinates and names saved to file',fich,'\n\n')
}

#- climatol2rclimdex.- Convert DAILY data from Climatol to RClimDex.
#Read homogenized data for the three daily variables, choose the reconstructions
#from the last homogeneous sub-period, and write them in an RClimDex file.
# varRR, varTX, varTN.- Name of the variables in the climatol files. If some
# variable is not available, name it as ''.
# yiRR, yfRR.- Initial and final years for the RR variable.
# yiTX, yfTX, yiTN, yfTN.- Initial and final years for the TX and TN variables.
# (By default they are the same for the three variables; otherwise, the output
# file will contain data for the common period only).
climatol2rclimdex <- function(varRR,varTX,varTN,yiRR,yfRR,yiTX=yiRR,
    yfTX=yfRR,yiTN=yiRR,yfTN=yfRR,prefix='hoclm',dir=NA,na='-99.9',
    nm=NA, dah=NA, nei=NA, est.c=NA) {
  anyi <- max(c(yiRR,yiTX,yiTN)) #initial year of the output
  anyf <- min(c(yfRR,yfTX,yfTN)) #final year of the output
  if(!is.na(dir)) if(!dir.exists(dir)) dir.create(dir) #output directory
  fech <- seq(as.Date(sprintf('%s-01-01',anyi)),as.Date(sprintf('%s-12-31',anyf)),by='1 day')
  ndd <- length(fech) #nr. of daily data per station
  avl <- rep(FALSE,3) #availability flags
  cod <- NULL
  #-------- read results for the three daily variables (if available):
  #precipitation:
  if(varRR != '') {
    load(sprintf('%s_%d-%d.rda',varRR,yiRR,yfRR))
    if(nm>0) stop(sprintf('Data in %s_%d-%d.rda does not seem to be DAILY!',
      varRR,yiRR,yfRR))
    #select series from last homogeneous fragments:
    fe <- seq(as.Date(sprintf('%s-01-01',yiRR)),as.Date(sprintf('%s-12-31',yfRR)),by='1 day')
    self <- match(fech,fe) #selected days
    dRR <- dah[self,1:nei] #selected data
    sRR <- unsufix(est.c[1:nei,4]) #selected stations
    if(length(sRR)>0) { avl[1] <- TRUE; cod <- sRR }
  }
  #maximum temperatures:
  if(varTX != '') {
    load(sprintf('%s_%d-%d.rda',varTX,yiTX,yfTX))
    if(nm>0) stop(sprintf('Data in %s_%d-%d.rda does not seem to be DAILY!',
      varTX,yiTX,yfTX))
    #select series from last homogeneous fragments:
    fe <- seq(as.Date(sprintf('%s-01-01',yiTX)),as.Date(sprintf('%s-12-31',yfTX)),by='1 day')
    self <- match(fech,fe) #selected days
    dTX <- dah[self,1:nei] #selected data
    sTX <- unsufix(est.c[1:nei,4]) #selected stations
    if(length(sTX)>0) { 
      avl[2] <- TRUE
      if(is.null(cod)) cod <- sTX else cod <- intersect(cod,sTX)
    }
  }
  #minimum temperatures:
  if(varTN != '') {
    load(sprintf('%s_%d-%d.rda',varTN,yiTN,yfTN))
    if(nm>0) stop(sprintf('Data in %s_%d-%d.rda does not seem to be DAILY!',
      varTN,yiTN,yfTN))
    #select series from last homogeneous fragments:
    fe <- seq(as.Date(sprintf('%s-01-01',yiTN)),as.Date(sprintf('%s-12-31',yfTN)),by='1 day')
    self <- match(fech,fe) #selected days
    dTN <- dah[self,1:nei] #selected data
    sTN <- unsufix(est.c[1:nei,4]) #selected stations
    if(length(sTN)>0) { 
      avl[3] <- TRUE
      if(is.null(cod)) cod <- sTN else cod <- intersect(cod,sTN)
    }
  }
  #-------- sort common station codes for the available variables:
  cod <- sort(cod)
  #-------- write RClimDex files (one per station):
  ne <- length(cod) #nr. of stations
  cat('\nGenerating',ne,'RClimDex files from Climatol homogenizations...:\n\n')
  for(i in 1:ne) { #for every station
    cat(' ',cod[i])
    if(is.na(dir)) stfile <- sprintf('%s%s.txt',prefix,cod[i])
    else stfile <- sprintf('%s/%s%s.txt',dir,prefix,cod[i])
    dat <- matrix(NA,ndd,3)
    if(avl[1]) dat[,1] <- dRR[,which(sRR==cod[i])]
    if(avl[2]) dat[,2] <- dTX[,which(sTX==cod[i])]
    if(avl[3]) dat[,3] <- dTN[,which(sTN==cod[i])]
    df <- data.frame(format(fech,'%Y'),format(fech,'%m'),format(fech,'%d'),dat)
    write.table(df,stfile,sep='\t',quote=FALSE,row.names=FALSE,
      col.names=FALSE,na=na)
  }
  cat('\n')
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
  fmeans <- sprintf('%s_%d-%d_msd.csv',varcli,anyip,anyfp)
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
dahstat <- function(varcli, anyi, anyf, anyip=anyi, anyfp=anyf, stat="me",
  ndc=NA, vala=2, cod=NULL, mnpd=0, mxsh=0, prob=.5, last=FALSE, long=FALSE,
  lsnh=FALSE, lerr=FALSE, relref=FALSE, mh=FALSE, pernys=100,
  estcol=c(1,2,4), sep=',', dec='.', eol="\n", nei=NA, x=NA) {
#anyip, anyfp= first and last year for statistics calculation.
#stat="me"(valores medios), "mdn"(medianas), "max"(máximos), "min"(mínimos),
#  "std"(desv.típ.), "q"(cuantiles), "tnd"(tendencias), "series"(series
#  secuenciales, en dos únicos ficheros), "mseries"(series mensuales, dos
#  ficheros por estación)
#  stat='tnd' genera también los p-valores, en *.pval
#ndc=no. de decimales (con prioridad sobre el ndec guardado en el *.rda)
#vala= 0(ninguno), 1(suma), 2(media), 3(máximo), 4(mínimo)
#cod=lista de códigos de las estaciones a listar
#mnpd=Mínimo porcentaje de datos originales
#mxsh=Máximo SNHT
#last=Listar solo las series del último periodo homogéneo
#long=Listar solo las series con mayor fragmento original
#lsnh=Listar solo las series con menor valor SNHT
#lerr=Listar solo las series con menor error cuadrático medio
#relref=Listar también las series de referencia confiables (*Cod)
#mh=Si TRUE, leer datos mensuales de la homogeneización diaria (*-mh_*.dat)
#pernys=No. de años sobre los que se expresan las tendencias (100 por defecto)
#estcol=Columnas de est.c seleccionadas para el listado
#Los parámetros sep, dec y eol permite personalizar los formatos de salida
#nei=No. inicial de estaciones. (Se obtiene del fichero *.rda).
#x=vector de fechas. (Se obtiene del fichero *.rda).
  #- inicializaciones
  if(anyi>anyf) stop ('First year of data greater than the last year!')
  if(anyip<anyi) stop("Asked initial year before first year of data!")
  if(anyfp>anyf) stop("Asked final year beyond last year of data!")
  #función elegida para el cálculo de los valores mensuales:
  fun <- c("mean","median","max","min","sd","quantile")[which(c("me","mdn","max","min","std","q","tnd")==stat)]
  #- si no se reconoce la opción stat, terminar aquí
  if(!length(fun) & stat!='series' & stat!='mseries')
    stop(sprintf("Option stat='%s' not recognized!",stat))
  estvar <- c('X','Y','Z','Code','Name','pod','ios','ope','snht')
  mes3 <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
  #- leer los datos de entrada
  load(sprintf('%s_%d-%d.rda',varcli,anyi,anyf))
  if(mh) { #leer los datos mensuales *-mh_*.dat (generados por dd2m)
    dah <- scan(sprintf('%s-mh_%d-%d.dat',varcli,anyi,anyf))
    nd <- length(dah)/ne
    nm <- 12 #solo datos mensuales si mh=TRUE
    na <- nd/nm #no. de años
    if(na!=round(na,9)) stop(sprintf('Number of data (%d) is not number of months (%d) times number of years (%d)',nd,nm,na))
    dim(dah) <- c(nm,na,ne)
    if(stat=='series' | stat=='mseries') {
      dat <- scan(sprintf('%s-m_%d-%d.dat',varcli,anyi,anyf))
      dim(dat) <- c(nm,na,nei)
      if(nm<1) stat <- 'series' #datos diarios no producen tablas mensuales
    }
  } else {
    if(nm>0) na <- nd/nm else na <- nd #no. de años o de datos (sub)diarios
    if(nm<2) dim(dah) <- c(1,na,ne)
  }
  if(nm<2 | stat=='series' | stat=='mseries') vala <- 0 #valor anual innecesario
  else {
    if(vala<0 | vala>4) vala <- 2 #valor medio en caso de vala erróneo
    funa <- c("sum","mean","max","min")[vala] #función para el valor anual
  }
  if(!is.na(ndc)) ndec <- ndc #cambiar el no. de decimales ndec del *.rda
  #- seleccionar estaciones solicitadas:
  if(!is.null(cod)) {
    ksel <- which(est.c[,4] %in% cod) #estaciones originales solicitadas
    esel <- est.c[,7] %in% ksel #id. madres e hijas
  } else esel <- rep(TRUE,ne)
  #seleccionar las estaciones con un mínimo de mnpd % de datos originales:
  if(mnpd>0) esel <- esel & est.c[,6]>=mnpd #vector de la selección
  #seleccionar las estaciones con un SNHT menor o igual a mxsh:
  if(mxsh>0) esel <- esel & est.c[,9]<=mxsh #vector de la selección
  if(last & ne>nei) esel[(nei+1):ne] <- FALSE #solo últimos fragmentos?
  else if(long | lsnh | lerr) {
    lsel <- rep(TRUE,length(esel)) #inicializar vector
    for(ko in 1:nei) { #para cada estación original
      kest <- which(est.c[,7]==ko) #series de la misma estación ko
      if(length(kest)>1) { #si hay más de un fragmento...
        if(long) ksel <- which.max(est.c[kest,6]) #mayor % de datos originales
        else if(lsnh) ksel <- which.min(est.c[kest,9]) #menor SNHT
        else ksel <- which.min(est.c[kest,10]) #menor error cuadrático medio
        lsel[kest[-ksel]] <- FALSE #selecciÓn deseada
      }
    }
    esel <- esel & lsel
  }
  #eliminar las estaciones de referencia confiables:
  if(!relref) esel <- esel & substr(est.c[,4],1,1)!='*'
  ne <- sum(esel) #no. de estaciones seleccionadas
  if(ne==0) stop("No station selected: No output")
  dah <- dah[,,esel]
  dim(dah) <- c(max(c(1,nm)),na,sum(esel))
  est.c <- est.c[esel,] #lista de estaciones seleccionadas
  iest <- est.c[,7] #índice de las correspondientes estaciones originales
  #- if(stat=="[m]series"), listar las series y sus flags en formato CSV
  if(stat=='series') { #series secuenciales, en dos únicos ficheros:
    dim(dah) <- c(nd,ne); dim(dat) <- c(nd,length(dat)/nd)
    #comparación datos homogeneizados con originales (no usar '=='!):
    df <- abs(dah-dat[,iest]) < 1e-9
    df <- as.numeric(df) #TRUE=1, FALSE=0
    df[df==0] <- 2 #datos distintos a los originales
    df[df==1] <- 0 #datos iguales a los originales
    df[is.na(df)] <- 1 #datos rellenados (originales ausentes)
    dim(df) <- dim(dah)
    #nombre de los ficheros de salida:
    ard <- sprintf('%s_%d-%d_series.csv',varcli,anyi,anyf)
    arf <- sprintf('%s_%d-%d_flags.csv',varcli,anyi,anyf)
    dah <- data.frame(cbind(format(x),dah))
    df  <- data.frame(cbind(format(x),df))
    colnames(dah) <- colnames(df) <- c('Date',est.c[,4])
    write.table(dah,ard,row.names=FALSE,quote=FALSE,sep=',',dec='.')
    write.table(df,arf,row.names=FALSE,quote=FALSE,sep=',',dec='.')
    cat(sprintf('Homogenized values written to %s,\nwith flags in %s:\n',ard,arf))
    cat('  0: Observed data\n')
    cat('  1: Missing data (filled in)\n')
    cat('  2: Corrected data\n')
    return(invisible())
  }
  if(stat=='mseries') { #series mensuales, cada estación dos ficheros:
    xa <- anyi:anyf #vector de años
    for(kest in 1:ne) { #para cada estación seleccionada
      dh <- dah[,,kest] #datos homogeneizados
      do <- dat[,,iest[kest]] #datos originales
      df <- abs(dh-do) < 1e-9
      df <- as.numeric(df) #TRUE=1, FALSE=0
      df[df==0] <- 2 #datos distintos a los originales
      df[df==1] <- 0 #datos iguales a los originales
      df[is.na(df)] <- 1 #datos rellenados (originales ausentes)
      dim(df) <- dim(dh)
      #nombre de los archivos de salida:
      ard <- sprintf('%s_%d-%d_%s_series.csv',varcli,anyi,anyf,est.c[kest,4])
      arf <- sprintf('%s_%d-%d_%s_flags.csv',varcli,anyi,anyf,est.c[kest,4])
      dh <- cbind(xa,t(dh))
      df  <- cbind(xa,t(df))
      if(nm==12) colnames(dh) <- c('Year',mes3)
      else colnames(dh) <- c('Year',as.character(1:nm))
      colnames(df) <- colnames(dh)
      write.table(dh,ard,row.names=FALSE,quote=FALSE,sep=',',dec='.')
      write.table(df,arf,row.names=FALSE,quote=FALSE,sep=',',dec='.')
    }
    cat(sprintf('Homogenized values written to %s_%d-%d_CODE_series.csv,\nwith flags in %s_%d-%d_CODE_flags.csv:\n',varcli,anyi,anyf,varcli,anyi,anyf))
    cat('  0: Observed data\n')
    cat('  1: Missing data (filled in)\n')
    cat('  2: Corrected data\n')
    return(invisible())
  }
  #- calcular la posición de los datos del periodo solicitado:
  if(nm>0) xk <- (anyip-anyi+1):(anyfp-anyi+1)
  else { #posición de las fechas solicitadas:
    datef <- as.Date(sprintf('%s-12-31',anyfp)) #fecha final (máxima)
    xs <- seq(as.Date(sprintf('%s-01-01',anyip)),datef,by='1 day')
    xk <- match(xs,x) #posiciones de xs dentro de x (vector de fechas completo)
  }
  #- if(vala), calcular los valores anuales
  if(vala & nm>0) { #calcular los valores anuales
    aval <- as.vector(apply(dah,2:3,funa))
    dim(dah) <- c(nm,na*ne)
    dah <- rbind(dah,aval)
    nm <- nm+1
    dim(dah) <- c(nm,na,ne)
  }
  #dimensionar valores a calcular:
  val <- matrix(NA,ne,max(1,nm))
  #- if(stat=="tnd"), calcular las tendencias del periodo escogido
  if(stat=="tnd") { #tendencias del periodo escogido
    if(nm>0) pernum <- pernys else pernum <- pernys*365.25
    pval=val #matriz para almacenar los p-valores
    for(i in 1:ne) {
      if(nm<2) { #una sola subserie
        aj <- lm(dah[1,xk,i]~xk) #regresión lineal
        val[i,] <- round(aj$coefficients[2]*pernum,ndec)
        pval[i,] <- round(summary(aj)$coefficients[2,4],3)
      }
      else {
        for(j in 1:nm) {
          aj <- lm(dah[j,xk,i]~xk) #regresión lineal
          val[i,j] <-  round(aj$coefficients[2]*pernum,ndec)
          pval[i,j] <- round(summary(aj)$coefficients[2,4],3)
        }
      }
    }
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
  if(stat=="q") ars <- sprintf('%s_%d-%d_%s%d.csv',varcli,anyip,anyfp,stat,round(100*prob))
  else ars <- sprintf('%s_%d-%d_%s.csv',varcli,anyip,anyfp,stat)
  write.table(dahs[order(est.c[,4]),],ars,row.names=FALSE,sep=',',dec='.')
  cat("\n  written to",ars,"\n")
  if(stat=="tnd") { #grabar los p-valores
    dahs2 <- data.frame(cbind(est.c[estcol],pval))
    names(dahs2) <- ndf
    ars <- sprintf('%s_%d-%d_pval.csv',varcli,anyip,anyfp)
    write.table(dahs2[order(est.c[,4]),],ars,row.names=FALSE,sep=',',dec='.')
    cat("P-values written to",ars,"\n")
  }
}

#- datsubset.- Subset data by selecting a subperiod and/or less missing data.
datsubset <- function(varcli,anyi,anyf,anyis=anyi,anyfs=anyf,minny=NA) {
#anyis, anyfs= first and last year for data subsetting.
#ninny= Minimum number of years with data to subset.
  if(anyis==anyi & anyfs==anyf & is.na(minny)) stop('No subsetting required!\n')
  if(anyis<anyi) stop("Asked initial selected year before first year of data!")
  if(anyfs>anyf) stop("Asked final selected year beyond last year of data!")
  na <- anyf-anyi+1 #nr. of years in original files
  nas <- anyfs-anyis+1 #nr. of years in selected subperiod
  fbas <- sprintf('%s_%d-%d',varcli,anyi,anyf) #raíz nombres de fichero
  fbas2 <- sprintf('%s_%d-%d',varcli,anyis,anyfs) #raíz nombres ficheros salida
  est.c <- read.table(sprintf('%s.est',fbas),colClasses=c("numeric","numeric","numeric","character","character"))
  ne <- nrow(est.c) #no. de estaciones
  dat <- scan(sprintf('%s.dat',fbas))
  numdat <- length(dat) #no. de datos leídos
  nd <- numdat/ne #no. de datos por estación
  dim(dat) <- c(nd,ne) #conversión de vector a matriz
  #calcular no. de datos por año y estación:
  z <- nd/na
  if(z>=1) nm <- ceiling(z)
  if(nm > 12) nm <- 0 #datos diarios
  else if(!nm%in%c(1,2,3,4,6,12)) {
    cat(sprintf('Computed nr. of data per year/station: %d.\n',nm))
    stop('  but it should be one of 1, 2, 3, 4, 6 or 12.\n')
  }
  #generar vector temporal (x):
  if(nm>0) tinc <- sprintf('%d months',12/nm) else tinc <- '1 day'
  x <- seq(as.Date(sprintf('%d-01-01',anyi)),length.out=nd,by=tinc)
  if(fbas==fbas2) { #renombrar ficheros de entrada para no pisarlos:
    file.rename(sprintf('%s.dat',fbas),sprintf('%s-bak.dat',fbas))
    file.rename(sprintf('%s.est',fbas),sprintf('%s-bak.est',fbas))
    cat(sprintf('Original files renamed to %s-bak.dat and %s-bak.est\n',fbas,fbas))
  }
  if(nas < na) { #subset a subperiod of data
    xa <- strftime(x,"%Y") #años de cada dato
    sel <- xa>=anyis & xa<=anyfs
    dat <- dat[sel,]
  }
  if(!is.na(minny)) { #subset data with nr. of years with data >= minny
    nad <- apply(!is.na(dat),2,sum) #nr. of available data per station
    if(nm>0) nyd <- nad/nm else nyd <- floor(nad/365.25)#nr. of years w data
    sel <- nyd >= minny
    if(sum(sel) < ne) {
      dat <- dat[,sel]
      est.c <- est.c[sel,]
    }
  }
  #write output files:
  if(nm>0) ncl <- nm else ncl <- 10
  write(dat,sprintf('%s.dat',fbas2),ncolumns=ncl)
  write.table(est.c,sprintf('%s.est',fbas2),col.names=FALSE,row.names=FALSE)
  cat(sprintf('Subset data written to files %s.dat and %s.est\n',fbas2,fbas2))
}

#- db2dat.- Get data from a database and build input files *.dat and *.est for
#the homogen() function. (ODBC must be intalled and properly configured.)
# ----------------------------------------------------------------------
#Example for a database called "climate", with user "USER" and password "PASS":
# R  #start R (version 3 or higher)
# library(RODBC)
# ch <- odbcConnect("climate",uid="USER",pwd="PASS") #connect to database
# db2dat('HRel',1961,2015,10,FALSE,ch,'%Y-%m-%d','monthly_relhum','Station',
# 'Date','Value','stations','Station','Name','Longitude','Latitude','Elevation')
# odbcClose(ch) #close connection to mcheng
# ----------------------------------------------------------------------
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
  #namax: Máximo no. permitido de datos diarios ausentes originalmente
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
  cat("\n\nMonthly",fun,"values saved to file",fichsal,"\n")
  if(namax>0 & !homog) cat('  (Months with more than',namax,'missing original daily data\n  have also been set to missing)\n')
}

#- homogen.- homogeneización automática de un conjunto de series de datos.
homogen <- function(varcli, anyi, anyf, suf=NA, nm=NA, nref=c(10,10,4), std=3,
swa=NA, ndec=1, dz.max=5, dz.min=-dz.max, wd=c(0,0,100), snht1=25, snht2=snht1,
tol=.02, maxdif=NA, mxdif=maxdif, maxite=999, force=FALSE, wz=.001, trf=0,
mndat=NA, gp=3, ini=NA, na.strings="NA", vmin=NA, vmax=NA, nclust=100,
cutlev=NA, grdcol=grey(.4), mapcol=grey(.4), hires=TRUE, expl=FALSE,
metad=FALSE, sufbrk='m', tinc=NA, tz='UTC', cex=1.2, verb=TRUE) {
  #varcli: variable climática (acrónimo usado)
  #anyi: año inicial
  #anyf: año final
  #suf: sufijo opcional a añadir al nombre de la variable para leer los datos.
  #nm: número de meses. (Si no se fija, se calcula por el no. de datos)
  #nref: no. (máximo) de estaciones de referencia
  #dz.max: límite superior de tolerancia de anomalías 
  #dz.min: límite inferior de tolerancia de anomalías
  #wd: Weight distance (km; distancia a la que el peso se reduce a la mitad;
  #    wd=0: todas las estaciones de referencia pesan lo mismo).
  #snht1: Umbral del SNHT en la primera pasada. (0 to skip stage 1)
  #snht2: Umbral del SNHT en la segunda pasada. (defaults to snht1. 0 to skip)
  #tol: factor de tolerancia (por referencia disponible) para cortes en cadena
  #swa: Semi-Window Amplitude (nr. of data; 5*nm by default, 365 if daily data).
  #trf: Transformar los datos? (0:no transformar; 1:log(x+1); >1:raíz trf)
  #mndat: Mínimo no. de datos para fragentar las series
  #gp: Parámetro de gráficos. 0=ninguno; 1=anomalías globales e histogramas;
  #    2=id+gráficos mensuales de diagnóstico; 3=id+gráficos de medias anuales
  #    móviles y correcciones; 4=id., pero con sumas anuales móviles
  #ini: fecha inicial (para datos diarios, en formato 'AAAA-MM-DD').
  #na.strings: strings marking missing data (NA by default).
  #maxdif: maximum data difference from previous iteration (ndec/2 by default).
  #mxdif: old maxdif parameter (maintained for compatilibility).
  #maxite: maximum number of iterations to compute means (999 by default).
  #vmin, vmax: rango de valores permitidos en la variable climática.
  #nclust: no. máximo de estaciones a usar en el análisis de agrupamiento
  #cutlev: level to cut dendrogram to define clusters (automatic by default).
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
  #tinc: time increment between data. Not set by default, but can be defined for subdaily data, as in e.g.: tinc='3 hour'
  #tz: Time zone. Only relevant for subdaily data. ('UTC' by default.)
  #cex: Character expansion factor para etiquetas y textos en los gráficos.
  #verb: Ver mensajes del proceso por pantalla (además de en el fichero *.txt).
  # -----------------------------------------------------------------  
  #- inicializaciones
  #funciones auxiliares:
  datmed.mean <- function(x) mean(datmed[x])
  datmed.sd <- function(x) sd(datmed[x])
  #en caso de error, cerrar archivos de salida:
  options(error=cerrar)
  #establecer maxdif en función de la precisión elegida:
  if(is.na(maxdif)) maxdif=10^(-ndec)/2 #0.05 para un decimal
  verde <- hsv(.33,1,.6) #color muy usado
  #skip detection stages if metad==TRUE or in exploratory mode:
  if(expl | metad) snht1 <- snht2 <- 0 #skip detection stages if metad==TRUE
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
  if(expl) nref[3] <- nref[2] <- nref[1] #keep nr. of references
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
  names(est.c) <- c('X','Y','Z','Code','Name')
  ne <- nrow(est.c) #no. de estaciones
  #comprobar si las coordenadas están en grados:
  if(max(abs(est.c[,1]))>180 | max(abs(est.c[,2]))>90) {
    deg <- FALSE
    if(mean(est.c[,1])>10000 | mean(est.c[,2])>10000) { #from m to km:
      est.c[,1] <- est.c[,1] / 1000.
      est.c[,2] <- est.c[,2] / 1000.
    }
  } else deg <- TRUE
# else {
#   deg <- TRUE
#   if(gp>0) {
#..................................
#Error in get(dbname) : object 'worldHiresMapEnv' not found
#     if(requireNamespace("maps",quietly=TRUE)) mapok<-TRUE else mapok<-FALSE
#     if(hires & requireNamespace("mapdata",quietly=TRUE))
#       maphr <- TRUE else maphr <- FALSE
#Workaround: put these packages in "Depends"
#..................................
#   }
# }
  fichd <- sprintf('%s.dat',fntr) #nombre del fichero de datos
  dat <- scan(fichd,na.strings=na.strings) #lectura de los datos
  numdat <- length(dat) #no. de datos leídos
  nd <- numdat/ne #no. de datos por estación
  if(nd-floor(nd)>1e-16) {
    cat(ne,"stations read from",fiche,"\n")
    cat(numdat,"data read from",fichd,"\n")
    stop("The number of data is not multiple of the number of stations!")
  } else cat(sprintf('Data matrix: %d data x %d stations\n',nd,ne))
  dim(dat) <- c(nd,ne) #conversión de vector a matriz
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
  #comprobar si los años están completos:
  if(nm>0 & nd%%nm==0) acomp <- TRUE else acomp <- FALSE
  if(is.na(swa)) { if(nm==0) swa=365 else swa=5*nm } #valores de swa por defecto
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
  else if(!is.na(tinc)) x <- seq(as.POSIXct(ini,tz=tz),length.out=nd,by=tinc)
  else x <- seq(as.Date(ini),length.out=nd,by='1 day')
  #- si hay algún término sin datos, emitir aviso y terminar
  numdat <- apply(!is.na(dat),1,sum) #no. de datos de cada término
  if(!min(numdat)) {
    fich=sprintf('%s_%d-%d-availability.pdf',varcli,anyi,anyf)
    pdf(fich,bg='white')
    plot(x,numdat,type='l',xlab='Time',ylab='Number of data',main='Number of data along time')
    grid(col=grey(.4))
    image(x,1:ne,dat,xlab='Time',ylab='Stations',main=paste(varcli,'data availability'),col=4)
    grid(col=grey(.4))
    graphics.off()
    zz <- which(numdat==0); z <- range(zz)
    cat('\n',sum(numdat==0),' time steps between terms ',z[1],' (',format(x[z[1]]),') and ',z[2],' (',format(x[z[2]]),')\n  have missing data in all stations!\n',sep='')
    if(length(zz)<=100) print(format(x[which(numdat==0)]))
    cat(sprintf('(See the figures in %s)',fich),'\n')
    stop("Cannot continue.\n(Shorten the study period, add series with data in those void terms.)\n\n")
  }
  nei <- ne  #no. inicial de estaciones
  est.i <- est.c #datos iniciales de las estaciones
  nsp <- rep(0,nei)  #no. de cortes de cada estación original
  iest <- 1:ne   #índice señalando la serie original de cada subserie
  outan <- matrix(NA,nd,ne) #anomalías de los outliers
  #- if(gp>0), generar gráficos iniciales
  if(gp>0) {
    #activar salida gráfica a documento pdf, con rótulo inicial:
    pdfname <- sprintf('%s_%d-%d.pdf',varcli,anyi,anyf)
    pdf(pdfname,title=pdfname,bg='white')
    old.par <- par(no.readonly=TRUE)
    plot(-1:1,-1:1,type="n",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
    text(0,0.4,sprintf("CLIMATOL %s",climatol.version),cex=4)
    text(0,-0.45,paste("Homogenization\ngraphic output of\n",varcli,"\n",anyi,"-",anyf,sep=""),cex=3)
    par(cex=cex)
    #datos disponibles en cada serie. (Solo si ne<=5*nclust, pues cuando el no.
    #de estaciones es muy grande, el gráfico ocupa mucho espacio):
    if(ne<=5*nclust) {
      image(x,1:ne,dat,xlab='Time',ylab='Stations',main=paste(varcli,'data availability'),col=4)
      grid(col=grdcol)
    } else cat(sprintf('\nMore than 5*%d stations: Per station data availability graph skipped\n\n',nclust))
    #no. de datos de cada término
    numdat <- apply(!is.na(dat),1,sum)
    plot(x,numdat,type="l",ylim=c(0,ne),col="blue",xlab='Time',ylab="Number of data",main=paste("Number of",varcli,"data in all stations"))
    grid(col=grdcol)
    abline(h=5,lty=2,col="green")
    abline(h=3,lty=2,col="red")
    #si hay algún término sin datos no podremos continuar:
    if(!min(numdat)) {
      cat("At least one term has missing data in all stations! (See the PDF graph)\n")
    }
    #boxplots de los datos de cada estación:
    if(nm>1 & acomp) { #boxplots para cada nm (mensuales, etc)
      dim(dat) <- c(nm,na,ne) #dimensiones provisionales
      for(me in 1:nm) { #para cada mes
        z <- data.frame(dat[me,,])
        names(z) <- 1:ne
        #etiqueta del mes (si nm!=12, poner solo el número):
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
    else if(ne<=nclust) { #un solo gráfico con los boxplots de cada estación
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
    if(maxdif>0.01) maxdif <- 0.01 #rebajar maxdif
  }
  if(gp>0) { #continuamos con los gráficos iniciales
    #histograma de todos los datos (distribución quasi-normal?)
    if(trf) main="Histogram of all (transformed) data"
    else main="Histogram of all data"
    hist(dat,xlab=varcli,main=main,col=hsv(.4,1,.8))
    #correlograma de series diferenciadas a nivel mensual (r <-> distancia)
    #(si hay más de nclust estaciones, solo de una muestra aleatoria de nclust)
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
        }
        dz <- (est.c[splc[i],3]-est.c[splc[j],3])*wz
        d2 <- dx*dx+dy*dy+dz*dz #distancia cuadrática
        est.d[i,j] <- sqrt(d2) #distancia
        est.d[j,i] <- est.d[i,j]  #matriz simétrica
      }
    }
    data <- dat[,splc] #copia de los datos
    if(nm>1 & acomp) { #calcular las series diferenciadas por meses
      dim(data) <- c(nm,na,nec) #dimensionamos por meses
      difd <- apply(data,c(1,3),diff)
      dim(difd) <- c(nd-nm,nec) #redimensionar
    }
    else difd <- diff(data) #series diferenciadas globalmente
    corm <- cor(difd,use="p") #matriz de correlaciones
    #eliminar |r|==1 (debidos a estaciones con solo 2 datos en común):
    corm[corm==1] <- NA; corm[corm==-1] <- NA
    if(ne>nclust) main <- sprintf('Correlogram of %d sampled series (first differences)',nclust)
    else if(nm>0) main <- "Correlogram of first difference series"
    else main <- "Correlogram of the daily series"
    if(trf) main <- paste(main," (transformed)",sep="")
    xd <- as.vector(est.d); y <- as.vector(corm)
    plot(xd,y,xlim=c(0,max(est.d,na.rm=TRUE)),col="blue",main=main,
      ylim=c(min(0,min(y,na.rm=TRUE)),1),xlab="Distance (km)",
      ylab="Correlation coefficient")
    grid(col=gray(.3)); abline(h=0,col=2)
    if(ne>2) {  #dendrograma de las estaciones
      dism <- dist(corm) #matriz de disimilaridad
      #si hay NA's en la matriz de disimilaridad, no intentar clustering
      if(!sum(is.na(dism))) {
        hc <- hclust(dism)
        if(ne>nclust) main <- paste("Dendrogram of",nclust,"sampled stations")
        else main <- "Dendrogram of station clusters"
        plot(hc,xlab="Stations",sub="",ylab="Dissimilarity",main=main)
        #clasificación de las estaciones hasta en un máximo de 9 grupos:
        if(is.na(cutlev)) cutlev <- mean(hc$height)+sd(hc$height)
        repeat {
          ct <- cutree(hc,h=cutlev)
          nc <- length(levels(factor(ct)))
          if(nc<10) break
          cutlev <- cutlev + .1
        }
        if(nc>1) {
          abline(h=cutlev,col="red",lty=2)
          if(ne<=nclust) { #list station clusters
            cat('\n-------------------------------------------\n')
            cat(sprintf('Stations in the %d clusters:\n\n',nc))
            print(split(est.c[,3:5],ct))
            cat('---------------------------------------------\n')
          }
        }
      } else { #un solo grupo si dism tiene NA's
        cat("\nNA's in similarity matrix: NO CLUSTER ANALYSIS\n\n")
        nc <- 1; ct <- 1
      }
      #mapa de las estaciones:
      if(nc==1) { col="blue"; main=paste(varcli,"station locations") }
      else {
        col=rainbow(nc,1,.55)[ct]
        main=paste(varcli," station locations (",nc," clusters)",sep="")
      }
      #compute map limits with a 5% margin over station coordinates:
      zxr <- range(est.c[,1]); zyr <- range(est.c[,2])
      zxm <- diff(zxr)*.05; zym <- diff(zyr)*.05
      xlim=c(min(est.c[,1])-zxm, max(est.c[,1])+zxm)
      ylim=c(min(est.c[,2])-zym, max(est.c[,2])+zym)
      #now draw the map
      if(deg) {
#..........................................
#Error in get(dbname) : object 'worldHiresMapEnv' not found
#Error in get(dbname) : object 'worldMapEnv' not found
  #     if(maphr) z<-try(maps::map('worldHires',col=mapcol,xlim=xlim,ylim=ylim))
  #     else if(mapok) z<-try(maps::map('world',col=mapcol,xlim=xlim,ylim=ylim))
#..........................................
#Workaround: Put these packages in 'Depends' to load them explicitely:
        if(hires) z<-try(map('worldHires',col=mapcol,xlim=xlim,ylim=ylim),TRUE)
        else z <- try(map('world',col=mapcol,xlim=xlim,ylim=ylim),TRUE)
#..........................................
        if (class(z)=="try-error") plot(0,0,xlim=xlim,ylim=ylim,xlab='',ylab='',asp=1/(cos(mean(zyr)*pi/180)),main=main,type='n')
#       else { maps::map.axes(); mtext(main,3,line=1,cex=1.7) }
        else { map.axes(); mtext(main,3,line=1,cex=1.7) }
      } else plot(0,0,xlim=xlim,ylim=ylim,xlab="X (km)",ylab="Y (km)",asp=1,main=main,type='n')
      if(ne>99) { #dibujar símbolos si hay más de 99 estaciones
        #dibujar primero en negro las estaciones que no están en la muestra:
        points(est.c[-splc,1:2],pch='+',cex=.5)
        #y ahora las estaciones de la muestra, con sus colores:
        points(est.c[splc,1:2],col=col,pch=ct)
      } else { #hasta nclust estaciones, poner el número
        text(est.c[,1:2],labels=1:ne,col=col)
      }
      grid(col=gray(.4))
    }
    rm(data,difd,corm) #borrar objetos temporales
  }
  #- if(gp==1), terminar (solo se deseaban los gráficos iniciales)
  if(gp==1) {
    par(old.par)
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
  write('"Code","Date","Observed","Suggested","Anomaly (std.devs.)"',Fout)
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
    else if(nchar(sufbrk>3)) fichbrk <- sprintf('%s_%d-%d_brk.csv',sufbrk,anyi,anyf)
    else fichbrk <- sprintf('%s-%s_%d-%d_brk.csv',varcli,sufbrk,anyi,anyf)
    brk <- read.csv(fichbrk,colClasses=c("character","character","character"))
    if(!is.na(tinc)) brk[,2] <- as.POSIXct(brk[,2],tz=tz)
      else brk[,2] <- as.Date(brk[,2])
    nbrk <- nrow(brk); nn <- 0
    if(nbrk<1) cat('No break-points in the metadata file.\n')
    else {
      for(kb in 1:nbrk) { #para cada break:
        i <- match(brk[kb,1],est.c[,4]) #estación a cortar
        if(is.na(i)) {
          cat(sprintf('\nCode %s not found in station list; break skipped',brk[kb,1]))
          next
        }
        kp <- match(brk[kb,2],x) #posición de corte
        if(is.na(tinc)) cat(sprintf('\n%s(%d) breaks at %s',est.c[i,4],i,format(x[kp])))
        else cat(sprintf('\n%s(%d) breaks at %s',est.c[i,4],i,format(x[kp],tz=tz,usetz=TRUE)))
        if(sum(!is.na(dat[1:(kp-1),i])) < mndat) {
          dat[1:(kp-1),i] <- NA
          cat(" Fragment with less than",mndat,"data DELETED\n")
        } else if(sum(!is.na(dat[kp:nd,i])) < mndat) {
          dat[kp:nd,i] <- NA
          cat(" Fragment with less than",mndat,"data DELETED\n")
        } else {
          nn <- nn+1 #incrementamos el no. de nuevas series
          iest <- c(iest,iest[i]) #añadir índice a la serie original
          nsp[iest[i]] <- nsp[iest[i]]+1 #y también su no. de saltos
          if(nm>0) { #contar no. de saltos por año
            z <- 1 + floor((kp-1)/nm) #término anual del salto
            nsy[z] <- nsy[z] + 1 #no. de saltos por año
          }
          dat <- cbind(dat,rep(NA,nd)) #nueva columna de datos
          #pasar los datos anteriores al corte a la nueva serie:
          dat[1:(kp-1),ne+nn] <- dat[1:(kp-1),i]
          dat[1:(kp-1),i] <- NA #borrar los datos pasados a la nueva serie
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
      refhom <- substr(est.c[,4],1,1)=='*' #actualizar referencias homogéneas
    }
    snht1=0 #pasar directamente a relleno de lagunas
  }
  #- for (ks in 1:3) #(snht en ventanas, snht total, y relleno de lagunas)
  for (ks in 1:3) { #para cada etapa:
    #- si snht1==0 en pasada 1 o  snht2==0 en pasada 2, saltar la pasada
    if((snht1==0 & ks==1) | (snht2==0 & ks==2)) next
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
      par(old.par)
      plot(-1:1,-1:1,type="n",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
      text(0,0.4,paste("Stage",ks),cex=4)
      if(ks==1) text(0,-0.3,paste("Binary splits on",swa,"term\nstepped windows\nwith SNHT >",snht1,"\nand wd =",wd[ks],"km"),cex=3)
      else if(ks==2) text(0,-0.3,paste("Binary splits on\nwhole series\nwith SNHT >",snht1,"\nand wd =",wd[ks],"km"),cex=3)
      else text(0,-0.3,paste("Anomalies after\nmissing data\nrecalculation\nwith wd =",wd[ks],"km\n( swa =",swa,")"),cex=2.5)
      par(cex=cex)
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
      oneref <- matrix(FALSE,nd,ne) #solo 1 referencia?
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
        cat("Stations :",formatC(which(numdat<mndat),0,4),"\n")
        cat("N of data:",formatC(numdat[numdat<mndat],0,4),"\n")
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
      #usar las medias y desviaciones típicas anteriores si existen:
      if(exists('dat.m0')) {
        dat.m <- dat.m0
        if(std==3) dat.s <- dat.s0
      }
      #- estandarizar los datos dat.d (obtener dat.z)
      switch(std,
        dat.z <- scale(dat.d,center=dat.m,scale=FALSE),     #std=1
        { dat.z <- dat.d; z <- which(dat.m > 1)
          dat.z[,z] <- scale(dat.d[,z],center=FALSE,scale=dat.m[z]) },
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
      cat("  Station(rank) Date: Observed -> Suggested (Anomaly, in std. devs.)\n")
      if(ks==3) cat("  Iteration Max.data.difference (Station_code)\n")
      repeat {
        ite <- ite+1
        #- ite+=1 y obtener las series estimadas (dat.e|c) con las vecinas
        #  actualizando used, nrefs y mindist:
        for(i in 1:ne) { #para cada estación
          if(!iest[i]) next #estación borrada
          if(refhom[i]) next #saltar estaciones confiables
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
              #si solo hay una referencia, marcar para no corregir la serie:
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
          { dat.c <- dat.e; z <- which(dat.m > 1)
            dat.c[,z] <- scale(dat.e[,z],center=FALSE,scale=1/dat.m[z]) },
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
        sanom <- scale(anom,center=anomm,scale=anoms)
        if(!expl & dz.maxk>.1) { #eliminar outliers
          elim <- sanom < dz.mink | sanom > dz.maxk #datos a eliminar
          elim[is.na(elim)] <- FALSE #eliminar los molestos NA
          elim[,refhom] <- FALSE #no modificar las series confiables
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
                #no eliminar si solo tenían una referencia!:
                if(oneref[j,i] & !force & nrefk>1) {
                  cat(" Only 1 reference! (Unchanged)")
                  elim[j,i] <- FALSE
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
          maxdif3 <- 3*maxddif #control de convergencia
          kmaxdif <- which.max(abs(dat.d-dat.d0)) #posición máx. dat. dif.
          kmaxest <- ceiling(kmaxdif/nd) #estación máx. dat. dif.
          maxsdif <- (dat.d-dat.d0)[kmaxdif %% nd,kmaxest]
        }
        dat.d0 <- dat.d #copia de los datos
        #- actualizar dat.m|s|z
        dat.m <- apply(dat.d,2,mean,na.rm=TRUE)
        if(std==3) dat.s <- apply(dat.d,2,sd,na.rm=TRUE)
        switch(std,
          dat.z <- scale(dat.d,center=dat.m,scale=FALSE), #std=1
          { dat.z <- dat.d; z <- which(dat.m > 1) #std=2
            dat.z[,z] <- scale(dat.d[,z],center=FALSE,scale=dat.m[z]) },
          dat.z <- scale(dat.d,center=dat.m,scale=dat.s), #std=3
          dat.z <- dat.d
        )
        #- if(!aref) break (no afinar los datos ausentes hasta el final)
        #  porque no parece necesario y alarga el tiempo de proceso 
        if(!aref) break
        #- if(ite>1), si los datos ya no varían, break
        if(ite>1) {
          if(maxddif>maxdif3) stop("Averages do not converge iteratively! (There must be something weird in the data). Cannot continue!")
          cat(ite,' ',round(maxsdif,ndec+2)," (",est.c[kmaxest,4],")\n",sep="")
          if(maxddif<=maxdif | ite==maxite) {
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
              #si solo hay una referencia, marcar para no corregir la serie:
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
          { dat.c <- dat.e; z <- which(dat.m > 1) #std=2
            dat.c[,z] <- scale(dat.e[,z],center=FALSE,scale=1/dat.m[z]) },
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
          { dat.z <- dat.d; z <- which(dat.m > 1) #std=2
            dat.z[,z] <- scale(dat.d[,z],center=FALSE,scale=dat.m[z]) },
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
      cat("\nPerforming shift analysis on the",ne,"series...\n")
      for(i in 1:ne) { #análisis de saltos en la media para cada estación
        if(refhom[i]) next #no analizar las series confiables
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
      #máximo tVx de todas las estaciones
      if(sum(!is.na(tVx))==0) tVxx <- 0 else tVxx <- max(tVx,na.rm=TRUE)
      while(tVxx > snht1) {
        i <- which.max(tVx) #estación con el máximo snht
        #si i usó referencias cortadas con un snht demasiado grande, iniciar
        #una nueva iteración:
        if(max(splt[used[i,]])>tVxx*(1+tol*min(nr,sum(used[i,])))) break
        kp <- kpx[i] #posición del tVx en la estación i
        if(oneref[kp,i] & !force & nrefk>1) { #no cortar con una sola referencia
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
          tit <- sprintf('%s   %d (%s)\n%s',varcli,i,est.c[i,4],est.c[i,5])
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
        #contar no. de saltos por año
        z <- as.integer(strftime(x[kp],"%Y"))-anyi+1 #término anual del salto
        nsy[z] <- nsy[z] + 1 #no. de saltos por año
        #dividir la serie por el punto de corte:
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
          dat <- cbind(dat,rep(NA,nd)) #nueva columna de datos
          #pasar los datos anteriores al corte a la nueva serie:
          dat[1:(kp-1),ne+nn] <- dat[1:(kp-1),i]
          dat[1:(kp-1),i] <- NA #borrar los datos pasados a la nueva serie
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
        refhom <- c(refhom,rep(FALSE,nn)) #actualizar vector de refhom
      }
      #- sin nuevos cortes? histogramas de snht y breaks
      if(!nn & !modif) {
        if(gp>1) {
          #histograma de máximos tV, globales (sin 0's, que no son reales):
          z <- tVx[!is.na(tVx) & tVx>0]
          main <- paste("Histogram of maximum SNHT (Stage ",ks,")",sep='')
          if(sum(!is.na(z))) hist(z,breaks=20,xlab='SNHT',col="purple",main=main)
          if(ks==2 | snht2<1) {
            #histograma de no. de cortes por estación:
            hist(nsp,breaks=0:max(9,max(nsp)+1)-.5,col="orange2",xlab="Number of splits",ylab="Number of stations",main="Number of splits per station")
            #frecuencias de fragmentación por años:
            w <- min(5,ceiling(400/na)) #anchura de las barras
            plot(anyi:anyf,nsy,type="h",lwd=w,col=2,ylim=c(0,max(10,max(nsy))),xlab="Years",ylab="Number of splits",main="Number of splits per year")
            grid(col=grdcol)
          }
        }
        #lista de posibles cortes que solo tienen una referencia:
        z <- which(tVx<0)
        if(length(z)>0) {
          cat('Series that could break but had only one reference:\n')
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
        tit <- sprintf('%s   %d (%s)\n%s',varcli,i,est.c[i,4],est.c[i,5])
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
    par(old.par)
    plot(-1:1,-1:1,type="n",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
    text(0,0.4,"Final graphics",cex=3.5)
    text(0,-0.3,"Adjusted series and\napplied corrections",cex=2.5)
    if(nm>0) xlab <- "Years" else xlab <- "Dates"
#   par(cex=cex)
    layout(matrix(1:2,2,1,byrow=TRUE))
    par(las=1,cex=.8*cex)
    #si hay pocos datos diarios (<700?), no filtrarlos:
    if(nm==0 & nd<700) { fltr=1; ylabd <- "Data" }
    else { #filtros para valores anuales:
      if(nm>0) fltr <- rep(1,nm) else fltr <- rep(1,365)
      if(gp>3) ylabd <- "Running annual totals"
      else {
        ylabd <- "Running annual means"
        if(nm>0) fltr <- fltr/nm else fltr <- fltr/365
      }
    }
    for(i in 1:nei) { #para cada estación original
      wi <- which(iest==i) #estaciones derivadas de la estación i
      lwi <- length(wi)
      if(!lwi) next #(estación totalmente borrada!)
      if(lwi>1) vi <- TRUE else vi <- FALSE
      #filtros para valores anuales 
      tit <- sprintf('%s   %d (%s)\n%s',varcli,i,est.c[i,4],est.c[i,5])
      yo <- as.vector(dat.o[,i]) #datos originales
      y <- dah[,wi] #datos homogeneizados
      par(mar=c(0,4,4,2),xaxt="n")
      matplot(x,filter(y,fltr),type="l",lty=1,col=2:20,ylab=ylabd,main=tit)
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
    par(old.par); par(cex=cex)
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
  rmse <- round(rmse,sedec) #redondear los RMSE
  pod <- floor(100*(nd-apply(dat.na,2,sum))/nd) #porcentaje de datos originales
  cat("\nPOD: Percentage of original data\n")
  print(summary(pod))
  #- imprimir resumen de resultados
  cat("\n")
  print(data.frame(ACmx=round(sac,2),SNHT=snhx,RMSE=rmse,POD=pod,Code=est.c[,4],Name=est.c[,5]),right=FALSE)
  #averiguar qué estaciones derivadas funcionan al final del periodo:
  cur <- apply(!is.na(dat[(nd-mndat+1):nd,]),2,sum) #últimos mndat términos
  cur[cur>0] <- 1
  #añadir cinco nuevas columnas a la tabla de estaciones (porcentaje de datos
  #originales, estación original, si funciona actualmente, SNHT y RMSE):
  est.c <- cbind(est.c,pod,iest,cur,snhx,rmse)
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
    if(sum(!is.na(z))) hist(z,breaks=20,xlab="SNHT",col=verde,main=main)
    #histograma de SNHT:
    z <- snhx; main <- "Histogram of maximum global SNHT"
    if(sum(!is.na(z))) hist(z,breaks=20,xlab="SNHT",col="purple",main=main)
    #gráfico de calidad/singularidad:
    plot(rmse,snhx,type="n",xlim=c(0,max(1,max(rmse,na.rm=TRUE))),ylim=c(0,max(50,max(snhx,na.rm=TRUE))),xlab="RMSE",ylab="SNHT",main="Station's quality/singularity")
    grid(col=grdcol)
    text(rmse,snhx,col=hsv(.7,1,.9))
  }
  if(gp>0) { par(old.par); graphics.off() } #cerrar la salida gráfica
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
  if(nm>0 & acomp) {
    dim(dat) <- c(nm,na,nei)
    dim(dah) <- c(nm,na,ne)
  }
  names(est.c) <- c('X','Y','Z','Code','Name','pod','ios','ope','snht','rmse')
  rownames(est.c) <- 1:ne
  save(dat,dah,est.c,nd,ne,nei,nm,x,ndec,std,ini, file=sprintf('%s.rda',fbas))
  #ordenar archivos de outliers y breaks:
  if(!metad) { 
    close(Fbrk)
    brk <- read.csv(sprintf('%s_brk.csv',fbas),colClasses=c("character","character","numeric"))
    brk <- brk[order(brk[,1],brk[,2]),]
    write.csv(brk,sprintf('%s_brk.csv',fbas),row.names=FALSE)
  }
  close(Fout)
  out <- read.csv(sprintf('%s_out.csv',fbas),colClasses=c("character","character","numeric","numeric","numeric"),check.names=FALSE)
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
  sink() #cerrar bitácora
}

#- homogsplit.- Apply homogen() on overlapping split areas.
homogsplit <- function(varcli, anyi, anyf, xc=NULL, yc=NULL, xo=.5, yo=.38,
  maponly=FALSE, suf=NA, nm=NA, nref=c(10,10,4), swa=NA, std=3, ndec=1,
  dz.max=5, dz.min=-dz.max, wd=c(0,0,100), snht1=25, snht2=snht1, tol=.02,
  maxdif=NA, mxdif=maxdif, force=FALSE, wz=.001, trf=0, mndat=NA, gp=3, ini=NA,
  na.strings="NA", maxite=999, vmin=NA, vmax=NA, nclust=100,
  grdcol=grey(.4), mapcol=grey(.4), hires=TRUE, expl=FALSE, metad=FALSE,
  sufbrk='m', tinc=NA, tz='UTC', cex=1.2, verb=TRUE, x=NA) {
  
  #adjust sufbrk if using metad=TRUE within this function:
  if(metad) sufbrk <- sprintf('%s-%s',varcli,sufbrk)
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
  if(is.na(nm)) { #calcular no. de datos por año y estación:
    z <- nd/na
    if(z>=1) nm <- ceiling(z)
    if(nm > 12) nm <- 0 #datos diarios
    else if(!nm%in%c(1,2,3,4,6,12)) {
      cat(sprintf('Computed nr. of data per year/station: %d.\n',nm))
      stop('Please set manually the right value of nm (one of 1,2,3,4,6,12)')
    }
  }
  #comprobar si los años están completos:
  if(nm>0 & nd%%nm==0) acomp <- TRUE else acomp <- FALSE
  #check whether coordinates are in degrees:
  if(max(abs(est.c[,1]))>180 | max(abs(est.c[,2]))>90) deg <- FALSE
  else deg <- TRUE
  if(is.null(xc) | is.null(yc)) { #plot sites to help choosing split borders:
    plot(est.c[,1:2],xlab='X',ylab='Y'); grid(col=grdcol)
    cat('mean(x)=',mean(est.c[,1]),';   mean(y)=',mean(est.c[,2]),'\n')
    if(deg) cat('Aspect ratio=',1/cos(mean(est.c[,2])*pi/180),'\n')
    cat('Choose the cut x and y values and call homogsplit() again specifying\nthem in vectors xc and yc.\n')
  }
  else { #apply homogen() on the selected areas:
    est0 <- est.c; dat0 <- dat #keep whole original data for selections
    pne <- 0 #previous number of stations (series)
    nxc <- length(xc) #nr. of x cut borders
    nyc <- length(yc) #nr. of y cut borders
    #save a map of available stations and split areas:
    f.map <- sprintf('%s-map.pdf',f.bas)
    main=paste('Split areas of the',nei,'available',varcli,'stations')
    pdf(f.map,bg='white')
    if(deg) {
      asp=1/(cos(mean(range(est.c[,2]))*pi/180)) #aspect ratio
      plot(est.c[,1:2],pch='+',col=hsv(.6,.7,1),asp=asp,xlab="Longitude (deg)",ylab="Latitude (deg)",main=main)
#     if(maphr) try(maps::map('worldHires',add=TRUE))
#     else if(mapok) try(maps::map('world',add=TRUE))
    } else plot(est.c[,1:2],pch='+',asp=1,xlab="X",ylab="Y",main=main)
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
    # -------------------- homogenize overlapping areas: -------------
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
        if(sum(est.ov)<5) { #too few stations in current selection?
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
        snht2=snht2, tol=tol, maxdif=maxdif, maxite=maxite, force=force, wz=wz,
        trf=trf, mndat=mndat, gp=gp, ini=ini, na.strings=na.strings, vmin=vmin,
        vmax=vmax, nclust=nclust, grdcol=grdcol, mapcol=mapcol, hires=hires,
        expl=expl, metad=metad, sufbrk=sufbrk, tinc=tinc, tz=tz, cex=cex,
        verb=verb)
      }
    }
    cat("\n======== End of homogenization of overlapping areas, after ")
    cat(format(round(Sys.time()-time1,2)),'\n')
    #joint areal homogenization files:
    f.out <- sprintf('%s_out.csv',f.bas) # out file name
    f.brk <- sprintf('%s_brk.csv',f.bas) # brk file name
    for(i in 1:noa) {
      f.baux <- sprintf('%s-%d_%d-%d',varcli,i,anyi,anyf)
      load(sprintf('%s.rda',f.baux))
      dim(dah) <- c(nd,ne)
      #keep only data from inner area (not from overlapping margins):
      qcod <- est0[ioa==i,4] #original station codes
      qest <- match(qcod,est.c[,4]) #indexes of original stations
      sel <- est.c[,7]%in%qest #selected stations and their derivatives
      dah <- dah[,sel]; est.c <- est.c[sel,]; ne <- sum(sel)
      nei <- length(qcod)
      if(i==1) {
        zne <- ne; znei <- nei
        write(dah,sprintf('%s_aux.dah',f.bas))
        write.table(est.c,sprintf('%s_aux.est',f.bas),row.names=FALSE,col.names=FALSE)
        zz <- read.csv(sprintf('%s_out.csv',f.baux),colClasses=c("character","character","numeric","numeric","numeric"))
        write.csv(subset(zz,zz[,1]%in%qcod),f.out,row.names=FALSE)
        if(!metad) {
          zz <- read.csv(sprintf('%s_brk.csv',f.baux),colClasses=c("character","character","numeric"))
          write.csv(subset(zz,zz[,1]%in%qcod),f.brk,row.names=FALSE)
        }
      } else {
        zne <- zne+ne; znei <- znei+nei
        write(dah,sprintf('%s_aux.dah',f.bas),append=TRUE)
        write.table(est.c,sprintf('%s_aux.est',f.bas),row.names=FALSE,col.names=FALSE,append=TRUE)
        zz <- read.csv(sprintf('%s_out.csv',f.baux),colClasses=c("character","character","numeric","numeric","numeric"))
        write.table(subset(zz,zz[,1]%in%qcod),f.out,sep=',',row.names=FALSE,col.names=FALSE,append=TRUE)
        if(!metad) {
          zz <- read.csv(sprintf('%s_brk.csv',f.baux),colClasses=c("character","character","numeric"))
          write.table(subset(zz,zz[,1]%in%qcod),f.brk,sep=',',row.names=FALSE,col.names=FALSE,append=TRUE)
        }
      }
    }
    dat <- dat0; ne <- zne; nei <- znei
    fich <- sprintf('%s_aux.dah',f.bas); dah <- scan(fich)
    dim(dah) <- c(nd,ne); file.remove(fich)
    fich <- sprintf('%s_aux.est',f.bas)
    est.c <- read.table(fich,colClasses=c('numeric','numeric','numeric','character','character','numeric','numeric','numeric','numeric'))
    file.remove(fich)
    #sort stations and their corresponding data:
    cod <- est.c[,4] #station codes
    zd <- regexpr('-',cod)>0 #daughter series (with a suffix)
    scod <- c(sort(cod[!zd]),sort(cod[zd]))#sorted codes (daughters at the end)
    ks <- match(scod,cod) #order of sorted codes in the original code list
    est.c <- est.c[ks,] #sorted station data frame
    dah <- dah[,ks] #data sorted as stations
    #sort original data also:
    z <- read.table(sprintf('%s.est',f.bas),colClasses=c("numeric","numeric","numeric","character","character"))
    dat <- dat[,match(z[,4],est.c[1:nei,4])] #original data also sorted
    scodi <- unsufix(scod) #codes without suffix
    est.c[,7] <- match(scodi,scod[1:nei]) #update references to mother series
    #set proper dimensions if not daily data, and save data and stations files:
    if(nm>0 & acomp) {
      dim(dat) <- c(nm,nd/nm,nei)
      dim(dah) <- c(nm,nd/nm,ne)
    }
    save(dat,dah,est.c,nd,ne,nei,nm,x,ndec,std,ini, file=f.rda)
    #sort lists of break-points and outliers:
    zout <- read.csv(f.out)
    out <- zout[order(zout[,1],zout[,2],zout[,3]),]
    write.csv(out,f.out,row.names=FALSE)
    if(!metad) {
      zbrk <- read.csv(f.brk)
      brk <- zbrk[order(zbrk[,1],zbrk[,2],zbrk[,3]),]
      write.csv(brk,f.brk,row.names=FALSE)
    }
    #display final messages to the user:
    cat("\n----------- Generated output files: -------------------------\n\n")
    cat(sprintf('%s.txt :  This text output\n',f.bas))
    cat(sprintf('%s :  List of corrected outliers\n',f.out))
    if(!metad) cat(sprintf('%s :  List of corrected breaks\n',f.brk))
    if(gp>0) cat(sprintf('%s-*_%d-%d.pdf :  Diagnostic graphics (one file per area)\n',varcli,anyi,anyf))
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
  fbn <- sprintf('%s_%d-%d',varcli,anyi,anyf) #original file base name
  #destination file base name:
  fbn2 <- sprintf('%s-%s_%d-%d',varcli,suffix,anyi,anyf)
  for(ext in c(".txt",".pdf")) {
    if(restore) file.rename(paste(fbn2,ext,sep=""),paste(fbn,ext,sep=""))
    else file.rename(paste(fbn,ext,sep=""),paste(fbn2,ext,sep=""))
  }
  if(restore) {
    name <- sprintf('%s.rda',fbn2)
    if(file.exists(name)) file.rename(name,sprintf('%s.rda',fbn))
    name <- sprintf('%s_out.csv',fbn2)
    if(file.exists(name)) file.rename(name,sprintf('%s_out.csv',fbn))
    name <- sprintf('%s_brk.csv',fbn2)
    if(file.exists(name)) file.rename(name,sprintf('%s_brk.csv',fbn))
  } else {
    name <- sprintf('%s.rda',fbn)
    if(file.exists(name)) file.rename(name,sprintf('%s.rda',fbn2))
    name <- sprintf('%s_out.csv',fbn)
    if(file.exists(name)) file.rename(name,sprintf('%s_out.csv',fbn2))
    name <- sprintf('%s_brk.csv',fbn)
    if(file.exists(name)) file.rename(name,sprintf('%s_brk.csv',fbn2))
  }
  return(invisible())
}

#- unsufix.- Remove numeric sufixes from the station codes.
unsufix <- function(str) {
  ns <- length(str) #nr. of strings
  for(i in 1:ns) {
    z <- strsplit(str[i],'')[[1]]
    w <- which(z=="-"); nw <- length(w)
    if(nw==0) next #no sufix in the string
    k <- w[length(w)] #last dash in the string
    if(k==nchar(str[i])) next #dash is the last character
    if(as.numeric(substring(str[i],k+1))<=0) next #sufix is no numeric
    str[i] <- substr(str[i],1,k-1)
  }
  return(str)
}

#- snht.- Standard Normal Homogeneity Test (Alexandersson)
snht <- function(x,nmt=3) {
#nmt: no. mínimo de términos de cada muestra
  n <- length(x)
  Tsnht <- rep(NA,n)
# Ssnht <- rep(NA,n)
  if(n<nmt*2) return(Tsnht) #insuficientes datos
  z <- (x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)
  for(i in (nmt+1):(n-nmt)) { #(despreciar los primeros y últimos nmt términos)
    if(is.na(x[i])) next
    n1 <- sum(!is.na(x[1:(i-1)])) #no. de términos de la muestra 1
    n2 <- sum(!is.na(x[i:n])) #no. de términos de la muestra 2
    if(n1<nmt | n2<nmt) next #al menos una muestra es demasiado pequeña
    z1 <- mean(z[1:(i-1)],na.rm=TRUE)
    z2 <- mean(z[i:n],na.rm=TRUE)
    Tsnht[i] <- n1*z1*z1 + n2*z2*z2
#   s1 <- sd(z[1:(i-1)],na.rm=TRUE)
#   s2 <- sd(z[i:n],na.rm=TRUE)
#   Ssnht[i] <- n1*z1*z1 + n2*z2*z2
##  if(is.na(x[i])) yc <- TRUE #marca de test ya calculado
##  else if(yc) yc <- FALSE    #quitar la marca
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

