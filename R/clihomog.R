#clihomog.R.- Homogenization functions for the Climatol package.
#Author: Jose A. Guijarro. Licence: GPL >= 3.0


climatol.version <- '4.0.0'
#- cerrar.- Close output files.
cerrar <- function(graphics=TRUE) {
  if(graphics) graphics.off()
  while (sink.number()>0) sink() #close logfile(s)
  closeAllConnections()
}

#- xls2csv.- Dump data from all *.xls* files into a csv file.
xls2csv <- function(tmpdir, archdir, var, datcols=1:4, codesep='-', dec='.',
sep=',') {
#tmpdir: temporal directory containing the files to read
#archdir: directory where to archive files after processing
#var: destination name of the variable
#datcols: data columns to be written to the output file
#codesep: character string separating the code from the rest of the file name ('-' by default)
# dec='.': character to use as decimal point in the output file
# sep=',': character separating data in the output file
  if(tmpdir==archdir) stop('Please, set archdir different from tmpdir')
  if(!requireNamespace("readxl", quietly=TRUE))
    stop('Please, install the package "readxl" and run this function again')    
  fdat <- sprintf('xls_%s_data.csv',var) #name of the data output file
  fsta <- sprintf('xls_%s_stations.csv',var) #name of the stations output file
  Fs <- file(fdat,'a') #open data output file
  Ft <- file(fsta,'a') #open stations output file
  fich <- dir(tmpdir,'*\\.xls*') #files in tmpdir
  nf <- 0 #counter of processed files
  for(f in fich) { #for every file
    cat(f,'\n')
    if(grepl(codesep,f)) {
      z <- strsplit(f,codesep)
      code <- z[[1]][1] #station code
      name <- z[[1]][2] #(likely) station name
    } else code <- name <- strsplit(f,'\\.')[[1]][1]
    z <- try(d <- as.data.frame(readxl::read_excel(sprintf('%s/%s',tmpdir,f)))[,datcols])
    if(inherits(z,"try-error")) next
    nalines <- apply(is.na(d),1,sum)>0 #lines with one or more missing data
    d <- d[!nalines,]
    d <- sapply (d,as.numeric)
    d <- data.frame(code,d)
    write.table(d,Fs,dec=dec,sep=sep,row.names=FALSE,col.names=FALSE)
    write(sprintf('%s%s%s',code,sep,name),Ft)
    file.rename(sprintf('%s/%s',tmpdir,f),sprintf('%s/%s',archdir,f))
    nf <- nf+1 #update no. of processed files
  }
  close(Fs); close(Ft)
  if(nf>0) {
    cat(sprintf('Data from %d %s/*.xls* files have been saved into %s\n',nf,tmpdir,fdat))
    cat(sprintf('Station codes and names have been written to %s\n',fsta))
    cat('   Coordinates and elevations should be added to this file\n')
    cat('   before running csv2climatol().\n')
    cat(sprintf('(Original files have been moved to the %s directory.)\n\n',archdir))
  } else cat(sprintf('No %s/*.xls* files found!\n\n',tmpdir))
}

#- csv2climatol.- Convert data in a single CSV file to Climatol input format.
#Station codes, names and coordinates can go in a separate CSV file.
csv2climatol <- function(csvfile, datacol=6:8, stnfile=csvfile, stncol=1:5,
varcli, anyi=NA, anyf=NA, mindat=NA, sep=',', dec='.', na.strings='NA',
dateformat='%Y-%m-%d', cf=1, ndec=1, header=TRUE) {
  #csvfile: name of the CSV file containing the data
  #datacol: column(s) holding station codes, dates and data. If 4 (5) values
  #  are provided, dates are expected to appear as year, month (and days) in
  #  separate columns. Otherwise, dates will be provided as character strings
  #  (see parameter dateformat below)
  #stnfile: name of the CSV file containing station codes, names and
  #  coordinates (if these data are not in the csvfile)
  #stncol: columns holding longitudes, latitudes, elevations and station
  #  codes and names. At least coordinates and station codes must be present
  #  in either csvfile or stnfile. Put a zero for any inexistent columns.
  #  Example when stnfile contains only, in this order, latitudes, longitudes
  #  and station names:   stncol=c(2,1,0,3,0)
  #varcli: (short) name of the climatic variable under study
  #anyi: first year to study
  #anyf: last year to study
  #mindat: minimum required number of data per station (by default, 60 monthly
  #  values or 365 daily values)
  #sep: data separator (',' by default: Comma Separated Values)
  #dec: decimal point ('.' by default)
  #na.strings: strings coding missing data ('NA' by default)
  #dateformat: format of dates (if not in separate columns. Default '%Y-%m-%d')
  #cf: conversion factor to apply if data units need to be changed
  #ndec: no. of decimals to round to
  #header: TRUE by default, set to FALSE if csvfile has no header
  #NOTE that if a stnfile is provided, then sep, dec, na.strings and header
  #  defined for csvfile will also be applied to stnfile.
  #----------------- Operation: -----------------------------------
  #read input table:
  d <- read.csv(csvfile,sep=sep,dec=dec,header=header,na.strings=na.strings,as.is=TRUE)
  #find out no. of stations and dates range:
  if(stnfile!=csvfile) stn <- read.csv(stnfile,sep=sep,dec=dec,header=header,
    na.strings=na.strings,as.is=TRUE)[,stncol]
  else stn <- unique(d[,stncol])
  ne <- nrow(stn) #no. of stations
  #if elevations are missing, set them to 99:
  if(stncol[3]==0) stn <- cbind(stn[,1:2],99,stn[,3])
  #if station names are missing, duplicate station codes:
  if(ncol(stn)==4) stn <- cbind(stn,stn[,4])
  stid <- as.character(stn[,4]) #station codes
  dupl <- duplicated(stid)
  if(sum(dupl)>0) { #remove duplicated codes:
    zz <- unique(stn[dupl,4])
    cat('Codes with different names or coordinates:\n')
    print(stn[stn[,4]%in%zz,])
    cat('Only one version of names and coordinates will be kept!\n')
    moda <- function(x) names(which.max(table(x))) #mode function
    for(qz in zz) {
      kz <- which(stn[,4]==qz)
      for(j in 1:5) stn[kz,j] <- moda(stn[kz,j])
    }
    stn <- unique(stn); stid <- stn[,4] #updated station list and codes
  }
  ne <- length(stid) #no. of stations
  ldc <- length(datacol); jd <- datacol[ldc] #data column
  if(ldc==3) {
    if(!inherits(d[,datacol[2]],'character')) d[,datacol[2]] <- 
      as.character(d[,datacol[2]])
    if(grepl('M',dateformat)) fech <- as.POSIXct(d[,datacol[2]],'utc',
      dateformat)
    else if(grepl('H',dateformat)) fech <- as.POSIXct(d[,datacol[2]],
      'utc','%Y%m%d%H')
    else fech <- as.Date(d[,datacol[2]],dateformat)
  } else if(length(datacol)==4) fech <- as.Date(sprintf('%d-%02d-01',
    d[,datacol[2]],d[,datacol[3]]))
  else  fech <- as.Date(sprintf('%d-%02d-%02d',d[,datacol[2]],d[,datacol[3]],
    d[,datacol[4]]))
  xstep <- min(diff(sort(unique(fech)))) #minimum time interval
  if(xstep==1) { nm <- 0; if(is.na(mindat)) mindat <- 365 } #daily values
  else { nm <- 12; if(is.na(mindat)) mindat <- 60 } #monthly values
  nas <- which(is.na(fech))
  z <- as.integer(strftime(range(fech,na.rm=TRUE),'%Y'))
  if(is.na(anyi)) anyi <- z[1] #initial year of data
  if(is.na(anyf)) anyf <- z[2] #final year of data
  #target dates vector:
  if(nm==0) {
    if(units(xstep)=='days') x <- seq(as.Date(sprintf('%d-01-01',anyi)),
      as.Date(sprintf('%d-12-31',anyf)),1)
    else x <- seq(as.POSIXct(sprintf('%d-01-01 00:00:00',anyi),'utc'),
      as.POSIXct(sprintf('%d-12-31 23:50:50',anyf),'utc'),
      by=sprintf('%d %s',xstep,units(xstep)))
  } else x <- seq(as.Date(sprintf('%d-01-01',anyi)),
    as.Date(sprintf('%d-12-01',anyf)),'1 month')
  nd <- length(x) #number of dates (=data per station)
  #initialize data matrix:
  dat <- matrix(NA,nd,ne)
  #populate data matrix:
  cat(sprintf('Creating %s input files for Climatol from %s ...\n',
    varcli,csvfile))
  for(i in 1:ne) { #for every station
    cat(sprintf(' %s',stid[i]))
    sel <- d[,datacol[1]]==stid[i] #select lines of current station
    ds <- d[sel,jd] #data 
    fe <- fech[sel] #dates
    kd <- match(fe,x) #match data dates with the dates vector
    #avoid "NAs are not allowed in subscripted assignments" error:
    z <- is.na(kd); if(sum(z>0)) { ds <- ds[!z]; kd <- kd[!z] }
    dat[kd,i] <- round(ds*cf,ndec)
  }
  cat('\n')
  #remove stations without mindat data:
  ndat <- apply(!is.na(dat),2,sum)
  sel <- ndat < mindat
  if(sum(sel)==ne) stop('Not enough data in any station. No files created!')
  if(sum(sel)>0) { dat <- dat[,!sel]; stn <- stn[!sel,] }
  #write data file:
  fich <- sprintf('%s_%s-%s.dat',varcli,anyi,anyf)
  write(dat,fich,ncolumns=ifelse(nm==0,10,12))
  cat('\nData saved to file',fich,':\n')
  print(summary(as.vector(dat)))
  #write stations file:
  fich <- sprintf('%s_%s-%s.est',varcli,anyi,anyf)
  stn[,1:3] <- sapply(stn[,1:3],as.numeric) #avoid coordinates as characters
  write.table(stn,fich,row.names=FALSE,col.names=FALSE)
  if(length(stncol==5)) {
    cat('\nStation coordinates and names saved to file',fich,':\n')
    names(stn) <- c('X (lon)','Y (lat)','Z (elev)','Code','Name')
    print(summary(stn))
  } else {
    cat('\nStation data saved to file',fich,'\n')
    cat('It should have columns: X (lon), Y (lat), Z (elev), Code, Name\n')
    cat('Please, edit the file to add the missing items in that order.\n\n')
  }
  if(length(nas)>0) {
    cat('Skipped data lines because of wrong dates:\n')
    print(d[nas,]); cat('\n')
  }
}

#- daily2climatol.- Convert daily data files to Climatol input format.
daily2climatol <- function(stfile, stcol=1:6, datcol=1:4, varcli='VRB',
anyi=NA, anyf=NA, mindat=365, sep=',', dec='.', na.strings='NA', header=TRUE) {
#stfile: file with file names and station coordinates, codes and names
#stcol: columns in stfile holding file names, longitudes, latitudes,
#  elevations and station codes and names. (Defaults to 1:6. Use 0 for codes
#  and/or names columns if they are missing, and numeric values will be
#  assigned.)
#datcol: columns in data files holding year,month,day,value (default to 1:4) 
#varcli: short name of the climatic variable under study
#anyi: first year to study (defaults to the first year available in data)
#anyf: last year to study (defaults to the last year available in data)
#mindat: minimum required number of data per station
#sep: Field separator in all files, whether data or stations. (',' by default.)
#dec: decimal point ('.' by default)
#na.strings: strings coding missing data ('NA' by default)
#header: TRUE by default, set to FALSE if files do not have a header
  #read stations file:
  st <- read.table(stfile,as.is=TRUE,sep=sep,dec=dec,header=header)
  ne <- nrow(st) #no. of stations
  if(is.na(anyi) | is.na(anyf)) { #check the time period of the data:
    cat('\nChecking the period covered by the data...\n')
    inidate <- as.Date('3000-12-31'); enddate <- as.Date('0001-01-01')
    for(i in 1:ne) { #for every station
      cat('',i)
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
  cat('\nCreating',varcli,'input files for Climatol from daily files...:\n')
  for(i in 1:ne) { #for every station
    cat(sprintf('%3d %s\n',i,st[i,stcol[1]]))
    d <- read.table(st[i,stcol[1]],sep=sep,dec=dec,header=header,na.strings=na.strings)
    ddates <- as.Date(sprintf('%d-%02d-%02d',d[,datcol[1]],d[,datcol[2]],d[,datcol[3]]))
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
  if(sum(sel)>0) { 
    cat('\nStations with less than',mindat,'data: ',which(sel),'\n')
    if(sum(sel)==ne) stop('No station has enough data!')
    dat <- dat[,!sel]; st <- st[!sel,]; ne <- nrow(st) }
  #write data file:
  anyi <- format(inidate,'%Y'); anyf <- format(enddate,'%Y')
  fich <- sprintf('%s_%s-%s.dat',varcli,anyi,anyf)
  write(dat,fich,ncolumns=10)
  cat('\nData saved to file',fich,':\n')
  print(summary(as.vector(dat)))
  #write stations file:
  nc <- ncol(st)
  #assign numeric codes and names if not provided:
  if(stcol[5]==0)  cod <- as.character(1:ne) else cod <- st[,stcol[5]]
  if(stcol[6]==0)  nam <- as.character(1:ne) else nam <- st[,stcol[6]]
  st <- cbind(st[,stcol[2:4]],cod,nam)
  fich <- sprintf('%s_%s-%s.est',varcli,anyi,anyf)
  write.table(st,fich,row.names=FALSE,col.names=FALSE)
  cat('\nStation coordinates and names saved to file',fich,':\n')
  names(st) <- c('X (lon)','Y (lat)','Z (elev)','Code','Name')
  print(summary(st))
}

#- rclimdex2climatol.- Convert RClimDex daily data files to Climatol format.
rclimdex2climatol <- function(stfile, kvar, chrcod=c(6,10), varcli='',
sep='\t', anyi=NA, anyf=NA, mis=-99.9, mindat=365, header=TRUE) {
#stfile: file with the data file names and station coordinates (HOMER format:
#   'dataFile latDeg latMin latSec lonDeg lonMin lonSec elev stationName')
#sep: column separator (tab by default)
#kvar: RClimDex variable to extract: 1(RR), 2(TX), 3(TN)
#chrcod: initial and final characters of data file names to use as codes
#anyi: initial year to study (defaults to the first year available in data)
#anyf: final year to study (defaults to the last year available in data)
#mindat: minimum number of data per station
  if(varcli=='') varcli=c('RR','TX','TN')[kvar] #(short) name of the variable
  cat('\nCreating',varcli,'Climatol input files from RClimDex files...:\n\n')
  st <- read.table(stfile,sep=sep,as.is=TRUE,header=header) #stations
  ne <- nrow(st)
  if(is.na(anyi) | is.na(anyf)) { #check the time period of the data:
    inidate <- as.Date('3000-12-31'); enddate <- as.Date('0001-01-01')
    for(i in 1:ne) { #for every station
      d <- read.table(st[i,1],header=header)
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
    d <- read.table(st[i,1],header=header) #data
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
  if(ncol(st)>8) df <- data.frame(X,Y,st[,8],cod,st[,9])
  else df <- data.frame(X,Y,st[,8],cod,cod)
  fich <- sprintf('%s_%s-%s.est',varcli,anyi,anyf)
  write.table(df,fich,row.names=FALSE,col.names=FALSE)
  cat('Station coordinates and names saved to file',fich,'\n\n')
}

#- climatol2rclimdex.- Convert DAILY data from Climatol to RClimDex.
climatol2rclimdex <- function(varRR, varTX, varTN, yiRR, yfRR, yiTX=yiRR,
yfTX=yfRR, yiTN=yiRR, yfTN=yfRR, header=TRUE, prefix='hoclm', dir=NA,
na='-99.9') {
#varRR, varTX, varTN: Name of the variables in the climatol files. If some
#  variable is not available, name it as ''.
#yiRR, yfRR: Initial and final years of the homogenized RR series.
#yiTX, yfTX, yiTN, yfTN: Initial and final years of the TX and TN series.
#  The same as yiRR and yfRR by default.
#header: include a header in the files? (TRUE by default)
#prefix: prefix to prepend to station codes to name the output RClimDex files.
#dir: Destination directory of the output RClimDex files.
#na: Missing data code to use in the output RClimDex files.
  nm <- x <- dah <- nei <- est.c <- NULL #(avoid invisible bidings)
  anyi <- max(c(yiRR,yiTX,yiTN)) #initial year of the output
  anyf <- min(c(yfRR,yfTX,yfTN)) #final year of the output
  if(!is.na(dir)) if(!dir.exists(dir)) dir.create(dir) #output directory
  fech <- seq(as.Date(sprintf('%s-01-01',anyi)),as.Date(sprintf('%s-12-31',anyf)),by='1 day')
  ndd <- length(fech) #no. of daily data per station
  avl <- rep(FALSE,3) #availability flags
  cod <- NULL
  #-------- read results for the three daily variables (if available):
  #precipitation:
  if(varRR != '') {
    load(sprintf('%s_%d-%d.rda',varRR,yiRR,yfRR))
    if(nm>0) stop(sprintf('Data in %s_%d-%d.rda does not seem to be DAILY!',
      varRR,yiRR,yfRR))
    self <- match(fech,x) #selected days
    dRR <- dah[self,1:nei] #selected data (from last homogeneous fragments)
    sRR <- est.c[1:nei,4] #selected stations
    if(length(sRR)>0) { avl[1] <- TRUE; cod <- sRR }
  }
  #maximum temperatures:
  if(varTX != '') {
    load(sprintf('%s_%d-%d.rda',varTX,yiTX,yfTX))
    if(nm>0) stop(sprintf('Data in %s_%d-%d.rda does not seem to be DAILY!',
      varTX,yiTX,yfTX))
    self <- match(fech,x) #selected days
    dTX <- dah[self,1:nei] #selected data (from last homogeneous fragments)
    sTX <- est.c[1:nei,4] #selected stations
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
    self <- match(fech,x) #selected days
    dTN <- dah[self,1:nei] #selected data (from last homogeneous fragments)
    sTN <- unsufix(est.c[1:nei,4]) #selected stations
    if(length(sTN)>0) { 
      avl[3] <- TRUE
      if(is.null(cod)) cod <- sTN else cod <- intersect(cod,sTN)
    }
  }
  #-------- sort common station codes for the available variables:
  cod <- sort(cod)
  #-------- write RClimDex files (one per station):
  ne <- length(cod) #no. of stations
  cat(sprintf('\nCreating %d RClimDex files from Climatol homogenizations for %d-%d ...:\n\n',ne,anyi,anyf))
  for(i in 1:ne) { #for every station
    if(is.na(dir)) stfile <- sprintf('%s%s.txt',prefix,cod[i])
    else stfile <- sprintf('%s/%s%s.txt',dir,prefix,cod[i])
    cat(' ',stfile)
    dat <- matrix(NA,ndd,3)
    if(avl[1]) dat[,1] <- dRR[,which(sRR==cod[i])]
    if(avl[2]) dat[,2] <- dTX[,which(sTX==cod[i])]
    if(avl[3]) dat[,3] <- dTN[,which(sTN==cod[i])]
    #exchange TX with TN when TX<TN:
    k <- which(dat[,2]<dat[,3])
    if(length(k)>0) {
      z <- dat[k,2]; dat[k,2] <- dat[k,3]; dat[k,3] <- z
      cat(sprintf('      %d days with TX < TN fixed\n',length(k)))
    } else cat('\n')
    #write the RClimDex file:
    df <- data.frame(format(fech,'%Y'),format(fech,'%m'),format(fech,'%d'),dat)
    names(df) <- c('Year','Month','Day','RR','TX','TN')
    write.table(df,stfile,sep='\t',quote=FALSE,row.names=FALSE,
      col.names=header,na=na)
  }
  cat('\n')
  kest=match(cod,est.c[,4])
  df <- data.frame(est.c[kest,c(4,5,2,1,3)],'XX')
  names(df) <- c('ID','STAT_NAME','LATITUDE','LONGITUDE','ALTITUDE','COUNTRY')
  write.table(df,'hoclm_stations.txt',sep='\t',row.names=FALSE)
  cat('Stations file "hoclm_stations.txt" has been saved.\n')
}

#- sef2climatol.- Convert SEF data files to CLIMATOL input files.
#SEF stands for Station Exchange Format. Visit:
#       https://datarescue.climate.copernicus.eu/node/80
#Missing elevations will be assigned the value 99
sef2climatol <- function(dr, Vbl, varcli=Vbl, ndec=1, na.strings="NA",
mindat=NA) {
  #dr: directory containing the SEF files
  #Vbl: name of the variable in the SEF files
  #varcli: name of the variable in the Climatol destination files
  #ndec: number of decimals to save
  #na.strings: missing data codes (specified as quoted strings)
  #mindat: minimum required number of data per station
  Fs <- file('SEFdata.csv','w') #open auxiliary file
  for(fich in dir(dr)) { #for every file in directory dr
    cat(fich,'\n')
    Fe <- file(sprintf('%s/%s',dr,fich),'r') #open for reading
    li <- readLines(Fe,1)
    if(substr(li,1,3)!='SEF') { cat(':  Not a SEF file'); next }
    li <- readLines(Fe,1)
    cod <- unlist(strsplit(li,'\t'))[2] #station code
    li <- readLines(Fe,1)
    nom <- unlist(strsplit(li,'\t'))[2] #station name
    li <- readLines(Fe,1)
    Y <- as.numeric(unlist(strsplit(li,'\t'))[2]) #Y (longitude)
    li <- readLines(Fe,1)
    X <- as.numeric(unlist(strsplit(li,'\t'))[2]) #X (latitude)
    li <- readLines(Fe,1)
    Z <- unlist(strsplit(li,'\t'))[2]
    if(is.na(Z)|Z==na.strings) Z <- 99 else  Z <- as.numeric(Z) #Z (elevation)
    li <- readLines(Fe,3)
    vrb <- unlist(strsplit(li[3],'\t'))[2] #Vbl
    if(vrb!=Vbl) { cat(':  Not variable',Vbl); next }
    li <- readLines(Fe,3)
    #read data table: (Some files may contain metadata like |DSFLAG="|, which
    #causes not reading the end of line until a pairing quoting is found in
    #the next line, hence skipping half of the data. Parameter quote='\\' 
    #has been set as a workaround.)
    d <- read.table(Fe,sep='\t',na.strings=na.strings,quote='\\',header=TRUE)
    nas <- is.na(d[,3])
    if(sum(nas)>0) d[nas,3] <- '01' #monthly values
    write(sprintf('%f,%f,%f,"%s","%s","%s",%s,%f',X,Y,Z,cod,nom,cod,
      sprintf('%s-%s-%s',d[,1],d[,2],d[,3]),d[,7]),Fs)
    close(Fe)
  }
  close(Fs)
  cat('\n')
  csv2climatol('SEFdata.csv',varcli=varcli,ndec=ndec,mindat=mindat,header=FALSE)
  file.remove('SEFdata.csv')
}

#- dahgrid.- Obtain grids of homogenized data.
dahgrid <- function(varcli, anyi, anyf, anyip=anyi, anyfp=anyf, grid, idp=2.0,
obsonly=TRUE, nmax=Inf) {
#anyip: first reference year for anomalies calculation.
#anyfp: final reference year for anomalies calculation.
#grid: base grid for interpolation, of class SpatialPoints.
#idp: Power of the inverse distance weights (2 by default).
#obsonly: Do not interpolate missing data estimated by homogen().
#nmax: Maximum number of nearest stations to use (all by default).
  nei <- nd <- nm <- std <- est.c <- dat <- NULL #(avoid invisible bidings)
  if(!requireNamespace("sp", quietly=TRUE)
   | !requireNamespace("gstat", quietly=TRUE)
   | !requireNamespace("raster", quietly=TRUE)
   | !requireNamespace("ncdf4", quietly=TRUE)
  ) stop('This function requires packages sp, gstat, raster and ncdf4.\nPlease, install the lacking packages and re-run the function')
  if(anyip<anyi) stop("Asked initial reference year before first year of data!")
  if(anyfp>anyf) stop("Asked final reference year beyond last year of data!")
  #- read original and homogenized data
  fbas <- sprintf('%s_%d-%d',varcli,anyi,anyf) #base file name
  load(sprintf('%s.rda',fbas))
  #- select series from the last fragments
  dah <- dah[,1:nei]
  #- calculate their means and standard deviations in the chosen period
  if(anyip==anyi & anyfp==anyf) { ki <- 1; kf <- nd } #initial and final pos.
  else if(nm==0) { #daily data
    ki <- which(x==as.Date(sprintf('%d-01-01',anyip)))
    kf <- which(x==as.Date(sprintf('%d-12-31',anyfp)))
  } else { ki <- (anyip-anyi)*nm+1; kf <- ki+(anyfp-anyip+1)*nm-1 }
  m <- apply(dah[ki:kf,],2,mean)
  if(std>2) s <- apply(dah[ki:kf,],2,sd)
  #- save the statistics with coordinates to allow their use with a GIS
  if(std<3) {
    df <- data.frame(est.c[1:nei,1:4],m)
    names(df) <- c('X','Y','Z','Code','Mean')
  } else {
    df <- data.frame(est.c[1:nei,1:4],m,s)
    names(df) <- c('X','Y','Z','Code','Mean','Std.Dev.')
  }
  fmeans <- sprintf('%s_msd.csv',fbas)
  write.csv(df,fmeans,row.names=FALSE)
  #- normalize the series
  switch(std,
    daz <- scale(dah,center=m,scale=FALSE), #std=1
    daz <- scale(dah,center=FALSE,scale=m),
    daz <- scale(dah,center=m,scale=s) #std=3 (default)
  )
  #- if(obsonly), blank data missing in the original series
  if(obsonly) daz[is.na(dat)] <- NA
  rg <- range(daz,na.rm=TRUE)
  #- interpolate means (and std. dev.), and save them in NetCDF
  df <- data.frame(est.c[1:nei,1:2],m)
  names(df) <- c('x','y','z')
  sp::coordinates(df) <- ~ x+y
  m <- gstat::idw(z~1,df,grid,nmax=nmax,idp=idp,debug.level=0) #means
  dimLon <- ncdf4::ncdim_def(name='lon', units='degrees_east', vals=unique(grid@coords[,1]))
  dimLat <- ncdf4::ncdim_def(name='lat', units='degrees_north', vals=rev(unique(grid@coords[,2])))
  varCli.m <- ncdf4::ncvar_def(name=sprintf('%s.m',varcli), units='', dim=list(dimLon,
    dimLat), missval=NA, longname=sprintf('%s %d-%d means',varcli,anyip,anyfp))
  listvar <- list(varCli.m)
  nc <- ncdf4::nc_create(sprintf('%s_m.nc',fbas), listvar) #open netCDF file
  zz <- raster::rasterFromXYZ(m)
  ncdf4::ncvar_put(nc,varCli.m,zz@data@values)
  ncdf4::nc_close(nc)
  if(std>2) {
    df <- data.frame(est.c[1:nei,1:2],s)
    names(df) <- c('x','y','z')
    sp::coordinates(df) <- ~x+y
    s <- gstat::idw(z~1,df,grid,nmax=nmax,idp=idp,debug.level=0) #std interp.
    varCli.s <- ncdf4::ncvar_def(name=sprintf('%s.s',varcli), units='',
      dim=list(dimLon, dimLat), missval=NA,
      longname=sprintf('%s %d-%d std. deviations',varcli,anyip,anyfp))
    listvar <- list(varCli.s)
    nc <- ncdf4::nc_create(sprintf('%s_s.nc',fbas), listvar) #open netCDF file
    zz=raster::rasterFromXYZ(s)
    ncdf4::ncvar_put(nc,varCli.s,zz@data@values)
    ncdf4::nc_close(nc)
  }
  #- === create a netCDF with grids interpolated at every time step
  if(is.na(ini)) ini <- sprintf('%d-01-01',anyi) #default initial date
  if(nm>0) x <- seq(as.Date(ini),length.out=nd,by=sprintf('%d months',12/nm))
  else x <- seq(as.Date(ini),length.out=nd,by='1 day')
  dimTime <- ncdf4::ncdim_def(name='Date', units='days since 1970-01-01',
    vals=as.numeric(x), calendar='standard')
  varCli <- ncdf4::ncvar_def(name=varcli, units='', dim=list(dimLon, dimLat,
    dimTime), missval=NA)
  listvar <- list(varCli)
  nc <- ncdf4::nc_create(sprintf('%s.nc',fbas), listvar) #open netCDF file
  #- for every time step:
  cat(sprintf('Interpolating %d grids...:      ',nd))
  kz <- max(10,round(nd/100))
  for(k in 1:nd) {
    if(!k%%kz) cat('\b\b\b\b\b',sprintf('%2s %%',round(k*100/nd)))
    #- interpolate (IDW) the normalized variable estandarizada at grid points
    df <- data.frame(est.c[1:nei,1:2],daz[k,])
    if(obsonly) df <- df[!is.na(df[,3]),]
    names(df) <- c('x','y','z')
    sp::coordinates(df) <- ~x+y
    z <- gstat::idw(z~1,df,grid,nmax=nmax,idp=idp,debug.level=0) #std interp.
    #from SpatialPointsDataFrame to RasterLayer:
    zz=raster::rasterFromXYZ(z)
    #- save values in the netCDF
    ncdf4::ncvar_put(nc,varCli,zz@data@values,start=c(1,1,k),count=c(-1,-1,1))
  }
  cat(' (done)\n\n')
  #- close the netCDF file and finish
  ncdf4::nc_close(nc)
  cat(sprintf('Normalized grids (%f to %f) saved to file %s.nc',rg[1],rg[2],fbas),'\n')
  cat('Means')
  if(std>2) cat(' and standard deviations')
  cat(' (of the whole series) saved to files\n')
  cat(sprintf('%s_m.nc',fbas))
  if(std>2) cat(',',sprintf('%s_s.nc',fbas))
  cat(' and',fmeans,'\n\n')
}

#- dahstat.- Extract series or statistics of the homogenized data.
dahstat <- function(varcli, anyi, anyf, anyip=anyi, anyfp=anyf, stat="me",
ndc=NA, vala=2, valm=vala, cod=NULL, prob=.5, all=FALSE, long=FALSE,
relref=FALSE, pernyr=10, estcol=c(1,2,4), sep=',', dec='.') {
# varcli: (Short) name of the homogenized climatic variable
# anyi: First year of the homogenized series
# anyf: Last year of the homogenized series
# anyip: First year for the statistical calculation
# anyfp: Last year for the statistical calculation
# stat: Statistic to calculate (one of "me"(means), "mdn"(medians), "max"(maxima), "min"(minima), "std"(standard deviations), "q"(quantiles), "tnd"(OLS trends and their p-values), "series"(none, just save homogenized series into a CSV file, plus another file with flags)
# ndc: No. of decimals (defaults to that used in the homogenization)
# vala: Annual value: 0(none), 1(sum), 2(mean), 3(maximum), 4(minimum)
# valm: Monthly value: 1(sum), 2(mean), 3(maximum), 4(minimum)
# cod: vector of requested station codes (all by default)
# prob: Probability to calculate quantiles (0.5 by default)
# all: If TRUE, all reconstructed series will be used. The default is FALSE, hence using only the series reconstructed from the last homogeneuos subperiod
# long: If TRUE (the default is FALSE), only series reconstructed from the longest homogeneuos subperiod will be used
# relref: Set to TRUE to use also any added reliable reference series
# pernyr: No. of years on which to express trend units (10 by default)
# estcol: Columns of est.c to include in the output tables (defaults to c(1,2,4): coordinates and station codes)
# sep: Field separator (',' by default)
# dec: Decimal point ('.' by default)
  #- inicializaciones
  if(valm==0) valm <- 2 #valm is needed for monthly aggregates
  na <- anyf-anyi+1 #no. of years
  if(anyi>anyf) stop ('First year of data greater than the last year!')
  if(anyip<anyi) stop("Asked initial year before first year of data!")
  if(anyfp>anyf) stop("Asked final year beyond last year of data!")
  #chosen function to calculate monthly aggregates:
  fun <- c("mean","median","max","min","sd","quantile")[which(c("me","mdn",
    "max","min","std","q","tnd")==stat)]
  lmcoef <- function(y,x) coef(summary(lm(y~x)))[2,c(1,4)] #regression coef.
  #- unrecognized stat option? finish
  if(!length(fun) & stat!='series') {
    cat('stat should be one of "mean","median","max","min","sd","quantile"\n')
    stop(sprintf("Option stat='%s' not recognized!",stat))
  }
  if(stat=='q') {
    if(length(prob)>1) stop('Please, provide a unique probability to calculate quantiles')
    if(prob<0 | prob>1) stop('prob must be a value between 0 and 1')
  }
  mes3 <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct",
    "Nov","Dec")
  #- read input data
  load(sprintf('%s_%d-%d.rda',varcli,anyi,anyf))
if(!is.na(ndc)) ndec <- ndc
  codo <- unsufix(est.c$Code) #original station codes of all series
  estvar <- names(est.c)
  if(nm==1 | stat=='series') vala <- 0 #annual value not needed
  else {
    if(nm==0) {
      z <- seq(as.Date(sprintf('%s-01-01',anyi)),as.Date(sprintf('%s-12-31',anyf)),1)
      if(length(z)!=nd) stop('Statistics need complete years to be calculated')
    }
    if(vala<1 | vala>4) vala <- 2 #if vala out of range, set it to mean
    funa <- c("sum","mean","max","min")[vala] #function for the annual value
  }
  #- locate data of the requested period:
  if(anyip!=anyi | anyfp!=anyf) {
    yy <- as.integer(strftime(x,'%Y'))
    xk <- min(which(yy==anyip)):max(which(yy==anyfp))
  } else xk <- 1:nd
  #- select requested stations:
  esel <- rep(TRUE,ne)
  if(!is.null(cod)) {
    ksel <- which(codo %in% cod) #requested stations
    esel[-ksel] <- FALSE
  } else cod <- est.c[1:nei,4]
  if(!all & ne>nei) esel[(nei+1):ne] <- FALSE #last fragments only
  else if(long) {
    lsel <- rep(TRUE,length(esel)) #initialize vector
    for(ko in 1:nei) { #for every original station
      kest <- which(codo==est.c$Code[ko]) #series of the same station ko
      if(length(kest)>1) { #if more than one fragment...
        ksel <- which.max(est.c$pod[kest]) #highest % of data
        lsel[kest[-ksel]] <- FALSE #chosen selection
      }
    }
    esel <- esel & lsel
  }
  #delete trusted reference series:
  if(!relref) esel <- esel & substr(est.c$Code,1,1)!='*'
  if(sum(esel)==0) stop("No station selected! (No output)")
  codo <- codo[esel] #original station codes
  est.c <- est.c[esel,] #selected stations
  #- select data from chosen period and estations
  dah <- dah[xk,esel]; dat <- dat[xk,esel[1:nei]] #selected data
  #update parameters:
  x <- x[xk]; ne <- sum(esel); nei <- length(cod); nd <- length(x)
  na <- anyfp-anyip+1 #no. of years
  iest <- match(codo,cod) #index of original stations
  #- if(stat=="series"), list series and flags in CSV format
  if(stat=='series') { #series, in two files (data and flags):
    #compare homogenized and original data (avoid '=='!):
    df <- abs(dah-dat[,iest]) < 1e-9
    df <- as.numeric(df) #TRUE=1, FALSE=0
    df[df==0] <- 2 #data different to the originals
    df[df==1] <- 0 #data equal to the originals
    df[is.na(df)] <- 1 #filled data (originally missing)
    dim(df) <- dim(dah)
    #name of output files:
    ard <- sprintf('%s_%d-%d_series.csv',varcli,anyi,anyf)
    arf <- sprintf('%s_%d-%d_flags.csv',varcli,anyi,anyf)
    dah <- data.frame(cbind(format(x),dah[,order(iest)]))
    df  <- data.frame(cbind(format(x),df[,order(iest)]))
    colnames(dah) <- colnames(df) <- c('Date',est.c[order(iest),4])
    write.table(dah,ard,row.names=FALSE,quote=FALSE,sep=',',dec='.')
    write.table(df,arf,row.names=FALSE,quote=FALSE,sep=',',dec='.')
    cat(sprintf('Homogenized values written to %s,\nwith flags in %s:\n',ard,arf))
    cat('  0: Observed data\n')
    cat('  1: Missing data (filled in)\n')
    cat('  2: Corrected data\n')
    return(invisible())
  }
  #- if(nm==0) calculate monthly aggregates:
  if(nm==0) {
    cat('Computing monthly aggregates... ')
    me <- strftime(x,"%m"); anyo <- strftime(x,"%Y")
    dm <- matrix(NA,na*12,ne) #monthly data
    funm <- c("sum","mean","max","min")[valm] #function for the monthly values
    for(ie in 1:ne) {
      z <- aggregate(dah[,ie],list(me,anyo),funm)
      dm[,ie] <- round(z[,3],ndec)
    }
    cat('Done.\n')
    nm=12; dah <- dm
  }
  dim(dah) <- c(nm,na,ne)
  #- if(vala), calculate annual values
  if(nm>1) { #calculate annual values
    aval <- as.vector(apply(dah,2:3,funa))
    dim(dah) <- c(nm,na*ne)
    dah <- rbind(dah,aval)
    nc <- nm+1 #no. of columns in the table
    dim(dah) <- c(nc,na,ne)
  } else nc <- 1
  #initialize matrix:
  val <- matrix(NA,ne,nc)
  #- if(stat=="tnd"), calculate trends
  if(stat=="tnd") {
    ndec <- ndec+1 #add a decimal place for trends
    pval <- val #matrix to allocate p-values
    if(ne==1) {
      z <- apply(dah,1,lmcoef,x=anyip:anyfp)
      val <- t(as.data.frame(round(z[1,]*pernyr,ndec)))
      pval <- t(as.data.frame(round(z[2,],3)))
    } else {
      z <- apply(dah,c(3,1),lmcoef,x=anyip:anyfp)
      val <- round(z[1,,]*pernyr,ndec)
      pval <- round(z[2,,],3)
    }
  }
  #- else, apply the requested function
  else {
    for(i in 1:ne) {
      if(nc==1) {
        if(stat=="q") val[i,] <- round(eval(call(fun,dah,prob)),ndec)
        else val[i,] <- round(eval(call(fun,dah[,,i])),ndec)
      }
      else { #monthly, bimonthly, quarterly or semester data:
        if(stat=="q") val[i,] <- round(apply(dah[,,i],1,fun,prob),ndec)
        else val[i,] <- round(apply(dah[,,i],1,fun),ndec)
      }
    }
  }
  #- issue message on the output files
  if(stat=="me") cat("Mean")
  else if(stat=="mdn") cat("Median")
  else if(stat=="max") cat("Maximum")
  else if(stat=="min") cat("Minimum")
  else if(stat=="std") cat("Standard deviation")
  else if(stat=="q") cat(prob,"prob. quantile")
  else if(stat=="tnd") cat("Trend")
  cat(" values of ",varcli," (",anyip,"-",anyfp,")",sep="")
  if(stat=="tnd") cat(", expressed in units per ",pernyr," years,",sep="")
  dahs <- data.frame(cbind(est.c[,estcol],val))
  if(nm==12) ndf <- c(estvar[estcol],mes3,"Annual")
  else if(nm==6) ndf <- c(estvar[estcol],sprintf('Bim%d',1:6),"Annual")
  else if(nm==4) ndf <- c(estvar[estcol],sprintf('Qrt%d',1:4),"Annual")
  else if(nm==3) ndf <- c(estvar[estcol],sprintf('4mn%d',1:3),"Annual")
  else if(nm==2) ndf <- c(estvar[estcol],sprintf('Sem%d',1:2),"Annual")
  else if(nm==1) ndf <- c(estvar[estcol],"Annual")
  else stop('Number of data per year is none of 1,2,3,4,6,12')
  names(dahs) <- ndf
  #- save values in output files
  #output file:
  if(stat=="q") ars <- sprintf('%s_%d-%d_%s%d.csv',varcli,anyip,anyfp,stat,round(100*prob))
  else ars <- sprintf('%s_%d-%d_%s.csv',varcli,anyip,anyfp,stat)
  write.table(dahs[order(est.c[,4]),],ars,row.names=FALSE,sep=',',dec='.')
  cat("\n  written to",ars,"\n")
  if(stat=="tnd") { #save p-values
    dahs2 <- data.frame(cbind(est.c[estcol],pval))
    names(dahs2) <- ndf
    ars <- sprintf('%s_%d-%d_pval.csv',varcli,anyip,anyfp)
    write.table(dahs2[order(est.c[,4]),],ars,row.names=FALSE,sep=',',dec='.')
    cat("P-values written to",ars,"\n")
  }
  cat("\n")
}

#- datsubset.- Subset data by subperiod, code list or no. of years with data.
datsubset <- function(varcli, anyi, anyf, anyis=anyi, anyfs=anyf, minny=NA, 
codes=NULL, na.strings=NA,  ini=NA) {
#varcli: (short) name of the climatic variable
#anyi, anyf: first and last year of the input data file
#anyis, anyfs: first and last year for data subsetting
#ninny: minimum number of years with data to subset
#codes: vector of chosen station codes. (Defaults to NULL, meaning all)
#na.strings: strings marking missing data (NA by default)
#ini:initial date (if not January 1st)
  if(anyis==anyi & anyfs==anyf & is.na(minny) & is.null(codes))
    stop('No subsetting required!\n')
  if(anyis<anyi) stop("Asked initial selected year before first year of data!")
  if(anyfs>anyf) stop("Asked final selected year beyond last year of data!")
  #- read input data
  z <- read.dat(varcli,anyi,anyf,ini=ini,na.strings=na.strings)
  est.c <- z$est.c; dat <- z$dat; na <- z$na; nd <- z$nd; ne <- z$ne; x <- z$x
  nm <- z$nm; rm(z) #free memory
  nas <- anyfs-anyis+1 #no. of years in selected subperiod
  fbas <- sprintf('%s_%d-%d',varcli,anyi,anyf) #base file name
  fbas2 <- sprintf('%s_%d-%d',varcli,anyis,anyfs) #base output file name
  if(fbas==fbas2) { #rename input data to avoid overwriting:
    file.rename(sprintf('%s.dat',fbas),sprintf('%s-ori.dat',fbas))
    file.rename(sprintf('%s.est',fbas),sprintf('%s-ori.est',fbas))
    cat(sprintf('Original files renamed to %s-ori.dat and %s-ori.est\n',fbas,fbas))
  }
  #subset a subperiod of data? :
  if(nas < na) {
    xa <- strftime(x,"%Y") #years of every data
    sel <- xa>=anyis & xa<=anyfs
    dat <- dat[sel,]
  }
  #subset series with at least minny years with data? :
  if(!is.na(minny)) {
    #minny must be an integer equal or greater than 1:
    if(minny!=round(minny)) minny <- round(minny); if(minny<1) minny <- 1
    nad <- apply(!is.na(dat),2,sum) #no. of available data per station
    if(nm>0) nyd <- nad/nm else nyd <- floor(nad/365.25)#no. of years w data
    sel <- nyd >= minny; nes <- sum(sel) #no. of selected stations
    if(nes==0) stop(sprintf('No series has >=%d years of data',minny))
    if(nes < ne) { dat <- dat[,sel]; est.c <- est.c[sel,] }
  }
  #subset series of selected stations? :
  if(!is.null(codes)) {
    ke <- match(codes,est.c[,4]); ke <- ke[!is.na(ke)]
    if(length(ke)==0) stop('Selected codes does not meet the other requirements')
    dat <- dat[,ke]; est.c <- est.c[ke,]
  }
  #write output files:
  if(nm>0) ncl <- nm else ncl <- 10
  write(dat,sprintf('%s.dat',fbas2),ncolumns=ncl)
  write.table(est.c,sprintf('%s.est',fbas2),col.names=FALSE,row.names=FALSE)
  cat(sprintf('Subset data written to files %s.dat and %s.est\n',fbas2,fbas2))
}

#- db2dat.- Get data from a database and build input files *.dat and *.est for
#  the homogen() function. (ODBC must be intalled and properly configured.)
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
db2dat <- function(varcli, anyi, anyf, minny=5, daily=TRUE,ch,
dformat='%Y-%m-%d', vtable, vcode, vdate, vval, stable, scode, sname, sx, sy,
sz) {
  #varcli: Short name of the climatic variable under study
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
  na <- anyf-anyi+1 #no. of years
  if(na<=0) stop('Last year must be greater than the first year')
  fini <- sprintf('%d-01-01',anyi)
  if(daily) {
    x <- seq(as.Date(sprintf('%d-01-01',anyi)),as.Date(sprintf('%d-12-31',anyf)),by='1 day')
    ndmin <- round(minny*365.25) #min. no. of daily data
    ffin <- sprintf('%d-12-31',anyf)
  } else {
    x <- seq(as.Date(sprintf('%d-01-01',anyi)),as.Date(sprintf('%d-12-01',anyf)),by='1 month')
    ndmin <- minny*12 #min. no. of monthly data
    ffin <- sprintf('%d-12-01',anyf)
  }
  nd <- length(x) #no. of data per station
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
    dd <- RODBC::sqlQuery(ch,sprintf("SELECT %s,%s FROM %s WHERE %s >= '%s' AND %s <= '%s' AND %s = '%s'",vdate,vval,vtable,vdate,fini,vdate,ffin,vcode,ds[i,4]))
    if(is.null(dim(dd))) next #no data for the variable at this station
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
    ne <- ne + 1 #count no. of saved series
    ndat <- ndat + !is.na(dat) #count no. of data at every time step
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

#- dd2m.- Calculate monthly values from daily or subdaily data.
dd2m <- function(varcli, anyi, anyf, ndec=1, valm=2, namax=30, x=NULL,
na.strings='NA', tz='utc') {
#varcli: Short name of the climatic variable
#anyi: Initial year
#anyf: Final year
#ndec: No. of decimals requested in the results (1 by default)
#valm: Monthly value (1=sum, 2=mean, 3=maximum, 4=minimum, 5=standart deviation)
#namax: Maximum allowed percentage of missing data in a month
#x: Time vector. Automatically set by default, but needed if data are taken
#  at irregular intervals.
#na.strings: Strings marking missing data (NA by default).
#tz: Time zone (if data are subdaily). ('utc' by default.)
  #- read input data
  z <- read.dat(varcli,anyi,anyf,x,na.strings=na.strings,tz=tz)
  est.c <- z$est.c; dat <- z$dat; na <- z$na; ne <- z$ne; x <- z$x
  me <- strftime(x,"%m")
  anyo <- strftime(x,"%Y")
  fun <- c("sum","mean","max","min","sd")[valm] #function for monthly values
  ndm <- na*12 #no. of monthly values
  dm <- matrix(NA,ndm,ne) #monthly values
  zp <- c(table(me,anyo)) #no. of possible data in every month
  for(ie in 1:ne) { #for every station
    cat(' ',ie)
    z <- aggregate(dat[,ie],list(me,anyo),fun,na.rm=TRUE)[,3] #monthly data
    z2 <- aggregate(is.na(dat[,ie]),list(me,anyo),sum)[,3] #no. of missings
    #with subdaily data some can be in year anyf+1:
    if(length(z)>ndm) { z <- z[1:ndm]; z2 <- z2[1:ndm]; zp <- zp[1:ndm] }
    nas <- z2/zp > namax/100. #months with not enough data
    z[nas] <- NA #set their values to missing
    dm[,ie] <- z #assign monthly data to the main matrix
  }
  dm[is.nan(dm)] <- NA #assign NA for missing data
  #save monthly data:
  fichsal <- sprintf("%s-m_%d-%d.dat",varcli,anyi,anyf)
  fichest <- sprintf("%s-m_%d-%d.est",varcli,anyi,anyf)
  write(round(dm,ndec),fichsal,ncolumns=12)
  write.table(est.c,fichest,row.names=FALSE,col.names=FALSE)
  cat("\n\nMonthly",fun,"values saved to file",fichsal,"\n")
  if(namax>0) cat('(Months with more than',namax,
    '% missing data have also been set to missing)\n\n')
}

#- exampleFiles.- Get the path to some example files.
#(Adapted from readxl::readxl_example)
exampleFiles <- function(file=NULL) {
#file: Name of the needed file. If NULL, all example files will be listed.
  if(is.null(file)) dir(system.file("files",package ="climatol"))
  else system.file("files", file, package = "climatol", mustWork=TRUE)
}

#- fix.sunshine.- Check homogenized daily sunshine hours and prune any excess.
fix.sunshine <- function(varcli, anyi, anyf) {
#varcli: Short name of the homogenized climatic variable
#anyi: First year of the homogenized series
#anyf: Last year of the homogenized series
  nm <- ndec <- x <- ne <- est.c <- nd <- NULL #(avoid invisible bidings)
  #- leer los datos de entrada
  frda <- sprintf('%s_%d-%d.rda',varcli,anyi,anyf) #file to load/save
  obj2save <- load(frda)
  if(nm!=0) stop('This function only applies to DAILY sunshine series')
  #- auxiliary functions to calculate maximum theoretical sunshine hours
  insolteor <- function(lat,fech) { #maximum theoretical sunshine hours
    nf <- length(fech) #no. of requested dates
    it <- rep(NA,nf) #maximum theoretical sunshine hours vector
    latr <- lat * 0.01745329 #latitude in radianes (0.01745329 = 2 * pi / 360)
    for(k in 1:nf) {
      dj <- as.numeric(strftime(fech[k],'%j'))
      dec <- declin(dj)
      c <- -tan(latr) * tan(dec)
      if(c <= -1.) it[k] = 24.
      else if(c >= 1.) it[k] = 0.
      else {
        b <- 1.570796 - atan(c / sqrt(-c*c+1.))
        r <- b * 24 / pi
        if(r > 12) d <- r else d <- 24-r
        it[k] <- r + 2 * d * 0.004627778 + .05 # 0.004627778 = .833 / 180
      }
    }
    return(it)
  }
  declin <- function(djul) { #sun declination
#Approximation from http://solardat.uoregon.edu/SolarRadiationBasics.html:
# declin = 23.45 * pi / 180 * sin(2 * pi * (284 + n) / 365) =
    return(0.4092797 * sin(4.888834 + 0.01721421 * djul))
  }
  #- prune exceeding values
  sink('fix.sunshine.txt',split=TRUE)
  cat('Checking sunshine durations of',frda,'and prunning any excess...\n')
  fixed <- FALSE; rmargin <- 1/10^ndec/2 #flag and rounding margin
  dec=declin(as.numeric(strftime(x,'%j')))
  for(j in 1:ne) {
    cat('------',est.c[j,4],est.c[j,5],'\n')
    c <- -tan(est.c[1,2]*0.01745329) * tan(dec)
    r <- (1.570796 - atan(c / sqrt(-c*c+1.)))*24/pi
    d <- r; z <- r>12; d[z] <- 24-d[z]
    it <- round(r + 2 * d * 0.004627778 + .05, ndec) #maximum possible value
    for(i in 1:nd) if(dah[i,j]>it[i]) {
      cat(est.c[j,4],format(x[i]),dah[i,j],'->',it[i],'\n')
      dah[i,j] <- it[i]; fixed <- TRUE
    }
  }
  if(fixed) {
    frda0 <- sprintf('%s.bak',frda)
    file.rename(frda,frda0)
    cat('\nOriginal file',frda,'renamed to',frda0,'\n')
    cat('Writing the new',frda,'file...\n')
    save(list=obj2save, file=frda)
    cat('List of fixed values saved to fix.sunshine.txt\n')
  } else cat('(No value has been modified)\n')
  sink()
}

#- homogen.- automatic homogenization of climate series.
homogen <- function(varcli, anyi, anyf, test='snht', nref=NULL, std=NA,
swa=NA, ndec=1, niqd=3, dz.max=.01, dz.min=-dz.max, cumc=NA, wd=NULL,
inht=25, sts=5, maxdif=NA, maxite=999, force=FALSE, wz=.001, mindat=NA,
onlyQC=FALSE, annual=c('mean','sum','total'), x=NULL, ini=NA, na.strings="NA",
vmin=NA, vmax=NA, hc.method='ward.D2', nclust=300, cutlev=NA, grdcol=grey(.4),
mapcol=grey(.4), expl=FALSE, metad=FALSE, sufbrk='m', tinc=NA, tz='utc',
rlemin=NA, rlemax=NA, cex=1.1, uni=NA, raway=TRUE, graphics=TRUE, verb=TRUE,
logf=TRUE, snht1=NA, snht2=NA, gp=NA) {
#varcli: Short name of the studied climatic variable
#anyi: Initial year
#anyf: Final year
#test: Break detection test to apply. One of 'snht' (the default) or 'cuct'
#nref: Maximum no. of reference stations at each stage
#std: Type of normalization. 1 (remove the mean), 2 (divide by the mean) or
#   3 (remove the mean and divide by the standard deviation).
#swa: Semi-Window Amplitude (no. of data; defaults to 60 months or 365 days).
#ndec: No. of required decimal places in the results (1 by default)
#niqd: No. of interquartilic distances to delete big outliers (3 by default)
#dz.max: If >1, upper tolerance limit for anomalies (if two values are given,
#  only those higher than the upper one will be rejected);
#  If <=1, percentage of anomalous data to reject (in each side of the
#  distribution (0.01 by default).
#dz.min: lower tolerance limit for anomalies (-dz.max by default).
#cumc: code of accumulated missing data.
#wd: Weight distance, in km. Distance at which the weight of a reference data
#  is halved. (If wd=0, all reference stations will have equal weight.)
#inht: Inhomogeneity threshold(s). (0 to skip the stage.)
#sts: Series tail size (defaults to 5).
#maxdif: maximum data difference from previous iteration (ndec/2 by default).
#maxite: maximum number of iterations to compute means (999 by default).
#force: force direct homogenization of (sub)daily series [FALSE].
#wz: Scale factor for elevation Z. The default value (0.001) equals vertical
#  differences in m to the horizontal differences in km. Can be used to
#  give more weight to Z, or to calculate horizontal distances only (wz=0).
#mindat: Minimum no. of data for a split fragment to become a new series
#  [swa/2 for daily series or 12 terms otherwise].
#onlyQC: Set to TRUE if only initial Quality Controls are requested [FALSE]
#annual: Running annual value to graph in the PDF output. One of 'mean' (the
#  default), 'sum' or 'total' (equivalent to 'sum').
#x: Time vector. Only needed if data are taken at irregular intervals.
#ini: Initial date, with format 'AAAA-MM-DD' (for daily data, if series does not begin on January first as recommended).
#na.strings: Strings marking missing data (NA by default).
#vmin, vmax: Range of allowed values for the climatic variable.
#hc.method: hierarchical clustering method ('ward.D2' by default).
#nclust: Maximum number of series for the cluster analysis [300].
#cutlev: Level to cut dendrogram to define clusters (automatic by default).
#grdcol: Color of the graphic background grids [grey(0.04)].
#mapcol: Color of coastlines and borders in the stations map [grey(0.04)].
#expl: Perform an exploratory analysis? [FALSE].
#metad: Use the breakpoints file as metadata? [FALSE].
#sufbrk: Suffix to add to varcli to form the name of the provided metadata file ['m'; set to '' if original data were monthly].
#tinc: Time increment between data. Not set by default. If defined
#  for subdaily data, should be in units 'hours', 'mins' or 'secs'.
#  E.g.: tinc='1 hours'. (Do not forget the last 's' in the units).
#tz: Time zone. Only relevant for subdaily data. ('utc' by default.)
#rlemin: Data run lengths will exclude values <= rlemin in quality control.
#rlemax: Data run lengths will exclude values >= rlemax in quality control.
#cex: Character expansion factor for graphic labels and titles [1.1].
#uni: Units to use in some axis labels. (None by default.)
#raway: Increase internal distances to reanalysis series to give more weight to observed series (TRUE by default).
#graphics: Output graphics in a PDF file [TRUE].
#verb: Verbosity [TRUE].
#logf: Save console messages to a log file?  [TRUE].
#snht1, snht2: Obsolete but kept for backwards compatibility.
#gp: Obsolete but kept for backwards compatibility.
#------------------------------------------------------------------
#backwards compatibility:
  if(!is.na(snht1)) inht <- snht1
  if(!is.na(snht2)) inht <- c(inht,snht2)
  if(!is.na(snht1)|!is.na(snht2)) cat('Please, note that parameters snht1 and snht2 are deprecated.\nUse inht in future applications of Climatol.\n')
  if(!is.na(gp)) {
    graphics <- TRUE
    switch(gp+1,
      graphics <- FALSE, #gp=0
      onlyQC <- TRUE,    #gp=1
      ,                  #gp=2
      ,                  #gp=3
      annual <- 'sum',   #gp=4
    )
    cat('Please, note that parameter gp is deprecated.\nUse graphics=FALSE for gp=0, onlyQC=TRUE for gp=1 or\nannual="total" for gp=4 in future applications of Climatol.\n')
  }
  if(onlyQC & !graphics) graphics <- TRUE #fix possible inconsistency
  #- initializations
  annual <- match.arg(annual); if(annual=='total') annual <- 'sum'
  verde <- hsv(.33,1,.6) #dark green
  #auxiliary functions:
  datmed.mean <- function(x) mean(datmed[x])
  datmed.sd <- function(x) sd(datmed[x])
  r3 <- function(x) sign(x)*abs(x)^(1/3) #cubic root
  #chosen inhomogeneity test:
  if(test=='cuct') { inhtest <- 'CucT'; test <- 'cuct' } #apply Cucconi
  else { inhtest <- 'SNHT'; test <- 'snht' } #apply SNHT
  #in case of error, close all output files:
  options(error=cerrar)
  if(!is.na(std)) {
    std <- as.integer(std) #std must be an integer between 1 and 3:
    if(std<1|std>3) cat('Warning:')
    if(std<1) { std <- 1; cat('std lower than 1 has been forced to 1\n') }
    if(std>3) { std <- 3; cat('std greater than 3 has been forced to 3\n') }
  }
  if(is.null(nref)) {
    nref <- c(10,10,4); if(!is.na(cumc)) nref[3] <- 2
  }
  if(is.null(wd)) {
    wd <- c(0,0,100); if(!is.na(cumc)) wd[3] <- 10
  }
  #if unset, set maxdif depending on ndec:
  if(is.na(maxdif)) maxdif=10^(-ndec)/2 #0.05 for one decimal
  #skip detection stages if metad==TRUE or in exploratory mode:
  if(expl | metad) inht <- c(0,0)
  #dz.min must be negative!:
  z=dz.min>0; if(sum(z)>0) dz.min[z] <- -dz.min[z]
  #- open log file and write header
  archlog <- paste(varcli,"_",anyi,"-",anyf,".txt",sep="")
  if(logf) sink(archlog,split=verb)
  cat("\nHOMOGEN() APPLICATION OUTPUT  (From R's contributed package 'climatol' ",climatol.version,")\n",sep='')
  cat("\n=========== Homogenization of ",varcli,", ",anyi,"-",anyf,". (",
    date(),")\n",sep="")
  time1 <- Sys.time() #time at the beginning of the process
  cat("\nParameters:")
  arg <- names(formals()) #list the function arguments
  nargs <- length(arg) #no. of arguments
  for(i in 1:nargs) {
    if(arg[i]=='x') next #do not print the time vector!
    cat(" ",arg[i],"=",sep="")
    cat(eval(as.symbol(arg[i])))
    if(i<nargs) cat(',')
  }
  cat("\n\n")
  #complete multivalue parameters:
  k <- length(inht); if(k<2) inht <- c(inht,inht) #same threshold in all stages
  k <- length(wd); if(k<3) wd <- c(rep(0,3-k),wd)
  k <- length(nref); if(k<3) nref <- c(nref,rep(nref[k],3-k))
  if(expl) nref[3] <- nref[1] #keep no. of references in exploratory mode
  #data anomaly tolerance (warning and delete):
  if(sum(dz.max) > 1) {
    dz.maxw <- min(dz.max); dz.maxd <- max(dz.max)
    dz.minw <- max(dz.min); dz.mind <- min(dz.min)
  }
  if(length(dz.max)==1) dz.maxw <- NA
  if(length(dz.min)==1) dz.minw <- NA
  fbas <- sprintf('%s_%d-%d',varcli,anyi,anyf) #file basename 
  #- read input data
  z <- read.dat(varcli,anyi,anyf,x,ini,tinc,na.strings,tz=tz)
  est.c <- z$est.c; dat <- z$dat; na <- z$na; nd <- z$nd; ne <- z$ne; x <- z$x
  nm <- z$nm; ini <- z$ini; tinc <- z$tinc; acomp <- z$acomp
  if(nm<1) {
    if(!onlyQC & !metad & !expl & !force & is.na(cumc)) stop('These series seem to be daily or sub-daily. Their direct homogenization\n  is not recommended. Please, use dd2m() and homogenize the monthly\n  series first or set force=TRUE to avoid this message.')
  } else if(!nm%in%c(1,2,3,4,6,12)) {
    cat(sprintf('Calculated no. of data per year and station: nm=%f\n',nm))
    stop('Complete years of monthly or seasonal data are required to avoid\n  missing data in the results. Please complete your series to have\n  a rounded number of data. (nm can be one of 1,2,3,4,6,12).')
  }
  #- disaggregate data if cumc code is set
  if(!is.na(cumc)) { #manage data coded as accumulated to the next:
    graphics <- FALSE #do not create graphics if using cumc
    cuml <- apply(dat==cumc,2,which) #list of accumulated terms
    if(length(cuml)==0) {
      cat('\nThere are no data coded as cumc =',cumc,'. Nothing done.\n')
      cerrar(graphics); return(invisible())
    }
    cumt <- function(z) { #first (cumA) and last (cumB) accumulated terms
      zdif <- diff(z)>1
      cumA <- z[c(TRUE,zdif)]
      cumB <- z[c(zdif,TRUE)]
      return(list(cumA,cumB))
    }
    z <- lapply(cuml,cumt)
    cumA <- sapply(z, function(x) x[1]) #first terms of accumulated runs
    cumB <- sapply(z, function(x) x[2]) #last terms of accumulated runs
    dat[dat==cumc] <- NA #delete accumulated data
    cumv <- vector('list',ne) #list of accumulation values:
    for(j in 1:ne) {
      cumv[[j]] <- dat[cumB[[j]]+1,j] #save accumulation values
      dat[cumB[[j]]+1,j] <- NA #delete them
      cumB[[j]] <- cumB[[j]]+1 #include them in the accumulated runs
    }
  } 
  nsy <- rep(0,na)   #no. of shifts per year
  if(is.na(swa)) { #swa default values:
    if(nm>0) swa <- 5*nm #5 years for monthly or lower frequency
    else { #1 year for daily or higher frequency
      z=which.max(table(unlist(rle(strftime(x,'%Y'))[1])))
      swa <- as.integer(names(z))
    }
  }
  if(swa>nd/4) swa <- ceiling(nd/4) #avoid too big a window
  #set mindat if unset by the user:
  if(is.na(mindat)) { 
    if(nm<=0) mindat <- nd/na/4 #three months
    else mindat <- max(5,nm)
  }
  if(graphics) { #activate graphic output to a pdf file:
    pdfname <- sprintf('%s_%d-%d.pdf',varcli,anyi,anyf)
    pdf(pdfname,title='homogen() graphics output',bg='white')
    old.par <- par(no.readonly=TRUE)
    plot(-1:1,-1:1,type="n",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
    text(0,0.4,sprintf("CLIMATOL %s",climatol.version),cex=4)
    text(0,-0.45,paste("Homogenization\ngraphic output of\n",varcli,"\n",anyi,"-",anyf,sep=""),cex=3)
    par(las=1,cex=cex,cex.axis=cex,cex.lab=cex)
    my.par <- par(no.readonly=TRUE)
  }
  #----------- Quality control of the series -----------------------
  Fout <- file(sprintf('%s_out.csv',fbas),'w') #open outliers file
  write('"Code","Date","Observed","Suggested","Anomaly (std.devs.)","Deleted"',Fout)
  #set vmin and std if unset by the user:
  mn1 <- apply(dat,2,quantile,prob=.1,na.rm=TRUE) #first deciles
  mn0 <- apply(dat,2,quantile,prob=.01,na.rm=TRUE) #first percentiles
  if(median(mn1)==0 & median(mn0)==0) { #series skewed and limited by zero:
    skewed <- TRUE
    if(is.na(vmin)) vmin <- 0 #minimum value
    if(is.na(rlemin)) rlemin <- 0.1 #minimum value for run lengths
    if(is.na(std) & vmin==0) std <- 2 #normal ratio normalization
  } else {
    skewed <- FALSE
    if(is.na(std)) std <- 3 #standardization
  }
  deld <- 0 #initialize deleted data count
  #check if there are values out of allowed range
  if(!is.na(vmin)) { #values lower than the minimum?
    zout <- dat<vmin; nout <- sum(zout,na.rm=TRUE)
    if(nout>0) {
      sout <- apply(zout,2,sum,na.rm=TRUE)
      for(k in 1:length(sout)) {
        if(sout[k]==0) next
        kz <- which(zout[,k])
        for(j in kz) write(c(est.c[k,4],format(x[j]),dat[j,k],NA,NA,1),
          Fout,ncolumns=6,sep=',')
      }
      dat[zout] <- NA #delete erroneous values
      cat('\nWarning: deleted',nout,'data lower than',vmin,'\n')
      deld <- deld + nout
    }
  }
  if(!is.na(vmax)) { #values higher than the maximum?
    zout <- dat>vmax; nout <- sum(zout,na.rm=TRUE)
    if(nout>0) {
      sout <- apply(zout,2,sum,na.rm=TRUE)
      for(k in 1:length(sout)) {
        if(sout[k]==0) next
        kz <- which(zout[,k])
        for(j in kz) write(c(est.c[k,4],format(x[j]),dat[j,k],NA,NA,1),
          Fout,ncolumns=6,sep=',')
      }
      dat[zout] <- NA #delete erroneous values
      cat('\nWarning: deleted',nout,'data greater than',vmax,'\n')
      deld <- deld + nout
    }
  }
  #check if there are series without any data:
  ksd <- which(apply(!is.na(dat),2,sum) == 0)
  if(length(ksd)>0) {
    cat("There are series with no data!!!:",'\n')
    print(est.c[ksd,])
    stop("Please, remove these series from the input files and run homogen() again")
  }
  #check if there are series with too few data:
  ksd <- which(apply(!is.na(dat),2,sum) < mindat)
  if(length(ksd)>0) {
    cat("There are series with too few data (less than",mindat,'):\n')
    print(est.c[ksd,])
    cat("Warning: Break-points cannot be corrected in these series",'\n')
  }
  #- if there are time steps without any data, issue a warning and finish
  numdat <- apply(!is.na(dat),1,sum) #no. of data in each time step
  if(!min(numdat)) {
    plot(x,numdat,type='l',xlab='Time',ylab='Number of data',
      main=sprintf('Number of %s data along time',varcli))
    z <- class(x)
    if(z[1]=='Date' | z[1]=='POSIXct') {
      grid(NA,NULL,col=grey(.4))
      if(z[1]=='Date') abline(v=axis.Date(1,x),lty=3,col=grey(.4))
      else if(z[1]=='POSIXct') abline(v=axis.POSIXct(1,x),lty=3,col=grey(.4))
    } else grid(col=grey(.4))
    zz <- which(numdat==0); z <- range(zz)
    cat('\n',sum(numdat==0),' time steps between terms ',z[1],' (',format(x[z[1]]),') and ',z[2],' (',format(x[z[2]]),')\n  have missing data in all stations!\n',sep='')
    if(length(zz)<=100) print(format(x[which(numdat==0)]))
    cat(sprintf('(See the figures in %s)',pdfname),'\n')
    cerrar(graphics)
    stop("Cannot continue.\n(Shorten the study period or add series with data in the void terms.)\n\n")
  }
  if(graphics) {
    #partition stations in groups for the boxplots:
    gs <- 40 #group size
    ge <- 1+((0:(ne-1))%/%gs); nge <- max(ge) #station groups
    #si the last group is very small, append it to the former group:
    if(nge>1 & sum(ge==nge)<=10) {
      ge[ge==nge] <- nge-1
      nge <- nge-1
    }
  }
  #detect too anomalous data (to be deleted to avoid their use as references):
  if(skewed) {
    da3 <- r3(dat); da3[dat==0] <- NA #cubic root of data without zeros
    bp <- boxplot(da3,range=niqd,plot=FALSE) #big outliers of the skewed series
  } else {
    bp <- boxplot(dat,range=niqd,plot=FALSE) #big outliers of the series
    da3 <- dat
  }
  nout <- length(bp$out)
  if(nout>0) {
    kj <- rep(NA,nout)
    for(k in 1:nout) {
      kj[k] <- which(da3[,bp$group[k]]==bp$out[k])
      bp$out[k] <- dat[kj[k],bp$group[k]]
    }
  }
  if(graphics) { #boxplots of the series, marking those too anomalous:
    for(ig in 1:nge) { #for every station group:
      ke <- which(ge==ig)
      boxplot(dat[,ke],ylab=ifelse(is.na(uni),'Values',uni),xaxt='n',col=5,
        xlab='Stations',main=paste(varcli,'data')); grid(col=grdcol)
      keg <- bp$group %in% ke; kex <- bp$group[keg]-(ig-1)*gs
      points(kex,bp$out[keg],pch=19,col=2,cex=1.3)
      gsx <- length(ke)
      if(max(ke)>100) axis(1,1:gsx,ke,las=2) else axis(1,1:gsx,ke)
    }
  }
  #delete too anomalous data:
  nout <- length(bp$out)
  if(nout>0) {
    cat('\nWarning:',nout,'big outliers deleted')
    if(!onlyQC) cat(' before homogenization:\n') else cat(':\n')
    grp <- unique(bp$group); ngrp <- length(grp)
    for(i in 1:ngrp) {
      j <- grp[i]; k <- which(bp$group==j); zout <- bp$out[k]
      cat(j,est.c[j,4],': ')
      if(length(zout)<=15) cat(zout,'\n') else cat(zout[1:15],'... etc\n')
      kk <- which(dat[,j]%in%zout); nkk <- length(kk)
      for(k in 1:nkk) write(c(est.c[j,4],format(x[kk[k]]),dat[kk[k],j],NA,NA,1),
        Fout,ncolumns=6,sep=',')
    }  
    #delete too anomalous data:
    for(k in 1:nout) dat[kk[k],bp$group[k]] <- NA 
    deld <- deld + nout #deleted data count
  }
  cat('\n')
  if(graphics) { #boxplots of increments between consecutive data:
    daf <- diff(dat); da3 <- r3(daf); da3[daf==0] <- NA
    bp <- boxplot(daf,range=niqd,plot=FALSE) #big outliers of the differences
    b0 <- boxplot(da3,range=niqd,plot=FALSE) #big outliers of their cubic roots
    if(length(b0$out)<length(bp$out)) { #choose the more conservative method
      bp <- b0 #change diff outliers by their cubic roots:
      nout <- length(bp$out); kj <- rep(NA,nout)
      if(nout>0) for(k in 1:nout) {
        kj[k] <- which(da3[,bp$group[k]]==bp$out[k])
        bp$out[k] <- daf[kj[k],bp$group[k]]
      }
      daf[daf==0] <- NA #remove zeros
    }
    ylab <- 'Increments'; if(!is.na(uni)) ylab <- sprintf('%s (%s)',ylab,uni)
    for(ig in 1:nge) { #for every station group:
      ke <- which(ge==ig)
      boxplot(daf[,ke],xlab='Stations',ylab=ylab,xaxt='n',col=hsv(.2,.5,1),
        main=paste(varcli,'data increments')); grid(col=grdcol)
      keg <- bp$group %in% ke; kex <- bp$group[keg]-(ig-1)*gs
      points(kex,bp$out[keg],pch=19,col=2,cex=1.3)
      gsx <- length(ke)
      if(max(ke)>100) axis(1,1:gsx,ke,las=2) else axis(1,1:gsx,ke)
    }
    #plot lengths of sequencies of constant values:
    dat.o <- dat; if(!is.na(rlemin)) dat.o[dat<rlemin] <- NA
    if(!is.na(rlemax)) dat.o[dat>rlemax] <- NA
    rlen <- sapply(apply(dat.o,2,rle), function(x) x[1]) #run lengths
    z <- sapply(rlen,unique); zx <- max(unlist(z)) #maximum run length
    ylab='No. of identical consecutive values'
    ylim=c(1,max(10,zx))
    for(ig in 1:nge) { #for every station group:
      ke <- which(ge==ig)
      if(max(ke)>100) las=2 else las=1
      plot(0,0,xlim=range(ke),ylim=ylim,xlab='Stations',ylab=ylab,type='n',
        las=las,main=paste(varcli,'run lengths')); grid(col=grdcol)
      for(k in ke) {
        if(inherits(z,'matrix')) uz <- z[,k] else uz <- z[[k]]; 
        s <- rep(k,length(uz))
        points(s,uz,pch=19,col=4)
      }
    }
  }
  #------------- End of initial quality control -----------------------
  dat.o <- dat #copy of original data
  if(nd<100) lw=3 #width of the bars of anomalies
  else if(nd<300) lw=2
  else lw=1
  nei <- ne  #initial no. of stations
  est.i <- est.c #initial stations metadata
  nsp <- rep(0,nei)  #no. of cuts in every original series
  iest <- 1:ne   #index of original series of every sub-series
  outan <- matrix(NA,nd,ne) #outlier anomalies
  #- if(graphics), go on with the initial graphics
  if(graphics) {
    #data availability in every series:
    z <- class(x)
    if(nm<0) cat("Per station data availability graphic skipped for subdaily data\n   (might be too heavy).\n")
    else { #per station data availability:
      if(sum(is.na(dat))==0) col=4 else col=c('white',4)
      image(!is.na(dat),col=col,xaxt='n',yaxt='n',useRaster=TRUE,
        xlab='Time',ylab='Series',main=paste(varcli,'data availability'))
      yy=as.integer(strftime(x,'%Y')); if(length(unique(yy))<3) yy <- 1:nd
      lb=pretty(yy); nt=length(lb)
      if(lb[nt] > max(yy)) { nt=nt-1; lb=lb[1:nt] }
      if(lb[1] < min(yy)) { lb=lb[2:nt]; nt=nt-1 }
      kp=integer(); for(k in 1:nt) kp=c(kp,min(which(yy==lb[k])))
      at=((1:nd-1)/(nd-1))[kp]
      axis(1,at,labels=lb,xlab='Time')
      abline(v=at,lty=3,col=grdcol)
      lb=pretty(1:ne); nt=length(lb)
      if(lb[nt] > ne) { nt=nt-1; lb=lb[1:nt] }
      if(lb[1] < 1) { lb=lb[2:nt]; nt=nt-1 }
      at=((1:ne-1)/(ne-1))[lb]
      axis(2,at,labels=lb)
      abline(h=at,lty=3,col=grdcol)
    }
    #no. of data in every time step:
    plot(x,numdat,type="l",col=4,ylab="Number of data",xlab='Time',
      ylim=c(0,ne),main=paste("Number of",varcli,"data in all stations"))
    if(z[1]=='Date' | z[1]=='POSIXct') {
      grid(NA,NULL,col=grdcol)
      if(z[1]=='Date') abline(v=axis.Date(1,x),lty=3,col=grdcol)
      else if(z[1]=='POSIXct') abline(v=axis.POSIXct(1,x),lty=3,col=grdcol)
    } else grid(col=grdcol)
    abline(h=5,lty=2,col=verde)
    abline(h=3,lty=2,col="red")
    #histogram of all data (near-normal distribution?)
    main="Histogram of all data"
    zh <- hist(dat,plot=FALSE)
    zx <- zh$breaks
    zy <- zh$counts; zy[zy==0] <- NA
    barplot(zy,log='y',space=0,ylim=c(.9,max(zy,na.rm=TRUE)*2),xlab=varcli,
      ylab='Frequency',main=main,col=hsv(.4,1,.8),names.arg=zh$mids)
    #correlogram of fisrt differences of the series (r <-> distance)
    #(if more than nclust series, use only a random sample of nclust series)
    if(ne>nclust) { splc <- sample(1:ne,nclust); nec <- nclust }
    else { splc <- 1:ne; nec <- ne }
    est.d <- matrix(NA,nec,nec) #distance matrix
    for(i in 1:(nec-1)) {
      for(j in (i+1):nec) {
        dx <- est.c[splc[i],1]-est.c[splc[j],1]
        dy <- est.c[splc[i],2]-est.c[splc[j],2]
        #distnaces in km (gross method as flat geometry):
        dx <- dx*111*cos((est.c[splc[i],2]+est.c[splc[j],2])*pi/360)
        dy <- dy*111
        dz <- (est.c[splc[i],3]-est.c[splc[j],3])*wz
        d2 <- dx*dx+dy*dy+dz*dz #quadratic distance
        est.d[i,j] <- sqrt(d2) #distance
        est.d[j,i] <- est.d[i,j]  #simmetric matrix
      }
    }
    data <- dat[,splc] #copy of the data
    difd <- diff(data) #first differences of the series
    corm <- cor(difd,use="p") #correlation matrix of the differences
    #change |r|==1 (due to cases with only 2 common data, apart from the
    #diagonal) by the mean correlation (to avoid NAs in the matrix):
    corm[abs(corm)==1] <- mean(corm,na.rm=TRUE)
    if(ne>2) {  #correlogram of the stations:
      if(ne>nclust) main <- sprintf('Correlogram of %d sampled %s series\n(First differences)',nclust,varcli)
      else main <- sprintf('Correlogram of first difference %s series',varcli)
      xd <- as.vector(est.d); y <- as.vector(corm)
      xmin <- floor(min(c(0,xd),na.rm=TRUE)); xmax <- ceiling(max(xd,na.rm=TRUE))
      ymin <- floor(min(c(0,y),na.rm=TRUE)); ymax <- ceiling(max(y,na.rm=TRUE))
      xbin <- seq(xmin,xmax,length=100)
      ybin <- seq(ymin,ymax,length=100)
      freq <- as.data.frame(table(findInterval(xd,xbin),findInterval(y,ybin)))
      freq[,1] <- as.integer(as.character(freq[,1]))
      freq[,2] <- as.integer(as.character(freq[,2]))
      freq2D <- matrix(0,100,100)
      freq2D[cbind(freq[,1], freq[,2])] <- freq[,3]
      freq2D[freq2D==0] <- NA
      col=rev(rainbow(16,start=0,end=2/3))
      nz=max(freq2D,na.rm=TRUE); if(nz<16) col=col[1:nz]
      image(xbin,ybin,freq2D,main=main,useRaster=TRUE,xlab="Distance (km)",
        ylab="Correlation coefficient",col=col)
      grid(col=gray(.3)); abline(h=0,col=2)
      #dendrogram of the stations:
      dism <- dist(corm) #dissimilarity matrix
      #if there are NAs in the dissimilarity matrix, set them to 1:
      kna=which(is.na(dism)); if(sum(kna)>0) dism[kna] <- 1
      hc <- hclust(dism,hc.method)
      if(ne>nclust) main <- paste("Dendrogram of",nclust,"sampled stations")
      else main <- "Dendrogram of station clusters"
      plot(hc,xlab="Stations",sub="",ylab="Dissimilarity",las=0,main=main)
      #station clusters up to a maximum of 9 groups:
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
          print(split(splc,ct))
          cat('---------------------------------------------\n')
        }
      }
      #stations map:
      if(nc==1) { col="blue"; main=paste(varcli,"station locations") }
      else {
        col=rainbow(nc,1,.55)[ct]
        main=sprintf('%s station locations (%d clusters)',varcli,nc)
      }
      #map limits with a 10% margin around station coordinates:
      zxr <- range(est.c[,1]); zyr <- range(est.c[,2])
      zxm <- diff(zxr)*.1; zym <- diff(zyr)*.1
      if(zxm==0 | zym==0) stop('X and/or Y coordinates are all the same!\nAllow some variation in the *.est input file')
      xlim <- zxr+c(-zxm,zxm); ylim <- zyr+c(-zym,zym)
      #now draw the map:
      par(mar=c(4,4,4,4))
      if(diff(zyr)<5. & requireNamespace('mapdata',quietly=TRUE))
        z <- try(maps::map('worldHires',col=mapcol,xlim=xlim,
        ylim=ylim,mar=par('mar')),silent=TRUE)
      else z <- try(maps::map('world',col=mapcol,xlim=xlim,ylim=ylim,
        mar=par('mar')),silent=TRUE)
      if(inherits(z,"try-error")) plot(0,0,xlim=xlim,ylim=ylim,xlab='',
        ylab='',asp=1/(cos(mean(zyr)*pi/180)),main=main,type='n')
      else { maps::map.axes(); title(main) }
      if(ne>99) { #if more than 99 stations, plot symbols
        #stations not in the sample are plotted first, in black:
        points(est.c[-splc,1:2],pch='+',col=2,cex=.5)
        #stations in the sample are plotted in color:
        points(est.c[splc,1:2],col=col,pch=ct)
      } else text(est.c[splc,1:2],labels=splc,col=col) #numbers
      grid(col=gray(.4))
    }
    par(my.par) #restore graphic parameters
    par(las=1,cex=cex,cex.axis=cex,cex.lab=cex)
  }
  #- if(onlyQC), finish (saving files if data were deleted)
  if(onlyQC) {
    if(deld>0) {
      fold <- sprintf('%s-old_%d-%d',varcli,anyi,anyf) #old file basename
      file.rename(sprintf('%s.dat',fbas),sprintf('%s.dat',fold))
      write(dat,sprintf('%s.dat',fbas))
      file.rename(sprintf('%s.est',fbas),sprintf('%s.est',fold))
      write.table(est.c,sprintf('%s.est',fbas),row.names=FALSE,col.names=FALSE)
      cat('New input data files free from the detected big outliers\n')
      cat('have been written, and original files have been renamed to\n')
      cat(sprintf('%s.dat and %s.est\n',fold,fold))
    }
    cat("\nOnly the initial exploratory graphics were demanded.\nSee them in ",varcli,"_",anyi,"-",anyf,".pdf\n\n",sep="")
    cerrar(graphics)
    if(exists('ct')) return(invisible(list(corm=corm,ct=ct)))
    else return()
  }
  #  Homogenization processs in three stages:
  #  1) splits in stepping windows
  #  2) splits in the whole series
  #  3) missing data filling
  #- open outliers and breaks files
  if(!metad) {
    Fbrk <- file(sprintf('%s_brk.csv',fbas),'w')
    write(sprintf('"Code","Date","%s"',inhtest),Fbrk)
  }
  #- compute distance and proximity rank matrices
  cat("Computing inter-station distances ...")
  refhom <- substr(est.c[,4],1,1)=='*' #trusted homogeneous references
  est.d <- matrix(0,ne,ne) #distance matrix
  for(i in 1:(ne-1)) {
    cat(" ",i)
    for(j in (i+1):ne) {
      dx <- est.c[i,1]-est.c[j,1]
      dy <- est.c[i,2]-est.c[j,2]
      #distances in km (gross method, using flat geometry):
      dx <- dx*111*cos((est.c[i,2]+est.c[j,2])*pi/360)
      dy <- dy*111
      dz <- (est.c[i,3]-est.c[j,3])*wz #vertical distance
      d2 <- dx*dx+dy*dy+dz*dz #3D quadratic distance
      est.d[i,j] <- sqrt(d2) #3D distance
      #prioritize observations by moving away series from reanalysis:
      if(raway & xor(refhom[i],refhom[j])) est.d[i,j] <- est.d[i,j]+1000
      est.d[j,i] <- est.d[i,j]  #simmetric matrix
    }
  }
  cat("\n")
  est.p <- t(apply(est.d,1,order)) #proximity ranks matrix
  #- Firts estimation of means and standard deviations
  datmed <- apply(dat,1,mean,na.rm=TRUE) #global mean series
  refmed <- mean(datmed) #reference global mean
  dat.m <- apply(dat,2,mean,na.rm=TRUE) #initial (raw) means of the series
  if(std==2 & min(dat.m)==0) { #avoid divisions by zero:
    j0 <- which(dat.m==0) #series with zero mean
    minN0 <- min(dat.m[dat.m>0]) #minimum non zero mean
    dat.m[j0] <- minN0 #assign minimum non zero mean to previous zeros
  }
  if(std==3) {
    refstd <- sd(datmed) #reference global standard deviation
    dat.s <- apply(dat,2,sd,na.rm=TRUE) #initial (raw) standard deviations
  }
  switch(std,
    dat.m <- dat.m + refmed - apply(!is.na(dat),2,datmed.mean),
    dat.m <- dat.m * refmed / apply(!is.na(dat),2,datmed.mean),
    {dat.m <- dat.m + refmed - apply(!is.na(dat),2,datmed.mean)
     dat.s <- dat.s + refstd - apply(!is.na(dat),2,datmed.sd)},
    dat.m <- dat.m + refmed - apply(!is.na(dat),2,datmed.mean)
  )
  #- metad==TRUE? Read *_brk.csv and split the series by the break-points
  if(metad) {
    cat('\nSplitting the series following the metadata file...:\n')
    if(sufbrk=='') fichbrk <- sprintf('%s_%d-%d_brk.csv',varcli,anyi,anyf)
    else fichbrk <- sprintf('%s-%s_%d-%d_brk.csv',varcli,sufbrk,anyi,anyf)
    brk <- read.csv(fichbrk,colClasses=c("character","character","character"))
    if(!is.na(tinc)) brk[,2] <- as.POSIXct(brk[,2],tz=tz)
      else brk[,2] <- as.Date(brk[,2])
    nbrk <- nrow(brk); nn <- 0
    if(nbrk<1) cat('No break-points in the metadata file.\n')
    else {
      for(kb in 1:nbrk) { #for every break:
        i <- match(brk[kb,1],est.c[,4]) #series to split
        if(is.na(i)) {
          cat(sprintf('\nCode %s not found in station list; break skipped',brk[kb,1]))
          next
        }
        kp <- match(brk[kb,2],x) #break-point location
        #check dates concordancies:
        if(is.na(kp)) kp <- which.max(x>brk[kb,2])
        if(is.na(tinc)) cat(sprintf('\n%s(%d) breaks at %s',est.c[i,4],i,format(x[kp])))
        else cat(sprintf('\n%s(%d) breaks at %s',est.c[i,4],i,format(x[kp],tz=tz,usetz=TRUE)))
        if(sum(!is.na(dat[1:(kp-1),i])) < mindat) {
          dat[1:(kp-1),i] <- NA
          cat(" Fragment with less than",mindat,"data DELETED\n")
        } else if(sum(!is.na(dat[kp:nd,i])) < mindat) {
          dat[kp:nd,i] <- NA
          cat(" Fragment with less than",mindat,"data DELETED\n")
        } else {
          nn <- nn+1 #increment no. of new series
          iest <- c(iest,iest[i]) #add original series index
          nsp[iest[i]] <- nsp[iest[i]]+1 #and its no. of breaks
          if(nm>0) { #count no. of breaks per year
            z <- 1 + floor((kp-1)/nm) #anual break location
            nsy[z] <- nsy[z] + 1 #no. of breaks per year
          }
          dat <- cbind(dat,rep(NA,nd)) #new data column
          #move pre-cut data to the new series:
          dat[1:(kp-1),ne+nn] <- dat[1:(kp-1),i]
          dat[1:(kp-1),i] <- NA #delete pre-cut data
          #copy coordinates and add a suffix to station code and name:
          z <- data.frame(est.i[iest[i],1:3],paste(est.i[iest[i],4],"-",1+nsp[iest[i]],sep=""),paste(est.i[iest[i],5],"-",1+nsp[iest[i]],sep=""))
          names(z) <- names(est.i)
          est.c <- rbind(est.c,z)
          #asign same means (and std. devs.?) to the new fragment:
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
      ne <- ne + nn  #update no. of stations
      cat(ne,"\n\n")
      refhom <- substr(est.c[,4],1,1)=='*' #update homogeneous references
    }
    inht <- c(0,0) #go to missing data filling
  }
  #- for (ks in 1:3) #(test in windows, in whole series and missing d. filling)
  for (ks in 1:3) { #for every stage:
    if(ks<3) {
      inh <- inht[ks] #inht threshold for stage ks
      if(inh==0) next #skip stage if inh==0
    }
    if(is.na(cumc)) { cat("\n\n========== STAGE",ks)
      switch(ks,
        cat(sprintf(" (%s on overlapping temporal windows) ===========\n\n",inhtest)),
        cat(sprintf(" (%s on the whole series) =======================\n\n",inhtest)),
        cat(" (Final calculation of all missing data) ==========\n\n")
      )
    } else cat(sprintf('\n====== Disaggregating daily accumulated precipitation coded as %d\n',cumc))
    #- compute weight matrix? (depends on the stage)
    if(ks==1) zz <- TRUE else if(wd[ks]!=wd[ks-1]) zz <- TRUE else zz <- FALSE
    if(zz) {
      est.w <- matrix(1,nei,nei) #weight matrix
      if(wd[ks]>0) { #weights different from 1
        cat("Computing inter-station weights...")
        wd2 <- wd[ks]*wd[ks]
        for(i in 1:(nei-1)) {
          for(j in (i+1):nei) {
            est.w[i,j] <- wd2/(wd2+est.d[i,j]*est.d[i,j])
            est.w[j,i] <- est.w[i,j]  #simmetric matrix
          }
        }
        cat(' (done)\n\n')
      }
    }
    #- if(graphics), issue page indicating the stage
    if(graphics) {
      plot(-1:1,-1:1,type="n",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
      text(0,0.4,paste("Stage",ks),cex=3)
      if(ks==1) text(0,-0.3,sprintf("Binary splits on %d term\nstepped windows with\nstd=%d, %s>%d\nand wd=%d km",round(swa),std,inhtest,round(inh),round(wd[ks])),cex=2.2)
      else if(ks==2) text(0,-0.3,sprintf("Binary splits on\nwhole series with\nstd=%d, %s>%d\nand wd=%d km",std,inhtest,round(inh),round(wd[ks])),cex=2.5)
      else text(0,-0.3,sprintf("Final anomalies of the\nhomogenized series with\nwd = %d km and nref = %d",round(wd[ks]),nref[ks]),cex=2)
    }
    #stage dependent values:
    if(ks==3) aref <- TRUE else aref <- FALSE
    nrefk <- nref[ks]
    #- repeat until no series is split
    repeat {
      #---------- Compute anomalies and delete anomalous data:
      #- initialize matrices dat.z|e|c oneref anom sanom mindist nrefs used
      dat.z <- matrix(NA,nd,ne) #observed data (normalized)
      dat.e <- matrix(NA,nd,ne) #estimated data (normalized)
      dat.c <- matrix(NA,nd,ne) #calculated data (estimated)
      oneref <- matrix(FALSE,nd,ne) # 1 reference only?
      anom <- matrix(NA,nd,ne) #anomalies
      sanom <- matrix(NA,nd,ne) #standardized anomalies
      mindist <- matrix(NA,nd,ne) #minimum distances
      nrefs <- matrix(NA,nd,ne) #no. of references
      used <- matrix(FALSE,ne,ne) #flags of used stations
      #working copy of the data:
      dat.d <- dat
      dat.na <- is.na(dat.d) #missing data index
      #- if there are time steps without any data, issue a warnig and finish
      numdat <- apply(!dat.na,1,sum)
      nmin=min(numdat)
      if(nmin==0) {
        cat("\nThere are terms with NO DATA!:\n")
        for(j in which(numdat==0)) cat(format(x[j]),"\n")
        stop("Cannot continue! Shorten the study period, add series with data in the empty terms, or be more tolerant to outliers.")
      }
      #use any existing former means and standard deviations:
      if(exists('dat.m0')) {
        dat.m <- dat.m0
        if(std==3) dat.s <- dat.s0
      }
      #- normalize dat.d (obtain dat.z)
      
      switch(std,
        for(ke in 1:ne) dat.z[,ke] <- dat.d[,ke]-dat.m[ke], #std=1
        for(ke in 1:ne) dat.z[,ke] <- dat.d[,ke]/dat.m[ke], #std=2
        for(ke in 1:ne) dat.z[,ke]<-(dat.d[,ke]-dat.m[ke])/dat.s[ke],#std=3
        dat.z <- dat.d
      )
      #- ite=0 and loop until estimated data converge
      #iterative process for estimating the means of the series:  
      ite <- 0
      if(expl) cat("\nCalculation of missing data\n")
      else if(is.na(cumc)) {
        cat("\nCalculation of missing data with outlier removal\n")
        cat('(Suggested data replacements are provisional)\n')
      }
      if(length(dat)>10000000) cat('This process may take a very long time (many days)\n')
      else if(length(dat)>1000000) cat('This process may take a long time (many hours)\n')
      if(is.na(cumc)) {
        if(ks==3 & !expl) cat("\nThe following lines will have one of these formats:\n")
        cat("  Station(rank) Date: Observed -> Suggested (Anomaly, in std. devs.)\n")
        if(ks==3) cat("  Iteration Max_data_difference (Station_code)\n")
      }
      maxddif0 <- 99999. #initialize the maximum data difference
      repeat {
        ite <- ite+1
        #- ite+=1 and estimate series (dat.e|c) from their neighbors
        #  update used, nrefs and mindist:
        for(i in 1:ne) { #for every station
          if(refhom[i]) next #skip trusted series
          ik <- iest[i] #original station index
          for(j in 1:nd) { #for every data
            se <- 0
            sw <- 0
            nr <- 0
            for(ir in 1:nei) { #for every station (possible reference)
              kr <- est.p[ik,ir]
              krf <- which(iest==kr) #fragments of the reference
              k <- which(!dat.na[j,krf]) #which one has observation?
              if(length(k)!=1) next #no fragment with observation
              k <- krf[k] #index of the fragment with observation
              if(i==k) next #it is the same station
              nr <- nr+1 #no. of references
              used[i,k] <- TRUE #flag used station
              #minimum distance to the nearest data:
              if(nr==1) mindist[j,i] <- max(est.d[ik,kr],1)
              w <- est.w[ik,kr]
              se <- se + w * dat.z[j,k]
              sw <- sw + w
              if(nr>=nrefk) break #if maximum no. of references, finish loop
            }
            if(!nr) { #no reference!
              dat.e[j,i] <- dat.z[j,i] #keep the original data
              nrefs[j,i] <- NA
            } else {
              nrefs[j,i] <- nr
              #with only one reference, flag to avoid series modification:
              if(nr==1 & !is.na(oneref[j,i])) oneref[j,i] <- TRUE
              #avoid negative estimations if std=2 (precipitation, etc):
              if(std==2 & se<0) se <- 0
              dat.e[j,i] <- se / sw #estimated data (normalized)
            }
          }
        }
        #change any NaN to NA. (It may happen with std=2):
        n <- sum(is.nan(dat.e))
        if(n>0) {
          cat(n,"NaN's in dat.e ! (changing them to NA's...)\n")
          dat.e[is.nan(dat.e)] <- NA
        }
        #calculate values unnormalizing dat.e:
        switch(std,
          for(ke in 1:ne) dat.c[,ke] <- dat.e[,ke]+dat.m[ke],     #std=1
          for(ke in 1:ne) dat.c[,ke] <- dat.e[,ke]*dat.m[ke],     #std=2
          for(ke in 1:ne) dat.c[,ke]<-dat.e[,ke]*dat.s[ke]+dat.m[ke], #std=3
          dat.c <- dat.e
        )
        #if cumc was set, dissaggregate data, save them in *.dat and finish:
        if(!is.na(cumc)) {
          for(ke in 1:ne) { #for every station
            cumvl <- length(cumv[[ke]]) #no. of accumulations
            if(cumvl==0) next #series without accumulations
            for(j in 1:cumvl) {
              #an embedded missing data? joint the accumulations:
              if(is.na(cumv[[ke]][j])) {
                if(j<cumvl & cumA[[ke]][j+1]==cumB[[ke]][j]+1)
                  cumA[[ke]][j+1] <- cumA[[ke]][j]
                next
              }
              sumc <- sum(dat.c[cumA[[ke]][j]:cumB[[ke]][j],ke])
              #if estimated data sum is 0, assign zeros and keep the
              #accumulated value in the reported day:
              if(sumc==0.0) { 
                dat[cumA[[ke]][j]:cumB[[ke]][j],ke] <- 0
                dat[cumB[[ke]][j],ke] <- cumv[[ke]][j]
                next #go to next accumulation
              }
              prop <- cumv[[ke]][j]/sumc #ratio accumulation/estimation_sum
              zc <- dat.c[cumA[[ke]][j]:cumB[[ke]][j],ke] * prop
              zk <- which.max(zc) #location of the maximum estimation
              zc <- round(zc,ndec)
              zd <- cumv[[ke]][j] - sum(zc) #rounding difference
              #add difference to the highest value:
              if(zd!=0.0) zc[zk] <- zc[zk]+zd
              dat[cumA[[ke]][j]:cumB[[ke]][j],ke] <- zc
            }
          }
          #save new input data and finish:
          fcum <- sprintf('%s-cum_%d-%d',varcli,anyi,anyf) #accum.file basename
          file.rename(sprintf('%s.dat',fbas),sprintf('%s.dat',fcum))
          write(dat,sprintf('%s.dat',fbas))
          file.rename(sprintf('%s.est',fbas),sprintf('%s.est',fcum))
          write.table(est.c,sprintf('%s.est',fbas),row.names=FALSE,col.names=FALSE)
          cat('\nAccumulated values have been distributed among the previous days and written\n')
          cat('as new input files. Original input files have been renamed to\n')
          cat(sprintf('%s.dat and %s.est\n\n',fcum,fcum))
          cerrar(graphics); return(invisible()) #finish
        }
        #- anomalies calculation (anom, sanom) and outtlier deletion
        anom <- dat.z-dat.e #anomalies
        anom[dat.na] <- NA  #forget anomalies of estimated data
        #normalize anomalies:
        anomm <- apply(anom,2,mean,na.rm=TRUE) #mean anomalies
        anoms <- apply(anom,2,sd,na.rm=TRUE) #std. dev. of anomalies
        for(i in 1:ne) sanom[,i] <- (anom[,i]-anomm[i])/anoms[i]
        if(!expl) { #delete outliers:
          if(ite==1 & sum(dz.max)<=1) { #set dz.* limits
            z <- c(sanom) #normalize anomalies vector
            dz.maxd <- quantile(z,probs=1-dz.max/100,na.rm=TRUE) #max.delete
            dz.maxw <- quantile(z,probs=1-dz.max/10,na.rm=TRUE) #max.warning
            dz.mind <- quantile(z,probs=dz.max/100,na.rm=TRUE) #min.delete
            dz.minw <- quantile(z,probs=dz.max/10,na.rm=TRUE) #min.warning
            dz.max <- 99 #do not repeat this in every stage
          }
          elim <- sanom < dz.mind | sanom > dz.maxd #anomalies to delete
          elim[is.na(elim)] <- FALSE #remove any NA
          elim[,refhom] <- FALSE #keep trusted series
          nelim <- sum(elim) #no. of data to delete
          if(nelim>0) { #delete the anomalous original data
            #list anomalous data to be deleted:
            for(i in 1:ne) {
              for(j in 1:nd) if(elim[j,i] & !is.na(oneref[j,i])) {
                outan[j,iest[i]] <- sanom[j,i] #save the outlier anomaly
                do <- dat.d[j,i] #original data
                dc <- dat.c[j,i] #calculated data
                cat(sprintf('%s(%d) %s',est.c[i,4],i,format(x[j])))
                cat(": ",do," -> ",round(dc,ndec)," (",round(sanom[j,i],2),
                  ")",sep="")
                #do not delete with only one reference!:
                if(oneref[j,i] & nrefk>1) {
                  cat(" Only 1 reference! (Unchanged)")
                  elim[j,i] <- FALSE
                }
                else { #write in Fout with flag 1
                  write(c(est.c[iest[i],4],format(x[j]),round(do,ndec),
                  round(dc,ndec),round(sanom[j,i],2),1),Fout,ncolumns=6,sep=',')
                }
                cat("\n")
              }
            }
            dat[elim] <- NA #delete anomalous data
            dat.na[elim] <- TRUE #update missing data flags
            maxddif0 <- 99999. #reset to avoid fake convergence break
          }
          else if(!aref) cat('(No detected outliers)\n')
        }
        #list suspect values? :
        if(ks==3 & !expl & (!is.na(dz.maxw) | !is.na(dz.minw))) {
          #suspect values but not big outliers:
          susp <- (sanom < dz.minw & sanom >= dz.mind) | 
                  (sanom > dz.maxw & sanom <= dz.maxd)
          susp[is.na(susp)] <- FALSE #remove any NA
          susp[,refhom] <- FALSE #unflag trusted series
          nsusp <- sum(susp) #no. of suspect data
          if(nsusp>0) { #list suspect observations:
            for(i in 1:ne) {
              for(j in 1:nd) if(susp[j,i] & !is.na(oneref[j,i])) {
                do <- dat.d[j,i] #original data
                dc <- dat.c[j,i] #calculated data
                #unflag suspect data with only one available reference!:
                if(oneref[j,i] & nrefk>1) susp[j,i] <- FALSE
                else { #write in Fout with flag 0
                  write(c(est.c[iest[i],4],format(x[j]),round(do,ndec),
                  round(dc,ndec),round(sanom[j,i],2),0),Fout,ncolumns=6,sep=',')
                }
              }
            }
          }
        }
        #- missing data filling
        dat.d[dat.na] <- dat.c[dat.na] 
        if(ite>1) {
          ddif <- dat.d-dat.d0 #data differences from previous iteration
          maxddif <- max(abs(ddif),na.rm=TRUE) #max. data difference
          kmaxdif <- which.max(abs(ddif)) #max. dat. dif. location
          kmaxest <- ceiling(kmaxdif/nd) #max. dat. dif. station
        }
        dat.d0 <- dat.d #data copy
        #- update dat.m|s|z
        dat.m <- apply(dat.d,2,mean,na.rm=TRUE)
        if(std==3) dat.s <- apply(dat.d,2,sd,na.rm=TRUE)
        switch(std,
          for(ke in 1:ne) dat.z[,ke] <- dat.d[,ke]-dat.m[ke], #std=1
          for(ke in 1:ne) dat.z[,ke] <- dat.d[,ke]/dat.m[ke], #std=2
          for(ke in 1:ne) dat.z[,ke]<-(dat.d[,ke]-dat.m[ke])/dat.s[ke], #std=3
          dat.z <- dat.d
        )
        #- if(!aref) break (no need to refine missing data until the end)
        if(!aref) break
        #- if(ite>1) and convergence reached or broken, break loop
        if(ite>1) {
          cat(ite,' ',round(maxddif,ndec+2)," (",est.c[kmaxest,4],")\n",sep="")
          if(maxddif>maxddif0) {
            cat("Data convergence broken\n\n")
            break
          }
          if(maxddif<=maxdif) {
            cat("Prescribed convergence reached\n\n")
            break
          }
          if(ite==maxite) {
            cat("\nAverage calculation skipped after",ite,"iterations\n")
            break
          }
          maxddif0 <- maxddif #save max. dat. dif. for next iteration
        }
      }
      #- save dat.m|s in dat0.m|s
      dat.m0 <- dat.m #copy of the means
      if(std==3) dat.s0 <- dat.s #copy of std. deviations
      oneref[is.na(oneref)] <- TRUE #reset oneref
      #- if(aref==TRUE) repeat missing data filling with self-correction
      if(aref==TRUE) {
        cat('Last series readjustment (please, be patient...)\n')
        #- obtain estimated series (dat.e, dat.c) from their neighbors
        #- and update used[ne,ne], nrefs[nd,ne] and mindist[ne,ne]:
        for(i in 1:ne) { #for every statin
          ik <- iest[i] #index of original station
          for(j in 1:nd) { #for every data
            se <- 0
            sw <- 0
            nr <- 0
            for(ir in 1:nei) { #for every possible reference
              kr <- est.p[ik,ir]
              krf <- which(iest==kr) #fragments of the reference
              k <- which(!dat.na[j,krf]) #which one has observation?
              if(length(k)!=1) next #no fragment with observation
              k <- krf[k] #index of the fragment with observation con dato
              if(i==k) next #same station
              nr <- nr+1 #no. of references
              used[i,k] <- TRUE #flag used station
              #minimum distance to the nearest data:
              if(nr==1) mindist[j,i] <- max(est.d[ik,kr],1)
              w <- est.w[ik,kr]
              se <- se + w * dat.z[j,k]
              sw <- sw + w
              #if maximum no. of references or self-reference, finish:
              if(nr>=nrefk | (aref & ir==1)) break
            }
            if(!nr) { #no reference!
              dat.e[j,i] <- dat.z[j,i] #keep observation
              nrefs[j,i] <- NA
            } else {
              nrefs[j,i] <- nr
              #if only one reference, flag to avoid any modification:
              if(nr==1 & !is.na(oneref[j,i])) oneref[j,i] <- TRUE
              #avoid negative estimations if std=2 (precipitation, etc):
              if(std==2 & se<0) se <- 0
              dat.e[j,i] <- se / sw #estimated data (normalized)
            }
          }
        }
        #change any NaN to NA. (It may happen with std=2):
        n <- sum(is.nan(dat.e))
        if(n>0) {
          cat(n,"NaN's in dat.e ! (changing them to NA's...)\n")
          dat.e[is.nan(dat.e)] <- NA
        }
        #calculate values unnormalizing dat.e:
        switch(std,
          for(ke in 1:ne) dat.c[,ke] <- dat.e[,ke]+dat.m[ke],     #std=1
          for(ke in 1:ne) dat.c[,ke] <- dat.e[,ke]*dat.m[ke],     #std=2
          for(ke in 1:ne) dat.c[,ke]<-dat.e[,ke]*dat.s[ke]+dat.m[ke],#std=3
          dat.c <- dat.e
        )
        #- missing data filling
        dat.d[dat.na] <- dat.c[dat.na]
        if(!is.na(vmax)) dat.d[dat.d > vmax] <- vmax
        if(!is.na(vmin)) dat.d[dat.d < vmin] <- vmin
        #- update dat.m|s|z
        dat.m <- apply(dat.d,2,mean,na.rm=TRUE)
        if(std==3) dat.s <- apply(dat.d,2,sd,na.rm=TRUE)
        switch(std,
          for(ke in 1:ne) dat.z[,ke] <- dat.d[,ke]-dat.m[ke], #std=1
          for(ke in 1:ne) dat.z[,ke] <- dat.d[,ke]/dat.m[ke], #std=2
          for(ke in 1:ne) dat.z[,ke]<-(dat.d[,ke]-dat.m[ke])/dat.s[ke], #std=3
          dat.z <- dat.d
        )
      }
      #- calculate final values of the anomalies (anom, sanom)
      anom <- dat.z-dat.e #anomalies
      anom[dat.na] <- NA  #forget anomalies of estimated data
      anomm <- apply(anom,2,mean,na.rm=TRUE) #mean anomalies
      anoms <- apply(anom,2,sd,na.rm=TRUE) #std. dev. of anomalies
      for(i in 1:ne) sanom[,i] <- (anom[,i]-anomm[i])/anoms[i]
      #- ----------- Shift detection/correction (binary split):
      #- if(ks>2) break (only missing data filling in the last stage)
      if(ks>2) break
      #split series when maximum test is higher than inh threshold:
      nn <- 0 #initialize no. of new series
      splt <- rep(0,ne)  #initialize tV of split series
      modif <- FALSE #initialize flag of series modification
      cat("\nPerforming shift analysis on the",ne,"series...\n")
      y <- sanom #normalized anomalies
      y[is.na(nrefs)] <- NA #remove anomalies without reference
      y[,refhom] <- NA #remove anomalies of trusted series
      if(ks==1) tkx <- apply(y,2,wtest,swa,sts,test)
      else tkx <- apply(y,2,test,sts)
      tVx <- tkx[1,]; kpx<- as.integer(tkx[2,])
      #- split series in decreasing order of tVx while not using as references
      #  recently split with tVx similar to the maximum tVx in all series:
      if(sum(!is.na(tVx))==0) tVxx <- 0 else tVxx <- max(tVx,na.rm=TRUE)
      while(tVxx > inh) {
        i <- which.max(tVx) #series with maximum tVx
        #if i used references split with a too high tVx, go to new iteration:
        if(max(splt[used[i,]])>tVxx*(1+.05*min(nr,sum(used[i,])))) break
        kp <- kpx[i] #location of tVx in series i
        if(oneref[kp,i] & nrefk>1) { #avoid split with only one reference
          tVx[i] <- -1 #set this series tVx to -1
          tVxx <- max(tVx,na.rm=TRUE) #maximum tVx of remaining series
          next
        }
        cat(sprintf('\n%s(%d) breaks at %s (%.1f)',est.c[i,4],i,
          format(x[kp]),tVx[i]))
        write(sprintf('%s,%s,%.1f',est.c[iest[i],4],format(x[kp]),tVx[i]),
          Fbrk,ncolumns=3)
        #graphic of anomalies with split location:
        if(graphics) {
          y <- sanom[,i] #vector of anomalies of this series
          ylab="Standardized anomalies (observed - computed)"
          tit <- sprintf('%s   %d (%s)\n%s',varcli,i,est.c[i,4],est.c[i,5])
          plot(x,y,type="h",lwd=lw,ylim=c(-5,5),main=tit,xlab='Time',
            ylab=ylab,col=hsv(.7,1,.9))
          z <- class(x)
          if(z[1]=='Date' | z[1]=='POSIXct') {
            grid(NA,NULL,col=grdcol)
            if(z[1]=='Date') abline(v=axis.Date(1,x),lty=3,col=grdcol)
            else if(z[1]=='POSIXct') abline(v=axis.POSIXct(1,x),lty=3,
              col=grdcol)
          } else grid(col=grdcol)
          abline(-3,0,lty=3,col=grdcol); abline(-5,0,lty=3,col=grdcol)
          lines(x,log10(nrefs[,i])-5,col='orange2')
          lines(x,log10(mindist[,i])-5,col=verde)
          mtext(" 1",4,las=1,adj=0,at=-5,col=verde)
          mtext(" 10",4,las=1,adj=0,at=-4,col=verde)
          mtext(" 100",4,las=1,adj=0,at=-3,col=verde)
          mtext("min.d.",4,las=1,adj=0,at=-5.4,col=verde)
          mtext(" (km)",4,las=1,adj=0,at=-2,col=verde)
          mtext("n.ref.",4,las=1,adj=0,at=-5.8,col='orange2')
          lines(rep(x[kp],2),c(-5,4.8),col="red",lty=2) #show split location
          text(x[kp],5,floor(tVxx))
        }
        #count no. of splits per year
        z <- as.integer(strftime(x[kp],"%Y"))-anyi+1 #year of the split
        nsy[z] <- nsy[z] + 1 #no. of splits per year
        #split series at the break point:
        nd1 <- sum(!is.na(dat[1:(kp-1),i])) #no. of data of fragment 1
        nd2 <- sum(!is.na(dat[kp:nd,i])) #no. of data of fragment 2
        if(nd1 < mindat & nd2 < mindat) stop(sprintf('\nBoth fragments have less than %d data. Please, remove this series\nfrom the input files or increase the mindat parameter.',mindat))
        del <- 0 #flag of deleted fragments
        if(std==2) {
          md1 <- mean(dat[1:(kp-1),i],na.rm=TRUE)
          md2 <- mean(dat[kp:nd,i],na.rm=TRUE)
        } else md1 <- md2 <- NA
        if(nd1 < mindat) {
          dat[1:(kp-1),i] <- NA; del <- del + 1
          cat(" Fragment with less than",mindat,"data DELETED")
        }
        if(nd2 < mindat) {
          dat[kp:nd,i] <- NA; del <- del + 1
          cat(" Fragment with less than",mindat,"data DELETED")
        }
        if(del==0) {
          nn <- nn+1 #update no. of new series
          iest <- c(iest,iest[i]) #add index of original series
          nsp[iest[i]] <- nsp[iest[i]]+1 #update no. of shifts
          dat <- cbind(dat,rep(NA,nd)) #new data column
          #move pre-cut data to the new series:
          dat[1:(kp-1),ne+nn] <- dat[1:(kp-1),i]
          dat[1:(kp-1),i] <- NA #delete pre-cut data
          #copy coordinates and add a suffix to station code and name:
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
        #update tVx and flags for next loop:
        modif <- TRUE #flag of modified series
        splt[i] <- tVx[i] #split tV of series i
        tVx[i] <- 0 #set its tVx to 0
        tVxx <- max(tVx,na.rm=TRUE) #maximum tVx of remaining series
      }
      if(nn) {
        cat("\n\nUpdate number of series: ",ne,"+",nn,"= ")
        ne <- ne+nn  #update no. of series
        cat(ne,"\n")
        refhom <- c(refhom,rep(FALSE,nn)) #update refhom vector
      }
      #- No new splits? histograms of tVx and breaks
      if(!nn & !modif) {
        if(graphics) {
          #histogram of global tVx (without unreallistic zeros):
          z <- tVx[!is.na(tVx) & tVx>0]
          main <- sprintf("Histogram of maximum %s (Stage %d)",inhtest,ks)
          if(sum(!is.na(z))) hist(z,breaks=20,xlab=inhtest,col="purple",main=main)
          if(ks==2 | inh<1) {
            #histogram of no. of splits per station
            #(too small fragments are not counted for, although they will be
            # reported in the final break list):
            hist(nsp,breaks=0:max(9,max(nsp)+1)-.5,col="orange2",xlab="Number of splits",ylab="Number of stations",main="Number of splits per station")
            #no. of splits per year:
            w <- min(5,ceiling(400/na)) #addapt bar width
            plot(anyi:anyf,nsy,type="h",lwd=w,col=2,ylim=c(0,max(10,max(nsy))),xlab="Years",ylab="Number of splits",main="Number of splits per year")
            grid(col=grdcol)
          }
        }
        #list of possible splits with only one referencia:
        z <- which(tVx<0)
        if(length(z)>0) {
          cat('Series that could break but had only one reference:\n')
          print(est.c[z,4])
        }
        break #break loop and go to next stage
      }
    }
  }
  #------------ End of the three homogeneization stages --------------
  #RMSE of estimated data:
  z <- dat.c; zo <- dat.o[,iest]; zo[dat.na] <- NA
  rmse <- apply((z-zo)^2,2,function(x) sqrt(mean(x,na.rm=TRUE)))
  #- graphics of anomalies of the homogenized series
  #  (with maximium tVx, grouped by original series):
  tVx <- rep(NA,ne) #(save final tVx in overlapping windows)
  inhx <- rep(NA,ne) #(save final tVx in complete series)
  for(io in 1:nei) { #for every original series
    wi <- which(iest==io) #series derived from station io
    lwi <- length(wi)
    for(i in wi) { #for every derived series
      y <- sanom[,i] #standardized anomalies of the series
      if(graphics) {
        ylab="Standardized anomalies (observed - computed)"
        tit <- sprintf('%s   %d (%s)\n%s',varcli,i,est.c[i,4],est.c[i,5])
        plot(x,y,type="h",lwd=lw,ylim=c(-5,5),main=tit,xlab='Time',ylab=ylab,col=hsv(.7,1,.9))
        z <- class(x)
        if(z[1]=='Date' | z[1]=='POSIXct') {
          grid(NA,NULL,col=grdcol)
          if(z[1]=='Date') abline(v=axis.Date(1,x),lty=3,col=grdcol)
          else if(z[1]=='POSIXct') abline(v=axis.POSIXct(1,x),lty=3,col=grdcol)
        } else grid(col=grdcol)
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
      #apply wtest and flag its maximum tV (if >=1):
      st <- wtest(y,swa,sts,test); tVx[i] <- st[1]; zz <- floor(st[1])
      if(zz) {
        kp <- as.integer(st[2])
        if(graphics) {
          lines(rep(x[kp],2),c(-5,4.8),col=verde,lty=2) #flag maximum wtest
          text(x[kp],5,zz,col=verde) #value
        }
      }
      #apply test and flag its maximum:
      st <- eval(call(test,y,sts))
      inhx[i] <- round(st[1],1); zz <- floor(st[1])
      if(!is.na(zz) & zz & graphics) {
        kp <- round(st[2])
        lines(rep(x[kp],2),c(-5,4.8),lty=4) #flag tVx (maximum test)
        text(x[kp],-5.2,zz) #value
      }
    }
  }
  #homogenized series:
  dah <- round(dat.d,ndec) #rount to the requested no. of decimals
  #- graphics of homogenized series and their corrections
  if(graphics) {
    plot(-1:1,-1:1,type="n",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
    text(0,0.4,"Final graphics",cex=3)
    text(0,-0.3,"Adjusted series and\napplied corrections",cex=2.2)
    if(nm>0) xlab <- "Years" else xlab <- "Dates"
    layout(matrix(1:2,2,1,byrow=TRUE))
    #filters to get annual values:
    ndpy <- round(nd/na) #no. of data per year
    if(nd<=120 | nd<ndpy*3) { fltr=1; ylabd <- "Data" }#few data? avoid filter
    else {
      fltr <- rep(1,ndpy)
      if(annual=='sum') ylabd <- "Running annual totals"
      else {
        ylabd <- "Running annual means"
        fltr <- fltr/ndpy
      }
    }
    if(!is.na(uni)) ylabd <- sprintf('%s (%s)',ylabd,uni)
    for(i in 1:nei) { #for every original station
      wi <- which(iest==i) #derived series of station i
      lwi <- length(wi)
      if(lwi>1) vi <- TRUE else vi <- FALSE
      #filtros para valores anuales 
      tit <- sprintf('%s   %d (%s)\n%s',varcli,i,est.c[i,4],est.c[i,5])
      yo <- as.vector(dat.o[,i]) #original observations
      y <- dah[,wi] #homogenized data
      par(mar=c(0,4,4,2),xaxt="n",cex=cex,cex.axis=cex,cex.lab=cex)
      yf <- stats::filter(y,fltr)
      ylim <- c(floor(min(yf,na.rm=TRUE)),ceiling(max(yf,na.rm=TRUE)))
      #(avoid matplot because of bad managing dates in X axes)
      plot(x,stats::filter(yo,fltr),type="n",ylim=ylim,ylab=ylabd,main=tit)
      matlines(x,yf,lty=1,col=2:20)
      lines(x,stats::filter(yo,fltr)) #redraw observations line
      par(xaxt="s")
      z <- class(x)
      if(z[1]=='Date' | z[1]=='POSIXct') {
        grid(NA,NULL,col=grdcol)
        if(z[1]=='Date') abline(v=axis.Date(1,x,tick=FALSE,labels=FALSE),
          lty=3,col=grdcol)
        else if(z[1]=='POSIXct') abline(v=axis.POSIXct(1,x,tick=FALSE,
          labels=FALSE),lty=3,col=grdcol)
      } else grid(col=grdcol)
      par(mar=c(5,4,0,2))
      #corrections:
      if(std==2) {
        yo[yo==0] <- NA; y[y==0] <- NA
        yd <- y/yo; ylab <- "Correction factors"
      } else {
        yd <- y-yo; ylab <- "Correction terms"
        if(!is.na(uni)) ylab <- sprintf('%s (%s)',ylab,uni)
      }
      ylim <- c(floor(min(yd,na.rm=TRUE)),ceiling(max(yd,na.rm=TRUE)))
      if(std==2 & ylim[2]<=2) ylim <- c(0,2)
      if(std!=2 & ylim[1]>=-1 & ylim[2]<=1) ylim <- c(-1,1)
      if(vi) plot(x,yd[,1],type="n",ylim=ylim,ylab=ylab,xlab='Time')
      else plot(x,yd,type="n",ylim=ylim,ylab=ylab,xlab='Time')
      matlines(x,yd,type="l",lty=1,col=2:20,xaxt='n')
      z <- class(x)
      if(z[1]=='Date' | z[1]=='POSIXct') {
        grid(NA,NULL,col=grdcol)
        if(z[1]=='Date') abline(v=axis.Date(1,x),lty=3,col=grdcol)
        else if(z[1]=='POSIXct') abline(v=axis.POSIXct(1,x),lty=3,col=grdcol)
      } else grid(col=grdcol)
    }
    par(my.par)
    par(las=1,cex=cex,cex.axis=cex,cex.lab=cex)
  }
  if(sum(inht)>0) cat("\n======== End of the homogenization process, after ")
  else cat("\n======== End of the missing data infilling process, after ")
  cat(format(round(Sys.time()-time1,2)),'\n')
  cat("\n----------- Final calculations:\n")
  #inhtest in each series
  if(inhtest=='SNHT') cat("\nSNHT: Standard normal homogeneity test (on anomaly series)\n")
  else cat("\nCucT: Cucconi test (on anomaly series)\n")
  print(summary(round(inhx,1)))
  #RMSE of estimated data (unnormalized):
  cat("\nRMSE: Root mean squared error of the estimated data\n")
  zz <- summary(rmse)
  print(zz)
  sedec <- max(1,2-ceiling(log10(zz[4]))) #no. of decimals of RMSE
  rmse <- round(rmse,sedec) #round RMSE
  pod <- floor(100*(nd-apply(dat.na,2,sum))/nd) #percentage of original data
  cat("\nPOD: Percentage of original data\n")
  print(summary(pod))
  #- print summary of results
  cat("\n")
  df <- data.frame(SNHT=inhx,RMSE=rmse,POD=pod,Code=est.c[,4],Name=est.c[,5])
  names(df)[1] <- inhtest
  print(df,right=FALSE)
  cat(sprintf('\nFrequency distribution tails of residual anomalies and %s:\n\n',inhtest))
  cat("Left tail of standardized anomalies:\n")
  print(round(quantile(sanom,probs=c(.001,.002,.005,.01,.02,.05,.1),na.rm=TRUE),1))
  cat("Right tail of standardized anomalies:\n")
  print(round(quantile(sanom,probs=c(.9,.95,.98,.99,.995,.998,.999),na.rm=TRUE),1))
  cat(sprintf('Right tail of %s on windows of %d terms with up to %d references:\n',inhtest,2*swa,nref[3]))
  print(round(quantile(tVx,probs=c(.9,.95,.98,.99,.995,.998,.999),na.rm=TRUE),1))
  cat(sprintf('Right tail of %s with up to %d references:\n',inhtest,nref[3]))
  print(round(quantile(inhx,probs=c(.9,.95,.98,.99,.995,.998,.999),na.rm=TRUE),1))
  #add new columns to the stations table (percentage of original data,
  #  inhx and RMSE):
  est.c <- cbind(est.c,pod,inhx,rmse)
  #round input data for graphics and final results
  dat <- round(dat,ndec)
  #- if(graphics), plot last graphics
  if(graphics) {
    #histogram of the anomalies (those of outliers in red):
    main <- "Histogram of standardized anomalies"
    z <- hist(c(sanom,outan),plot=FALSE)
    zx <- z$breaks
    zy <- z$counts; zy[zy==0] <- NA
    barplot(zy,log='y',space=0,ylim=c(.9,max(zy,na.rm=TRUE)*2),ylab='Frequency',col='green',main=main,xlab='Anomalies (standard deviations)')
    axis(1,1:length(zx)-1,labels=as.character(zx))
    if(sum(!is.na(outan))) { #redraw outliers in red
      zy <- hist(outan,breaks=zx,plot=FALSE)$counts; zy[zy==0] <- NA
      barplot(zy,log='y',space=0,ylim=c(.9,max(zy,na.rm=TRUE)*2),col=hsv(0,.75),add=TRUE)
    }
    #histogram of tVx by overlapping windows (withoud unrealistic zeros):
    z <- tVx[!is.na(tVx) & tVx>0]
    main <- sprintf("Histogram of maximum windowed %s",inhtest)
    if(sum(!is.na(z))) hist(z,breaks=20,xlab=inhtest,col=verde,main=main)
    #histogram of tVx of the complete series:
    z <- inhx; main <- sprintf("Histogram of maximum global %s",inhtest)
    if(sum(!is.na(z))) hist(z,breaks=20,xlab=inhtest,col="purple",main=main)
    #quality/singularity graphic:
    if(is.na(uni)) xlab <- 'RMSE' else xlab <- sprintf('RMSE (%s)',uni)
    plot(rmse,inhx,type="n",xlim=c(0,max(1,max(rmse,na.rm=TRUE))),ylim=c(0,max(50,max(inhx,na.rm=TRUE))),xlab=xlab,ylab=inhtest,main="Station's quality/singularity")
    grid(col=grdcol)
    text(rmse,inhx,col=hsv(.7,1,.9))
  }
  if(graphics) { par(old.par); graphics.off() } #close graphic output
  #- save results in a rda file
  dat <- dat.o
  names(est.c) <- c('X','Y','Z','Code','Name','pod',test,'rmse')
  rownames(est.c) <- 1:ne
  if(graphics & exists('ct')) save(dat,dah,nd,ndec,uni,est.c,corm,ct,nei,ne,nm,std,x,ini, file=sprintf('%s.rda',fbas))
  else save(dat,dah,nd,ndec,uni,est.c,nei,ne,nm,std,x,ini, file=sprintf('%s.rda',fbas))
  #sort fies of outliers and breaks:
  if(!metad) { 
    close(Fbrk)
    brk <- read.csv(sprintf('%s_brk.csv',fbas),colClasses=c("character","character","numeric"))
    brk <- brk[order(brk[,1],brk[,2]),]
    write.csv(brk,sprintf('%s_brk.csv',fbas),row.names=FALSE)
  }
  close(Fout)
  out <- read.csv(sprintf('%s_out.csv',fbas),colClasses=c("character","character","numeric","numeric","numeric","numeric"),check.names=FALSE)
  #remove duplicates keeping the last rows:
  out <- out[as.integer(row.names(unique(out[,1:2],fromLast=TRUE))),]
  out <- out[order(out[,1],out[,2]),] #sort by stations and dates
  write.csv(out,sprintf('%s_out.csv',fbas),row.names=FALSE)
  cat("\n----------- Generated output files: -------------------------\n\n")
  if(logf) cat(sprintf('%s.txt :  Text output of the whole process',fbas),'\n')
  cat(sprintf('%s_out.csv :  List of corrected outliers',fbas),'\n')
  cat(sprintf('%s_brk.csv :  List of corrected breaks',fbas),'\n')
  if(graphics) cat(sprintf('%s.pdf :  Diagnostic graphics',fbas),'\n')
  cat(sprintf('%s.rda :  Homogenization results.',fbas))
  cat(' Postprocess with (examples):\n')
  cat(sprintf('   dahstat(\'%s\',%d,%d) #averages',varcli,anyi,anyf,fbas),'\n')
  cat(sprintf('   dahstat(\'%s\',%d,%d,stat=\'tnd\') #OLS trends and p-values',varcli,anyi,anyf),'\n')
  cat(sprintf('   dahgrid(\'%s\',%d,%d,grid=YOURGRID) #homogenized grids',varcli,anyi,anyf),'\n')
  cat('   ... (See other options in the package documentation)\n\n')
  while (sink.number()>0) sink() #close log file(s)
}

#- outrename.- Append a suffix to the output files, to avoid overwrites.
outrename <- function(varcli, anyi, anyf, suffix, restore=FALSE) {
#varcli: Short name of the studied climatic variable.
#anyi: Initial year of the study period.
#anyf: Final year of the study period.
#suffix: Suffix to be inserted (or removed) in the output file names.
#restore: If TRUE, the suffix will be removed.
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

#- QCthresholds.- Obtain monthly thresholds for Quality Control alerts.
QCthresholds <- function(dat, ndec=1, probs=c(0.,.001,.01,.99,.999,1.),
minval=NA, maxval=NA, homog=TRUE, verb=TRUE) {
#dat: either the name of a *.rda file of Climatol homogenization results or a
#  data.frame of daily (or subdaily) data in columns, dates or date/times (of
#  class Date or POSIXct) in the first column and station codes in the header
#ndec: number of decimals of output values [1] (defaults shown between brackets)
#probs: probabilities of the quantiles to be computed [0.,.001,.01,.99,.999,1.]
#minval: minimum value to compute runs of constant values (e.g., set to .1 with
#  daily precipitation to avoid the calculation of long runs of zeros)
#maxval: maximum value to compute runs of constant values (e.g., set to 97 with relative humidity to avoid reporting long runs of values greater than 97% in episodes of persistent fog).
#homog: use homogenized data if a *.rda file is used as input [TRUE]
#verb: list all calculated values? [TRUE]
  nei <- est.c <- dah <- NULL #(avoid invisible bidings)
  if(is.data.frame(dat)) { #check the data frame:
    nd <- nrow(dat) #no. of data rows
    x <- dat[,1] #time vector
    if(!inherits(x,'Date') & !inherits(x,'POSIXct')) stop('The first column should be of class Date (or POSIXct in the case of subdaily data')
    ne <- ncol(dat) #no. of data columns
    codes <- names(dat[,2:ne]) #station codes
    dat <- as.matrix(dat[,2:ne]); ne <- ne-1 #no. of stations
  } else if(is.character(dat)) { #check the input file:
    n <- nchar(dat)
    if(substr(dat,n-3,n)!='.rda') stop('Files of Climatol homogenization results have extension .rda')
    if(!file.exists(dat)) stop(sprintf('File %s not found',dat))
    load(dat)
    ne <- nei; codes <- est.c[1:nei,4]
    if(homog) dat <- dah[,1:nei]
  }
  np <- length(probs) #no. of probability levels
  probs2 <- probs[probs>0.5]; np2 <- length(probs2) #right tail probabilities
  st <- as.factor(t(matrix(rep(codes,nd),ne,nd))) #factor of station codes
  mm <- as.factor(matrix(rep(strftime(x,'%m'),ne),nd,ne)) #factor of months
  #quantiles of the data:
  thr1 <- round(unlist(tapply(dat,list(st,mm),quantile,probs=probs,
    na.rm=TRUE)),ndec)
  dim(thr1) <- c(np,ne,12); dimnames(thr1) <- list(probs,codes,1:12)
  #quantiles of first differences:
  st <- as.factor(t(matrix(rep(codes,nd-1),ne,nd-1))) #station codes
  thr2 <- round(unlist(tapply(abs(diff(dat)),list(st),quantile,probs=probs2,
    na.rm=TRUE)),ndec)
  dim(thr2) <- c(np2,ne); thr2 <- t(thr2); dimnames(thr2) <- list(codes,probs2)
  #quantiles of longest runs of constant values
  if(!is.na(minval)) dat[dat<minval] <- NA #delete too low values
  if(!is.na(maxval)) dat[dat>maxval] <- NA #delete too high values
  z <- sapply(apply(dat,2,rle), function(x) x[1]) #run lengths
  thr3 <- t(sapply(z,quantile,probs=probs2,na.rm=TRUE))
  dimnames(thr3) <- list(codes,probs2)
  if(verb) { #list all calculated thresholds:
    cat('\n=========== thr1: Monthly quantiles of the data\n\n')
    for(i in 1:ne) {
      cat('--------- Station',codes[i],'\n')
      print(thr1[,i,])
    }
    cat('\n=========== thr2: Quantiles of the first differences\n\n')
    print(thr2)
    cat('\n=========== thr3: Quantiles of run lengths of constant values')
    if(!is.na(minval)) cat(' >=',minval)
    cat('\n\n')
    print(thr3)
  }
  #save the results in R binary format:
  save(thr1,thr2,thr3, file='QCthresholds.Rdat')
  cat('\nThresholds thr1,thr2,thr3 saved into QCthresholds.Rdat\n')
  cat('(Rename this file to avoid overwriting it in the next run.)\n')
}

#- read.dat.- Read input data (for several climatol functions).
read.dat <- function(varcli, anyi, anyf, x=NULL, ini=NA, tinc=NA, tz='utc',
na.strings=NA) {
#varcli: Short name of the studied climatic variable.
#anyi: Initial year of the study period.
#anyf: Final year of the study period.
#x: Time vector. Automatically set by default. 
#ini: Initial date, with format 'AAAA-MM-DD' (for daily data).
#tinc: Time increment between data. Not set by default, but must be defined
#   for subdaily data, with units 'hours', 'mins' or 'secs'.
#tz: Time zone (if data are subdaily). ('utc' by by default)
#na.strings: strings marking missing data (NA by default).
  fbas <- sprintf('%s_%d-%d',varcli,anyi,anyf) #base file name 
  fiche <- sprintf('%s.est',fbas) #stations file name
  #read coordinates, codes and names of the stations:
  est.c <- read.table(fiche,colClasses=c("numeric","numeric","numeric","character","character"))
  names(est.c) <- c('X','Y','Z','Code','Name')
  ne <- nrow(est.c) #no. of stations
  #check for any missing X, Y, Code:
  xnas <- sum(is.na(est.c[,1])); ynas <- sum(is.na(est.c[,2]))
  cnas <- sum(is.na(est.c[,4]))
  if(xnas>0 | ynas>0 | cnas>0) stop('There are missing coordinates and/or station codes!\nPlease, complete the *.est file and run this function again.')
  #set any missing elevation to 99 m:
  knas <- is.na(est.c[,3])
  if(sum(knas>0)) {
    cat('Warning: Missing elevations are replaced by 99\n')
    est.c[knas,3] <- 99
  }
  #change any hyphen in codes to underscores:
  z <- gsub('-','_',est.c[,4])
  zn <- sum(z != est.c[,4])
  if(zn>0) {
    cat('In',zn,'codes containing the "-" character, they have been changed to "_"\n\n')
    est.c[,4] <- z
  }
  #check for duplicate codes:
  z <- duplicated(est.c[,4])
  if(sum(z)>0) {
    zz <- unique(est.c[z,4])
    cat('Duplicated codes detected:\n')
    print(est.c[est.c[,4]%in%zz,4:5])
    stop('The station file *.est must contain unique codes.')
  }
  #use station codes for any missing station names:
  knas<- is.na(est.c[,5])
  if(sum(knas>0)) {
    cat('Warning: Missing station names are replaced by their codes\n')
    est.c[knas,5] <- est.c[knas,4]
  }
  #check if coordinates are in degrees:
  if(max(abs(est.c[,1]))>180 | max(abs(est.c[,2]))>90) 
    stop('Station coordinates must be in geographic degrees with decimals:\n-180 to 180 for X (longitude) and -90 to 90 for Y (latitude)')
  fichd <- sprintf('%s.dat',fbas) #data file name
  dat <- scan(fichd,na.strings=na.strings) #read data
  numdat <- length(dat) #no. of read data
  nd <- numdat/ne #no. of data per station
  if(nd-floor(nd)>1e-16) {
    cat(ne,"stations read from",fiche,"\n")
    cat(numdat,"data read from",fichd,"\n")
    stop("The length of data is not a multiple of the number of stations!")
  } else cat(sprintf('Data matrix: %d data x %d stations\n',nd,ne))
  dim(dat) <- c(nd,ne) #convert vector to matrix
  na <- anyf-anyi+1 #no. of years
  if(!is.null(x)) nm <- 0
  else { #calculate nm (no. of data per year and station):
    z <- nd/na
    if(z>=1) nm <- floor(z)
    if(nm > 366) nm <- -1 #flag of subdaily data
    else if(nm > 12) nm <- 0 #flag of daily data
  }
  #check if years are complete:
  ndp <- length(seq(as.Date(sprintf('%s-01-01',anyi)),
    as.Date(sprintf('%s-12-31',anyf)),1)) #no. of days in the period
  if(nm>0 & nd%%nm==0) acomp <- TRUE
  else if(nd%%ndp==0) acomp <- TRUE
  else acomp <- FALSE
  #set tinc?
  if(is.na(tinc) & nd>ndp) {
    nh <- 24*ndp/nd #no. of hours
    if(nh>=1 & 24%%nh==0) tinc <- sprintf('%d hours',nh)
    else if((nh*60)==floor(nh*60)) tinc <- sprintf('%d mins',nh*60)
  }
  #- build time vector x if not supplied by the user:
  if(is.null(x)) {
    if(is.na(ini)) ini <- sprintf('%d-01-01',anyi) #default initial date
    if(nm>0) x <- seq(as.Date(ini),length.out=nd,by=sprintf('%d months',round(12/nm)))
    else if(nm==0) x <- seq(as.Date(ini),length.out=nd,by='1 day')
    else x <- seq(as.POSIXct(ini,tz=tz),length.out=nd,by=tinc)
  } else ini <- x[1]
  if(length(x) != nd) stop(sprintf('The provided vector x has %d dates, but there appear to be %d data per station. Please, fix this inconsistency.',length(x),nd))
  return(list(est.c=est.c,dat=dat,na=na,nd=nd,ne=ne,nm=nm,x=x,ini=ini,
    tinc=tinc,acomp=acomp))
}

#- unsufix.- Remove numeric sufixes from station codes.
unsufix <- function(str) sapply(strsplit(str,'-'), function(x) x[1])

#- snht.- Maximum Standard Normal Homogeneity Test (allowing missing data).
snht <- function(y, mints=3, allT=FALSE) {
#mints: minimum tail size (minimum number of terms in the tails of the series)
#allT: return all T values? (By default, only the maximum and its position)
  yav <- which(!is.na(y)) #available data
  x <- y[yav] #series without missing data
  n <- length(x)
  if(n<mints*2) return(c(0,0)) #insuficientes datos
  if(sd(x)==0) return(c(0,0)) #serie constante
  T <- rep(NA,n)
# z <- scale(x) #a bit slower than:
  z <- (x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)
  for(i in mints:(n-mints)+1) { #(skip tails with less than mints terms)
    if(is.na(x[i])) next
    n1 <- sum(!is.na(x[1:(i-1)])) #no. of terms of sample 1
    n2 <- sum(!is.na(x[i:n])) #no. of terms of sample 2
    if(n1<mints | n2<mints) next #at least one sample is too small
    z1 <- mean(z[1:(i-1)],na.rm=TRUE)
    z2 <- mean(z[i:n],na.rm=TRUE)
    T[i] <- n1*z1*z1 + n2*z2*z2
  }
  if(allT) return(T) else return(c(max(T,na.rm=TRUE),yav[which.max(T)]))
}

#- wtest.- Inhomogeneity test on overlapping 2*nt term windows.
wtest <- function(x, nt=48, sts=3, test) {
  ntt <- length(x) #no. of total terms of the series
  ntv <- sum(!is.na(x)) #no. of valid (non-missing) data of the series
  if(2*nt>ntv) return(c(0,0)) #not enough valid data for the test
  tV <- 0 #initialization of maximum tV (test value) to return
  pk <- 0 #initialization of the tV location to return
  #initialization of sample limits (a1-b1, a2-b2):
  k <- 1; while(k<ntt & is.na(x[k])) k <- k+1; a1 <- k
  n<-1; while(n<nt & k<ntt) { k <- k+1; if(!is.na(x[k])) n <- n+1; }
  b1 <- k
  k <- k+1; while(k<ntt & is.na(x[k])) k <- k+1; a2 <- k
  n<-1; while(n<nt & k<ntt) { k <- k+1; if(!is.na(x[k])) n <- n+1; }
  b2 <- k
  #apply the test to overlapping windows:
  repeat {
    st <- eval(call(test,x[a1:b2],sts))
    stx <- st[1]
    if(!is.na(tV) & stx>tV) { tV <- stx; pk <- round(st[2])+a1-1 }
    if(b2==ntt) return(c(tV,pk))
    #shift windows forwards:
    a1 <- a2; b1 <- b2
    k <- b2+1; while(k<ntt & is.na(x[k])) k <- k+1
    if(is.na(x[k])) return(c(tV,pk)) else a2 <- k
    n<-1; while(n<nt & k<ntt) { k <- k+1; if(!is.na(x[k])) n <- n+1; }
    b2 <- k
  }
}

#- cuct.- Maximum Cucconi test in a series (allowing missing data and zeros).
cuct <- function(y, mints=3) {
#mints: minimum tail size ((minimum number of terms in the tails of the series)
  yav <- which(!is.na(y)) #available data
  x <- y[yav] #series without missing data
  n <- length(x) #series length
  if(n<mints*2) return(c(NA,NA)) #not enough data
  Tx <- 0; kx <- 0 #initialize maximum test value and location
  rkx <- rank(x,ties.method="first") #ranks of the series
  rkz <- n+1-rkx #reverse ranks of the series
  #(use as.numeric() to avoid "exceeded integers" in long series:)
  rkx <- as.numeric(rkx)*rkx; rkz <- as.numeric(rkz)*rkz #squared ranks
  nz1 <- n+1; nz2 <- 2*n+1; nz8 <- 8*n+11 #auxiliary variables
  r <- 2*(as.numeric(n)*n-4) / (nz2*nz8) - 1
  #calculate test along the series:
  for(i in mints:(n-mints)+1) { #(skip too short tails)
    n1 <- i-1; n2 <- n-n1 #no. of terms of the samples
    S1 <- sum(rkx[i:n]); S2 <- sum(rkz[i:n]) #summations
    den <- sqrt(n1*n2*nz1*nz2*nz8/5) #denominator
    U <- (6*S1-n2*nz1*nz2)/den
    V <- (6*S2-n2*nz1*nz2)/den
    T <- (U*U+V*V-2*r*U*V) / (2*(1-r*r)) #value of the test
    if(T>Tx) { Tx <- T; kx <- i }
  }
  return(c(Tx,yav[kx]))
}

