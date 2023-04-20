#clitools.R.- Graphic tools for the climatol package.
#Author: Jose A. Guijarro. Licence: GPL >= 3.0


#- dens2Dplot.- Two dimensional density plot.
#Adapted from Elzizi's answer at http://stackoverflow.com/questions/
#  18089752/r-generate-2d-histogram-from-raw-data
dens2Dplot <- function(x, y, nbins=100, pal=NULL, xlab='', ylab='',
    xlim=c(NA,NA), ylim=c(NA,NA), ...) {
#x, y: Variables for the scatterplot
#nbins: Number of bins in X and Y coordinates of the scatterplot.
#pal: Color palette
#xlab, ylab: Labels for X and Y axis
#xlim, ylim: Limits for X and Y axis
#...: Other graphic parameters
  if(is.null(pal)) pal=rev(rainbow(16,start=0,end=.65))
  xmin <- floor(min(c(x,xlim[1]),na.rm=TRUE))
  ymin <- floor(min(c(y,ylim[1]),na.rm=TRUE))
  xmax <- ceiling(max(c(x,xlim[2]),na.rm=TRUE))
  ymax <- ceiling(max(c(y,ylim[2]),na.rm=TRUE))
  xbin <- seq(xmin,xmax,length=nbins)
  ybin <- seq(ymin,ymax,length=nbins)
  freq <- as.data.frame(table(findInterval(x,xbin),findInterval(y,ybin)))
  freq[,1] <- as.integer(as.character(freq[,1]))
  freq[,2] <- as.integer(as.character(freq[,2]))
  freq2D <- matrix(0,nbins,nbins)
  freq2D[cbind(freq[,1], freq[,2])] <- freq[,3]
  freq2D[freq2D==0] <- NA
  image(xbin,ybin,freq2D,col=pal,xlab=xlab,ylab=ylab,...)
}

#- IDFcurves.- Obtain Intensity-Duration-Frequency curves.
IDFcurves <- function(prdat, stname, clmn=1:2, tz='utc', na.code=NA,
prunits='mm', mindpy=0.8, gumbel=TRUE, timeaggr=c(10,20,30,60,120,180,360,720),
retper=c(5,10,20,30,50,75,100),...) {
#prdat: Data frame with Time (as POSIXct) and sub-hourly precipitation data.
#stname: Station name.
#clmn: Columns where Time and precipitation data are located in prdat.
#tz: Time zone ['utc' by default].
#na.code: Numeric missing data code.
#prunits: Precipitation units [mm].
#mindpy: Minimum available data proportion to process data in any year.
#gumbel: Adjust a Gumbel distribution? [TRUE].
#timeaggr: Time intervals (in minutes) on which to aggregate precipitation.
#retper: Return periods (in years) for extreme precipitation estimation.
#...: Additional graphic parameters.
  if(!requireNamespace("evd", quietly=TRUE))
    stop('Please, install the package "evd" and run this function again')
  z <- eval(substitute(alist(...)))
  if(!is.null(z$lty)) lty <- z$lty else lty <- 1:5
  if(!is.null(z$col)) col <- z$col else col <- 1:5
  tdif <- diff(prdat[,clmn[1]])
  tint <- quantile(tdif,probs=c(0,.05,.5,.95,1))
  if(length(unique(tint))>1) { #check if time intervals are constant:
    cat('Time intervals between observations are not constant:\n')
    print(tint)
    cat('\nPlease, provide a continuous series of sub-hourly precipitation\n')
    cat('data. (Observations with missing data are allowed.)\n')
    stop()
  } else tint <- tint[3] #constant time interval
  #check if data are sub-hourly:
  tunit <- attr(tint,'units')
  if(tunit=='secs') tint <- tint/60 #tint in minutes
  else if(tunit!='mins') stop(sprintf('Time interval between observations is %s %s.\nPlease, provide sub-hourly data.',tint,tunit))
  year <- strftime(prdat[,clmn[1]],'%Y',tz=tz) #year of the observations
  pr <- prdat[,clmn[2]] #precipitation data
  if(!is.na(na.code)) pr[pr==na.code] <- NA #apply missing data code
  ndy <- tapply(!is.na(pr),year,sum) #no. of non-missing data per year
  mdy <- median(table(year)) #median no. of data per year (missing included)
  ny <- length(ndy) #no. of years
  #maximum precipitation in each aggregation time:
  ndaggr <- timeaggr/as.numeric(tint) #no. of data in each aggregation time
  nt <- length(timeaggr) #no. of aggregation times
  prmax <- matrix(NA,ny,nt)
  for(i in 1:nt) {
    z <- filter(pr,rep(1,ndaggr[i]))
    prmax[,i] <- tapply(z,year,max,na.rm=TRUE)
  }
  #disregard annual maximums of years with less than mindpy data proportion:
  dis <- ndy < mdy*mindpy
  ny <- ny-sum(dis) #no. of years with data
  period <- paste(range(as.integer(names(dis)[dis==FALSE])),collapse='-')
  prmax <- prmax[!dis,] #delete years without enough data
  #calculate maximum expected precipitation for every return period:
  prob <- 1/retper; nrp <- length(retper)
  pmax=matrix(NA,nt,nrp)
  if(gumbel) {
    for(j in 1:nt) {
      aj <- evd::fgev(prmax[,j],shape=0); z <- aj$estimate
      pmax[j,] <- round(evd::qgumbel(prob,z[1],z[2],lower.tail=FALSE),1)
    }
  } else {
    for(j in 1:nt) {
      aj <- evd::fgev(prmax[,j]); z <- aj$estimate
      pmax[j,] <- round(evd::qgev(prob,z[1],z[2],z[3],lower.tail=FALSE),1)
    }
  }
  #plot IDF curves:
  tagghours <- timeaggr/60 #aggregation times in hours
  pxh=scale(t(pmax),center=FALSE,scale=tagghours) #pmax in mm/h
  matplot(tagghours,t(pxh),type='l',lwd=2,xlab='Hours',las=1,
    ylab=sprintf('Intensity (%s/h)',prunits),...)
  grid(col=grey(.4))
  legend(9,120,lwd=2,lty=lty,col=col,legend=retper,title='Years',bg='white')
  title(sprintf('IDF at %s (%s)',stname, period))
  #return maximum precipitation accumulations:
  rownames(pmax) <- timeaggr
  colnames(pmax) <- retper
  return(invisible(as.data.frame(pmax)))
}

#- MHisopleths.- Isopleths on a months-hours plot.
MHisopleths <- function(dat,vrb,fun='mean',xlab='Months',ylab='Hours',cex=1.2,
col4RP=c('cyan','yellow','red'),title='') {
#dat: dataframe with POSIX times and data in columns
#vrb: name of the column containing the chosen data
#fun: function to aggregate subhourly data into hourly
#xlab, ylab: labels for the X and Y axis
#cex: character expansion parameter for the size of labels
#col4RP: vector of colors for colorRampPalette()
#title: main title
  if(!requireNamespace('fields', quietly=TRUE)) stop('This function requires the package fields.\nPlease, install it and re-run the function')
  old.par=par(no.readonly=TRUE)
  par(cex=cex,cex.lab=cex)
  Dtime <- dat[,1]
  dat <- dat[,which(names(dat)==vrb)]
  mm=as.integer(strftime(Dtime,'%m'))
  hh=as.integer(strftime(Dtime,'%H'))
  df=aggregate(dat,list(mm,hh),fun,na.rm=T)
  df[,1] <- df[,1]-0.5; df[,2] <- df[,2]+0.5
  z <- df[,3]; dim(z) <- c(12,24) #data matrix
  #interpolate matrix margins at 0.5 resolution:
  df.list <- list(x=.5:11.5,y=.5:23.5,z=z)
  grid.list <- list(x=seq(0,12,.5),y=seq(0,24,.5))
  dat.list <- fields::interp.surface.grid(df.list, grid.list)
  z=dat.list[[3]] #dim: 25 49
  z[,1] <- z[,49] <- (z[,2]+z[,48])/2
  z[1,] <- z[25,] <- (z[2,]+z[24,])/2
  #plot the diagram:
  x=seq(0,12,.5); y=seq(0,24,.5)
  levels=pretty(z,10); nl=length(levels)
  col <- colorRampPalette(col4RP)(nl+1)
  plot(NA,xlim=c(0,12),ylim=c(0,24),xaxs='i',yaxs='i',axes=FALSE,frame=TRUE,
    xlab='Months',ylab='Hours')
  .filled.contour(x=x, y=y, z=z, levels=levels,col=col)
  contour(x,y,z,add=TRUE,levels=levels,labcex=1)
  for(i in 0:12) segments(i,0,i,-1.5)
  axis(1,0:12,labels=FALSE); mtext(1:12,1,0.5,at=.5:11.5,cex=cex)
  axis(2,seq(0,24,2),las=1)
  grid(ny=6,col=grey(.4))
  if(title!='') title(title)
  par(old.par)
}

#- diagwl.- Walter & Lieth climatic diagram.
diagwl <- function(dat, cols=1:6, format='%Y-%m-%d', yeari=NA, yearf=NA,
stname="", alt=NA, per="", mlab="", shem=FALSE, p3line=FALSE, ...) {
#dat: Data frame with the required climatic data (see details).
#cols: Columns containing the needed data [1:6].
#format: Format of the dates if data are provided in 4 columns ['%Y-%m-%d'].)
#yeari, yearf: If dat is a file name, initial and final years of the period to
#  use (defaults to the period contained in the data file).
#stname: Name of the climatological station.
#alt: Elevation (altitude) of the climatological station.
#per: If data is a data frame with already calculated climate averages, the
#  original period of the data. 
#mlab: Vector of 12 monthly labels for the X axis (see the details).
#shem: Southern hemisphere? (FALSE by default.)
#p3line: Draw a supplementary precipitation line referenced to three times the
# temperature? (FALSE by default; this parameter was suggested by Bogdan Rosca)
#...: Other optional graphic parameters.
  old.par <- par(no.readonly = TRUE)
  on.exit(par(old.par))
  par(mar=c(4,4,5,4), las=1, new=FALSE)
  pcol="#005ac8"; tcol="#e81800"; pfcol="#79e6e8"; sfcol="#09a0d1" #used colors
  #etiquetas de los meses
  if(length(mlab)!=12) {
    if(mlab=="es") mlab=c("E","F","M","A","M","J","J","A","S","O","N","D")
    else if(mlab=="en") mlab=c("J","F","M","A","M","J","J","A","S","O","N","D")
    else mlab=c(1:12) #numeric labels
  }
  if(is.null(cols)) { #read file and calculate climate summary
    if(ncol(dat)<12) stop('Data frame has less than 12 columns!')
    if(ncol(dat)==13) dat <- dat[,1:12]
    nr <- nrow(dat) #no. of raws of monthly data
    switch(nr,
      stop('At least two monthly rows (average precipitation and temperature)\n  must be supplied'),
      cat("Warning: When only monthly precipitation and mean temperature\n         are provided, no frost risk will be drawn.\n"),
      cat("Warning: When no absolute minimum temperatures are provided,\n         likely frost will not be drawn.\n")
    )
  } else {
    if(length(cols)==4) { #disaggregate dates into year, month, day:
      dates <- as.Date(dat[,cols[1]],format=format)
      dat <- data.frame(YY=as.integer(strftime(dates,'%Y')),
        MM=as.integer(strftime(dates,'%m')),
        DD=as.integer(strftime(dates,'%d')),dat[,cols[2:4]])
    } else dat <- dat[,cols]
    z <- range(dat[,1])
    if(is.na(yeari)) yeari <- z[1]
      else if(yeari < z[1]) yeari <- z[1]
      else if(yeari > z[1]) d <- d[d[,1] >= yeari,]
    if(is.na(yearf)) yearf <- z[2]
      else if(yearf > z[2]) yearf <- z[2]
      else if(yearf < z[1]) d <- d[d[,1] <= yearf,]
    per=sprintf('%d-%d',yeari,yearf) #period
    ny <- yearf-yeari+1 #no. of years
    datcli <- matrix(NA,4,12)
    datcli[1,] <- round(aggregate(dat[,4],list(dat[,2]),sum)$x / ny , 1)
    datcli[2,] <- round(aggregate(dat[,5],list(dat[,2]),mean)$x , 1)
    datcli[3,] <- round(aggregate(dat[,6],list(dat[,2]),mean)$x , 1)
    datcli[4,] <- round(aggregate(dat[,6],list(dat[,2]),min)$x , 1)
    dat <- datcli; nr <- nrow(dat)
  }
  dat <- as.matrix(dat)
  if(shem) { #Southern hemisphere: shift data six months forward
    m1 <- dat[,1:6]
    m2 <- dat[,7:12]
    dat <- cbind(m2,m1)
    mlab <- c(mlab[7:12],mlab[1:6])
  }
  p <- dat[1,] #monthly average precipitations
  if(nr==2) tm <- dat[2,]
  else tm <- apply(dat[2:3,],2,mean)  #monthly average temperatures
  pmax <- max(p) #maximum precipitation
  ymax <- 60  #default maximum Y-axis ordinate
  if(pmax > 300) ymax <- 50 + 10*floor((pmax+100)/200)
  ymin <- min(-1.5,min(tm)) #minimum Y-axis ordinate
  #ejes:
  if(ymin < -1.5) {
    ymin=floor(ymin/10)*10 #rounded minimum Y-axis ordinate
    labT <- paste(ymin)
    labP <- ""
    if(ymin < -10) {
      for(i in (ymin/10+1):-1) {
        labT <- c(labT,i*10)
        labP <- c(labP,"")
      }
    }
    labT <- c(labT,"0","10","20","30","40","50","")
    labP <- c(labP,"0","20","40","60","80","100","300")
  }
  else {
    labT <- c("0","10","20","30","40","50","")
    labP <- c("0","20","40","60","80","100","300")
  }
  if(ymax > 60) {
    for(i in 6:(ymax/10-1)) {
      labT <- c(labT,"")
      labP <- c(labP,100*(2*i-7))
    }
  }
  plot(0:13-0.5,c(tm[12],tm[1:12],tm[1]),xlim=c(0,12),ylim=c(ymin,ymax),type="n",xaxs="i",yaxs="i",xaxp=c(0,12,12),xlab="",ylab="",xaxt="n",yaxt="n",bty="n")
  lmin <- ymin #minimum Y-axis ordinate to label
  if(lmin==-1.5) lmin=0
  axis(2,((lmin/10):(ymax/10))*10,labels=labT,col.axis=tcol)
  axis(4,((lmin/10):(ymax/10))*10,labels=labP,col.axis=pcol)
  mtext("C",2,col=tcol,las=1,line=3,adj=0,at=55)
  mtext("mm",4,col=pcol,las=1,line=3,adj=1,at=55)
  abline(0,0)
  abline(50,0)
  #labels:
  if(is.na(alt)) mtext(stname,line=2,adj=0)
  else mtext(paste(stname," (",alt," m)",sep=""),line=2,adj=0)
  mtext(per,line=1,adj=0)
  mtext(paste(round(mean(tm),1)," C        ",round(sum(p)),
    " mm",sep=""),line=1,adj=1)
  x <- 0:13-0.5
  p2 <- c(p[12],p[1:12],p[1])
  if(p3line) { #additional precipitation line at 1:3 scale
    yl3 <- c(p[12],p[1:12],p[1])/3
    yl3[yl3>50] <- 50
  }
  if(pmax<=100) {
    xl <- x
    yl <- c(p[12],p[1:12],p[1])/2
    n2 <- 14
  }
  else { #scale change when precipitation > 100 mm
    xp <- numeric(30)
    yp <- numeric(30)
    xl <- numeric(25)
    yl <- numeric(25)
    n <- 0
    n2 <- 0
    gr <- FALSE
    if(p2[1]>100) { #first point
      n <- n+1
      xp[n] <- x[1]
      yp[n] <- 50
      n <- n+1
      xp[n] <- x[1]
      yp[n] <- 50 + (p2[1]-100)/20
      n2 <- n2+1
      xl[n2] <- x[1]
      yl[n2] <- 50
      gr <- TRUE
    }
    else {
      n2 <- n2+1
      xl[n2] <- x[1]
      yl[n2] <- p2[1]/2
    }
    for(i in 2:14) {  #remaining points
      if(gr) {  #if previous p > 100
        n <- n+1
        if(p2[i]>100) {
          xp[n] <- x[i]
          yp[n] <- 50 + (p2[i]-100)/20
        }
        else {
          xp[n] <- x[i-1] + (100-p2[i-1])/(p2[i]-p2[i-1])
          yp[n] <- 50
          n2 <- n2+1
          xl[n2] <- xp[n]
          yl[n2] <- 50
          n <- n+1
          xp[n] <- NA
          yp[n] <- NA
          n2 <- n2+1
          xl[n2] <- x[i]
          yl[n2] <- p2[i]/2
          gr <- FALSE
        }
      }
      else {  # if previos p <=100
        if(p2[i]>100) { #if p > 100
          n <- n+1
          xp[n] <- x[i-1] + (100-p2[i-1])/(p2[i]-p2[i-1])
          yp[n] <- 50
          if(xl[n2]!=x[i-1]) {  #avoid repeating points!
            n2 <- n2+1
            xl[n2] <- x[i-1]
            yl[n2] <- p2[i-1]/2
          }
          n2 <- n2+1
          xl[n2] <- xp[n]
          yl[n2] <- 50
          n <- n+1
          xp[n] <- x[i]
          yp[n] <- 50 + (p2[i]-100)/20
          gr <- TRUE
        }
        else { # p <=100
          n2 <- n2+1
          xl[n2] <- x[i]
          yl[n2] <- p2[i]/2
        }
      }
    }
    if(!is.na(yp[n])) {  #close last poligon
      n <- n+1
      xp[n] <- xp[n-1]
      yp[n] <- 50
      n2 <- n2+1
      xl[n2] <- 12.5
      yl[n2] <- 50
    }
    polygon(xp[1:n],yp[1:n],col=pcol,border=pcol)
  }
  #patterns:
  pi <- approx(xl[1:n2],yl[1:n2],n=66)$y
  ti <- approx(x,c(tm[12],tm[1:12],tm[1]),n=66)$y
  ti[ti<0] <- 0 #avoid patterns below zero
  d <- pi - ti
  xi <- (1:66)/5-0.7
  xw <- subset(xi,d>0) #humid period
  y1 <- subset(pi,d>0)
  y2 <- subset(ti,d>0)
  if(length(xw)>0) segments(xw,y1,xw,y2,col=pcol,lty=1,lwd=1)
  xw <- subset(xi,d<0) #dry period
  y1 <- subset(pi,d<0)
  y2 <- subset(ti,d<0)
  if(length(xw)>0) segments(xw,y1,xw,y2,col=tcol,lty=3,lwd=2)
  if(nr>2) {
    #sure frosts
    for(i in 1:12) if(dat[3,i]<=0) rect(i-1,-1.5,i,0,col=sfcol)
    if(nr>3) #likely frosts
    for(i in 1:12) { if(dat[4,i]<=0 && dat[3,i]>0) rect(i-1,-1.5,i,0,col=pfcol)}
    else mtext('(Likely frost months not provided)',1,line=1.5)
  } else mtext('(No monthly frost risk provided)',1,line=1.5)
  #curvas de P y T:
  lines(xl[1:n2],yl[1:n2],col=pcol,lwd=2)
  if(p3line) lines(x,yl3)
  lines(x,c(tm[12],tm[1:12],tm[1]),col=tcol,lwd=2)
  if(nr>2) {
    #mean maximum temperatures of the warmest month
    mtext(formatC(max(as.matrix(dat[2,])),digits=1,format="f"),2,las=1,
      line=2,at=35)
    #mean minimum temperatures of the coldest month
    mtext(formatC(min(as.matrix(dat[3,])),digits=1,format="f"),2,las=1,
      line=2,at=15)
  }
  #tick month limits:
  for(i in 0:13) segments(i,0,i,-1.5)
  #label months:
  mtext(mlab,1,las=1,line=0.5,adj=0.5,at=x[2:13])
  #reset old.par (reset former graphic parameters):
  invisible()
}

#- meteogram.R.- Daily meteogram of 8 meteorological variables.
meteogram <- function(df, code='', name='', cols=1:9, tz='utc', hlab='Hours',
datefm='%Y-%m-%d', vlab=c('Wind direction (deg)','Wind speed (m/s)',NA,NA,
'Temperature (C)','Rel. humidity (%)','Precip. (mm)','Pressure (hPa)'),
vcol=c(hsv(.1,1,.9),hsv(.1,1,.9),2,2,2,hsv(.4,1,.7),4,'brown'),
llim=c(0,0,NA,NA,0,0,0,NA), ulim=c(360,20,NA,NA,20,100,4,NA)) {
#df: Data frame with (around) one day of data
#code: Code of the station
#name: Name of the station
#cols: Column order of the expected variables (see details)
#tz: Time zone of the supplied time vector ('utc' by default)
#hlab: Label for hours (default='Hours')
#datefm: Date format for the title of the meteogram (the default is
#  '%Y-%m-%d', the ISO 8601 date format)
#vlab: Variable lables
#vcol: Colors for every variable
#llim: Lower graphic limits (if fixed)
#ulim: Upper graphic limits (if fixed)
  old.par <- par(no.readonly=TRUE); on.exit(par(old.par)) #reset par on exit
  tv <- df[,cols[1]]; dt <- df[,cols[2:9]] #time and data
  date <- names(which.max(table(as.Date(tv)))) #most common date of data
  date <- strftime(as.Date(date),datefm)
  z1 <- as.integer(strftime(tv,'%H',tz=tz))
  z2 <- as.integer(strftime(tv,'%M',tz=tz))
  hv <- z1+z2/60 #decimal hours vector
  ht <- pretty(hv,12) #time ticks
  if(sum(!is.na(dt[,4]))==0) vx <- ulim[2]
  else vx <- max(ulim[2],max(dt[,4],na.rm=TRUE))
  if(sum(!is.na(dt[,5]))==0) { tn <- llim[5]; tx <- ulim[5] }
  else {
    tn <- floor(min(dt[,5],na.rm=TRUE)/5)*5
    tx <- max(tn+20,max(dt[,5],na.rm=TRUE))
  }
  if(sum(!is.na(dt[,7]))==0) px <- ulim[7]
  else px <- max(ulim[7],max(dt[,7],na.rm=TRUE))
  if(sum(!is.na(dt[,8]))==0) psn <- 1000
  else psn <- floor(min(dt[,8],na.rm=TRUE)/5)*5
  layout(matrix(1:4,4,1,byrow=TRUE),heights=c(4,3,3,4))
  par(mar=c(0,4,4,4),las=1,xaxt="n",cex=.8,xaxs="i")
  plot(hv,dt[,1],ylim=c(llim[1],ulim[1]),col=vcol[1],ylab=vlab[1],yaxt="n",
    cex.lab=1.1,main=sprintf('%s - %s  (%s)',code,name,date)) #10'av.wind dir.
  points(hv,dt[,3],col=vcol[3],pch='+') #maximum wind gust direction
  axis(2,0:4*90)
  abline(v=ht,lty=3,col=grey(.5))
  abline(h=0:4*90,lty=3,col=grey(.5))
  par(mar=c(0,4,0,4))
  plot(hv,dt[,2],type="l",ylim=c(0,vx),col=vcol[2],ylab=vlab[2],
    cex.lab=1.1) #10' average wind speed
  lines(hv,dt[,4],type="l",col=vcol[4]) #maximum wind gust speed
  grid(nx=NA,ny=NULL,col=grey(.5))
  abline(v=ht,lty=3,col=grey(.5))
  plot(hv,dt[,5],type="l",ylim=c(tn,tx),col=vcol[5],ylab=vlab[5],
    col.lab=vcol[5],cex.lab=1.1) #temperature
  grid(nx=NA,ny=NULL,col=grey(.5))
  abline(v=ht,lty=3,col=grey(.5))
  lines(hv,dt[,6]*(tx-tn)/100+tn,col=vcol[6]) #relative humidity
  axis(4,tn+0:4*(tx-tn)/4,0:4*25)
  mtext(vlab[6],4,3,las=0,col=vcol[6],cex=.9)
  par(mar=c(4,4,0,4))
  plot(hv,dt[,7],type="S",ylim=c(0,px),col=vcol[7],xlab="",ylab=vlab[7],
    col.lab=vcol[7],cex.lab=1.1) #precipitation
  lines(hv,dt[,7],type="h",col=vcol[7])
  grid(nx=NA,ny=NULL,col=grey(.5))
  abline(v=ht,lty=3,col=grey(.5))
  lines(hv,(dt[,8]-psn)/5,col=vcol[8]) #barometric pressure
  axis(4,0:4*px/4,labels=psn+0:4*5)
  mtext(vlab[8],4,3.1,las=0,col=vcol[8],cex=.9)
  axis(1,ht,xaxt='s') #label hours
  mtext(sprintf('%s %s',hlab,tz),1,3)
}


#- runtnd.- Running trends on time windows of different lengths.
runtnd <- function(d, anyi, minyr=10, units='Units', pernyr=10, stname=NA,
k=NULL, palneg=c('blue','white'), palpos=c('white','red'), ...) {
#d: Series of annual values (without missing data)
#anyi: Initial year of the series
#units: Units label for the legend
#minyr: Minimum no. of years to compute trends (10 by default)
#pernyr: Factor for trend units (per 10 years by default)
#stname: Station name
#k: Vector of trend intervals (automatically set by default)
#palneg: Color gradation for negative trends
#palpos: Color gradation for positive trends
#...: Optional graphic parameters
  if(sum(is.na(d))>0) stop('Missing data detected! The series must be complete')
  nd <- length(d)
  if(minyr<0) { #make a graphic of running trends of minyr years:
    minyr <- -minyr; nt <- nd-minyr+1
    tn <- pv <- rep(NA,nt); x <- 1:minyr
    for(i in 1:nt) {
      z <- coef(summary(lm(d[i:(minyr+i-1)]~x)))[2,c(1,4)]
      tn[i] <- z[1] #trend
      pv[i] <- z[2] #p-value
    }
    main <- sprintf('%d years running trends',minyr)
    if(!is.na(stname)) main <- sprintf('%s at %s',main,stname)
    z <- tn*pernyr; x <- (anyi+minyr-1):(anyi+nd-1)
    plot(x,z,type='l',las=1,col=4,
      xlab='Last year of the running window',
      ylab=sprintf('Trend (%s/%dyears)',units,pernyr),main=main,...)
    z[pv>.1] <- NA; lines(x,z,col=4,lwd=2)
    z[pv>.05] <- NA; lines(x,z,col=4,lwd=3)
    grid(col=grey(.4)); abline(h=0,col=2)
    return(invisible(data.frame(tnd=tn,pvl=pv)))
  } else { #make a graphic of running trends of different period lengths:
    xw <- minyr:nd; nc <- length(xw) #window lengths and their number
    x <- 1:nd; if(!is.na(anyi)) x <- x+anyi-1 #time vector
    xi <- x[1:nc]; xf <- x[xw] #inicial and final year vectors
    #--- compute trends and p-values:
    tnd <- pvl <- matrix(NA,nc,nc) #trend and p-val matrices
    for(i in 1:nc) {
      for(j in i:nc) {
        y <- xi[j-i+1]:xf[j]; dx <- d[y-anyi+1]
        lmc <- coef(summary(lm(dx~y)))[2,c(1,4)]
        tnd[j,i] <- lmc[1] #trend
        pvl[j,i] <- lmc[2] #p-value
      }
    }
    tnd <- tnd*pernyr #adjust trend units per pernyr years
    if(is.null(k)) k <- pretty(tnd) #set bin intervals
    nk <- length(k); k0 <- which(k==0)
    np <- 2*sum(k>0); nn <- 2*sum(k<0) #no. of positive and negative intervals
    nx <- max(np,nn) #maximum no. of intervals positive or negative
    col <- character()
    if(nn>0) col <- c(col,colorRampPalette(palneg)(nx+1)[(nx-nn+1):nx])
    if(np>0) col <- c(col,colorRampPalette(palpos)(nx+1)[2:(np+1)])
    fields::image.plot(xf,xw,tnd,xlab='Last year',ylab='Window length (years)',
      col=col,las=1,zlim=range(k))
    grid(col=grey(.4))
    lxw <- length(xw)
    mtext(sprintf('%s/%dyr',units,pernyr),4,1,las=1,at=xw[lxw])
    if(lxw<=101) { #add white points to non-significant values:
      cex <- 1.3-lxw/100
      zp=pvl; zp[!is.na(zp)]=0
      zp[pvl>=.05] <- cex*.6; zp[pvl>=.1] <- cex
      points(row(zp)+anyi+minyr-2,col(zp)+minyr-1,col='white',pch=19,cex=zp)
      legend("topleft",pch=15,col='grey50',bg='white',pt.cex=2,
        legend=c(as.expression(bquote(~'0.05 <='~alpha~'< 0.10')),
        as.expression(bquote(~'0.10 <='~alpha))))
      legend("topleft",pch=19,col='white',pt.cex=c(.6,1),
        legend=c('',''),bty='n')
      if(is.na(stname)) title('Running trends')
      else title(sprintf('Running trends at %s',stname))
    }
    return(invisible(list(tnd=tnd,pvl=pvl)))
  }
}

#- windrose.- Plot a windrose from a table with columns 'DateTime, Dir, Speed'.
windrose <- function(dat, cols=1:3, code='', name='', uni='m/s', maxnsc=8,
fnum=4, fint=5, flab=2, pal=c('cyan','yellow','orange','red','brown'),
ang=-3*pi/16, margin=c(0,0,4,0), ...) {
#dat: Data frame with columns DateTime, Wind direction and Wind speed.
#cols: Columns containing DateTime, Wind direction and Wind speed.
#code: Station code.
#name: Station name.
#uni: Speed units for the legend header.
#maxnsc: Maximum number of wind speed classes.
#fnum: Number of reference circles to plot.
#fint: Frequency interval (in %) between reference circles.
#flab: Parameter indicating which circles must be labelled (1: only outer circle; 2: all circles, the default; any other value will not lable any circle).
#pal: Color gradation to fill the frequency polygons.
#ang: Angle along which circles will be labelled.
#margin: Margins vector for the plot (to be passed to \code{par}).
#...: Other graphic parameters.
  #----------- compute frequencies by directions ---------------
  z <- range(dat[,cols[1]]); startdate <- z[1]; enddate <- z[2]
  dd <- dat[,cols[2]]; vv <- dat[,cols[3]]; rm(dat)
  dd[dd<0 | dd>360] <- NA #skip bad or variable direction values
  vv[is.na(dd)] <- NA #skip speeds with missing direction
  nd <- sum(!is.na(vv)) #nr. of available data
  if(nd==0) {
    cat('No wind data available in columns',cols,
    'of the provided data frame!\n'); stop()
  }
  cal <- sum(vv==0,na.rm=TRUE) #nr. of calm wind observations
  dd <- round(dd/360*16+1); dd[dd==17] <- 1 #direction classes
  vm <- tapply(vv,as.factor(dd),mean,na.rm=TRUE) #mean speed by direction
  vx <- tapply(vv,as.factor(dd),max,na.rm=TRUE) #max. speed by direction
  vmt <- mean(vv,na.rm=TRUE) #overall mean speed
  vmx <- max(vv,na.rm=TRUE) #overall max. speed
  vv[vv==0] <- NA #remove calm observations
  dd[is.na(vv)] <- NA #avoid directions with missing speed
  vvc <- pretty(vv) #speed classes
  nvvc <- length(vvc)-1 #number of speed classes
  #convert wind speeds to factor classes:
  if(nvvc>maxnsc) { #if too many speed classes, limit their number:
    nvvc <- maxnsc
    vvf <- cut(vv,c(vvc[1:maxnsc],999))
    spclasses <- paste(vvc[1:nvvc],vvc[2:(nvvc+1)],sep='-')
    spclasses[nvvc] <- paste('>=',vvc[maxnsc])
    geflag <- TRUE
  } else {
    vvf <- cut(vv,vvc)
    spclasses <- paste(vvc[1:nvvc],vvc[2:(nvvc+1)],sep='-')
    geflag <- FALSE
  }
  #compute the frequency table:
  fr <- table(vvf,dd)
  #if there are void direction classes, fill them with zeroes:
  z <- as.integer(colnames(fr))
  if(length(z)<16) {
    zrn <- rownames(fr)
    zm <- matrix(rep(0,16*nrow(fr)),dim(fr))
    zm[,z] <- fr
    fr <- zm
    colnames(fr) <- as.character(1:16) 
    rownames(fr) <- zrn
  }
  #distribute calms in the first speed class:
  nd1 <- sum(fr[1,]); fr[1,] <- fr[1,] * (nd1+cal)/nd1
  if(nd>0) fr <- fr*100./nd #compute frequencies as percentages
  frt <- apply(fr,1,sum) #sum frequencies by speed classes
  frtd <- apply(fr,2,sum) #sum frequencies by direction classes
  #build the frequency table:
  tab <- cbind(fr,frt)
  tab <- round(rbind(tab,c(frtd,sum(frtd))),1)
  tab <- data.frame(rbind(tab,round(c(vm,vmt),1),round(c(vx,vmx),1)))
  names(tab) <- c('N','NNE','NE','ENE','E','ESE','SE','SSE','S','SSW','SW','WSW','W','WNW','NW','NNW','Total')
  row.names(tab) <- c(spclasses,'Total','Mean Sp.','Mx.M.Sp.')
  #----------- plot the windrose -------------------------------
  old.par <- par(no.readonly=TRUE)
  on.exit(par(old.par))
  fr <- tab[1:nvvc,1:16]
  if(geflag) row.names(fr)[nvvc] <- paste('>=',vvc[nvvc])
  ndir <- 16 #nr. of direction classes
  nr <- nvvc #nr. of speed classes
  fmax <- fnum*fint #max. frequency to be circled
  key <- (nr>1) #legend if more than one speed class
  #make room for the legend at the left side:
  if(key) mlf <- 3 else mlf <- 1  #left margin factor
  par(mar=margin, new=FALSE)    #windrose margin
  # x,y components for every wind direction and plot settings:
  fx <- cos(pi/2-(2*pi/ndir*0:(ndir-1)))
  fy <- sin(pi/2-(2*pi/ndir*0:(ndir-1)))
  plot(fx,fy,xlim=c(-fmax-mlf*fint,fmax+fint),ylim=c(-fmax-fint,fmax+fint),
    xaxt="n",yaxt="n",xlab="",ylab="",bty="n",asp=1,type="n")
  if(nr==1) {  #only one speed class
    cx <- fx*fr
    cy <- fy*fr
  }
  else {  #more than one speed classes
    f <- apply(fr,2,sum)
    cx <- fx*f
    cy <- fy*f
    for(i in nr:2) {
      f <- f-fr[i,]
      cx <- c(cx,NA,fx*f)
      cy <- c(cy,NA,fy*f)
    }
  }
  col <- colorRampPalette(pal)(maxnsc)
  polygon(cx,cy,col=col[nr:1])
  symbols(c(0*1:fnum),c(0*1:fnum),circles=c(fint*1:fnum),inches=FALSE,add=TRUE)
  segments(0*1:ndir,0*1:ndir,fmax*fx,fmax*fy)
  fmaxi <- fmax+fint/4
  text(0,fmaxi,"N")
  text(0,-fmaxi,"S")
  text(fmaxi,0,"E")
  text(-fmaxi,0,"W")
  if(flab==2)
    for(i in 1:fnum) text(i*fint*cos(ang),i*fint*sin(ang),paste(i*fint,"%"))
  else if(flab==1)
    text(fmax*cos(ang),fmax*sin(ang),paste(fmax,"%"))
  if(key) { #legend
    legend(-fmaxi-2.3*fint,fmaxi+2,fill=col,legend=spclasses)
    text(-fmaxi-1.5*fint,fmaxi+.9*fint,uni)
  }
  if(code=='') title(sprintf('%s windrose\n%d obs. from %s to %s',name,nd,startdate,enddate))
  else title(sprintf('%s-%s windrose\n%d obs. from %s to %s',code,name,nd,startdate,enddate))
  return(invisible(tab))
}

