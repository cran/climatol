dat.z <- dat.e <- dat.c <- oneref <- anom <- sanom <- mindist <- NA 
used <- dat.d <- dat.na <- dat.m <- dat.s <- datmed <- refmed <- NA
refstd <- dat.m0 <- dat.s0 <- outan <- dat <- est.c <- ne <- NA
nd <- na <- anyi <- anyf <- est.d <- est.w <- est.p <- NA
dah <- dahs <- iest <- refhom <- NA
verde <- hsv(.33,1,.6)

depudm <- function(varcli, anyi, anyf, nm=12, a=0, b=1,  wz=.001, deg=FALSE,
  std=3, nref=10, wa=10000, dz.max=4, mxdif=.005, xf=0, nmd=5, rtrans=0, ndec=1,
  graf=2, aref=FALSE, maxite=50, vmin=NA, vmax=NA) {

  nmd <- max(5,nmd)

  if(class(xf)=="Date") fechas <- TRUE else fechas <- FALSE

  matpesos(wa=wa, wz=wz, deg=deg)

  dat.z <<- matrix(NA,nd,ne) 
  dat.e <<- matrix(NA,nd,ne) 
  dat.c <<- matrix(NA,nd,ne) 
  oneref <<- matrix(FALSE,nd,ne) 
  anom <<- matrix(NA,nd,ne) 
  sanom <<- matrix(NA,nd,ne) 
  mindist <<- matrix(NA,nd,ne) 
  used <<- matrix(FALSE,ne,ne) 

  repeat {
    dat.d <<- dat

    numdat <- apply(!is.na(dat.d),1,sum)
    nmin=min(numdat)
    if(!nmin) {
      cat("\nThere are terms with NO DATA!:\n")
      for(j in which(numdat==0)) {
        if(nm>1) {
          me <- j%%nm; if(me==0) me <- 12
          cat(anyi+j%/%nm,me,"\n")
        }
        else if(fechas) cat(format(xf[j]),"\n")
        else cat(j,"\n")
      }
      stop("Cannot continue! Shorten the study period, add series with data in the empty terms, or loose the outlier rejection.")
    }

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

      estelim(which(numdat<nmd),std=std,wa=wa,wz=wz,deg=deg)
      next 
    }
    dat.na <<- is.na(dat.d) 

    if(exists("dat.m0")) { 
      dat.m <<- dat.m0
      if(std==3) dat.s <<- dat.s0
    }
    else { 
      datmed <<- apply(dat.d,1,mean,na.rm=TRUE) 
      refmed <<- mean(datmed) 
      dat.m <<- apply(dat.d,2,mean,na.rm=TRUE) 
      if(std==3) {
        refstd <<- sd(datmed) 
        dat.s <<- apply(dat.d,2,sd,na.rm=TRUE) 
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
      dat.m0 <<- dat.m 
    }

    stdrz(std=std,mscalc=FALSE)

    ite <- 0
    cat("Iterative computation of missing data, with (optional) outlier removal\n")
    cat('(Max. difference <=',mxdif,'or until',maxite,'iterations):\n')
    cat("Max. data change (station)\n")
    repeat { 
      ite <- ite+1

      datest(std=std,nr=nref,aref=FALSE)

      anom <<- dat.z-dat.e 
      anom[dat.na] <<- NA  

      anomm <- apply(anom,2,mean,na.rm=TRUE) 
      anoms <- apply(anom,2,sd,na.rm=TRUE) 
      for(j in 1:ne) sanom[,j] <<- (anom[,j]-anomm[j])/anoms[j] 
      elim <- abs(sanom)>dz.max & !refhom 
      elim[is.na(elim)] <- FALSE 
      nelim <- sum(elim) 
      if(nelim>0) { 

        for(i in 1:ne) {
          for(j in 1:nd) if(elim[j,i] & !is.na(oneref[j,i])) {
            outan[j,iest[i]] <<- sanom[j,i] 
            do <- dat.d[j,i] 
            dc <- dat.c[j,i] 
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
            if(oneref[j,i]) { 
              cat(" Only 1 reference! (Unchanged)")
              elim[j,i] <- FALSE
              oneref[j,i] <<- NA 
            }
            cat("\n")
          }
        }
        dat[elim] <<- NA 
        dat.na[elim] <<- TRUE 
      }
      dat.d[dat.na] <<- dat.c[dat.na] 
      if(ite>1) {
        maxddif <- max(abs(dat.d-dat.d0),na.rm=TRUE) 
        kmaxdif <- ceiling(which.max(abs(dat.d-dat.d0))/nd) 
      }
      dat.d0 <- dat.d 
      stdrz(std=std) 
      if(aref==FALSE) break 
      if(ite>1) {
        cat(round(maxddif,ndec+2)," (",est.c[kmaxdif,4],")\n",sep="")
        if(maxddif<=mxdif || ite==maxite) {
          if(ite==maxite) cat("\nAverage computation skipped after",ite,"iterations\n")
          else cat("\n")
          break
        }
      }
    }
    dat.m0 <<- dat.m 
    if(std==3) dat.s0 <<- dat.s 

    oneref[is.na(oneref)] <<- TRUE

    if(aref==TRUE) {
      datest(std=std,nr=nref,aref=TRUE)
      dat.d[dat.na] <<- dat.c[dat.na] 
      stdrz(std=std) 
    }

    arch <- paste(varcli,"_",anyi,"-",anyf,".dah",sep="") 
    z <- dat.d
    if(rtrans>0) z[z>1] <- z[z>1]^rtrans

    if(!is.na(vmin)) z[z<vmin] <- vmin
    if(!is.na(vmax)) z[z>vmax] <- vmax
    if(nm>0) ncols <- nm else ncols <- 10
    write(round(z,ndec),arch,ncolumns=ncols)

    anom <<- dat.z-dat.e 
    anom[dat.na] <<- NA  
    anomm <- apply(anom,2,mean,na.rm=TRUE) 
    anoms <- apply(anom,2,sd,na.rm=TRUE) 
    for(j in 1:ne) sanom[,j] <<- (anom[,j]-anomm[j])/anoms[j] 
    break
  }
}

leerdm <- function(varcli,anyi,anyf,b=1,a=0,na.strings="NA") {
  arch <- paste(varcli,"_",anyi,"-",anyf,sep="") 
  arche <- paste(arch,".est",sep="") 

  est.c <<- read.table(arche,colClasses=c("numeric","numeric","numeric","character","character"))
  ne <<- nrow(est.c) 
  archd <- paste(arch,".dat",sep="") 
  dat <<- scan(archd,na.strings=na.strings) 
  dat <<- a+dat*b 
  numdat <- length(dat) 
  nd <<- numdat/ne 
  if(nd-floor(nd)>1e-16) {
    cat(ne,"stations read from",arche,"\n")
    cat(numdat,"data read from",arch,"\n")
    stop("The number of data is not multiple of the number of stations!")
  }
  dim(dat) <<- c(nd,ne) 
  anyi <<- anyi 
  anyf <<- anyf 
}

matpesos <- function(wa=100, wz=.001, deg=FALSE) {

  est.d <<- matrix(NA,ne,ne) 
  est.w <<- matrix(NA,ne,ne) 
  est.p <<- matrix(NA,ne,ne) 
  cat("Computing inter-station weights")
  if(ne>100) cat(" (",date(),")",sep="")
  else cat(":")
  cat("\n")
  for(i in 1:(ne-1)) {
    cat(" ",i)
    for(j in (i+1):ne) {
      dx <- est.c[i,1]-est.c[j,1]
      dy <- est.c[i,2]-est.c[j,2]
      if(deg) {  
        dx <- dx*111*cos((est.c[i,2]+est.c[j,2])*pi/360)
        dy <- dy*111
      }
      dz <- (est.c[i,3]-est.c[j,3])*wz
      d2 <- dx*dx+dy*dy+dz*dz 

      if(d2==0 & iest[i]!=iest[j]) est.d[i,j] <<- 0.01 
      else est.d[i,j] <<- sqrt(d2) 
      est.d[j,i] <<- est.d[i,j]  
      if(wa<=0) est.w[i,j] <<- 1 
      else est.w[i,j] <<- wa/(wa+d2) 
      est.w[j,i] <<- est.w[i,j]  
    }
  }
  cat("\nComputing proximity ranks...\n")

  for(i in 1:ne) est.p[i,] <<- order(est.d[i,])
  if(ne>100) cat("Done! (",date(),")\n",sep="")
}

stdrz <- function(std=3,mscalc=TRUE) {
  if(mscalc) { 
    dat.m <<- apply(dat.d,2,mean,na.rm=TRUE) 
    if(std==3) dat.s <<- apply(dat.d,2,sd,na.rm=TRUE) 
  }

  switch(std,
    for(i in 1:ne) dat.z[,i] <<- dat.d[,i] - dat.m[i],

    for(i in 1:ne) { if(dat.m[i]>1) dat.z[,i] <<- dat.d[,i] / dat.m[i]
      else dat.z[,i] <<- dat.d[,i] },
    for(i in 1:ne) dat.z[,i] <<- (dat.d[,i]-dat.m[i]) / dat.s[i],
    dat.z <<- dat.d
  )
}

datest <- function(std=3,nr=10,aref=FALSE) {
  for(i in 1:ne) { 
    for(j in 1:nd) { 
      se <- 0
      sw <- 0
      nref <- 0
      for(ir in 1:ne) { 
        kr <- est.p[i,ir]
        if(kr==i) next 
        if(dat.na[j,kr]) next 

        if(refhom[kr]) nref <- nref+2 else nref <- nref+1
        used[i,kr] <<- TRUE 

        if(nref==1) mindist[j,i] <<- max(est.d[i,kr],1)
        se <- se + est.w[i,kr] * dat.z[j,kr]
        sw <- sw + est.w[i,kr]

        if(nref>=nr | (aref & !est.d[i,kr])) break
      }
      if(!nref) dat.e[j,i] <<- dat.z[j,i] 
      else {

        if(nref==1 & !is.na(oneref[j,i])) oneref[j,i] <<- TRUE

        if(std==2 & se<0) se <- 0
        dat.e[j,i] <<- se / sw 
      }
    }
  }

  n <- sum(is.nan(dat.e))
  if(n>0) {
    cat(n,"NaN's in dat.e ! (changing them to NA's...)\n")
    dat.e[is.nan(dat.e)] <<- NA
  }

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
}

dahstat <- function(varcli, anyi,anyf, anyip=anyi, anyfp=anyf, nm=12, ndec=1,
  vala=2, cod=NULL, mnpd=0, mxsh=0, out="med", prob=.5, func=FALSE, long=FALSE,
  pernum=100, eshcol=4, sep=" ", eol="\n",leer=TRUE,a=0,b=1,na.strings='NA') {

  if(anyip<anyi) stop("Asked initial year before first year of data!")
  if(anyfp>anyf) stop("Asked final year beyond last year of data!")
  estvar <- c("X","Y","Z","Code","Name","PD","io","op","SNHT")
  mes3 <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
  na <- anyf-anyi+1 
  if(leer & (out=='csv' | out=='csv2')) { 
    leerdm(varcli,anyi,anyf,b=b,a=a,na.strings=na.strings)
    dim(dat) <<- c(nm,na,ne)
    rm(est.c,pos=1)
  }
  if(nm==1 | out=='csv' | out=='csv2') vala <- 0 
  else if(vala<0 || vala>4) vala <- 2 
  funa <- c("sum","mean","max","min")[vala] 

  fun <- c("mean","median","max","min","sd","quantile")[which(c("med","mdn","max","min","std","q","tnd")==out)]
  arch <- paste(varcli,"_",anyi,"-",anyf,sep="") 
  arche <- paste(arch,".esh",sep="") 
  est.c <- read.table(arche,as.is=TRUE) 
  ne <- nrow(est.c) 
  archd <- paste(arch,".dah",sep="") 
  dah <- scan(archd) 
  dim(dah) <- c(nm,na,ne) 

  if(!length(fun) & out!='csv' & out!='csv2') {
    cat("Data available in dah[",nm,",",na,",",ne,"] for the ",ne," stations in est.c\n",sep="")
    cat("(No further output, since option out='",out,"' is not recognized)\n",sep="")
    dah <<- dah; est.c <<- est.c
    return(invisible(ne))
  }

  if(!is.null(cod)) esel <- est.c[,4] %in% cod 
  else { 

    esel <- est.c[,6]>=mnpd 

    if(mxsh) esel <- esel & est.c[,9]<=mxsh 
    if(func) esel <- esel & as.logical(est.c[,8]) 
    else if(long) {
      lsel <- rep(TRUE,ne) 
      for(ko in 1:max(est.c[,7])) { 
        kest <- which(est.c[,7]==ko) 
        if(length(kest)>1) { 
          kmax <- which.max(est.c[kest,6]) 
          lsel[kest[-kmax]] <- FALSE 
        }
      }
      esel <- esel & lsel
    }
  }
  ne <- sum(esel) 
  if(ne==0) stop("No station selected: No output")
  dah <- dah[,,esel]
  if(nm==1 | ne==1) dim(dah) <- c(nm,na,ne)
  est.c <- est.c[esel,]
  if(vala) { 
    aval <- as.vector(apply(dah,2:3,funa))
    dim(dah) <- c(nm,na*ne)
    dah <- rbind(dah,aval)
    nm <- nm+1
    dim(dah) <- c(nm,na,ne)
  }

  if(out!="csv" & out!="csv2") val <- matrix(NA,ne,nm)

  x <- anyip:anyfp 
  xk <- x-anyi+1 
  if(out=="tnd") { 
    for(i in 1:ne) {
      if(nm==1) { 
        aj <- lm(dah[1,xk,i]~x) 
        val[i,] <- round(aj$coefficients[2]*pernum,ndec)
      }
      else {
        for(j in 1:nm) {
          aj <- lm(dah[j,xk,i]~x) 
          val[i,j] <-  round(aj$coefficients[2]*pernum,ndec)
        }
      }
    }
  }
  else if(out=="csv" | out=="csv2") { 

    nas <- length(x) 
    for(kest in 1:ne) { 
      kori <- est.c[kest,7] 
      dh <- dah[,xk,kest] 
      do <- dat[,xk,kori] 

      df <- abs(dh-do) < 1e-9
      df <- as.numeric(df) 
      df[df==0] <- 2 
      df[df==1] <- 0 
      df[is.na(df)] <- 1 
      dim(df) <- dim(dh)

      ard <- paste(varcli,"_",anyip,"-",anyfp,"_",est.c[kest,4],".csv",sep="")
      arf <- paste(varcli,"_",anyip,"-",anyfp,"_",est.c[kest,4],".flg",sep="")
      dh <- cbind(x,t(dh))
      df <- cbind(x,t(df))
      if(nm==12) colnames(dh) <- c('Year',mes3)
      else colnames(dh) <- c('Year',as.character(1:nm))
      colnames(df) <- colnames(dh)
      if(out=="csv") {
        write.csv(dh,ard,row.names=F)
        write.csv(df,arf,row.names=F)
      }
      else {
        write.csv2(dh,ard,row.names=F)
        write.csv2(df,arf,row.names=F)
      }
    }
    cat('Homogenized values written in *.csv, with flags in *.flg\n')
    return(invisible())
  }
  else { 
    for(i in 1:ne) {
      if(nm==1) {
        if(out=="q") val[i,]<-round(eval(call(fun,dah[,xk,i]),prob),ndec)
        else val[i,]<-round(eval(call(fun,dah[,xk,i])),ndec)
      }
      else { 
        if(out=="q") val[i,]<-round(apply(dah[,xk,i],1,fun,prob),ndec)
        else val[i,]<-round(apply(dah[,xk,i],1,fun),ndec)
      }
    }
  }

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

  ars <- paste(varcli,"_",anyip,"-",anyfp,".",out,sep="") 

  if(out=="q") ars <- paste(ars,formatC(100*prob,width=2,flag="0"),sep="")
  write.table(dahs,ars,row.names=FALSE,sep=sep,eol=eol)
  cat("\n  written to",ars,"(and remain in memory as 'dahs')\n")
}

snhtw <- function(x,nt=48) {
  ntt <- length(x) 
  ntv <- sum(!is.na(x)) 
  if(2*nt>ntv) return(c(0,0)) 
  tV <- 0 
  pk <- 0 

  k <- 1; while(k<ntt && is.na(x[k])) k <- k+1; a1 <- k
  n<-1; while(n<nt && k<ntt) { k <- k+1; if(!is.na(x[k])) n <- n+1; }
  b1 <- k
  k <- k+1; while(k<ntt && is.na(x[k])) k <- k+1; a2 <- k
  n<-1; while(n<nt && k<ntt) { k <- k+1; if(!is.na(x[k])) n <- n+1; }
  b2 <- k

  repeat {
    st <- snht(x[a1:b2])
    stx <- max(st,na.rm=TRUE)
    if(stx>tV) { tV <- stx; pk <- which.max(st)+a1-1 }
    if(b2==ntt) return(c(tV,pk))

    a1 <- a2; b1 <- b2
    k <- b2+1; while(k<ntt && is.na(x[k])) k <- k+1
    if(is.na(x[k])) return(c(tV,pk)) else a2 <- k
    n<-1; while(n<nt && k<ntt) { k <- k+1; if(!is.na(x[k])) n <- n+1; }
    b2 <- k
  }
}

snht <- function(x,nmt=3) {

  n <- length(x)
  Tsnht <- rep(NA,n)
  if(n<nmt*2) return(Tsnht) 
  z <- (x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)
  yc <- FALSE 
  for(i in (nmt+1):(n-nmt)) { 
    if(is.na(x[i]) && yc) next 
    n1 <- sum(!is.na(x[1:(i-1)])) 
    n2 <- sum(!is.na(x[i:n])) 
    if(n1<nmt | n2<nmt) next 
    z1 <- mean(z[1:(i-1)],na.rm=TRUE)
    z2 <- mean(z[i:n],na.rm=TRUE)
    Tsnht[i] <- n1*z1*z1 + n2*z2*z2
    if(is.na(x[i])) yc <- TRUE 
    else if(yc) yc <- FALSE    
  }
  return(Tsnht)
}

estelim <- function(x,std,wa=100,wz=.001,deg=FALSE) {
  est <- setdiff(1:ne,x) 
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
  if(exists("refhom")) refhom <<- refhom[est]

  ne <<- ne-length(x)

  matpesos(wa=wa, wz=wz, deg=deg)
}

plotstan <- function(i,varcli,x,xlab,swa,lw) {
  y <- sanom[,i] 
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

outrename <- function(varcli, anyi, anyf, suffix, varsuf=TRUE, restore=FALSE) {

  fbn <- paste(varcli,"_",anyi,"-",anyf,sep="") 

  if(varsuf) fbn2 <- paste(varcli,"-",suffix,"_",anyi,"-",anyf,sep="")
  else fbn2 <- paste(fbn,"-",suffix,sep="")
  for(ext in c(".txt",".dah",".esh",".pdf")) {
    if(restore) file.rename(paste(fbn2,ext,sep=""),paste(fbn,ext,sep=""))
    else file.rename(paste(fbn,ext,sep=""),paste(fbn2,ext,sep=""))
  }
}

cerrar <- function() {
  sink()
  graphics.off()
}

dd2m <- function(varcli, anyi, anyf, ini, anyip=anyi, anyfp=anyf, ndec=1,
  valm=2, nmin=15, na.strings="NA", dah=TRUE) {
  arch <- paste(varcli,"-d_",anyi,"-",anyf,sep="") 
  do <- scan(paste(arch,".dat",sep=""),na.strings=na.strings) 
  est <- read.table(paste(arch,".est",sep="")) 
  neo <- nrow(est) 
  nd <- length(do)/neo 
  dim(do) <- c(nd,neo)
  if(dah) {
    dh <- scan(paste(arch,".dah",sep="")) 
    esh <- read.table(paste(arch,".esh",sep="")) 
    ne <- nrow(esh) 
    dim(dh) <- c(nd,ne)
    iest <- esh[,7] 
  }
  else { dh <- do; ne <- neo; iest <- 1:ne }
  fech <- as.Date(0:(nd-1),origin=ini) 
  fun <- c("sum","mean","max","min")[valm] 
  na <- anyfp-anyip+1 
  dm <- array(NA,c(12,na,ne)) 
  for(ie in 1:ne) { 
    cat(' ',ie)
    for(aa in anyip:anyfp) { 
      ka <- aa-anyip+1 
      for(me in 1:12) { 
        aamm <- sprintf('%04d-%02d',aa,me)
        sel <- substring(fech,1,7)==aamm 
        if(sum(!is.na(do[sel,iest[ie]])) < nmin) next 
        dm[me,ka,ie] <- eval(call(fun,dh[sel,ie]))
      }
    }
  }
  dm[is.nan(dm)] <- NA 

  if(dah) {
    fichsal <- paste(varcli,"-m_",anyi,"-",anyf,".dah",sep="")
    fichest <- paste(varcli,"-m_",anyi,"-",anyf,".esh",sep="")
  }
  else {
    fichsal <- paste(varcli,"-m_",anyi,"-",anyf,".dat",sep="")
    fichest <- paste(varcli,"-m_",anyi,"-",anyf,".est",sep="")
  }
  write(round(dm,ndec),fichsal,ncolumns=12)
  write.table(est[iest,],fichest,row.names=FALSE,col.names=FALSE)
  cat("\nMonthly",fun,"values output to file",fichsal,"\n")
}

homogen <- function(varcli, anyi, anyf, nm=12, nref=10, dz.max=5,
  wd=c(0,0,100), swa=60, snht1=25, tVt=snht1, snht2=25, snhtt=snht2, tol=.02,
  tVf=tol, mxdif=.05, force=FALSE, a=0, b=1, wz=.001, deg=TRUE, rtrans=0,
  std=3, ndec=1, mndat=0, leer=TRUE, gp=3, na.strings="NA", nclust=100,
  maxite=50, ini="", vmin=NA, vmax=NA, verb=TRUE) {

  options(error=cerrar)
  if(deg==TRUE) require(maps)
  verde <<- hsv(.33,1,.6) 

  if(nm<1) varcli <- paste(varcli,"-d",sep="")

  archlog <- paste(varcli,"_",anyi,"-",anyf,".txt",sep="")
  sink(archlog,split=verb) 
  cat("\nHOMOGEN() APPLICATION OUTPUT  (From R's contributed package 'climatol' CR 2.2)\n\n")
  cat("\n=========== Homogenization of ",varcli,", ",anyi,"-",anyf,". (",
        date(),")\n",sep="")
  cat("homogen:")
  arg <- names(formals()) 
  for(i in 1:length(arg)) {
    cat(" ",arg[i],"=",sep="")
    cat(eval(as.symbol(arg[i])),sep=",")
  }
  cat("\n\n")

  if(mndat<=0) { 
    if(nm>5) mndat <- nm
    else if(nm<=0) mndat <- swa/2
    else mndat <- 5
  }

  k <- length(wd); if(k<3) wd <- c(wd,rep(wd[k],3-k))
  k <- length(nref); if(k<3) nref <- c(nref,rep(nref[k],3-k))
  k <- length(dz.max); if(k<3) dz.max <- c(dz.max,rep(dz.max[k],3-k))

  mes3 <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")

  if(leer) leerdm(varcli,anyi,anyf,b=b,a=a,na.strings=na.strings)
  else {
    ne <<- nrow(est.c) 
    nd <<- length(dat[,1]) 
  }
  if(!is.na(vmin)) { 
    n <- sum(dat<vmin,na.rm=TRUE)
    if(n) {
      dat[dat<vmin] <<- vmin 
      cat("Number of data forced to",vmin,":",n,"\n")
    }
  }
  if(!is.na(vmax)) { 
    n <- sum(dat>vmax,na.rm=TRUE)
    if(n) {
      dat[dat>vmax] <<- vmax 
      cat("Number of data forced to",vmax,":",n,"\n")
    }
  }
  if(gp>2) dat.o <- dat 
  if(nd<100) lw=3 
  else if(nd<300) lw=2
  else lw=1

  if(swa>=nd) swa <- floor(swa/2) 
  if(nm>0) {
    x <- (anyi*nm):(anyf*nm+nm-1)/nm 
    xlab <- "Years"
    na <<- anyf-anyi+1 
    nsy <- rep(0,na)   
    if(nd != nm*na) {
      cat(nd,"data *",ne,"stations =",nd*ne,"data, but\n")
      cat(nd,"data /",na,"years = nm =",nd/na,"data per year!\n")
      stop("Inconsistent number of data!")
    }
  }
  else {
    if(ini>0) {
      x <- as.Date(0:(nd-1),origin=ini) 
      xlab <- "Dates"
    }
    else {
      x <- 1:nd
      xlab <- "Index"
    }
  }
  nei <- ne  
  est.i <- est.c 
  nsp <- rep(0,nei)  
  iest <<- 1:ne   
  outan <<- matrix(NA,nd,ne) 

  if(gp>0) {

    pdfname <- paste(varcli,"_",anyi,"-",anyf,".pdf",sep="")
    pdf(pdfname,title=pdfname)
    plot(-1:1,-1:1,type="n",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
    text(0,0.4,"CLIMATOL",cex=4)
    text(0,-0.45,paste("Homogenization\ngraphic output of\n",varcli,"\n",anyi,"-",anyf,sep=""),cex=3)

    numdat <- apply(!is.na(dat),1,sum)
    plot(x,numdat,type="l",ylim=c(0,ne),col="blue",xlab=xlab,ylab="Nr. of data",main=paste("Nr. of",varcli,"data in all stations"))
    abline(h=5,lty=2,col="green")
    abline(h=3,lty=2,col="red")
    grid()

    if(!min(numdat)) {
      stop("At least one term has missing data in all stations! (See the PDF graph)\nCannot continue. (Shorten the study period, or add series with data in those void terms)")
    }

    if(nm>1) { 
      dim(dat) <<- c(nm,na,ne) 
      for(me in 1:nm) { 
        z <- data.frame(dat[me,,])
        names(z) <- 1:ne

        if(nm==12) labm <- mes3[me] else labm <- me
        labm <- paste(" (",labm,")",sep="")
        if(ne>nclust) hist(as.matrix(z),xlab=varcli,main=paste("Data values of ",varcli,labm,sep=""),col="wheat")
        else {
          boxplot(z,xlab="Stations",ylab="Values",main=paste("Data values of ",varcli,labm,sep=""),col="wheat",border=hsv(.7,1,.9))
          grid(col=grey(.4))
          abline(h=0)
        }
      }
      dim(dat) <<- c(nd,ne) 
    }
    else { 
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
  if(rtrans>0) { 
    dat[!is.na(dat)&dat>1] <<- dat[!is.na(dat)&dat>1]^(1/rtrans)

    z <- (100+mxdif)^(1/rtrans)-100^(1/rtrans)
    if(mxdif > z) mxdif <- z
  }
  if(gp>0) { 

    if(rtrans>0) main="Histogram of all (transformed) data"
    else main="Histogram of all data"
    hist(dat,xlab=varcli,main=main,col=hsv(.4,1,.8))

    if(ne>nclust) { splc <- sample(1:ne,nclust); nec <- nclust }
    else { splc <- 1:ne; nec <- ne }
    est.d <<- matrix(NA,nec,nec) 
    cat("Computing inter-station distances:")
    for(i in 1:(nec-1)) {
      cat(" ",i)
      for(j in (i+1):nec) {
        dx <- est.c[splc[i],1]-est.c[splc[j],1]
        dy <- est.c[splc[i],2]-est.c[splc[j],2]
        if(deg) {  
          dx <- dx*111*cos((est.c[splc[i],2]+est.c[splc[j],2])*pi/360)
          dy <- dy*111
        }
        dz <- (est.c[splc[i],3]-est.c[splc[j],3])*wz
        d2 <- dx*dx+dy*dy+dz*dz 
        est.d[i,j] <<- sqrt(d2) 
        est.d[j,i] <<- est.d[i,j]  
      }
    }
    cat("\n")
    data <- dat[,splc] 
    if(nm>1) { 
      dim(data) <- c(nm,na,nec) 
      difd <- apply(data,c(1,3),diff)
      dim(difd) <- c(nd-nm,nec) 
    }
    else difd <- diff(data) 
    corm <- cor(difd,use="p") 

    corm[corm==1] <- NA; corm[corm==-1] <- NA
    if(ne>nclust) main <- paste("Correlogram of first difference",nclust,"sampled series")
    else main <- "Correlogram of first difference series"
    if(rtrans>0) main <- paste(main," (root ",rtrans," transformed)",sep="")
    xd <- as.vector(est.d); y <- as.vector(corm)
    plot(xd,y,xlim=c(0,max(est.d,na.rm=TRUE)),xlab="Distance (km)",ylab="Correlation coefficient",main=main,col="blue")
    grid(col=gray(.4))
    if(ne>2) {  
      dism <- dist(corm) 

      if(!sum(is.na(dism))) {
        hc <- hclust(dism,method="ward")
        if(ne>nclust) main <- paste("Dendrogram of",nclust,"sampled stations")
        else main <- "Dendrogram of station clusters"
        plot(hc,xlab="Stations",sub="",ylab="Dissimilarity",main=main)

        cutlev <- mean(hc$height)+sd(hc$height)
        repeat {
          ct <- cutree(hc,h=cutlev)
          nc <- length(levels(factor(ct)))
          if(nc<10) break
          cutlev <- cutlev + .1
        }
        if(nc>1) abline(h=cutlev,col="red",lty=2)
      } else { nc <- 1; ct <- 1 } 

      if(nc==1) { col="blue"; main=paste(varcli,"station locations") }
      else {
        col=rainbow(nc,1,.55)[ct]
        main=paste(varcli," station locations (",nc," clusters)",sep="")
      }
      if(deg) asp=1/(cos(mean(est.c[,2])*pi/180)) 
      if(ne>nclust) { 

        if(deg) plot(est.c[,1:2],asp=asp,xlab="Longitude (deg)",ylab="Latitude (deg)",main=main,type='n')
        else plot(est.c[,1:2],asp=1,xlab="X (km)",ylab="Y (km)",main=main,type='n')

        points(est.c[-splc,1:2],pch='+',cex=.5)

        points(est.c[splc,1:2],col=col,pch=ct)
      }
      else { 
        if(deg) plot(est.c[,1:2],type="n",asp=1/(cos(mean(est.c[,2])*pi/180)),xlab="Longitude (deg)",ylab="Latitude (deg)",main=main)
        else plot(est.c[,1:2],type="n",asp=1,xlab="X (km)",ylab="Y (km)",main=main)
        text(est.c[,1:2],labels=1:ne,col=col)
      }
      grid(col=gray(.4))
      if(deg==TRUE) try(map('world',col=grey(.65),add=TRUE))
    }
    rm(data,difd,corm) 
  }

  if(gp==1) { 
    graphics.off() 
    cat("\nOnly the initial exploratory graphics were demanded.\nSee them in ",varcli,"_",anyi,"-",anyf,".pdf\n",sep="")
  }
  else { 
    if(exists("dat.m0")) rm(dat.m0,pos=1) 
    refhom <<- substr(est.c[,4],1,1)=='*' 
    for (k in 1:3) { 

      if(snht1==0 & k<3) next 
      if(k==2) {
        if(snhtt>0) tVt <- snhtt 
        else next 
      }
      cat("\n\n============ STAGE",k,"=======================\n\n")

      if(gp>0) {
        plot(-1:1,-1:1,type="n",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
        text(0,0.4,paste("Stage",k),cex=4)
        if(k==1) text(0,-0.3,paste("Binary splits on",swa,"term\nstepped windows\nwith tV >",tVt,"\nand wd =",wd[k],"km"),cex=3)
        else if(k==2) text(0,-0.3,paste("Binary splits on\nwhole series\nwith SNHT >",snhtt,"\nand wd =",wd[k],"km"),cex=3)
        else text(0,-0.3,paste("Anomalies after\nmissing data\nrecalculation\nwith wd =",wd[k],"km\n( swa =",swa,")"),cex=2.5)
      }
      repeat { 

        depudm(varcli, anyi, anyf, nm=nm, a=a, b=b, wz=wz, deg=deg, std=std,
          nref=nref[k], wa=wd[k]*wd[k], dz.max=dz.max[k], mxdif=mxdif, xf=x,
          nmd=mndat, rtrans=rtrans, ndec=ndec, aref=(k==3 & tVt>0),
          maxite=maxite, vmin=vmin, vmax=vmax)

        if(k>2) break

        nn <- 0 
        tVx <- rep(0,ne) 
        kpx <- rep(NA,ne) 
        splt <- rep(0,ne)  
        modif <- FALSE 
        cat("Performing shift analysis for the",ne,"stations.")
        if(k==1) cat("  tV values:\n") else cat("  SNHT values:\n")
        for(i in 1:ne) { 
          cat(" ",i,":",sep="")
          if(refhom[i]) cat('*') 
          else {
            y <- sanom[,i] 
            if(k==1) { 
              st <- snhtw(y,swa) 
              tVx[i] <- st[1]
              kpx[i] <- st[2]
              cat(round(tVx[i],1)) 
            }
            else { 
              st <- snht(y)
              if(sum(!is.na(st))>0) {
                tVx[i] <- max(st,na.rm=TRUE)
                kpx[i] <- which.max(st)
              }
              cat(round(tVx[i],1)) 
            }
          }
          if(!i%%5) cat("\n") 
        }
        if(i%%5) cat("\n")

        tVxx <- max(tVx,na.rm=TRUE) 
        while(tVxx > tVt) {
          i <- which.max(tVx) 

          if(max(splt[used[i,]])>tVxx*(1+tVf*min(nref,sum(used[i,])))) break
          kp <- kpx[i] 
          cat("\n",as.character(est.c[i,4]),"(",i,")",sep="")
          if(oneref[kp,i] & !force) cat(" could break at ") 
          else cat(" breaks at ")
          if(nm>1) cat(anyi+floor((kp-1)/nm)," ",(kp-1)%%nm+1,sep="")
          else if(nm==1) cat(anyi+kp-1)
          else cat(paste(x[kp]))
          cat(" (",round(tVx[i],1),")",sep="")
          if(oneref[kp,i] & !force) { 

            cat(", but it has only one reference\n")
            tVx[i] <- -1 
            tVxx <- max(tVx,na.rm=TRUE) 
            next
          }

          if(gp>1) {
            plotstan(i,varcli,x,xlab,swa,lw)
            lines(rep(x[kp],2),c(-5,4.8),col="red",lty=2) 
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
            nn <- nn+1 
            iest <<- c(iest,iest[i]) 
            nsp[iest[i]] <- nsp[iest[i]]+1 
            if(nm>0) { 
              z <- 1 + floor((kp-1)/nm) 
              nsy[z] <- nsy[z] + 1 
            }
            dat <<- cbind(dat,rep(NA,nd)) 

            dat[kp:nd,ne+nn] <<- dat[kp:nd,i]
            dat[kp:nd,i] <<- NA 

            z <- data.frame(est.i[iest[i],1:3],V4=paste(est.i[iest[i],4],"-",1+nsp[iest[i]],sep=""),V5=paste(est.i[iest[i],5],"-",1+nsp[iest[i]],sep=""))
            est.c <<- rbind(est.c,z)
            attr(est.c,"row.names")[ne+nn] <<- as.integer(ne+nn)

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

          modif <- TRUE 
          splt[i] <- tVx[i] 
          tVx[i] <- 0 
          tVxx <- max(tVx,na.rm=TRUE) 
        }
        if(nn) {
          cat("\n\nUpdate number of stations: ",ne,"+",nn,"new = ")
          ne <<- ne+nn  
          cat(ne,"\n\n")
          refhom <<- c(refhom,rep(FALSE,nn)) 
        }
        if(!nn && !modif) {
          if(gp>1) {

            w <- ceiling(200/ne) 
            if(w>10) w <- 10
            ylim <- c(0,80)
            fcol <- 4 
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

            z <- tVx[!is.na(tVx) & tVx!=0]
            if(k==1) main <- "Histogram of maximum tV"
            else main <- "Histogram of maximum SNHT"
            if(sum(!is.na(z))) hist(z,breaks=min(0,floor(min(z))):max(20,ceiling(max(z))),xlab=xlab2,col="purple",main=main)

            if(min(z,na.rm=TRUE)<0) hist(z[z<0],breaks=1,col=2,add=TRUE)
            if(k==2 | snhtt<1) {

              hist(nsp,breaks=0:max(9,max(nsp)+1)-.5,col="orange",xlab="Number of splits",ylab="Number of stations",main="Number of splits per station")
              if(nm>0) { 
                w <- min(5,ceiling(400/na)) 
                plot(anyi:anyf,nsy,type="h",lwd=w,col=2,ylim=c(0,max(10,max(nsy))),xlab="Years",ylab="Number of splits",main="Number of splits per year")
                grid(col=grey(.4))
              }
            }
          }
          break 
        }
      }
    }

    tVx <- rep(NA,ne) 
    snhx <- rep(NA,ne) 
    if(gp>1) for(io in 1:nei) { 
      wi <- which(iest==io) 
      lwi <- length(wi)
      if(!lwi) next 
      for(i in wi) { 
        y <- sanom[,i] 
        plotstan(i,varcli,x,xlab,swa,lw)

        st <- snhtw(y,swa); zz <- floor(st[1]); tVx[i] <- st[1]
        if(zz) {
          kp <- st[2]
          lines(rep(x[kp],2),c(-5,4.8),col=verde,lty=2) 
          text(x[kp],5,zz,col=verde) 
        }

        st <- snht(y)
        if(sum(!is.na(st))>0) {
          kp <- which.max(st)
          snhx[i] <- floor(max(st,na.rm=TRUE))
          zz <- floor(snhx[i])
          lines(rep(x[kp],2),c(-5,4.8),lty=4) 
          text(x[kp],-5.2,zz) 
        }

        if(cor.test(y,1:nd)$p.value <= 0.05) {
          aj <- lm(y~x)
          abline(aj,col=4)
        }
      }
    }

    tVx[tVx==0] <- NA

    if(gp>2) {
      plot(-1:1,-1:1,type="n",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
      text(0,0.4,"Final graphics",cex=3.5)
      text(0,-0.3,"Recalculated series\nand applied corrections",cex=2.5)
      if(nm>0) xlab <- "Years" else xlab <- "Dates"

      dh <- scan(paste(varcli,"_",anyi,"-",anyf,".dah",sep="")) 
      dim(dh) <- c(nd,ne) 
      old.par <- par(no.readonly=TRUE)
      layout(matrix(1:2,2,1,byrow=TRUE))
      par(las=1,cex=.8)
      for(i in 1:nei) { 
        wi <- which(iest==i) 
        lwi <- length(wi)
        if(!lwi) next 
        if(lwi>1) vi <- TRUE else vi <- FALSE
        if(nm>1) { 
          fltr <- rep(1,nm) 
          if(gp>3) ylab <- "Running annual totals"
          else { ylab <- "Running annual means"; fltr <- fltr/nm }
        }
        else { fltr <- 1; ylab <- "Values" }
        tit <- paste(varcli," at ",est.i[i,4],"(",i,"), ",est.i[i,5],sep="")
        yo <- as.vector(dat.o[,i]) 
        y <- dh[,wi] 
        par(mar=c(0,4,4,2),xaxt="n")
        matplot(x,filter(y,fltr),type="l",lty=1,col=2:20,ylab=ylab,main=tit)
        lines(x,filter(yo,fltr))
        grid(col=grey(.4))
        par(mar=c(5,4,0,2),xaxt="s")

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

    cat("\nACmx: Station maximum absolute autocorrelations of anomalies\n")
    sac <- rep(NA,ne) 
    for(i in 1:ne) {
      zz <- acf(anom[,i],plot=FALSE,na.action=na.pass)$acf
      zz[1] <- 0 
      sac[i] <- max(abs(zz)) 
    }
    print(summary(round(sac,2)))

    cat("\nSNHT: Standard normal homogeneity test (on anomaly series)\n")
    print(summary(round(snhx,1)))

    cat("\nRMSE: Root mean squared error of the estimated data\n")

    if(std==2 && dat.m[i]>1) { 
      for(i in 1:ne) anom[,i] <- anom[,i] * dat.m[i]
    }
    else if(std==3) {
      for(i in 1:ne) anom[,i] <- anom[,i] * dat.s[i]
    }
    if(rtrans) { 
      rmse <- rep(NA,ne) 
      sdz <- apply(anom,2,sd,na.rm=TRUE)
      for(i in 1:ne) {
        l1 <- (dat.m[i]-sdz[i])^rtrans
        l2 <- (dat.m[i]+sdz[i])^rtrans
        rmse[i] <- (l2-l1)/2 
      }
    }
    else rmse <- apply(anom,2,sd,na.rm=TRUE)
    zz <- summary(rmse)
    print(zz)
    sedec <- max(1,2-ceiling(log10(zz[4]))) 
    pod <- floor(100*(nd-apply(dat.na,2,sum))/nd) 
    cat("\nPD: Percentage of original data\n")
    print(summary(pod))

    cat("\n")
    print(data.frame(ACmx=round(sac,2),SNHT=snhx,RMSE=round(rmse,sedec),PD=pod,Code=est.c[,4],Name=est.c[,5]),right=FALSE)

    cur <- apply(!is.na(dat[(nd-mndat+1):nd,]),2,sum) 
    cur[cur>0] <- 1

    arch <- paste(varcli,"_",anyi,"-",anyf,".esh",sep="") 
    write.table(cbind(est.c,pod,iest,cur,snhx),arch,row.names=FALSE,col.names=FALSE)
    cat("\n----------- Generated outputs: --------------------------------\n\n")
    cat(varcli,"_",anyi,"-",anyf,".txt : This text output\n",sep="")
    cat(varcli,"_",anyi,"-",anyf,".dah : Homogenized data (postprocess with 'dahstat()')\n",sep="")
    cat(varcli,"_",anyi,"-",anyf,".esh : List of homogenized stations (original and split)\n",sep="")
    if(gp>1) { 

      main <- "Histogram of normalized anomalies"
      if(sum(!is.na(outan))) {
        xlim <- range(sanom,na.rm=TRUE)
        xlim[1] <- min(xlim[1], min(outan,na.rm=TRUE))
        xlim[2] <- max(xlim[2], max(outan,na.rm=TRUE))
        zz <- (xlim[2]-xlim[1])/100 
        xlim[1] <- xlim[1] - zz; xlim[2] <-  xlim[2] + zz
        zz <- (xlim[2]-xlim[1])/30 
        breaks=seq(xlim[1],xlim[2],zz)
        hist(sanom,breaks=breaks,xlim=xlim,xlab="Anomaly",col="green",main=main)
        hist(outan,breaks=breaks,col="red",add=TRUE)
      }
      else hist(sanom,breaks=30,xlab="Anomaly",col="green",main=main)
      abline(v=-dz.max,col=2); abline(v=dz.max,col=2)

      z <- snhx; main <- "Histogram of maximum SNHT"
      if(sum(!is.na(z))) hist(z,breaks=20,xlab="SNHT",col="purple",main=main)

      plot(rmse,snhx,type="n",xlim=c(0,max(1,max(rmse,na.rm=TRUE))),ylim=c(0,max(50,max(snhx,na.rm=TRUE))),xlab="RMSE",ylab="SNHT",main="Station's quality/singularity")
      grid(col=grey(.4))
      text(rmse,snhx,col=hsv(.7,1,.9))
    }
    if(gp>0) { 
      graphics.off() 
      cat(varcli,"_",anyi,"-",anyf,".pdf : Diagnostic graphics\n",sep="")
    }
    cat("\n")
  }
  sink() 
}

