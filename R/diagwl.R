#Diagrama de Walter y Lieth
diagwl <- function( dat, est="", alt=NA, per="", margen=c(4,4,5,4),
    mlab="", pcol="#005ac8", tcol="#e81800", pfcol="cyan",
    sfcol="#0eb6d7", shem=FALSE, ...) {
  plot.new()
  old.par <- par(no.readonly = TRUE)
  on.exit(par(old.par))
  par(mar=margen, pty="s", las=1, new=FALSE)
  nr <- nrow(dat) #no. de filas de datos mensuales
  if(nrow(dat) < 4) { cat("4 rows of monthly data must be provided\n"); break; }
  #etiquetas de los meses
  if(mlab=="es") mlab=c("E","F","M","A","M","J","J","A","S","O","N","D")
  else if(mlab=="en") mlab=c("J","F","M","A","M","J","J","A","S","O","N","D")
  else mlab=c(1:12) #etiquetas numéricas
  dat <- as.matrix(dat)
  if(shem) { #Hemisferio sur: desplazar los datos medio año
    m1 <- dat[,1:6]
    m2 <- dat[,7:12]
    dat <- cbind(m2,m1)
    mlab <- c(mlab[7:12],mlab[1:6])
  }
  p <- dat[1,] #precipitaciones medias mensuales
  if(nr==2) tm <- dat[2,]
  else tm <- apply(dat[2:3,],2,mean)  #temperaturas medias mensuales
  pmax <- max(p) #precipitación máxima
  ymax <- 60  #máxima ordenada por defecto
cat("ymax=",ymax,"\n");
  if(pmax > 300) ymax <- 50 + 10*floor((pmax+100)/200)
cat("pmax=",pmax,"ymax=",ymax,"\n");
  ymin <- min(-1.5,min(tm)) #mínima ordenada sin redondear
  #ejes:
  if(ymin < -1.5) {
    ymin=floor(ymin/10)*10 #mínima ordenada redondeada
cat("ymin=",ymin,"\n")
    labT <- paste(ymin)
    labP <- ""
    if(ymin < -10) {
      for(i in (ymin/10+1):-1) {
        labT <- c(labT,i*10)
        labP <- c(labP,"")
cat(i,i*10,"\n",labT,"\n",labP,"\n")
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
  plot(0:13-0.5,c(tm[12],tm[1:12],tm[1]),xlim=c(0,12),ylim=c(ymin,ymax),type="l",col=tcol,xaxs="i",yaxs="i",xaxp=c(0,12,12),xlab="",ylab="",xaxt="n",yaxt="n",bty="n",...)
  lmin <- ymin #mínima ordenada a rotular
  if(lmin==-1.5) lmin=0
  axis(2,((lmin/10):(ymax/10))*10,labels=labT,col.axis=tcol)
  axis(4,((lmin/10):(ymax/10))*10,labels=labP,col.axis=pcol)
  mtext("C",2,col=tcol,las=1,line=3,adj=0,at=55)
  mtext("mm",4,col=pcol,las=1,line=3,adj=1,at=55)
  abline(0,0)
  abline(50,0)
  #rótulos:
  if(is.na(alt)) mtext(est,line=2,adj=0)
  else mtext(paste(est," (",alt," m)",sep=""),line=2,adj=0)
  mtext(per,line=1,adj=0)
  mtext(paste(round(mean(tm*10))/10,"C        ",round(sum(p))," mm",sep=""),line=1,adj=1)
  x <- 0:13-0.5
  p2 <- c(p[12],p[1:12],p[1])
  if(pmax<=100) {
    xl <- x
    yl <- c(p[12],p[1:12],p[1])/2
    n2 <- 14
  }
  else { #cambio de escala en prec. > 100
    xp <- numeric(30)
    yp <- numeric(30)
    xl <- numeric(25)
    yl <- numeric(25)
    n <- 0
    n2 <- 0
    gr <- FALSE
    if(p2[1]>100) { #primer punto
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
    for(i in 2:14) {  #demás puntos
      if(gr) {  #si p anterior > 100
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
      else {  # p anterior <=100
        if(p2[i]>100) { #si p > 100
          n <- n+1
          xp[n] <- x[i-1] + (100-p2[i-1])/(p2[i]-p2[i-1])
          yp[n] <- 50
          if(xl[n2]!=x[i-1]) {  #no repetir puntos!
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
    if(!is.na(yp[n])) {  #cerrar último polígono
      n <- n+1
      xp[n] <- xp[n-1]
      yp[n] <- 50
      n2 <- n2+1
      xl[n2] <- 12.5
      yl[n2] <- 50
    }
    polygon(xp[1:n],yp[1:n],col=pcol,border=pcol)
  }
  #tramas:
  pi <- approx(xl[1:n2],yl[1:n2],n=66)$y
  ti <- approx(x,c(tm[12],tm[1:12],tm[1]),n=66)$y
  d <- pi - ti
  xi <- (1:66)/5-0.7
  xw <- subset(xi,d>0) #periodo húmedo
  y1 <- subset(pi,d>0)
  y2 <- subset(ti,d>0)
  if(length(xw)>0) segments(xw,y1,xw,y2,col=pcol,lty=1,lwd=1)
  xw <- subset(xi,d<0) #periodo seco
  y1 <- subset(pi,d<0)
  y2 <- subset(ti,d<0)
  if(length(xw)>0) segments(xw,y1,xw,y2,col=tcol,lty=3,lwd=2)
  #curvas de P y T:
  lines(xl[1:n2],yl[1:n2],col=pcol,lwd=2)
  lines(x,c(tm[12],tm[1:12],tm[1]),col=tcol,lwd=2)
  #media de las máximas del mes más cálido
  mtext(formatC(max(as.matrix(dat[2,])),digits=1,format="f"),2,las=1,
    line=2,at=35)
  #media de las mínimas del mes más frío
  mtext(formatC(min(as.matrix(dat[3,])),digits=1,format="f"),2,las=1,
    line=2,at=15)
  #marcar límites de los meses:
  for(i in 0:13) segments(i,0,i,-1.5)
  #heladas seguras
  for(i in 1:12) if(dat[3,i]<=0) rect(i-1,-1.5,i,0,col=sfcol)
  #heladas probables
  for(i in 1:12) if(dat[4,i]<=0 && dat[3,i]>0) rect(i-1,-1.5,i,0,col=pfcol)
  #rótulos meses:
  mtext(mlab,1,las=1,line=0.5,adj=0.5,at=x[2:13])
  #reset old.par (restablecemos parámetros gráficos anteriores):
  invisible()
}
