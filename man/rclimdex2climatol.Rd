\name{rclimdex2climatol}
\alias{rclimdex2climatol}
\title{Convert RClimDex daily data files to \code{climatol} input format}
\description{
  This function can be useful to prepare the \code{climatol} input files when the user have their daily data in RClimDex format.
}
\usage{
rclimdex2climatol(stfile, stcol=1:5, kvar, chrcod=c(6,10), sep='', anyi=NA,
anyf=NA, mis=-99.9, mindat=365, header=TRUE)
}
\arguments{
  \item{stfile}{Name of the file with the list of data file names and station coordinates, codes and names.}  
  \item{stcol}{Columns in \code{stfile} holding data file names, longitudes,
  latitudes, elevations and station codes and names. (Defaults to 1:5.)}
  \item{kvar}{RClimDex variable to extract: 1 (RR), 2 (TX), 3 (TN).}
  \item{chrcod}{Initial and final characters of data file names to be used as station codes. c(6,10) by default.}
  \item{sep}{Field separator in \code{stfile} and data files (space or tab by default).}
  \item{anyi}{First year to study. (Defaults to the first year of available data.)}
  \item{anyf}{Last year to study. (Defaults to the last year of available data.)}
  \item{mis}{Missing data code. (Defaults to -99.9.)}
  \item{mindat}{Minimum required number of data per station. (Defaults to 365 daily data.)}
  \item{header}{Do files have a header line? \code{TRUE} by default).}
}
\details{
Users of the RClimDex program can convert their daily data files to the
\code{climatol} format. All files listed in \code{stfile} will be read, and the
selected variable (precipitation, maximum or minimum temperature) will be
stored in a unique \code{*.dat} file, with its companion \code{*.est} station
file.  Therefore, if you want to convert all three variables, you must run this
function three times.

Coordinates must be given in degrees with decimals, using the minus sign for sourthern latitudes and western longitudes.
}
\seealso{\code{\link{homogen}}, \code{\link{climatol2rclimdex}}}
\examples{
## Set a temporal working directory:
wd <- tempdir()
wd0 <- setwd(wd)

## Prepare a few files in RClimDex format:
data(climatol_data)
gY=c(46,46,46); mY=c(06,15,14); sY=c(42,25,53)
gX=c(14,15,14); mX=c(03,09,50); sX=c(05,06,05)
df=data.frame(File=c('p064.txt','p084.txt','p082.txt'),
LatDeg=gY,LatMin=mY,LatSec=sY,LonDeg=gX,LonMin=mX,LonSec=sX,
elev=SIstations[,3],name=SIstations[,5])
write.table(df,'stations.txt',sep='\t',row.names=FALSE)
write.table(p064.df,'p064.txt',sep='\t',row.names=FALSE,quote=FALSE)
write.table(p084.df,'p084.txt',sep='\t',row.names=FALSE,quote=FALSE)
write.table(p082.df,'p082.txt',sep='\t',row.names=FALSE,quote=FALSE)

## Now run the example:
rclimdex2climatol('stations.txt',3,chrcod=c(1,4))

## Return to user's working directory:
setwd(wd0)

## Input and output files can be found in directory:
print(wd)
}
\keyword{manip}
