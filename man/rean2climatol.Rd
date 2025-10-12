\name{rean2climatol}
\alias{rean2climatol}
\title{Append a reanalysis series in CSV to \code{climatol} input files.}
\description{
  This function reads a reanalysis series in CSV format and appends their corresponging data to the climatol input file \code{*.dat}. A line is also appended to the stations file \code{*.est}.
}
\usage{
  rean2climatol(varcli, anyi, anyf, csvfile, X, Y, skip=0, datacol=1:2,
  dateformat='\%Y-\%m-\%d', a=0, b=1, ndec=1, rname='rean', sep=',', dec='.')
}
\arguments{
  \item{varcli}{Short name of the climatic variable under study.}
  \item{anyi, anyf}{First and last years to study.}
  \item{csvfile}{Name of the CSV (text) file containing the reanalysis series.}
  \item{X, Y}{Longitude and Latitude of the series in degrees.}
  \item{skip}{Number of header lines to be skipped. [0]}
  \item{datacol}{Columns holding dates and data. [1:2]}
  \item{dateformat}{Format of dates if not in separate columns. ['\%Y-\%m-\%d']}
  \item{a, b}{Parameters of the optional transformation \code{a+b*dat} to be applied to the reanalyis series to adapt them to the units used in our data files. [0, 1]}
  \item{ndec}{No. of decimals to round. [1]}
  \item{rname}{Name of the reanalysis. ['rean']}
  \item{sep}{Data separator (',' by default: Comma Separated Values).}
  \item{dec}{Decimal point ('.' by default).}
}
\details{
  If 3 (4) values are provided in the \code{datacol} parameter, dates are expected to appear as year, month (and day) in separate columns. Otherwise, dates will be provided as character strings with \code{dateformat} format.

  The reanalysis series will be assigned an elevation of 99 m. (The user can edit the \code{*.est} file to change it.)
}
\seealso{\code{\link{homogen}}}
\examples{
## Set a temporal working directory and write input files:
wd <- tempdir()
wd0 <- setwd(wd)

## Save example climatol input files:
load('~/cc/climatol/data/climatol_data.rda')
dat <- as.matrix(RR3st[,2:4])
write(dat,'Prec_1981-1995.dat')
write.table(SIstations,'Prec_1981-1995.est',row.names=FALSE,col.names=FALSE)

## Copy the provided simulated reanalysis series:
file.copy(exampleFiles('bogusERA5.txt'),'.')

## Now run the function. (Reanalysis precipitation rates in kg/m2/s
##  are converted to mm with parameter b):
rean2climatol('Prec', 1981, 1995, 'bogusERA5.txt', 15, 46.2, datacol=c(2:4,1),
  rname='bogusERA5', b=86400)

## Return to user's working directory:
setwd(wd0)

## New output files can be found in directory:
print(wd)
}
\keyword{manip}
