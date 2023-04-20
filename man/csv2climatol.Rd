\name{csv2climatol}
\alias{csv2climatol}
\title{Convert data in a single CSV file to \code{climatol} input format}
\description{
  This function helps to prepare the \code{climatol} input files when the users
 have their data in a single CSV file, as the output of xls2csv().
}
\usage{
  csv2climatol(csvfile, datacol=6:8, stnfile=csvfile, stncol=1:5, varcli,
  anyi=NA, anyf=NA, mindat=NA, sep=',', dec='.', na.strings='NA',
  dateformat='\%Y-\%m-\%d', cf=1, ndec=1, header=TRUE)
}
\arguments{
  \item{csvfile}{name of the CSV file containing the data.}
  \item{datacol}{column(s) holding station codes, dates and data.}
  \item{stnfile}{name of the CSV file containing station codes, names and
     coordinates (if these data are not in the csvfile).}
  \item{stncol}{columns holding longitudes, latitudes, elevations and station
     codes and names.}
  \item{varcli}{short name of the climatic variable under study.}
  \item{anyi}{first year to study.}
  \item{anyf}{last year to study.}
  \item{mindat}{minimum required number of data per station (by default,
    60 monthly data or 365 daily data).}
  \item{sep}{data separator (',' by default: Comma Separated Values).}
  \item{dec}{decimal point ('.' by default).}
  \item{na.strings}{strings coding missing data ('NA' by default).}
  \item{dateformat}{format of dates (if not in separate columns.
    Default: \code{'\%Y-\%m-\%d')}}
  \item{cf}{conversion factor to apply if data units need to be changed.}
  \item{ndec}{no. of decimals to round to.}
  \item{header}{\code{TRUE} by default, set to \code{FALSE} if \code{csvfile} has no header.}
}
\details{
  If \code{datacol} holds  4 (or 5) values, dates are expected to appear as
  year, month (and days) in separate columns. Otherwise, dates will be provided
  as character strings (see parameter \code{dateformat}).
  Station codes, names and coordinates can go in a separate file
  \code{stnfile}. At least coordinates and station codes must be present in
  either \code{csvfile} or \code{stnfile}. Put a zero for any inexistent
  columns. Example when \code{stnfile} contains only, in this order, latitudes,
  longitudes and station names: \code{stncol=c(2,1,0,3,0)}.
  Note that if a stnfile is provided, then sep, dec, na.strings and header
  defined for csvfile will also be applied to stnfile.
}
\seealso{\code{\link{xls2csv}}, \code{\link{homogen}}}
\examples{
## Set a temporal working directory:
wd <- tempdir()
wd0 <- setwd(wd)

## Create origin and destination directories and copy example input files:
dir.create('dir1'); dir.create('dir2')
file.copy(exampleFiles('p064.xlsx'),'dir1')
file.copy(exampleFiles('p082.xlsx'),'dir1')
file.copy(exampleFiles('p084.xlsx'),'dir1')

## Create input files for csv2climatol with the function xls2csv:
xls2csv('dir1','dir2','RR')

## Add bogus coordinates and elevations to the station file:
est=read.table('xls_RR_stations.csv',sep=',')
est=data.frame(1:3,21:23,101:103,est)
write.table(est,'xls_RR_stations.csv',sep=',',row.names=FALSE,col.names=FALSE)

## Now run the example of csv2climatol:
csv2climatol('xls_RR_data.csv', datacol=1:5, stnfile='xls_RR_stations.csv',
  varcli='RR',header=FALSE)

## Return to user's working directory:
setwd(wd0)

## Input and output files can be found in directory:
print(wd)
}
\keyword{manip}
