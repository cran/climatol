\name{datsubset}
\alias{datsubset}
\title{Subset data by subperiod, code list or no. of years with data}
\description{
  This function allows saving a subset of \code{climatol} input data into new input
files by selecting a subperiod, a minimum number of years with data and/or a group of stations.
}
\usage{
datsubset(varcli, anyi, anyf, anyis=anyi, anyfs=anyf, minny=NA, codes=NULL,
na.strings=NA, ini=NA)
}
\arguments{
  \item{varcli}{Short name of the studied climatic variable.}
  \item{anyi}{Initial year of the data present in the file.}
  \item{anyf}{Final year of the data present in the file.}
  \item{anyis}{First year of the output subperiod. (Defaults to \code{anyi}).}
  \item{anyfs}{Last year of the output subperiod. (Defaults to \code{anyf}).}
  \item{minny}{Minimum number of years with data to retain the series.}
  \item{codes}{Vector of chosen station codes. (Defaults to \code{NULL}, meaning all).}
  \item{na.strings}{Strings marking missing data (\code{NA} by default).}
  \item{ini}{Initial date (if not January 1st).}
}
\details{
Homogenization by \code{climatol} requires that no time step be totally void of data in all stations simultaneously. This function allows subsetting already existing \code{climatol} input files by selecting a subperiod and/or stations with a minimum number of years with data (may contain gaps).

  Another possibility is to choose a group of stations, useful when the initial cluster analysis reveals areas with different climate regimes that should be homogenized independently.
}
\examples{
## Set a temporal working directory and write input files:
wd <- tempdir()
wd0 <- setwd(wd)
data(climatol_data)
write.table(Temp.est,'Temp_1961-2005.est',row.names=FALSE,col.names=FALSE)
write(Temp.dat,'Temp_1961-2005.dat',ncolumns=12)

## Now run the examples:
datsubset('Temp',1961,2005,1971,2000,minny=20)
datsubset('Temp',1971,2000,codes=c('st02','st03'))

## Return to user's working directory:
setwd(wd0)

## Input and output files can be found in directory:
print(wd)
}
\keyword{manip}
