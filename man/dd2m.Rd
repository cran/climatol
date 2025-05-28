\name{dd2m}
\alias{dd2m}
\title{Compute monthly data from daily (or subdaily) series}
\description{
  Daily or sub-daily series are aggregated into total, mean, maximum, or
  minimum monthly values, and saved to files in \code{climatol} input format.
}
\usage{
dd2m(varcli, anyi, anyf, ndec=1, valm=2, namax=30, x=NULL, na.strings="NA",
tz='utc')
}
\arguments{
  \item{varcli}{Short name of the studied climatic variable,
    as in the data file name.}
  \item{anyi}{Initial year of the data present in the file.}
  \item{anyf}{Final year of the data present in the file.}
  \item{ndec}{Number of decimal places to be saved in the output file.}
  \item{valm}{Monthly value to compute:
    \describe{
      \item{1:}{Sum,}
      \item{2:}{Mean,}
      \item{3:}{Maximum,}
      \item{4:}{Minimum,}
      \item{5:}{Standard deviation.}
    }
  }
  \item{namax}{Maximum percentage of missing data in any month to compute its
  monthly value. (30 by default)}
  \item{x}{Time vector. If not provided, it will be built as dates (or date-time for sub-daily data) beginning January 1st of the initial year. The user must provide if data are taken at irregular intervals or they don't start on January 1st.}
  \item{na.strings}{Missing data code in the original daily data. (\code{NA} by default.)}
  \item{tz}{Time zone (\code{'utc'} by default). Only relevant for subdaily data.}
}
\details{
Data are read from files \file{VRB_YEAR-YEAR.dat} and
\file{VRB_YEAR-YEAR.est}, and output monthly data will be saved to
files with the same names but with the suffix \code{-m} appended to the name of the variable.
Series with less than 12 monthly data will be removed from the output monthly series.
}
\seealso{\code{\link{homogen}}, \code{\link{dahstat}}, \code{\link{dahgrid}}}
\examples{
## Set a temporal working directory and write input files:
wd <- tempdir()
wd0 <- setwd(wd)
data(climatol_data)
write.table(SIstations,'RR_1981-1995.est',row.names=FALSE,col.names=FALSE)
write(as.matrix(RR3st[,2:4]),'RR_1981-1995.dat')

## Now run the example:
dd2m('RR',1981,1995,valm=1)

## Return to user's working directory:
setwd(wd0)

## Input and output files can be found in directory:
print(wd)
}
\keyword{datagen}
\keyword{manip}

