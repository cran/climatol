\name{daily2climatol}
\alias{daily2climatol}
\title{Convert daily data files to \code{climatol} input format}
\description{
  This function can be useful to prepare the \code{climatol} input files when the users have their daily data in per station individual files.
}
\usage{
daily2climatol(stfile, stcol=1:6, datcol=1:4, varcli='VRB', anyi=NA, anyf=NA,
  mindat=365, sep=',', dec='.', na.strings='NA', header=TRUE)
}
\arguments{
  \item{stfile}{File with file names and station coordinates, codes and names.}
  \item{stcol}{Columns in \code{stfile} holding data file names, longitudes,
  latitudes, elevations and station codes and names. (Defaults to 1:6. Use 0
  for codes and/or names columns if they are missing, and numeric values will
  be assigned.)}
  \item{datcol}{Columns in data files holding year, month, day, value.}
  \item{varcli}{(Short) name of the studied climatic variable.}
  \item{anyi}{First year to study (defaults to the first year of available data).}
  \item{anyf}{Last year to study (defaults to the last year of available data).}
  \item{mindat}{Minimum required number of data per station. (Defaults to 365
  daily data.)}
  \item{sep}{Field separator in all files, whether data or stations. (Defaults
  to white space.)}
  \item{dec}{Decimal point. ('.' by default.)}
  \item{na.strings}{Strings coding missing data (\code{'NA'} by default).}
  \item{header}{Logical value indicating whether input files have a header line
  or not. (\code{TRUE} by default.)}
}
\details{
Many users have their daily series in separate files (one per station). This
function can be used to read these daily data files and write the input files
needed by the \code{homogen} function of this \code{climatol} package.

When either station codes or names are missing in the stations file, its
corresponding column must be set to 0. In this case, codes and/or names will be
assigned with numeric values.

Field separator, decimal point and the presence of a header line must be consistent in all files (data files and stations file).

If your files follow the RClimDex convention, you can use the \code{rclimdex2climatol} function instead.
}
\seealso{\code{\link{rclimdex2climatol}}, \code{\link{homogen}}}
\examples{
## Set a temporal working directory and write example input files:
wd <- tempdir()
wd0 <- setwd(wd)
data(climatol_data)
df=cbind(File=c('p064.csv','p084.csv','p082.csv'),SIstations)
write.csv(df,'stations.csv',row.names=FALSE,quote=FALSE)
write.csv(p064.df,'p064.csv',row.names=FALSE,quote=FALSE)
write.csv(p084.df,'p084.csv',row.names=FALSE,quote=FALSE)
write.csv(p082.df,'p082.csv',row.names=FALSE,quote=FALSE)

## Now run the example:
daily2climatol('stations.csv',varcli='RR')

## Return to user's working directory:
setwd(wd0)

## Input and output files can be found in directory:
print(wd)
}
\keyword{manip}
