\name{daily2climatol}
\alias{daily2climatol}
\title{Convert daily data files to Climatol input format}
\description{
  This function can be useful to prepare the Climatol input files when the user have their daily data in per station individual files.
}
\usage{
daily2climatol(stfile, stcol=1:6, datcol=1:4, varcli, anyi=NA, anyf=NA,
  mindat=365, sep='', dec='.', na.strings='NA', header=FALSE)
}
\arguments{
  \item{stfile}{Name of the file with the list of data file names and station coordinates, codes and names.}  
  \item{stcol}{Columns in \code{stfile} holding data file names, longitudes, latitudes, elevations and station codes and names. (Defaults to 1:6. Use 0 for codes and/or names columns if they are missing, and numeric values will be assigned.)}
  \item{datcol}{Columns in data files (named as in the first column of \code{stfile}) holding year, month, day, value.}
  \item{varcli}{Acronym of the name of the studied climatic variable.}
  \item{anyi}{First year to study (defaults to the first year of available data).}
  \item{anyf}{Last year to study (defaults to the last year of available data).}
  \item{mindat}{Minimum required number of data per station. (Defaults to 365 daily data.)}
  \item{sep}{Field separator in all files, whether data or stations. (Defaults to white space.)}
  \item{dec}{Symbol of the decimal point. (Defaults to '.'.)}
  \item{na.strings}{Missing data code in the original daily data. (Defaults to 'NA'.)}
  \item{header}{Logical value indicating whether input files have a header line or not. (Defaults to \code{FALSE}, but this parameter does not apply to CSV files, which are expected to have a header.)}
}
\details{
Many users have their daily series in separate files (one per station), either
in text format (as the files used by RClimDex) or in spread-sheets (from which
they can be exported to CSV text files). This function can be used to read
these daily data files and write the input files needed by the \code{homogen}
function of this Climatol package.

When either station codes or names are missing in the stations file, its
corresponding column must be set to 0. In this case, codes and/or names will be
assigned with numeric values.

Field separator, decimal point and the presence of a header line must be consistent in all files (data files and station list).
}
\value{
  This function does not return any value.
}
\seealso{\code{\link{homogen}}}
\examples{
wd <- tempdir() #temporal working directory for the examples
wd0 <- setwd(wd)
data(dailies)
write.table(stations,'stations.txt',sep=',',row.names=FALSE,col.names=FALSE)
write.table(WY003,'WY003.txt',sep=',',row.names=FALSE,col.names=FALSE)
write.table(WY018,'WY018.txt',sep=',',row.names=FALSE,col.names=FALSE)
write.table(WY020,'WY020.txt',sep=',',row.names=FALSE,col.names=FALSE)
daily2climatol('stations.txt',c(1,3,2,4,0,0),2:5,'TestVar',sep=',',na.strings='-99.9')
#Return to user's working directory:
setwd(wd0)
#Input and output files can be found in directory:
print(wd)
}

