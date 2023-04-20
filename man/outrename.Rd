\name{outrename}
\alias{outrename}
\title{Rename homogen's output files}
\description{
  This function inserts a suffix to the output file names of \code{homogen},
  to prevent them from being overwritten by any further run.
}
\usage{
outrename(varcli, anyi, anyf, suffix, restore=FALSE)
}
\arguments{
  \item{varcli}{Short name of the studied climatic variable,
    as in the data file name.}
  \item{anyi}{Initial year of the study period.}
  \item{anyf}{Final year of the study period.}
  \item{suffix}{Suffix to be inserted (or removed) in the output file names.}
  \item{restore}{Set this parameter to \code{TRUE} to remove the \code{suffix}
    previously inserted by this function. (\code{FALSE} by default.}
}
\details{
  The suffix is appended to the \code{varcli} after a hyphen. The purpose of
  this function is to allow a new application of \code{homogen} to the same
  data with different parameters without overwriting the previous results.
}
\seealso{\code{\link{homogen}}}
\examples{
## Set a temporal working directory, write input files and homogenize them:
wd <- tempdir()
wd0 <- setwd(wd)
data(climatol_data)
write.table(Temp.est,'Temp_1961-2005.est',row.names=FALSE,col.names=FALSE)
write(Temp.dat,'Temp_1961-2005.dat',ncolumns=12)
datsubset('Temp',1961,2005,1991) #subset data to shorten example run time
homogen('Temp',1991,2005) #obtain homogenization output files

## Now run the example:
outrename('Temp',1991,2005,'bak') #rename them to avoid being overwritten

## Return to user's working directory:
setwd(wd0)

## Input and output files can be found in directory:
print(wd)
}
\keyword{utilities}
