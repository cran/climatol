\name{datrestore}
\alias{datrestore}
\title{Rename homogen's output files}
\description{
  This function restores some deleted outliers into the dah matrix of the \code{*.rda} output file.
}
\usage{
datrestore(varcli, anyi, anyf, QCout=FALSE)
}
\arguments{
  \item{varcli}{Short name of the studied climatic variable,
    as in the data file name.}
  \item{anyi}{Initial year of the study period.}
  \item{anyf}{Final year of the study period.}
  \item{QCout}{Set this parameter to \code{TRUE} to read the selected outliers to be restored from the \code{*-QC*_out.csv} output by a previous run of \code{homogen(..., onlyQC=TRUE)}.}
}
\details{
  When the user checks the list of outliers in the output file \code{*_out.csv}, true extreme values (or sequences of identical values) that have been deleted can be restored by changing their \code{deleted} field to negative.
  This accepted values will be restored in the \code{dah} matrix of homogenized series contained in the \code{*.rda} file output by the \code{homogen} function, but only in the series reconstructed from the last homogeneous subperiod.
}
\seealso{\code{\link{homogen}}}
\examples{
## Set a temporal working directory, write input files and homogenize them:
wd <- tempdir()
wd0 <- setwd(wd)
data(climatol_data)
write.table(Temp.est,'Temp_1961-2005.est',row.names=FALSE,col.names=FALSE)
write(Temp.dat,'Temp_1961-2005.dat')
homogen('Temp',1961,2005) #obtain homogenization output files
out <- read.csv('Temp_1961-2005_out.csv') #read list of outliers
## Change the sign of a couple of deleted values to be restored:
out[2,6] <- -1; out[6,6] <- -9
write.csv(out,'Temp_1961-2005_out.csv',row.names=FALSE)

## Now run the example:
datrestore('Temp',1961,2005) #restore the selected values

## Return to user's working directory:
setwd(wd0)

## Input and output files can be found in directory:
print(wd)
}
\keyword{utilities}
