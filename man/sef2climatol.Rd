\name{sef2climatol}
\alias{sef2climatol}
\title{Convert SEF data files to \code{climatol} input files.}
\description{
  This function reads all SEF files contained in a directory and writes their data in \code{*.dat} and \code{*.est} \code{climatol} input files.
}
\usage{
  sef2climatol(dr,Vbl,varcli=Vbl,ndec=1,na.strings="NA",mindat=NA)
}
\arguments{
  \item{dr}{directory containing the SEF files}
  \item{Vbl}{name of the variable in the SEF files}
  \item{varcli}{name of the variable in the \code{climatol} destination files}
  \item{ndec}{number of decimals to save}
  \item{na.strings}{missing data codes (specified as quoted strings)}
  \item{mindat}{minimum required number of data per station}
}
\details{
  SEF (Station Exchange Format) is the Copernicus Climate Change Service
format for Data Rescue projects. Visit https://datarescue.climate.copernicus.eu/node/80

  Missing elevations will be assigned the value 99

  Some files may contain a single quotation mark in the metadata field,
  causing not reading the end of line until a pairing quoting is found in
  the following line, hence skipping half of the data. Parameter quote='\\'
  has been set in the reading command as a workaround.

  All data are dumped into a temporary file named \code{SEFdata.csv}, which is used by the function \code{csv2climatol} to write the input files for \code{climatol}.
}
\seealso{\code{\link{csv2climatol}}, \code{\link{homogen}}}
\examples{
## Set a temporal working directory and write input files:
wd <- tempdir()
wd0 <- setwd(wd)

## Create a directory and copy all SEF files to be processed:
dir.create('dir1')
file.copy(exampleFiles('GHCN_v4_Bhamo.tsv'),'dir1')
file.copy(exampleFiles('GHCN_v4_Diamond_Island.tsv'),'dir1')

## Now run the function:
sef2climatol('dir1','ta')

## Return to user's working directory:
setwd(wd0)

## Input and output files can be found in directory:
print(wd)
}
\keyword{manip}
