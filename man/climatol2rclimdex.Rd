\name{climatol2rclimdex}
\alias{climatol2rclimdex}
\title{Convert DAILY data from \code{climatol} to RClimDEX/RClimPACT input
format}
\description{
  This function reads homogenized daily series of precipitation (RR) and
  extreme temperatures (TX, TN), adjusted from the last homogeneous sub-period,
  and writes them in files (one per station) with RClimDEX format.
}
\usage{
climatol2rclimdex(varRR, varTX, varTN, yiRR, yfRR, yiTX=yiRR, yfTX=yfRR,
  yiTN=yiRR, yfTN=yfRR, header=TRUE, prefix='hoclm', dir=NA, na='-99.9')
}
\arguments{
  \item{varRR, varTX, varTN}{Name of the variables in the \code{climatol} destination files. If some variable is not available, name it as ''.}  
  \item{yiRR, yfRR}{Initial and final years of the homogenized RR series.}
  \item{yiTX, yfTX, yiTN, yfTN}{Initial and final years of the TX and TN series. (The same as yiRR and yfRR by default.)}
  \item{header}{include a header in the files? (\code{TRUE} by default)}
  \item{prefix}{Prefix to prepend to station codes to name the output RClimDEX files.}
  \item{dir}{Destination directory of the output RClimDEX files. (If not set, they will be saved into the current R working directory).}
  \item{na}{Missing data code to use in the ouput files. (\code{'-99.9'} by
  default.)}
}
\details{
After homogenizing daily series with \code{climatol}, the user may be
interested in applying the RClimDEX or RClimPACT/RClimPACTs programs to the
homogenized series. This function automatizes the conversion between both file
formats.

Note that if there are some days with TX<TN (can happen because of the
independent homogenization of extreme temperatures), a trivial fix will be
applied by just exchanging the problem values.
}
\seealso{\code{\link{homogen}}}
\examples{
## Set a temporal working directory and generate input files:
wd <- tempdir()
wd0 <- setwd(wd)

## copy example daily RR, TX and TN homogenization results:
file.copy(exampleFiles('RR_1981-1995.rda'),'.')
file.copy(exampleFiles('TX_1981-1995.rda'),'.')
file.copy(exampleFiles('TN_1981-1995.rda'),'.')

## Now run the example:
climatol2rclimdex('RR','TX','TN',1981,1995)

## Return to user's working directory:
setwd(wd0)

## Input and output files can be found in directory:
print(wd)
}
\keyword{manip}
