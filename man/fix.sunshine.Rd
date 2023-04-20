\name{fix.sunshine}
\alias{fix.sunshine}
\title{Check homogenized daily sunshine hours and prune any excess}
\description{
  This function loads homogenization results of daily sunshine series and
  prunes any excess over maximum theoretical sunshine duration.
}
\usage{ fix.sunshine(varcli, anyi, anyf) }
\arguments{
\item{varcli}{Short name of the homogenized climatic variable.}
\item{anyi}{First year of the homogenized series.}
\item{anyf}{Last year of the homogenized series.}
}
\details{
  Any modified value is listed to the console and written to \code{fix.sunshine.txt}

The original \code{*.rda} file is saved as \code{*.rda.bak} and a new
  \code{*.rda} file is written with the fixed sunshine values.
}
\seealso{\code{\link{homogen}}}
\examples{
## Set a temporal working directory:
wd <- tempdir()
wd0 <- setwd(wd)

## copy example daily sunshine homogenization results:
file.copy(exampleFiles('SS_1991-2000.rda'),'.')

## Now run the example:
fix.sunshine('SS',1991,2000)

## Return to user's working directory:
setwd(wd0)

## Input and output files can be found in the directory:
print(wd)
}
\keyword{manip}
\keyword{utilities}
