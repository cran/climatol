\name{leerdm}
\alias{leerdm}
\title{Read the climatological series to homogenize}
\description{
  Climatological series, and coordinates and names of the stations, are
  read from \file{VAR_ANYI-ANYF.dat} and \file{VAR_ANYI-ANYF.est}
  files respectively.
}
\usage{
  leerdm(varcli, anyi, anyf, b=1, a=0, na.strings="NA")
}
\arguments{
  \item{varcli}{Acronym of the name of the studied climatic variable,
    as in the data file name.}
  \item{anyi}{Initial year of the data present in the file.}
  \item{anyf}{Final year of the data present in the file.}
  \item{a, b}{Parameters of the optional transformation \code{a+b*dat} to be
    applied to data when read from the files.}
}
\details{
  This function was written to be called by \code{homogen}, but may be
  invoqued by the user to read the climatic series.
}
\value{
  No value is returned by this function, but it generates the following
  objects in the space memory:
  \item{ne}{Number of stations.}
  \item{nd}{Number of data in each station.}
  \item{dat}{Matrix data read from the input files, with
  dimensions \code{nd,ne}.}
  \item{est.c}{Data frame with the coordinates, codes and names of
  the stations.}
}
\seealso{\code{\link{homogen}}}
\examples{
  #As this function only works with files, you must uncompress the
  #example data in your working directory first. (The example files
  #can be obtained from http://webs.ono.com/climatol/climatol-dat.zip).
  #After uncompressing them, we can run:
  leerdm("Tmin", 1956, 2005)
  #data are now in the memory space, in objects dat, est.c, ne and nd.
}
\keyword{internal}
\keyword{IO}
