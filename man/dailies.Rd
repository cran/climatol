\name{dailies}
\alias{dailies}
\alias{stations}
\alias{WY003}
\alias{WY018}
\alias{WY020}
\docType{data}
\title{Daily temperature data}
\description{
  Daily temperatures (1981-1990) for three stations, in separate files.
}
\usage{data(dailies)}
\format{This data set contains:
\describe{
\item{stations}{Data frame with file names, latitudes, longitudes and
elevations (invented)}
\item{WY003, WY018, WY020}{Data frames containing station codes, years, months,
days and daily mean temperatures}
}
}
\source{Subset of synthetic temperature data from Killick, R. E. (2016):
Benchmarking the Performance of Homogenisation Algorithms on Daily Temperature
Data, PhD Thesis, University of Exeter.}
\seealso{\code{\link{daily2climatol}}}
\examples{
data(dailies)
print(stations)
head(WY003)
#clean memory space:
rm(stations, WY003, WY018, WY020)
}
\keyword{datasets}

