\name{Ttest}
\alias{Ttest}
\docType{data}
\title{Daily temperature data}
\description{
  Example of a network of 12 daily temperature series for 1981-2000.
}
\usage{data(Ttest)}
\format{This data set contains:
\describe{
\item{dat}{Daily mean temperature data arranged in a 7305x12 matrix}
\item{est.c}{Data frame with station's coordinates (X,Y,Z), codes and names (invented)}
}
}
\source{Subset of synthetic temperature data from Killick, R. E. (2016):
Benchmarking the Performance of Homogenisation Algorithms on Daily Temperature Data, PhD Thesis, University of Exeter.}
\seealso{\code{\link{homogen}}, \code{\link{dd2m}}, \code{\link{homogsplit}}}
\examples{
data(Ttest)
print(est.c)
str(dat)
rm(est.c,dat)
}
\keyword{datasets}

