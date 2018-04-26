\name{tmax}
\alias{tmax}
\docType{data}
\title{Daily temperature data}
\description{
  Maximum daily temperatures for three stations in 2001-2003.
}
\usage{data(tmax)}
\format{This data set contains:
\describe{
\item{dat}{Daily maximum temperature data arranged in a 1095x3 matrix}
\item{est.c}{Data frame with station's coordinates (X,Y,Z), codes and names
(invented)}
}
}
\source{Subset of synthetic temperature data from Killick, R. E. (2016):
Benchmarking the Performance of Homogenisation Algorithms on Daily Temperature
Data, PhD Thesis, University of Exeter.}
\seealso{\code{\link{climatol2rclimdex}}}
\examples{
data(tmax)
print(est.c)
head(dat)
#clean memory space:
rm(est.c,dat)
}
\keyword{datasets}

