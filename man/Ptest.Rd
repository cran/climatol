\name{Ptest}
\alias{Ptest}
\alias{dat}
\alias{est.c}
\docType{data}
\title{Monthly precipitation data}
\description{
  Example of a network of 20 monthly precipitation series for 1951-2010.
}
\usage{data(Ptest)}
\format{This data set contains:
\describe{
\item{dat}{Monthly precipitation data arranged in a 720x20 matrix}
\item{est.c}{Data frame with station's coordinates (X,Y,Z), codes and names}
}
}
\source{Subset of synthetic precipitation data from the Spanish research project MULTITEST (2015-2017)}
\seealso{\code{\link{homogen}}, \code{\link{homogsplit}}}
\examples{
data(Ptest)
print(est.c)
str(dat)
rm(est.c,dat)
}
\keyword{datasets}

