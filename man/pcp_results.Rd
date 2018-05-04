\name{pcp_results}
\alias{pcp_results}
\alias{dah}
\alias{ini}
\alias{nd}
\alias{ndec}
\alias{ne}
\alias{nei}
\alias{nm}
\alias{std}
\alias{x}
\docType{data}
\title{Results of a monthly precipitation homogenization}
\description{
  Results of a homogenization of monthly precipitation totals in 2001-2010.
}
\format{This data set contains:
\describe{
\item{dah}{Homogenized precipitations arranged in a 12x10x6 array}
\item{dat}{Row precipitations arranged in a 12x10x5 array}
\item{est.c}{Data frame with station's coordinates (X,Y,Z), codes and names
(invented)}
\item{ini}{Initial date of the time vector}
\item{nd}{Number of data per series}
\item{ndec}{Number of decimals of the data}
\item{ne}{Number of series after the homogenization}
\item{nei}{Number of series before the homogenization}
\item{nm}{Number of months or number of data per year and station}
\item{std}{Type of data normalization used during the homogenization}
\item{x}{Time vector associated to the data}
}
}
\source{Synthetic precipitation data generated for the \code{dahgrid} example.}
\seealso{\code{\link{dahgrid}}}
\examples{
data(pcp_results)
#station information after the homogenization:
print(est.c)
#raw monthly precipitations of the third station:
print(t(dat[,,3]))
#homogenized monthly precipitations of the third station:
print(t(dah[,,3]))
#range of the dates of the data:
range(x)
#other parameters of the homogenization:
cat('First date:',ini,'\nNumber of data per station:',nd,
'\nNumber of decimals:',ndec,
'\nNumber of series before the homogenization:',nei,
'\nNumber of series after the homogenization:',ne,
'\nNumber of "months" (number of data per year and station):',nm,
'\nType of data normalization:',std,'(normal ratio)\n')
#clean memory space:
rm(est.c,dat,dah,ini,nd,ndec,ne,nei,nm,std,x)
}
\keyword{datasets}

