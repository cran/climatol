\name{leerdm}
\alias{leerdm}
\title{Read the (monthly) climatological series}
\description{
  Climatological series, and coordinates and names of the stations, are
  read from \file{VAR\_AI-AF.dat} and \file{VAR\_AI-AF.est} files respectively.
}
\usage{
  leerdm(varcli, anyi, anyf, nm=12, b=1, a=0)
}
\arguments{
  \item{varcli}{Acronym of the name of the studied climatic variable,
    as in the data file name.}
  \item{anyi}{Initial year of the data present in the file (four digits)}
  \item{anyf}{Final year of the data present in the file (four digits)}
  \item{nm}{Number of series in each station. (Default=12, for monthly data).}
  \item{a, b}{Parameters of the optional transformation a+b*dat to be applied
    to data when read from the files.}
}
\details{
  \code{leerdm} is an internal function used by \code{depudm}, not intended
  to be called directly by the user.
}
\value{
  Objects \code{dat} (series data) and \code{est.c} (coordinates and names of
  the stations) are created, and number of years (\code{na}), number of
  stations (\code{ne}) and per station number of data (\code{nd}) are
  also computed.
}
\seealso{\code{\link{depudm}}}
\keyword{internal}
