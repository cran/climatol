\name{destipif}
\alias{destipif}
\title{Series unstandardization}
\description{
  Undo the standardization performed by \code{tipif}.
}
\usage{
  destipif(ttip = 3)
}
\arguments{
  \item{ttip}{Type of standardization:
    \describe{
      \item{0:}{none,}
      \item{1:}{deviations from the mean,}
      \item{2:}{proportions of the mean (only for means greater than 1),}
      \item{3:}{full standardization (subtract the mean and divide by the
        standard deviation).}
    }
  }
}
\details{
  \code{destipif} is an internal function used by \code{depudm}, not intended
  to be called directly by the user.
  It takes the data from object \code{dat.d} previously created.
}
\value{
  Previously created object \code{dat.e} is unstandardized, using mean values
  (\code{dat.m}) and, eventually, standard deviations (\code{dat.s}).
}
\seealso{\code{\link{depudm}}}
\keyword{internal}
