\name{tipif}
\alias{tipif}
\title{Series standardization}
\description{
  Standardization of the series by the chosen method.
}
\usage{
  tipif(ttip = 3)
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
  \code{tipif} is an internal function used by \code{depudm}, not intended
  to be called directly by the user.
  It takes the data from object \code{dat.d} previously created.
}
\value{
  Standardized series \code{dat.z}, mean values \code{dat.m} and,
  optionally, standard deviations \code{dat.s} are created.
}
\seealso{\code{\link{depudm}}}
\keyword{internal}
