\name{datest}
\alias{datest}
\title{Estimate series from their neighbors}
\description{
  Data of each series are estimated from all the other stations,
  through weighting averages of their data, standardized by the chosen method.
}
\usage{
  datest(verb = FALSE)
}
\arguments{
  \item{verb}{Verbose information about the progress of the operation}
}
\details{
  \code{datest} is an internal function used by \code{depudm}, not intended
  to be called directly by the user.
  It takes the data of the normalized series from object \code{dat.d}
  previously generated, and computes estimated data (dat.e) as averages,
  weighted by the weights est.w, of all the other stations.
}
\value{
  An object \code{dat.e}, with the estimated (normalized) data, is generated.
}
\seealso{\code{\link{depudm}}}
\keyword{internal}
