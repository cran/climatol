\name{grabeps}
\alias{grabeps}
\title{Save the displayed graphic as an EPS file}
\description{
  The last displayed graphic is saved in a file named \file{RAAMMDDHHMM-NN.eps}
  (in EPS format).
}
\usage{
  grabeps()
}
\details{
  \code{grabeps} is a function called by \code{grafanom}, but it can
  be called directly by the user.
  The name of the file saved begins by R, followed by digits indicating
  the year, month, day, hour and minute (two digits each) of the first
  call of this function in the present R session, plus a rank number
  (two digits also).
}
\value{
  No value is returned by this function.
}
\note{The last graphic will not be completely dumped to the file until the
  session ends.}
\seealso{\code{\link{grafanom}}}
\keyword{utilities}
