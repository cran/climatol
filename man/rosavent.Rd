\name{rosavent}
\alias{rosavent}
\title{Wind-rose plot}
\description{
  Plot of a wind-rose. Wind frequencies may be supplied by speed intervals,
  and can be absolute or relative.
}
\usage{
  rosavent(frec, fnum = 4, fint = 5, flab = 2, ang = 3 * pi/16,
    col = rainbow(10, 0.5, 0.92, start = 0.33, end = 0.2),
    margen = c(0, 0, 4, 0), key = TRUE, uni = "m/s", ...)
}
\arguments{
  \item{frec}{Wind frequencies.}
  \item{fnum}{Number of reference circumferences to plot.}
  \item{fint}{Frequency steps (in \%) between reference circumferences.}
  \item{flab}{
    Parameter indicating which circumferences must be labelled:
    \describe{
      \item{1:}{Label outer circumference only,}
      \item{2:}{Label all circumferences,}
      \item{Other value:}{Do not label any circumference.}
    }
  }
  \item{ang}{Angle along which circumferences will be labelled.}
  \item{col}{Colors to fill the frequency polygons.}
  \item{margen}{Margins vector for the plot (to be passed to \code{par}).}
  \item{key}{Set to FALSE if you do not want a legend of the wind-rose, that
    will otherwise be plotted if frequencies are supplied by speed intervals.}
  \item{uni}{Speed units for the legend header.}
  \item{\dots}{Other graphic parameters.}
}
\details{
  The number of direction classes is indifferent: As many as elements of the
  frequency vector (or columns of the matrix) will be considered, but the
  first element (column) must correspond to North.
  Frequencies will be converted to percentages, no matter how they were
  supplied.
}
\value{
  No value is returned by this function.
}
\seealso{\code{\link{plot}}, \code{\link{par}}}
\examples{
  data(windfreq.dat)
  rosavent(windfreq.dat,4,4,ang=-3*pi/16,main="Annual windrose")
}
\keyword{hplot}
