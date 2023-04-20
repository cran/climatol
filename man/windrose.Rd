\name{windrose}
\alias{windrose}
\title{Wind-rose plot}
\description{
  This function plots a wind-rose from a data frame with columns DateTime,
  Wind direction and Wind speed.
}
\usage{
windrose(dat, cols=1:3, code='', name='', uni='m/s', maxnsc=8, fnum=4, fint=5,
flab=2, pal=c('cyan','yellow','orange','red','brown'), ang=-3*pi/16,
margin=c(0,0,4,0), ...)
}
\arguments{
  \item{dat}{Data frame with columns DateTime, Wind direction and Wind speed.}
  \item{cols}{Columns containing DateTime, Wind direction and Wind speed [1:3].}
  \item{code}{Station code.}
  \item{name}{Station name.}
  \item{uni}{Speed units for the legend header ['m/s'].}
  \item{maxnsc}{Maximum number of wind speed classes [8].}
  \item{fnum}{Number of reference circles to plot [4].}
  \item{fint}{Frequency interval (in \%) between reference circles [5].}
  \item{flab}{Parameter indicating which circles must be labelled:
    \describe{
      \item{1:}{Label outer circle only,}
      \item{2:}{Label all circles (the default),}
      \item{Other value:}{Do not label any circle.}
    }
  }
  \item{pal}{Color gradation to fill the frequency polygons.}
  \item{ang}{Angle along which circles will be labeled, in radians [\code{-3*pi/16}].}
  \item{margin}{Margins vector for the plot (to be passed to \code{par}) [\code{c(0,0,4,0)}].}
  \item{\dots}{Other graphic parameters.}
}
\details{
  After reading the data, a frequency table is calculated in 16 wind directions
  and a variable number of wind speed classes. The wind direction data must be
  provided in degrees.

  This table, which covers all available pairs of wind direction and speed present in the data frame, is the basis of the wind-rose plot.
}
\value{
  The table of wind frequencies by direction and speed classes is returned invisibly.
}
\examples{
data(climatol_data) #load example data
windtable <- windrose(AWS_1year, 1:3, 'st123', 'My airport') #plot windrose
print(windtable) #display the table of calculated wind frequencies
}
\keyword{datagen}
\keyword{hplot}
