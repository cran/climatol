\name{meteogram}
\alias{meteogram}
\title{Daily meteogram of eight meteorological variables}
\description{
  This function plots a meteogram from hourly or sub-hourly data of eight
  meteorological variables available in a data frame spanning one day.
}
\usage{
  meteogram(df, code='', name='', cols=1:9, tz='utc', hlab='Hours', 
  datefm='\%Y-\%m-\%d', vlab=c('Wind direction (deg)','Wind speed (m/s)',NA,NA,
  'Temperature (C)','Rel. humidity (\%)','Precip. (mm)','Pressure (hPa)'),
  vcol=c(hsv(.1,1,.9),hsv(.1,1,.9),2,2,2,hsv(.4,1,.7),4,'brown'),
  llim=c(0,0,NA,NA,0,0,0,NA), ulim=c(360,20,NA,NA,20,100,4,NA))
}
\arguments{
\item{df}{Data frame with (around) one day of data.}
\item{code}{Code of the station.}
\item{name}{Name of the station.}
\item{cols}{Column order of the expected variables (see details).}
\item{tz}{Time zone of the supplied time vector (\code{'utc'} by default).}
\item{hlab}{Label for hours (\code{'Hours'} by default).}
\item{datefm}{Date format for the title of the meteogram (the default is
  \code{'\%Y-\%m-\%d'}, the ISO 8601 date format).}
\item{vlab}{Variable labels.}
\item{vcol}{Colors for every variable.}
\item{llim}{Lower graphic limits (if fixed).}
\item{ulim}{Upper graphic limits (if fixed).}
}
\details{
  This function expects a data frame containing observation time
  and eight meteorological variables in this column order:
  \enumerate{
    \item{Time of the observation (as POSIXct)}
    \item{10 minutes average wind direction in degrees}
    \item{10 minutes average wind speed in m/s}
    \item{3 sec. maximum gust direction in degrees}
    \item{3 sec. maximum gust speed in m/s}
    \item{Air temperature in degrees Celsius}
    \item{Relative humidity in \%}
    \item{Precipitation in mm}
    \item{Barometric pressure in hPa}
  }

  However, if the data frame has these variables in a different order, it can
  be specified with the parameter \code{cols}.

  See \code{\link{strftime}} for ways to specify date formats.
}
\seealso{\code{\link{strftime}}}
\examples{
data(climatol_data)
meteogram(AWS_1day, 'S123', 'My airport')
}
\keyword{hplot}
