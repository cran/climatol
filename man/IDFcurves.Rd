\name{IDFcurves}
\alias{IDFcurves}
\title{Obtain Intensity-Duration-Frequency curves}
\description{
  Intensity-Duration-Frequency curves are obtained from a sub-hourly time series of precipitation by adjusting Generalized Extreme Value distribution to annual maxima of different time intervals.
}
\usage{
IDFcurves(prdat, stname, clmn=1:2, tz='utc', na.code=NA,
prunits='mm', mindpy=0.8, gumbel=TRUE, timeaggr=c(10,20,30,60,120,180,360,720),
retper=c(5,10,20,30,50,75,100),...)
}
\arguments{
\item{prdat}{Data frame with Time (as POSIXct) and sub-hourly precipitation data.}
\item{stname}{Station name.}
\item{clmn}{Columns where Time and precipitation data are located in \code{prdat}.}
\item{tz}{Time zone [\code{'utc'} by default].}
\item{na.code}{Numeric missing data code.}
\item{prunits}{Precipitation units [mm].}
\item{mindpy}{Minimum available data proportion to process data in any year.}
\item{gumbel}{Adjust a Gumbel distribution? [\code{TRUE}].}
\item{timeaggr}{Time intervals (in minutes) on which to aggregate precipitation.}
\item{retper}{Return periods (in years) for extreme precipitation estimation.}
\item{\dots}{Additional graphic parameters.}
}
\details{
  The precipitation time series must be provided as a data frame with POSIXct
  times in the first column and precipitation in the second. However, these
  data can be in other columns of a wider data frame if the columns containing
  these variables are defined in the parameter \code{clmn}.

  When setting \code{gumbel=FALSE} a Generalized Extreme Value distribution will be adjusted instead of the particular case of a Gumbel distribution.
}
\value{
  A table of maximum precipitation accumulations (totals, not mm/h as in the graphic) is returned invisibly.
}
\examples{
\dontrun{
data(climatol_data)
tab <- IDFcurves(prec10min,'My airport',cex.axis=1.2,cex.lab=1.2) #IDF plot

## See the maximum precipitation accumulations in the different time intervals:
tab
}
}
\keyword{datagen}
\keyword{hplot}
