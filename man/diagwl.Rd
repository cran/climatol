\name{diagwl}
\alias{diagwl}
\title{Walter & Lieth climatic diagram}
\description{
  Plot a Walter & Lieth climatic diagram of a station.
}
\usage{
diagwl(dat, cols=1:6, format='\%Y-\%m-\%d', yeari=NA, yearf=NA,
stname='', alt=NA, per='', mlab='', shem=FALSE, p3line=FALSE, ...)
}
\arguments{
  \item{dat}{Data frame with the required climatic data (see details).}
  \item{cols}{Columns containing the needed data.}
  \item{format}{Format of the dates if data are provided in 4 columns ['\%Y-\%m-\%d'].}
  \item{yeari, yearf}{If \code{dat} is a file name, initial and final years of
  the period to use (defaults to the period contained in the data file).}
  \item{stname}{Name of the climatological station.}
  \item{alt}{Elevation (altitude) of the climatological station.}
  \item{per}{If data is a data frame with already calculated climate averages,
  the original period of the data.}
  \item{mlab}{Vector of 12 monthly labels for the X axis (see the details).}
  \item{shem}{Southern hemisphere? (\code{FALSE} by default.)}
  \item{p3line}{Draw a supplementary precipitation line referenced to three
  times the temperature? (\code{FALSE} by default.)}
  \item{\dots}{Other optional graphic parameters.}
}
\details{
  The data frame can contain daily data of precipitation and extreme temperatures or 12 columns with pre-calculated monthly climate parameters.

  In the latter case the monthly values from January to December occupy the 12
  first columns (additional columns may be present, although they will be
  disregarded) and four rows, in the following order:
  \describe{
    \item{Row 1:}{Mean precipitation}
    \item{Row 2:}{Mean maximum daily temperature}
    \item{Row 3:}{Mean minimum daily temperature}
    \item{Row 4:}{Absolute monthly minimum temperature}
  }
  This last row is only used to determine the probable frost months
  (when absolute monthly minimums are equal or lower than 0 C).
  It is very important to set \code{cols=NULL} in this case.

  Alternatively, if daily data of precipitation and extreme temperatures are provided, dates can be given in three separate columns (year, month, day) or in a single column with the specified \code{format} (\code{'\%Y-\%m-\%d'} by default).

  \code{cols} indicate in which columns are located the dates and climatic
  data. By default they are expected in columns 1 to 3 for year, month and day,
  and columns 4 to 6 for precipitation, maximum and minimum temperature
  respectively. (Remember to set \code{cols=NULL} if you provide pre-calculated monthly climatic parameters.)

  \code{mlab} is the vector for the 12 monthly labels, but it may be set to just \code{'en'} or \code{'es'} to use the first letter of the month names in English or Spanish respectively.

  For stations located in the southern hemisphere it is useful to set
  \code{shem=TRUE}, in order to keep the summer period in the central zone of
  the graphic (the diagram will begin the plot with the July data).

  As described by Walter and Lieth, when monthly precipitation is greater than
  100 mm, the scale is increased from 2 mm/C to 20 mm/C to avoid too high
  diagrams in very wet locations. This change is indicated by a black
  horizontal line, and the graph over it is filled in solid blue.

  When the precipitation graph lies under the temperature graph (P < 2T) we
  have an arid period (filled in dotted red vertical lines). Otherwise the
  period is considered wet (filled in blue lines), unless \code{p3line=TRUE},
  that draws a precipitation black line with a scale P = 3T; in this case
  the period in which 3T > P > 2T is considered semi-arid. (Parameter
  \code{p3line} was suggested by Bogdan Rosca.)

  Daily maximum average temperature of the hottest month and daily minimum
  average temperature of the coldest month are frequently used in vegetation
  studies, and are labeled in black at the left margin of the diagram.
}
\references{Walter H & Lieth H (1960): Klimadiagramm Weltatlas. G. Fischer,
  Jena.}
\examples{
data(climatol_data)

## from pre-calculated monthly climatic data:
diagwl(datcli,cols=NULL,est="My Airport",alt=100,per="1961-90",mlab="en")

## from daily series of precipitation and extreme temperatures:
diagwl(p064.df, stname="Cold Place", alt=100, per="1961-1990", mlab="en")
}
\keyword{hplot}
