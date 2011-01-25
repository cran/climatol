\name{diagwl}
\alias{diagwl}
\title{Walter \& Lieth climatic diagram}
\description{
  Plot of a Walter \& Lieth climatic diagram of a station.
}
\usage{
diagwl(dat, est = "", alt = NA, per = "", margen = c(4, 4, 5, 4), mlab = "",
  pcol = "#005ac8", tcol = "#e81800", pfcol = "#79e6e8", sfcol = "#09a0d1",
  shem = FALSE, p3line=FALSE, ...)
}
\arguments{
  \item{dat}{Monthly climatic data for which the diagram will be plotted.}
  \item{est}{Name of the climatological station}
  \item{alt}{Altitude of the climatological station}
  \item{per}{Period on which the averages have been computed}
  \item{margen}{Margins vector for the plot (to be passed to \code{par}).}
  \item{mlab}{Month labels for the X axis:
    \describe{
      \item{"en":}{Month initials in English.}
      \item{"es":}{Month initials in Spanish.}
      \item{Other:}{Numeric labels (1-12).}
    }
  }
  \item{pcol}{Color pen for precipitation.}
  \item{tcol}{Color pen for temperature.}
  \item{pfcol}{Fill color for probable frosts.}
  \item{sfcol}{Fill color for sure frosts.}
  \item{shem}{Set to \code{TRUE} for southern hemisphere stations.}
  \item{p3line}{Set to \code{TRUE} to draw a suplementary precipitation line referenced to three times the temperature (as suggested by Bogdan Rosca).}
  \item{\dots}{Other graphic parameters}
}
\details{
  Climatic data must be passed as a 4x12 matrix of monthly (January to
  December) data, in the following order:
  \describe{
    \item{Row 1:}{Mean precipitation}
    \item{Row 2:}{Mean maximum daily temperature}
    \item{Row 3:}{Mean minimum daily temperature}
    \item{Row 4:}{Absolute monthly minimum temperature}
  }
  This last row is only used to determine the probable frost months
  (when absolute monthly minimums are equal or lower than 0°C).

  If \code{shem=TRUE}, the diagram will begin the plot with the July data,
  to keep the summer period in the central zone of the graphic.
}
\value{
  No value is returned by this function.
}
\references{WALTER H & LIETH H (1960): Klimadiagramm Weltatlas. G. Fischer,
  Jena.}
\seealso{\code{\link{plot}}, \code{\link{par}}}
\examples{
  data(datcli)
  diagwl(datcli,est="Example station",alt=100,per="1961-90",mlab="en")
}
\keyword{hplot}
