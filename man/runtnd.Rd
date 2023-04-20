\name{runtnd}
\alias{runtnd}
\title{Running trends on time windows of different lengths}
\description{
  This function plots running trends on time windows of different lengths in a
  colored grid with axis 'Last year' and 'Window length'.
}
\usage{
  runtnd(d, anyi, minyr=10, units='Units', pernyr=10, stname=NA, k=NULL,
  palneg=c('blue','white'), palpos=c('white','red'), ...)
}
\arguments{
\item{d}{Series of annual values (without missing data).}
\item{anyi}{Initial year of the series.}
\item{units}{Units label for the legend.}
\item{minyr}{Minimum no. of years to compute trends (10 by default).}
\item{pernyr}{Factor for trend units (per 10 years by default).}
\item{stname}{Station name (for the title).}
\item{k}{Vector of breaks for the trend scale colors (automatically set by
default).}
\item{palneg}{Color gradation for negative trends [\code{c('blue','white')}].}
\item{palpos}{Color gradation for positive trends [\code{c('white','red')}].}
\item{\dots}{Additional graphic parameters.}
}
\details{
  The input must be a complete (no missing data) series of annual values.

  If \code{minyr} is negative, running trends calculated on \code{-minyr} years
  will be plotted, with increasing line widths when significance reaches 0.10
  and 0.05 levels. Otherwise, a colored graphic of running trends calculated on
  different window widths will be displayed, masking low significance values
  with white dots.
}
\value{
  A data frame or a list with \code{tnd} (trends) and \code{pvl} (p-values) is
  returned invisibly when \code{minyr} is negative or positive, respectively.
}
\examples{
data(climatol_data)
runtnd(Tav, 1901, -30, units='C', stname='Oslo', cex.axis=1.2, cex.lab=1.2)
runtnd(Tav[31:120], 1931, 30, units='C', stname='Oslo')
}
\keyword{hplot}
