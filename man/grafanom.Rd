\name{grafanom}
\alias{grafanom}
\title{Plot of anomalies, and original vs. estimated data}
\description{
  Plot of the standardized anomalies (original minus estimated data) and,
  optionally, original and estimated data (unstandardized).
}
\usage{
  grafanom(me=me, labm="", ttip=3, ndec=1, sqrtrans=FALSE, pval=.05)
}
\arguments{
  \item{me}{month (or sub-series number, if not monthly) to plot.}
  \item{labm}{month (or sub-series) label.}
  \item{ttip}{Type of standardization (to undo it for plotting the series).}
  \item{ndec}{Number of decimal places of the standard deviation label.}
  \item{sqrtrans}{Set to \code{TRUE} if a square root transformation was
    applied to the data (to undo it for plotting the series).}
  \item{pval}{If greater than 0, t-test of mean difference will be applied
    to moving windows of 10 and 20 terms (between samples of 5 and 10 terms),
    and an overall trend test. If \code{graf=TRUE}, p-values of the t-test
    will be plotted and, provided that the trend is significant (lower than
    \code{pval}), the regression line with time will be plotted as well.}
}
\details{
  \code{grafanom} is an internal function used by \code{depudm}, not intended
  to be called directly by the user.
  Anomalies are plotted centered on their own average.
  After displaying the plot, the user is prompted to enter the number of
  the station to jump to (RETURN will jump to the next), or to view the
  plot of the original and estimated series, save the graphic, or end.
  Graphics are saved in EPS format, in files named \file{RAAMMDDHHMM-NN.eps}:
  an initial R followed by digits indicating the year, month, day, hour
  and minute (two digits each) of the first call of this function in the
  present R session, plus a rank number (two digits also).
}
\value{
  No value is returned by this function.
}
\note{
  The last graphic will not be completely dumped to the file until the session
  ends.
}
\seealso{\code{\link{depudm}}, \code{\link{grabeps}}}
\keyword{internal}
\keyword{hplot}
