\name{movttest}
\alias{movttest}
\title{t-test applied to a running window}
\description{
  A t-test of mean sample differences is successively applied to a window
  of \code{N} terms that runs along a series.
}
\usage{
  movttest(x, indme, nterm=10, pval=.05, verb=TRUE)
}
\arguments{
  \item{x}{Time series to be tested}
  \item{indme}{Label to identify the station and month (or sub-series)}
  \item{nterm}{Number of terms of the running window}
  \item{pval}{Significance threshold under which, provided that
    \code{verb=TRUE}, a line indicating the p-value and its position in
    the series will be shown.}
  \item{verb}{If \code{TRUE}, p-values under the threshold \code{pval} will be
    shown.}
}
\details{
  \code{movttest} is an internal function used by \code{depudm}, not intended
  to be called directly by the user.
  It is applied to a running window of \code{N} terms; therefore, the t-test is
  performed on the difference of the means of the first and last \code{N/2}
  terms of the window.
}
\value{
  An object \code{pv} is created, holding the p-values of all the t-test
  applied to the studied series.
}
\seealso{\code{\link{depudm}}, \code{\link{grafanom}}}
\keyword{internal}
\keyword{htest}
\keyword{ts}
