\name{dd2m}
\alias{dd2m}
\title{Compute monthly data from daily series}
\description{
  Monthly data are computed from daily series stored in
  \file{VAR-d_ANYI-ANYF.dah} and saved in \file{VAR-m_ANYI-ANYF.dah}.
}
\usage{
dd2m(varcli, anyi, anyf, ini, anyip=anyi, anyfp=anyf, ndec=1,
valm=2, nmin=15, na.strings="NA", dah=TRUE)
}
\arguments{
  \item{varcli}{Acronym of the name of the studied climatic variable,
    as in the data file name.}
  \item{anyi}{Initial year of the data present in the file.}
  \item{anyf}{Final year of the data present in the file.}
  \item{ini}{Initial date of the daily data, with format \code{"YYYY-MM-DD"}}.
  \item{anyip}{First year of the desired period.}
  \item{anyfp}{Last year of the desired period.}
  \item{ndec}{Number of decimal places to be saved in the output file.}
  \item{valm}{Monthly value to compute:
    \describe{
      \item{1:}{Sum,}
      \item{2:}{Mean,}
      \item{3:}{Maximum,}
      \item{4:}{Minimum.}
    }
  }
  \item{nmin}{Minimum number of daily available data in a month to compute its
  monthly value.}
  \item{na.strings}{Missing data code in the original daily data.}
  \item{dah}{Logical value indicating if input daily data are homogenized.
    (Defaults to \code{TRUE}.}
}
\details{
  The daily input data may not begin by the first of January, hence the need
  to indicate the initial date to allow a proper attribution to the correct
  month. Output monthly data has always 12 data per year, using the \code{NA}
  code to provide for any missing data.

  Input daily data are read from files with a base name
  \file{varcli-d_anyi-anyf} (note the \code{-d} appended to the name of the
  variable to indicate that they are daily data), and extensions \code{.dah}
  and \code{.esh} for data and station coordinates if you want to compute
  monthly data from homogenized daily data, or \code{.dat} and \code{.est} if
  they are raw data. (In this latter case, \code{dah=FALSE} must be specified).
}
\value{
  This function does not return any value.
}
\seealso{\code{\link{homogen}}}
\examples{
#This function only works with files. In order to run the example,
#you should run the examples of the homogen function first.
#Once the results from homogen have been generated, you can run:
#
#   dd2m("Tmax", 2007, 2010, ini='2007-01-01')
#
#See the resulting monthly mean data in the generated file
#'Tmin-m_1971-2000.dah'.
#To take advantage of all the missing daily data estimations, run:
#
#   dd2m("Tmax", 2007, 2010, ini='2007-01-01',nmin=0)
#
#(The last quarter of 2010 is still NA (not available), because
#there were no original data it that period).
}
\keyword{datagen}
