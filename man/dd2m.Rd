\name{dd2m}
\alias{dd2m}
\title{Compute monthly data from daily series}
\description{
  Monthly data are computed from daily series stored in
  \file{VAR-d_ANYI-ANYF.dah} and saved in \file{VAR-m_ANYI-ANYF.dah}.
}
\usage{
  dd2m(varcli, anyi, anyf, ini, anyip=anyi, anyfp=anyf, ndec=1,
  valm=2, nmin=15, na.strings="NA")
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
  \item{nmin}{Minimum number of available data in a month to compute the
  monthly value.}
  \item{na.strings}{Missing data code in the original daily data.}
}
\details{
  The daily input data may not begin by the first of January, hence the need
  to indicate the initial date to allow a proper attribution to the correct
  month. Output monthly data has always 12 data per year, using the \code{NA}
  code to provide for any missing data.

  No stations file is output, but the monthly data are save in the same
  station ordering as the daily input file.
}
\value{
  This function does not return any value.
}
\seealso{\code{\link{homogen}}}
\examples{
  #This function only works with files. In order to run the example,
  #you should run the examples of the homogen function first.
  #Once the results from homogen have been generated, you can run:
  \dontrun{dd2m("Tmax", 2007, 2010, ini='2007-01-01')}
  #See the resulting monthly mean data in the generated file
  #'Tmin-m_1971-2000.dah'.
  #To take advantage of all the missing daily data estimations, run:
  \dontrun{dd2m("Tmax", 2007, 2010, ini='2007-01-01',nmin=0)}
  #(The last quarter of 2010 is still NA (not available), because
  #there were no original data it that period).
}
\keyword{datagen}
