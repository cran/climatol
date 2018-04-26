\name{dd2m}
\alias{dd2m}
\title{Compute monthly data from daily series}
\description{
  Daily series are aggregated into total, mean, maximum, or minimum monthly
  values, and saved to files for further processing.
}
\usage{
dd2m(varcli, anyi, anyf, anyip=anyi, anyfp=anyf, ndec=1, suf=NA, valm=2,
namax=10, na.strings="NA", homog=FALSE, ini=NA)
}
\arguments{
  \item{varcli}{Acronym of the name of the studied climatic variable,
    as in the data file name.}
  \item{anyi}{Initial year of the data present in the file.}
  \item{anyf}{Final year of the data present in the file.}
  \item{anyip}{First year of the output period. (Defaults to \code{anyi}).}
  \item{anyfp}{Last year of the output period. (Defaults to \code{anyf}).}
  \item{ndec}{Number of decimal places to be saved in the output file.}
  \item{suf}{Optional suffix to be added (by a hyphen) to the acronym of the
  variable to form the name of the input data.}
  \item{valm}{Monthly value to compute:
    \describe{
      \item{1:}{Sum,}
      \item{2:}{Mean,}
      \item{3:}{Maximum,}
      \item{4:}{Minimum.}
    }
  }
  \item{namax}{Maximum number of missing data in any month to compute its
  monthly value. (10 by default)}
  \item{na.strings}{Missing data code in the original daily data.}
  \item{homog}{If \code{TRUE}, monthly aggregation will be performed on daily
  series adjusted from a previous monthly homogenization.}
  \item{ini}{Initial date of the daily data, with format \code{"YYYY-MM-DD"}.
  Its default \code{NA} value assume it to be the first day of the period of
  study.}
}
\details{
This function can be applied to either raw or homogenized daily data.

Raw data are read from files \file{VAR_YEAR-YEAR.dat} and
\file{VAR_YEAR-YEAR.est}, and monthly data will be saved to
files \file{VAR-m_YEAR-YEAR.dat} and \file{VAR-m_YEAR-YEAR.est}.

Homogenized data are read from \file{VAR_YEAR-YEAR.rda}. In this case
monthly data will be saved to files named \file{VAR-mh_YEAR-YEAR.dat} and
\file{VAR-mh_YEAR-YEAR.est}.
}
\value{
  This function does not return any value.
}
\seealso{\code{\link{homogen}}, \code{\link{dahstat}}, \code{\link{dahgrid}}}
\examples{
#Set a temporal working directory and write input files:
wd <- tempdir()
wd0 <- setwd(wd)
data(Ttest)
write(dat,'Ttest_1981-2000.dat')
write.table(est.c,'Ttest_1981-2000.est',row.names=FALSE,col.names=FALSE)
rm(dat,est.c) #remove loaded data from memory space
#Now run the example:
dd2m('Ttest',1981,2000)
#Return to user's working directory:
setwd(wd0)
#Input and output files can be found in directory:
print(wd)
}
\keyword{datagen}
