\name{Datasets}
\alias{Datasets}
\alias{climatol_data}
\alias{AWS_1day}
\alias{AWS_1year}
\alias{datcli}
\alias{p064.df}
\alias{p082.df}
\alias{p084.df}
\alias{prec10min}
\alias{RR3st}
\alias{SIstations}
\alias{Tav}
\alias{Temp.dat}
\alias{Temp.est}
\alias{TN3st}
\alias{TX3st}
\title{Data sets to run examples of the functions in the \code{climatol} package.}
\description{This object contains several small datasets needed to run the examples of most functions.}
\usage{
  data(climatol_data)
}
\details{
This data set holds the following collection of data objects:
\describe{
  \item{TX3st}{Maximum daily temperature of 3 stations during 1981-1995.}
  \item{TN3st}{Minimum daily temperature of 3 stations during 1981-1995.}
  \item{RR3st}{Daily precipitation of 3 stations during 1981-1995.}
  \item{SIstations}{Stations coordinates, codes and names of the \code{*3st} data.}
  \item{p064.df}{Data frame with RR, TX and TN data for station p064.}
  \item{p084.df}{Data frame with RR, TX and TN data for station p084.}
  \item{p082.df}{Data frame with RR, TX and TN data for station p082.}
  \item{AWS_1year}{Hourly data from an Automatic Weather Station during one year.}
  \item{AWS_1day}{10 minutes data from an Automatic Weather Station during one day.}
  \item{datcli}{Monthly climatic parameters to plot a Walter&Lieth diagram.}
  \item{Temp.dat}{Monthly temperature of five stations during 1961-2005.}
  \item{Temp.est}{Stations coordinates, codes and names of the \code{Temp.dat} data.}
  \item{Tav}{Annual average temperature at Oslo (Norway) during 1901-2020.}
  \item{prec10min}{Ten minutes precipitation data during 1991-2020.}
}

Some examples need the use of files rather than these data objects. In that case they are provided in a special folder of the installation directory tree and are made accessible through the function \code{\link{exampleFiles}}.
}
\source{
\code{RR3st, TX3st, TN3st, p064.df, p082.df and p084.df} data were obtained
from the historical run (1950-2005) of the Regional Atmospheric Climate Model
version 2 of the Royal Netherlands Meteorological Institute (KNMI) in the frame of the INDECIS project <https://indecis.eu>.

Oslo annual average temperatures \code{Tav} were downloaded from the HCLIM
database.

The other objects contain real data, but anonimized to avoid data policy
restrictions.
}
\seealso{\code{\link{exampleFiles}}}
\examples{
data(climatol_data)
datcli
head(p064.df)
head(AWS_1year)
}
\references{
Lundstad, Elin; Brugnara, Yuri; Broennimann, Stefan (2022): Global Early Instrumental Monthly Meteorological Multivariable Database (HCLIM). PANGAEA, https://doi.org/10.1594/PANGAEA.940724
}
\keyword{datasets}
