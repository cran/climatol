\name{MHisopleths}
\alias{MHisopleths}
\title{Isopleths on a months-hours diagram}
\description{
  This function takes hourly or subhourly data (spanning at least one
  year) and plots isopleths of the chosen variable in a colored two-dimensional
  (months, hours) diagram.
}
\usage{
MHisopleths(dat, vrb, fun='mean', xlab='Months', ylab='Hours', cex=1.2,
  col4RP=c('cyan','yellow','red'), title='')
}
\arguments{
  \item{dat}{ dataframe containing the data in columns with date/time
         of class POSIX in the first column.}
  \item{vrb}{name of the column containing the chosen data.}
  \item{fun}{function to aggregate subhourly data into hourly.}
  \item{xlab, ylab}{labels for the X and Y axis.}
  \item{cex}{character expansion parameter for the size of labels.}
  \item{col4RP}{vector of colors for the \code{colorRampPalette} function.}
  \item{title}{main title.}
}
\details{
  The user can choose any column of data present in \code{dat}. (Depending on
  the variable the default colors may not be the most appropriate.)
}
\examples{
data(climatol_data)
MHisopleths(AWS_1year,'Temp',title='Mean temperature (C) -- My airport, 2002')
MHisopleths(AWS_1year,'WSpd',title='Wind speed (m/s) -- My airport, 2002')
}
\keyword{hplot}
