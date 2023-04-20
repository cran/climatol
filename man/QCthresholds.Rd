\name{QCthresholds}
\alias{QCthresholds}
\title{Obtain monthly thresholds for Quality Control alerts}
\description{
  This function calculate monthly quantiles of daily or subdaily series that can be used as thresholds for Quality Control alerts.
}
\usage{
  QCthresholds(dat, ndec=1, probs=c(0.,.001,.01,.99,.999,1.), minval=NA,
  maxval=NA, homog=TRUE, verb=TRUE)
}
\arguments{
\item{dat}{Either the name of a *.rda file of \code{climatol} homogenization results or a data.frame of daily (or subdaily) data in columns, dates or date/times (of class Date or POSIXct) in the first column and station codes in the header}
\item{ndec}{number of decimals of output values [1] (defaults shown between brackets)}
\item{probs}{probabilities of the quantiles to be computed [0., .001, .01, .99, .999, 1.]}
\item{minval}{minimum value to compute runs of constant values [\code{NA}].}
\item{maxval}{maximum value to compute runs of constant values [\code{NA}].}
\item{homog}{use homogenized data if a *.rda file is used as input [\code{TRUE}].}
\item{verb}{list all calculated values? [\code{TRUE}].}
}
\details{
\code{minval} and \code{maxval} allow to exclude frequent values that would
result in the report of long runs of identical data. Examples: set
\code{minval=0.1} in daily precipitation to avoid long runs of zeros or set
\code{maxval=97} in relative humidity to avoid long runs of near saturation
values in episodes of persistent fog.

Calculated thresholds are shown in the text output and are also saved in a
binary R file named \code{QCthresholds.Rdat}, which contains the matrices
\code{thr1}, \code{thr2} and \code{thr3}. Load this file and write the
thresholds in the required format for importation into a Climate Data
Management System.
}
\seealso{\code{\link{homogen}}}
\examples{
## Set a temporal working directory and write input files:
wd <- tempdir()
wd0 <- setwd(wd)
data(climatol_data)

## Now run the examples:
QCthresholds(RR3st,minval=0.1) #daily precipitation of three stations
QCthresholds(TX3st) #daily maximum temperatures of three stations
load('QCthresholds.Rdat') #load last calculated thresholds
thr1[1,,] #thresholds with 0% probability to find lower values
thr1[,3,] #monthly thresholds of the third station
thr2 #thresholds of absolute increments between consecutive data
thr3 #thresholds for equal data run lengths

## Return to user's working directory:
setwd(wd0)

## Input and output files can be found in directory:
print(wd)
}
\keyword{datagen}
