\name{weekendaccum}
\alias{weekendaccum}
\title{Check for possible precipitation weekend accumulations.}
\description{
  This function is intended to check for possible undeclared daily
  precipitation weekend accumulations, assigning an accumulation code if
  necessary to the false zeroes.}
\usage{
  weekendaccum(varcli,anyi,anyf,na.strings='NA',cumc=-1,wdsl=1)
}
\arguments{
  \item{varcli}{Short name of the climatic variable.}
  \item{anyi}{Initial year of the data.}
  \item{anyf}{Final year of the data.}
  \item{na.strings}{Character strings to be treated as missing data [\code{'NA'}].}
  \item{cumc}{Code of accumulated missing data [-1].}
  \item{wdsl}{Weekday significance level (in \%) to detect 1 to 3 false consecutive zeroes followed by an accumulation of precipitation. Only relevant for daily precipitation [1].}
}
\details{
  This function is intended to be applied to DAILY PRECIPITATION only, whenever
  there is a suspicion that observations of precipitation on an accumulation
  rain gauge (as the Hellmann type) has been skipped on weekends (although the
  test is applied to any day of the week).

  Parameter \code{'wdsl'} needs a value between 0.1 and 10 to enable this test, equivalent to significant levels between 0.001 and 0.1. (Default value: 1, meaning a significant level of 0.01.)

  The test is applied station by station and year by year. When the difference
  between the frequency of zeroes of any week day and the mean of 1, 2 or 3
  preceding days is significant, those preceding days will be assigned the
  chosen cumulative code \code{cumc} (-1 by default) if and only if all the
  precipitation of those preceding days is zero.

  In this way, ulterior application of the \code{homogen} function with the
  same cumulative code will distribute the accumulated precipitation into
  their corresponding days.
}
\seealso{\code{\link{homogen}}}
\examples{
## Set a temporal working directory and write input files:
wd <- tempdir()
wd0 <- setwd(wd)

## Load example daily precipitation data:
data(climatol_data)
dat <- as.matrix(RR3st[,2:4])

## Set zero precipitation in some days and years of two of the three stations:
x <- seq(as.Date('1981-01-01'),as.Date('1995-12-31'),1)
year <- as.integer(strftime(x,'\%Y'))
wkd  <- as.integer(strftime(x,'\%w')) #weekday: 0(Sunday) to 6
dat[year\%in\%1986:1990 & wkd==6,1] <- 0
dat[year\%in\%1989:1992 & wkd\%in\%5:6,2] <- 0
write(dat,'Prec_1981-1995.dat') #write the input data file
## write stations input file:
write.table(SIstations,'Prec_1981-1995.est',row.names=FALSE,col.names=FALSE)

## Now run the function:
weekendaccum('Prec',1981,1995)

## Return to user's working directory:
setwd(wd0)

## Input and output files can be found in directory:
print(wd)
}
\keyword{manip}
