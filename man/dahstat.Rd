\name{dahstat}
\alias{dahstat}
\title{Extract series or statistics of the homogenized data}
\description{
  Lists series, means, medians, standard deviations, quantiles or trends, for a
  specified period, from series homogenized by \code{\link{homogen}}.
}
\usage{
dahstat(varcli, anyi, anyf, anyip=anyi, anyfp=anyf, stat="me", ndc=NA, vala=2,
valm=vala, cod=NULL, prob=.5, all=FALSE, long=FALSE, relref=FALSE, pernyr=10,
estcol=c(1,2,4), sep=',', dec='.')
}
\arguments{
  \item{varcli}{Short name of the studied climatic variable,
    as in the data file name.}
  \item{anyi}{Initial year of the homogenized period.}
  \item{anyf}{Final year of the homogenized period.}
  \item{anyip}{First year of the period to analyze. (Defaults to \code{anyi}).}
  \item{anyfp}{Last year of the period to analyze. (Defaults to \code{anyf}).}
  \item{stat}{Statistical parameter to compute for the selected period:
    \describe{
      \item{"me":}{Means (default),}
      \item{"mdn"}{Medians,}
      \item{"max"}{Maxima,}
      \item{"min"}{Minima,}
      \item{"std"}{Standard deviations,}
      \item{"q"}{Quantiles (see the \code{prob} parameter),}
      \item{"tnd"}{OLS trends and their p-values,}
      \item{"series"}{Do not compute any statistic; only write homogenized
        series and flags into two CSV files.}
      \item{"mseries"}{As before, but output series of homogenized monthly aggregates into a CSV file. (Only applicable to daily or subdaily series.)}
    }
  }
  \item{ndc}{Number of decimal places to be saved in the output file (defaults to that used in the homogenization).}
  \item{vala}{Annual values to compute from the sub-annual data:
    \describe{
      \item{0:}{None,}
      \item{1:}{Sum,}
      \item{2:}{Mean (default),}
      \item{3:}{Maximum,}
      \item{4:}{Minimum.}
    }
  }
  \item{valm}{Monthly values to calculate from sub-monthly data (defaults to
  \code{vala}):
    \describe{
      \item{1:}{Sum,}
      \item{2:}{Mean,}
      \item{3:}{Maximum,}
      \item{4:}{Minimum.}
    }
  }
  \item{cod}{Vector of requested station codes (all by default).}
  \item{prob}{Probability for the computation of quantiles (0.5 by default,
    i.e., medians). You can set probabilities with more than 2 decimals, but
    the name of the output file will be identified with the rounded percentile.}
  \item{all}{If \code{TRUE}, all reconstructed series will be used. The default
    is \code{FALSE}, hence using only the series reconstructed from the last
    homogeneous subperiod.}
  \item{long}{If \code{TRUE} (the default is \code{FALSE}), only series
    reconstructed from the longest homogeneous subperiod will be used.}
  \item{relref}{If \code{TRUE}, statistics from reliable reference series will
  also be listed. (\code{FALSE} by default).}
  \item{pernyr}{Number of years on which to express trend units (10 by default).}
  \item{estcol}{Columns of the homogenized stations file to be included in the
  output file. (Defaults to c(1,2,4), the columns of station coordinates and codes).}
  \item{sep}{Field separator (',' by default).}
  \item{dec}{Decimal point ('.' by default).}
}
\details{
  Homogenized data are read from the file \file{VRB_ANYI-ANYF.rda}
  saved by \code{\link{homogen}}, while this function saves the
  computed data for the specified period in \file{VRB_ANYIP-ANYFP.STAT},
  where \code{STAT} is substituted by the \code{stat} requested
  statistic. An exception is when \code{stat="q"}, since then the
  extension of the output file will be \code{qPP}, where \code{PP}
  stands for the specified \code{prob} probability (in percent).

  The output period \code{ANYIP-ANYFP} must of course be comprised
  within the period of the input data, \code{ANYI-ANYF}.

  \code{stat='tnd'} computes trends by Ordinary Least Squares linear regression
  on time, listing them in a CSV file \file{*_tnd.csv} and their p-values in
  \file{*_pval.csv}

  If \code{stat='series'} is chosen, two text files in CSV format will be
  produced for every station, one with the data and another with their flags: 0
  for original, 1 for infilled and 2 for corrected data.
}
\seealso{\code{\link{homogen}}, \code{\link{dahgrid}}.}
\examples{
## Set a temporal working directory:
wd <- tempdir()
wd0 <- setwd(wd)

## Copy an example file of homogenization results:
file.copy(exampleFiles('Temp_1991-2000.rda'),'.')

## Now run the examples:
dahstat('Temp', 1991, 2000)
dahstat('Temp', 1991, 2000, stat='q', prob=0.4)
dahstat('Temp', 1991, 2000, stat='tnd')
dahstat('Temp', 1991, 2000, stat='series')

## Return to user's working directory:
setwd(wd0)

## Input and output files can be found in directory:
print(wd)
}
\keyword{datagen}
\keyword{utilities}
