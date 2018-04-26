\name{dahstat}
\alias{dahstat}
\title{Statistical summaries of the homogenized data}
\description{
  Lists means, standard deviations, quantiles or trends, for a
  specified period, from series homogenized by \code{\link{homogen}}.
}
\usage{
dahstat(varcli, anyi, anyf, anyip=anyi, anyfp=anyf, stat="me", ndc=NA, vala=2,
cod=NULL, mnpd=0, mxsh=0, prob=.5, last=FALSE, long=FALSE, lsnh=FALSE,
lerr=FALSE, relref=FALSE, mh=FALSE, pernys=100, estcol=c(1,2,4), sep=',',
dec='.', eol='\n', nei=NA, x=NA)
}
\arguments{
  \item{varcli}{Acronym of the name of the studied climatic variable,
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
      \item{"tnd"}{Trends,}
      \item{"series"}{Do not compute any statistics; only output all homogenized
        series in individual *.csv files.}
    }
  }
  \item{ndc}{Number of decimal places to be saved in the output file (1 by default).}
  \item{vala}{Annual values to compute from the sub-annual data:
    \describe{
      \item{0:}{None,}
      \item{1:}{Sums,}
      \item{2:}{Means (default),}
      \item{3:}{Maxima,}
      \item{4:}{Minima.}
    }
  }
  \item{cod}{Optional vector of codes of the stations to be processed.}
  \item{mnpd}{Minimum percentage of original data. (0 = no limit).}
  \item{mxsh}{Maximum SNHT. (0 = no limit).}
  \item{prob}{Probability for the computation of quantiles (0.5 by default,
    i.e., medians). You can set probabilities with more than 2 decimals, but
    the name of the output file will be identified with the rounded percentile.}
  \item{last}{Logical value to compute statistics only for stations
    working at the end of the period of study. (\code{FALSE} by default).}
  \item{long}{Logical value to compute statistics only for series built from
    the longest homogeneous sub-period. (\code{FALSE} by default).}
  \item{lsnh}{Logical value to compute statistics from series built from
    the homogeneous sub-period with lowest SNHT. (\code{FALSE} by default).}
  \item{lerr}{Logical value to compute statistics only for series built from
    the homogeneous sub-period with lowest RMSE. (\code{FALSE} by default).}
  \item{relref}{If \code{TRUE}, statistics from reliable reference series will be also listed. (\code{FALSE} by default).}
  \item{mh}{If \code{TRUE}, read monthly data computed from daily adjusted series. (\code{FALSE} by default).}
  \item{pernys}{Number of years on which to compute trends. (Defaults to 100).}
  \item{estcol}{Columns of the homogenized stations file to be included in the
  output file. (Defaults to c(1,2,4), the columns of station coordinates and codes).}
  \item{sep}{String to use for separating the output data. (',').}
  \item{dec}{Character to use as decimal point in the output data. ('.').}
  \item{eol}{Line termination style. ('\\n').}
  \item{nei}{Number of stations in the input files. (To be read from the *.rda file.)}
  \item{x}{Vector of dates. (To be read from the *.rda file.)}
}
\details{
  Homogenized data are read from the file \file{VAR_ANYI-ANYF.rda}
  saved by \code{\link{homogen}}, while this function saves the
  computed data for the specified period in \file{VAR_ANYIP-ANYFP.STAT},
  where \code{STAT} is substituted by the \code{stat} requested
  statistic. An exception is when \code{stat="q"}, since then the
  extension of the output file will be \code{qPP}, where \code{PP}
  stands for the specified \code{prob} probability (in percent).
  The output period \code{ANYIP-ANYFP} must of course be comprised
  within the period of the input data, \code{ANYI-ANYF}.

  Parameters \code{mnpd} and \code{mxsh} act as filters
  to produce results only for series that have those minimum percentages
  of original data and maximum SNHT values. Alternatively, \code{long},
  \code{last}, \code{lsnh} and \code{lerr} allow the selection of series reconstructed from the preferred homogeneous sub-period, depending on the parameter set to \code{TRUE}. However, note that in many cases the shorter sub-periods may have lower SNHT and RMSE values, and therefore parameters \code{lsnh} and \code{lerr} should be used with caution. The most advisable paramenters to select most suitable reconstructions are \code{long} for computing normal values and \code{last} for climate monitoring of new incoming data.

 to select only those stations
  working at the end of the period studied. No selection is performed
  by default, listing the desired statistic for all the reconstructed
  series (from every homogeneous sub-period).

  \code{stat='tnd'} computes trends by OLS linear regression on time, listing
  them in a CSV file \file{*_tnd.csv} and their p-values in \file{*_pval.csv}

  If \code{stat='series'} is chosen, two text files in CSV format will be
  produced for every station, one with the data and another with their flags: 0
  for original, 1 for infilled and 2 for corrected data. (Not useful for daily series.)
}
\value{
  This function does not return any value, since outputs are saved to files.
}
\seealso{\code{\link{homogen}}, \code{\link{dahgrid}}.}
\examples{
#Set a temporal working directory and write input files:
wd <- tempdir()
wd0 <- setwd(wd)
data(Ptest)
dim(dat) <- c(720,20)
dat[601:720,5] <- dat[601:720,5]*1.8
write(dat[481:720,1:5],'pcp_1991-2010.dat')
write.table(est.c[1:5,1:5],'pcp_1991-2010.est',row.names=FALSE,col.names=FALSE)
homogen('pcp',1991,2010,std=2)
#Now run the examples:
dahstat('pcp',1991,2010)
dahstat('pcp',1991,2010,stat='tnd')
#Return to user's working directory:
setwd(wd0)
#Input and output files can be found in directory:
print(wd)
}
\keyword{datagen}
