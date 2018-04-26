\name{homogsplit}
\alias{homogsplit}
\title{Apply homogen() on overlapping rectangular areas.}
\description{
  If the number of series is too big to be homogenized all at a time (normally several thousands, depending on the computer resources), this function can homogenize them by splitting the geographical domain in overlapping rectangular areas.
}
\usage{
homogsplit(varcli, anyi, anyf, xc=NULL, yc=NULL, xo=.5, yo=.38, maponly=FALSE,
suf=NA, nm=NA, nref=c(10,10,4), swa=NA, std=3, ndec=1, dz.max=5,
dz.min=-dz.max, wd=c(0,0,100), snht1=25, snht2=snht1, tol=.02, maxdif=NA,
mxdif=maxdif, force=FALSE, wz=.001, trf=0, mndat=NA, gp=3, ini=NA,
na.strings="NA", maxite=999, vmin=NA, vmax=NA, nclust=100, grdcol=grey(.4),
mapcol=grey(.4), hires=TRUE, expl=FALSE, metad=FALSE, sufbrk='m', tinc=NA,
tz='UTC', cex=1.2, verb=TRUE, x=NA)
}
\arguments{
\item{varcli}{Acronym of the name of the studied climatic variable,
  as in the data file name.}
\item{anyi}{Initial year of the data present in the file.}
\item{anyf}{Final year of the data present in the file.}
\item{xc}{Vector of X axis coordinates setting the domain splitting meridians.}
\item{yc}{Vector of Y axis coordinates setting the domain splitting parallels.}
\item{xo}{Overlapping width in the East-West direction.}
\item{yo}{Overlapping width in the North-South direction.}
\item{maponly}{Do not homogenize. Only draw a map with stations locations and domain partitioning.}
\item{suf}{Optional suffix appended with a '-' to the name of the variable
in the input files.}
\item{nm}{Number of data per year in each station. (Defaults to NA, and then it will be computed from the total number of data).}
\item{nref}{Maximum number of references for data estimation. (Defaults to 10 in the detection stages, and to 4 in the final series adjustments).}
\item{swa}{Size of the step forward to be applied to the staggered window
application of SNHT. If not set (the default), 365 terms (one year) will be
used for daily data, and 60 otherwise.}
\item{std}{Type of normalization:
  \describe{
    \item{1:}{deviations from the mean,}
    \item{2:}{rates to the mean (only for means greater than 1),}
    \item{3:}{standardization (subtract the mean and divide by the sample
      standard deviation).}
  }
}
\item{ndec}{Number of decimal digits to which the homogenized data must be
rounded.}
\item{dz.max}{Threshold of outlier tolerance, in standard deviations.  (5 by
default).}
\item{dz.min}{Lower threshold of outlier tolerance if different from the higher one. (By default, they will be the same, with opposite signs).}
\item{wd}{Distance (in km) at which reference data will weigh half of that
of another located at zero distance. (Defaults to \code{c(0,0,100)}, meaning
that no weighting will be applied in the first two stages, and 100 km in the
third).}
\item{snht1}{Threshold value for the stepped SNHT window test applied in stage
1. (25 by default. No SNHT analysis will be performed if \code{snht1=0}).}
\item{snht2}{Threshold value for the SNHT test when applied to the complete
series in stage 2 (same value as snht1 by default).}
\item{tol}{Tolerance factor to split several series at a time. The default is 0.02, meaning that a 2\% will be allowed for every reference data. (E.g.: if the maximum SNHT test value in a series is 30 and 10 references were used to compute the anomalies, the series will be split if the maximum test of the reference series is lower than 30*(1+0.02*10)=36. Set \code{tol=0} to disable further splits when any reference series has already been split at the same iteration).}
\item{maxdif}{Maximum difference of any data item in consecutive iterations. If not set, defaults to half of the data precision (defined by the number of decimals).}
\item{mxdif}{Old maxdif parameter (maintained for compatilibility).}
\item{force}{Force break even when only one reference is available.
(\code{FALSE} by default).}
\item{wz}{Scale parameter of the vertical coordinate \code{Z}. 0.001 by default, which gives the vertical coordinate (in m) the same weight as the horizontal coordinates (internally managed in km).}
\item{trf}{By default, data are not transformed (\code{trf=0}), but if the data
frequency distribution is very skewed, the user can choose to apply a log(x+1)
transformation (\code{trf=1}) or any root of index \code{trf>1} (2 for square
root, 3 for cubic root, etc. Fractional numbers are allowed).}
\item{mndat}{Minimum number of data for a split fragment to become a new
series. It defaults to half of the \code{swa} value for daily data, or to \code{nm} otherwise, with a minimum of 5 terms.}
\item{gp}{Graphic parameter:
  \describe{
    \item{0:}{no graphic output,}
    \item{1:}{only descriptive graphics of the input data,}
    \item{2:}{as with 1, plus diagnostic graphics of anomalies,}
    \item{3:}{as with 2, plus graphics of running annual means and
    applied corrections,}
    \item{4:}{as with 3, but running annual totals (instead of means)
    will be plotted in the last set of graphics. (Better when working with
    precipitation data).}
  }
}
\item{ini}{Initial date, with format \code{'YYYY-MM-DD'}. If not set, it will be assumed that the series begin the first of January of the initial year \code{anyi}.}
\item{na.strings}{Character string to be treated as a missing value.  (It can
be a vector of strings, if more than one is needed). Defaults to 'NA', the
standard missing data code in R.}
\item{maxite}{Maximum number of iterations when computing the means of the
series. (999 by default).}
\item{vmin}{Minimum possible value (lower limit) of the studied variable.
Unset by default, but note that \code{vmin=0} will be applied if \code{std} is
set to 2.}
\item{vmax}{Maximum possible value (upper limit) of the studied variable.
(E.g., for relative humidity or relative sunshine hours it is advisable to set
\code{vmax=100}).}
\item{nclust}{Maximum number of stations for the cluster analysis. (If much
greater than 100, the default value, the process may be too long and the
graphic too dense).}
\item{grdcol}{Color of the graphic background grids. (Gray by default.)}
\item{mapcol}{Color of the background map. (Gray by default).}
\item{hires}{By default, the background map will be drawn in high resolution. Set this parameter to \code{FALSE} if you are studying a big geographical area (>1000 km).}
\item{expl}{Set this to \code{TRUE} to perform an exploratory analysis.}
\item{metad}{Set this to \code{TRUE} if a metadata file is provided (see the details).}
\item{sufbrk}{Suffix to add to \code{varcli} to form the name of the provided
metadata file. This parameter is only relevant when \code{metad=TRUE}. Its
default value \code{'m'} is meant to read the file of break-points detected at
the monthly scale.}
\item{tinc}{Time increment between data. Not set by default, but can be defined
for subdaily data, as in e.g.: \code{tinc='3 hour'}.}
\item{tz}{Time zone. Only relevant for subdaily data. (\code{'UTC'} by default.)}
\item{cex}{Character expansion factor for graphic labels and titles. (Defaults to 1.2. Note that if station names are long, they will not fit in titles when increasing this parameter too much.)}
\item{verb}{Verbosity. Set to \code{FALSE} to avoid messages being output to
the console. (They will be in the output log file anyway).}
  \item{x}{Vector of dates. (To be read from the *.rda file.)}
}
\details{
First of all take into account that this is an experimental function, and will
fail if there are time steps completely void of data in any sub-area.

If you have not decided the splitting meridians and parallels, do not set them, and the function will provide a map to help in selecting the areas.

If you set the \code{xc} and \code{yc} splitting borders, setting \code{maponly=TRUE} will also produce a map with the stations, plus the required overlapping areas, without doing any homogenization. In this way you can review the limits of the areas to choose new ones if you are not happy with the current partitioning.

All parameters except \code{xc}, \code{yc}, \code{xo}, \code{yo} and \code{maponly} are the same as in the \code{homogen} function, and will be passed to it to perform the homogenization. 

If a rectangular area include less than 10 stations, these will be added to the next area. Warning! If it is the last area, they will not be processed, and the homogenization results will have inconsistent number of stations. In this case the user should try a new set of cutting limits.

One graphic output file will be produce for every area containing stations, but the rest of the output will be merged into single files.
}
\value{
  This function does not return any value.
}
\seealso{\code{\link{homogen}}, \code{\link{dahstat}}, \code{\link{dahgrid}}}
\examples{
#Set a temporal working directory and write input files:
wd <- tempdir()
wd0 <- setwd(wd)
data(Ptest)
dim(dat) <- c(720,20)
dat[601:720,5] <- dat[601:720,5]*1.8
write(dat[481:720,1:12],'pcp_1991-2010.dat')
write.table(est.c[1:12,1:5],'pcp_1991-2010.est',row.names=FALSE,col.names=FALSE)
#Now run the example:
homogsplit('pcp',1991,2010,2.9,39.6,0,0,std=2)
#Return to user's working directory:
setwd(wd0)
#Input and output files can be found in directory:
print(wd)
#
# Note that this is just a trivial example; this function is intended to
# homogenize a network of thousands of long series which might overload
# computer resources, and is still experimental.
}
\keyword{datagen}
