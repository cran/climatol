\name{homogen}
\alias{homogen}
\title{Automatic homogenization of climatological series}
\description{
  Automatic homogenization of climatological series, including missing
  data filling and detection and correction of outliers and shifts in the
  mean of the series.
}
\usage{
homogen(varcli, anyi, anyf, test='snht', nref=NULL, std=NA, swa=NA,
ndec=1, niqd=3, dz.max=.01, dz.min=-dz.max, cumc=NA, wd=NULL, inht=25,
sts=5, maxdif=NA, maxite=999, force=FALSE, wz=.001, mindat=NA, onlyQC=FALSE,
annual=c('mean','sum','total'), x=NULL, ini=NA, na.strings="NA", vmin=NA,
vmax=NA, hc.method='ward.D2', nclust=300, cutlev=NA, grdcol=grey(.4),
mapcol=grey(.4), expl=FALSE, metad=FALSE, sufbrk='m', tinc=NA, tz='utc',
rlemin=NA, rlemax=NA, cex=1.1, uni=NA, raway=TRUE, graphics=TRUE, verb=TRUE,
logf=TRUE, snht1=NA, snht2=NA, gp=NA)
}
\arguments{
\item{varcli}{Short name of the studied climatic variable.}
\item{anyi}{Initial year of the data.}
\item{anyf}{Final year of the data.}
\item{test}{Inhomogeneity test to apply: 'snht' (the default) or 'cuct'
(Cucconi test, experimental).}
\item{nref}{Maximum number of references for data estimation [defaults to 10
in the detection stages, and to 4 in the final series adjustments].}
\item{std}{Type of normalization:
  \describe{
    \item{1:}{deviations from the mean,}
    \item{2:}{rates to the mean (only for means greater than 1),}
    \item{3:}{standardization (subtract the mean and divide by the sample
      standard deviation).}
  }
}
\item{swa}{Size of the step forward to be applied to the overlapping window
application of the detection test [365 terms (one year) for daily data, and 60
otherwise].}
\item{ndec}{Number of decimal digits to round the homogenized data [1].}
\item{niqd}{Number of interquartilic distances to delete big outliers [3].}
\item{dz.max}{Threshold of outlier tolerance, in standard deviations if greater than one, or as a percentage of data to reject otherwise [0.01].}
\item{dz.min}{Lower threshold of outlier tolerance if different from \code{dz.max}.}
\item{cumc}{Code of accumulated missing data.}
\item{wd}{Distance (in km) at which reference data will weight half of those
located at zero distance [\code{c(0,0,100)}].}
\item{inht}{Thresholds for the change in the mean detection tests [25].}
\item{sts}{Series tail size (no. of terms) not tested for inhomogeneities [5].}
\item{maxdif}{Maximum data difference from previous iteration [\code{ndec/2}].}
\item{maxite}{Maximum number of iterations to compute means of the series
[999].}
\item{force}{Force direct homogenization of daily or sub-daily series [\code{FALSE}].}
\item{wz}{Scale parameter of the vertical coordinate \code{Z} [0.001].}
\item{mindat}{Minimum number of data for a split fragment to become a new
series [\code{swa/2 for daily series or 12 terms otherwise}].}
\item{onlyQC}{Set to \code{TRUE} if only initial Quality Controls are requested [\code{FALSE}]}
\item{annual}{Running annual value to graph in the PDF output. One of 'mean' (the default), 'sum' or 'total' (equivalent to 'sum').}
\item{x}{Time vector. Only needed if data are taken at irregular intervals.}
\item{ini}{Initial date, with format \code{'YYYY-MM-DD'}, if series does not begin on January first (as recommended).}
\item{na.strings}{Character strings to be treated as missing data [\code{'NA'}].}
\item{vmin}{Minimum possible value (lower limit) of the studied variable.}
\item{vmax}{Maximum possible value (upper limit) of the studied variable.}
\item{hc.method}{Hierarchical clustering method ['ward.D2'].}
\item{nclust}{Maximum number of series for the cluster analysis [300].}
\item{cutlev}{Level to cut the dendrogram to define clusters [\code{NA}].}
\item{grdcol}{Color of the graphic background grids [\code{grey(0.04)}].}
\item{mapcol}{Color of coastlines and borders in the stations map [\code{grey(0.04)}].}
\item{expl}{Perform an exploratory analysis? [\code{FALSE}].}
\item{metad}{Use the breakpoints file as metadata? [\code{FALSE}].}
\item{sufbrk}{Suffix to add to \code{varcli} to form the name of the provided
metadata file ['m'].}
\item{tinc}{Time increment between data [\code{NA}].}
\item{tz}{Time zone [\code{'utc'}]. Only relevant for subdaily data.}
\item{rlemin}{Data run lengths will exclude values \code{<= rlemin} in quality control [\code{NA}].}
\item{rlemax}{Data run lengths will exclude values \code{>= rlemax} in quality control [\code{NA}].}
\item{cex}{Character expansion factor for graphic labels and titles [1.1].}
\item{uni}{Units to use in some axis labels [''].}
\item{raway}{Increase internal distances to reanalysis series to give more
weight to observed series [\code{TRUE}].}
\item{graphics}{Output graphics in a PDF file [\code{TRUE}].}
\item{verb}{Verbosity [\code{TRUE}].}
\item{logf}{Save console messages to a log file? [\code{TRUE}].}
\item{snht1, snht2}{Obsolete (use \code{inht} instead), but kept for backwards compatibility.}
\item{gp}{Obsolete (use \code{graphics=FALSE} for \code{gp=0}, \code{onlyQC=TRUE} for \code{gp=1} or \code{annual="total"} for \code{gp=4}), but kept for backwards compatibility.}
}
\details{
Input data must be provided in two text files, one with the data (with
extension \code{dat}) and another with the station coordinates (with extension
\code{est}). Both have as base name, \file{VAR_YEAR-YEAR}, composed by
the short name of the climatological variable, and the initial and final years
of the data, as set in the first three parameters of the call, \code{varcli},
\code{anyi} and \code{anyf}.

Data are stored in a free blank separated format (any number of data items per
line is allowed), in chronological order, station by station (all data from
station 1 go first, then all data from station 2, and so on). As dates are not
stored in this file, all data must be present in the file, using a code for any
missing data in the records (\code{NA} by default, but any other code can be
used, provided that they are specified in the parameter \code{na.strings}).

The stations file, with extension \code{est}, is also a blank separated text
file where each line identifies a single station, with structure \code{'X Y Z
CODE NAME'}. Coordinates \code{X} and \code{Y} are expected in geographical
degrees (longitude and latitude, in this order and in decimal form). Otherwise
they will be assumed to be in km, or in m if the mean of either \code{X} and
\code{Y} is greater than 10000; elevation \code{Z} must be supplied in m; and
the identification \code{CODE} and the full \code{NAME} of the station must be
quoted if they contains blanks). Fully reliable series may be marked by putting
an asterisk (*) at the beginning of their \code{CODE} to skip their outlier and
break-point analysis. This is not recommended with observed series, but can be
useful when using reanalysis series as references in data sparse regions.

This function will stop with an error condition if any time step becomes void
of data in all stations at the same time. One or more series with data in the
void time steps must be added to successfully run \code{homogen} again. If no
other series are available in the area, reanalysis series of the closer
grid-points can be used, adding their coordinates to the *.est file and
prepending an asterisk (*) to the codes assigned to the series as mentioned
above.

\code{dz.max} (and \code{dz.min} if different from \code{dz.max}) can be a
vector of two values, one for suspect data and the other for probable errors.
Only the latter will be deleted, but all will be listed in the \file{*_out.csv}
output file. By default, the more extreme 0.01\% in each tail of the
distribution will be considered errors, and values exceeding 0.1\% will be
suspect data. Inspection of the anomalies histogram near the end of the PDF
output file will help in tuning these parameters by setting number of standard
deviations to be used as rejection thresholds.

\code{inht} has a default value of 25, which is a suitable conservative value
for monthly values of temperature, but not so much for precipitation or for
daily series. Therefore it is advisable to adjust it empirically with the help
of the histograms available by the end of the graphic output. Anyway,
inhomogeneities in daily or subdaily series should be detected on their monthly
aggregates, which can be easily obtained by means of the function \code{dd2m}.
Two values can be given to this parameter (one for each of the two detection
stages), as in e.g.  \code{inht=c(30,25)}. When only one value is provided, it
will be used for both stages. If any or both values are zeros, the
corresponding homogenization stage will be skipped.

The default value \code{wz=0.001} gives to the vertical coordinate (in m) the
same weight as the horizontal coordinates (internally managed in km). Other
values can be set to overweight elevation differences (wz>0.001) or to
calculate only horizontal distances (wz=0).

\code{vmin} and \code{vmax} are unset by default, but if the variable is found
to have a skewed probability distribution with a minimum value of zero,
\code{vmin} will be set to zero. The same will happen if the user sets
\code{std=2}.

\code{sufbrk} is only relevant when \code{metad=TRUE}. Its default value
\code{'m'} is meant to read the file of break-points detected at the monthly
scale, but if the data were originally monthly, \code{sufbrk=''} should be set.

\code{tinc}, unset by default, can be defined for subdaily data, as in e.g.:
\code{tinc='3 hours'}, especially if first and/or last years are incomplete.
Units can be 'hours', 'mins' or 'secs'.

The default \code{cex=1.1} increase by a 10\% the size of labels in the graphic
output. Note that if station names are long, they will not fit in titles when
increasing this parameter too much.

The graphic output file (in PDF format) begins with a first quality control of
the series, providing box-plots for every series showing (1) the range of their
values, (2) their increments between successive terms and (3) the length of
segments with constant data. Too big outliers are deleted at this stage because
they would compromise the quality of the homogenization and missing data
filling. During the rest of the process outlier detection is based on spatial
differences between neighboring normalized data.

The following pages offer: (a) a summary of the data availability and frequency
distribution; (b) a correlogram of the first differences of the series, (c) a
dendrogram based on these correlations and a map with the station locations
(marked with numbers if less than 100, and with symbols otherwise; (d) graphics
of normalized spatial anomalies showing the detected breaks, the minimum
distance to a reference data and the number of references used; (e) a histogram
of maximum inht values found in overlapping window analysis; (d) and (e) are
repeated for the analysis on the whole series; (f) histograms of number of
splits per station and per year; (g) graphics of final anomalies of the series;
(h) graphics of the reconstructed series and applied corrections; (i) a
histogram of the normalized anomalies of all data (useful to set rejection
thresholds for the outliers); (j) final histograms of inht values; and (k) a
plot of quality/singularity of the stations (a bad score may be due to a bad
quality of the series, but also to a singular siting with a peculiar
micro-climate).

Note that every time that a significant shift in the mean of the series is
detected, it will be split into two (potentially) homogeneous sub-periods, and
hence the final number of homogenized series will be increased, as complete
homogeneous series will be reconstructed from all of them. When several
homogeneous series have been yielded for the same location, the user can choose
to use that reconstructed from the last sub-period (the usual behavior of
other homogenization packages), which is perfect for climate monitoring of
newly incoming data. However, statistics derived from all of them can be useful
for climatic mapping, when no a priori knowledge can indicate which of the
sub-periods will be more representative at the studied spatial scale).

The processing time can be from seconds (a few monthly series) to many hours
(hundreds of daily series). If you must process a huge amount of data, you
should consider splitting your study region into smaller areas to be
homogenized independently.
}
\value{
This function does not return any value, its results being saved to files
with the same base name as the input files, and extensions:  
\describe{
\item{*.txt:}{A text file that logs all the processing output,} 
\item{*_out.csv:}{List of corrected outliers,} 
\item{*_brk.csv:}{List of corrected breaks,} 
\item{*.pdf:}{PDF file with a collection of diagnostic graphics,}
\item{*.rda:}{Homogenization results in R binary format, used by the
\code{dahstat} and other post-processing functions, but can be loaded
by the user for further data manipulation with the function \code{load}. This
file contains the following objects:
  \describe{
  \item{dat}{matrix of the original series,}
  \item{dah}{matrix of the homogenized series,}
  \item{nd}{number of data (time steps) in every series,}
  \item{ndec}{number of decimals in the data,}
  \item{uni}{data units,}
  \item{est.c}{data frame with columns:
    \describe{
    \item{X}{longitude,}
    \item{Y}{latitude,}
    \item{Z}{elevation,}
    \item{Code}{station code,}
    \item{Name}{station name,}
    \item{pod}{percentage of original data,}
    \item{snht}{(or \code{cuct} when \code{test='cuct'}): Remaining
      inhomogeneity test values in the homogenized series. Can be greater
      than the set \code{inht} threshold because of a lower number of reference
      stations,}
    \item{rmse}{estimated root mean squared errors of the homogenized series}
    }}
  \item{ct}{Cluster Analysis series groups,}
  \item{nei}{number of input series,}
  \item{ne}{number of series after the homogenization,}
  \item{nm}{number of "months" (data items) in a year (0=daily data),}
  \item{std}{type of normalization applied to the data,}
  \item{x}{vector of the time dimension,}
  \item{ini}{initial date of the period under study.}
  }
} }
}
\references{
Guijarro JA (2014): Quality Control and Homogenization of Climatological
Series. In Eslamian S (Ed.), Handbook of Engineering Hydrology, Vol. 1:
Fundamentals and Applications. Francis and Taylor, CRC Group, USA, ISBN
9781466552357, 636 pp.

Azorin-Molina C, Guijarro JA, McVicar TR, Trewin BC, Frost AJ, Chen D (2019):
An approach to homogenize daily peak wind gusts: An application to the
Australian series. Int. J. Climatol., 39:2260-2277. doi: 10.1002/joc.5949

Dumitrescu A, Cheval S, Guijarro JA (2019): Homogenization of a combined hourly
air temperature dataset over Romania. Int. J. Climatol., 40:2599-2608, DOI:
10.1002/joc.6353

Visit <https://climatol.eu/> for updates of code and documentation (user's
guide, links to videos, etc).
}
\seealso{\code{\link{dahstat}}, \code{\link{dahgrid}}, \code{\link{outrename}},
\code{\link{dd2m}}}
\examples{
## Set a temporal working directory and write input files:
wd <- tempdir()
wd0 <- setwd(wd)
data(climatol_data)
write.table(Temp.est,'Temp_1961-2005.est',row.names=FALSE,col.names=FALSE)
write(Temp.dat,'Temp_1961-2005.dat',ncolumns=12)

## Now run the example:
homogen('Temp',1961,2005)

## Return to user's working directory:
setwd(wd0)

## Input and output files can be found in directory:
print(wd)
}
\keyword{cluster}
\keyword{datagen}
\keyword{graphs}
\keyword{hplot}
\keyword{ts}
