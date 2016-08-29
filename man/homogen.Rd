\name{homogen}
\alias{homogen}
\title{Automatic homogenization of climatological series}
\description{
  Automatic homogenization of climatological series, including missing
  data filling and detection and correction of outliers and shifts in the
  mean of the series.
}
\usage{
homogen(varcli, anyi, anyf, suf=NA, nm=NA,  nref=c(10,10,4), std=3, swa=NA,
ndec=1, dz.max=5, dz.min=-dz.max, wd=c(0,0,100), snht1=25, snht2=snht1,
tol=.02, mxdif=NA, force=FALSE, wz=.001, trf=0, mndat=NA, gp=3, ini=NA,
na.strings="NA", maxite=50, vmin=NA, vmax=NA, nclust=100,
clustmethod='ward.D2', grdcol=grey(.5), mapcol=grey(.65), hires=TRUE,
expl=FALSE, metad=FALSE, sufbrk='m', verb=TRUE)
}
\arguments{
\item{varcli}{Acronym of the name of the studied climatic variable, as in the
input data files.}
\item{anyi}{Initial year of the data present in the file.}
\item{anyf}{Final year of the data present in the file.}
\item{suf}{Optional suffix appended with a '-' to the name of the variable
in the input files.}
\item{nm}{Number of data per year in each station. (Defaults to NA, and then it will be computed from the total number of data).}
\item{nref}{Maximum number of references for data estimation. (Defaults to 10 in the detection stages, and to 4 in the final series adjustments).}
\item{std}{Type of normalization:
  \describe{
    \item{1:}{deviations from the mean,}
    \item{2:}{rates to the mean (only for means greater than 1),}
    \item{3:}{standardization (subtract the mean and divide by the sample
      standard deviation).}
  }
}
\item{swa}{Size of the step forward to be applied to the staggered window
application of SNHT. If not set (the default), 365 terms (one year) will be
used for daily data, and 60 otherwise.}
\item{ndec}{Number of decimal digits to which the homogenized data must be
rounded.}
\item{dz.max}{Threshold of outlier tolerance, in standard deviations.  (5 by
default).}
\item{dz.min}{Lower threshold of outlier tolerance may be different from the higher one. (By default, they are the same, with opposite signs).}
\item{wd}{Distance (in km) at which reference data will weigh half of that
of another located at zero distance. (Defaults to \code{c(0,0,100)}, meaning
that no weighting will be applied in the first two stages, and 100 km in the
third).}
\item{snht1}{Threshold value for the stepped SNHT window test applied in stage
1. (25 by default).}
\item{snht2}{Threshold value for the SNHT test when applied to the complete
series in stage 2 (same value as snht1 by default).}
\item{tol}{Tolerance factor to split several series at a time. The default is 0.02, meaning that a 2\% will be allowed for every reference data. (E.g.: if the maximum SNHT test value in a series is 30 and 10 references were used to compute the anomalies, the series will be split if the maximum test of the reference series is lower than 30*(1+0.02*10)=36. Set \code{tol=0} to disable further splits when any reference series has already been split at the same iteration).}
\item{mxdif}{Maximum difference of any data item in consecutive iterations. If not set, defaults to half of the data precision (defined by the number of decimals).}
\item{force}{Force break even when only one reference is available.
(\code{FALSE} by default).}
\item{wz}{Scale parameter of the vertical coordinate \code{Z}. 0.001 by default, which gives the vertical coordinate (in m) the same weight as the horizontal coordinates (internally managed in km).}
\item{trf}{By default, data are not transformed (\code{trf=0}), but if the data
frequency distribution is very biased, the user can choose to apply a log(x+1)
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
standard R missing data code.}
\item{maxite}{Maximum number of iterations when computing the means of the
series. (50 by default).}
\item{vmin}{Minimum possible value (lower limit) of the studied variable.
Unset by default, but note that \code{vmin=0} will be applied if \code{std} is
set to 2.}
\item{vmax}{Maximum possible value (upper limit) of the studied variable.
(E.g., for relative humidity or relative sunshine hours it is advisable to set
\code{vmax=100}).}
\item{nclust}{Maximum number of stations for the cluster analysis. (If much
greater than 100, the default value, the process may be too long and the
graphic too dense).}
\item{clustmethod}{Clustering method. Defaults to \code{'ward.D2'}, but the user can choose any hierarchical method available in the \code{hclust} function.}
\item{grdcol}{Color of the graphic background grids. (Middle gray by default).}
\item{mapcol}{Color of the background map (gray by default).}
\item{hires}{By default, the background map will be drawn in high resolution. Set this parameter to \code{FALSE} if you are studying a big geographical area (>1000 km).}
\item{expl}{Set this to \code{TRUE} to perform an exploratory analysis.}
\item{metad}{Set this to \code{TRUE} if a metadata file is provided (see the details).}
\item{sufbrk}{Suffix to add to \code{varcli} to form the name of
the provided metadata file. \code{'m'} by default, to read the
breaks detected at the monthly scale. (Only relevant when \code{metad=TRUE}).}
\item{verb}{Verbosity. Set to \code{FALSE} to avoid messages being output to
the console. (They will be in the output log file anyway).}
}

\details{
Input data must be provided in two text files, one with the data (with
extension \code{dat}) and another with the station coordinates (with extension
\code{est}). Both have as base name, \file{VAR_YEAR-YEAR}, composed by
the acronym of the climatological variable, and the initial and final years of
the data, as set in the first three parameters of the call, \code{varcli},
\code{anyi} and \code{anyf}.

Data are stored in a free blank separated format (any number of data items per
line is allowed), in chronological order, station by station (all data from
station 1 go first, then all data from station 2, and so on). As dates are not
stored in this file, all data must be present in the file, using a code for any
missing data in the records (\code{NA} by default, but any other code can be
used, provided that they are specified in the parameter \code{na.strings}).

The stations file, with extension \code{est}, is also a blank separated text
file where each line identifies a single station, with structure \code{X Y Z
CODE NAME}: Coordinates \code{X} and \code{Y} may be in geographical degrees
(longitude and latitude, in decimal form), or in m if you are using a UTM
projection; elevation \code{Z} must be supplied in m; and the identification
\code{CODE} and the full \code{NAME} of the station must be quoted if they
contains blanks).

The transformation of the input data may be very useful to normalize highly
biased L-shape distributed variables (as is often the case with precipitation,
wind speed, ...), but use it with caution. Alternatively, you can use the rate
normalization (\code{std=2}) on the raw data if the variable has a natural zero
lower limit. (This alternative has yielded better results than transformation
in some applications, provided that no homogeneous sub-period has
means lower than 1. If this is the case, a workaround may be to multiply all
data values by a constant prior to their homogenization).

The default values of \code{dz.max}, \code{snht1} and \code{snht2} can be
appropriate for monthly values of temperature, but not so much for
precipitation or for daily series. Therefore it is advisable to adjust them
empirically with the help of the histograms in the graphic output of a first
exploratory application, using \code{expl=TRUE} in the call to the function.

This graphic output includes: a) a summary of the data availability and frequency distribution; b) a correlogram of the first differences of the series;
c) a dendrogram based on these correlations and a map with the station locations; d) graphics of anomalies showing the detected breaks, the minimum distance to a reference data and the number of references used; e) a histogram of maximum SNHT values found in overlapping window analysis; d) and e) are repeated for the analysis on the whole series; f) histograms of number of splits per station and per year; g) graphics of final anomalies of the series h) graphics of the reconstructed series and applied corrections; i) a histogram of the normalized anomalies of all data (may help to set rejection thresholds for the outliers); final histograms of SNHT values; and j) a plot of quality/singularity of the stations (a bad score may be due to a bad quality of the series, but also to a singular siting with a peculiar micro-climate).

Note that every time that a significant shift in the mean of the series is
detected, it will be split into two (potentially) homogeneous sub-periods, and
hence the final number of homogenized series will be increased, as complete
homogeneous series will be reconstructed from them. When several homogeneous
series have been yielded for the same location, the user can choose to use that
reconstructed from the last sub-period (for climate monitoring), from the
period with the highest percentage of original data or the highest final SNTH
(the more robust reconstruction) or all of them (e.g., for climatic mapping,
when no a priori knowledge can indicate which of the sub-periods will be more
representative at the studied spatial scale). Aditional columns \code{pod},
\code{ios}, \code{ope} and \code{snht} in the stations data frame can help the
user to choose the most appropriate series for each location when using the
post-processing functions \code{dahstat} and \code{dahgrid} (see below).
}
\value{
This function does not return any value, its results being saved to files
with the same base name as the input files, and extensions:  
\describe{
\item{*.txt:}{A text file that logs all the processing output,} 
\item{*_out.csv:}{List of corrected outliers,} 
\item{*_brk.csv:}{List of corrected breaks,} 
\item{*.pdf:}{PDF file with a collection of diagnostic graphics,}
\item{*.rda:}{Homogenization results in R binary format, used by \code{dahstat}
and \code{dahgrid} post-processing functions, but can be loaded by the user for
further with the function \code{load}). This file contains the following
objects:
  \describe{
  \item{dat}{matrix or array of original data,}
  \item{dah}{matrix or array of homogenized data,}
  \item{est.c}{data frame with columns:}
    \describe{
    \item{X}{X coordinate,}
    \item{Y}{Y coordinate,}
    \item{Z}{elevation,}
    \item{Code}{code of the station,}
    \item{Name}{name of the station,}
    \item{pod}{percentage of original data,}
    \item{ios}{index of original station,}
    \item{ope}{operating at the end of the period? (0=no, 1=yes),}
    \item{snht}{relative SNTH of the homogenized series}
    }
  \item{nd}{number of time steps in every series,}
  \item{ne}{number of series after the homogenization,}
  \item{nei}{initial number of series,}
  \item{nm}{number of "months" in a year (0=daily data),}
  \item{ndec}{number of decimals in the data,}
  \item{std}{type of standardization used (as explained in the details),}
  \item{ini}{initial date of the period under study}
  }
}
}
}
\references{
Guijarro JA (2014): Quality Control and Homogenization of Climatological
Series. In Eslamian S (Ed.), Handbook of Engineering Hydrology, Vol. 1:
Fundamentals and Applications. Francis and Taylor, CRC Group, USA, ISBN
9781466552357, 636 pp.

Guijarro JA (2015): Homogenization of Spanish mean wind speed monthly series.
In 8th Seminar for Homogenization and Quality Control in Climatological
Databases and Third Conference on Spatial Interpolation Techniques in
Climatology and Meteorology, Budapest, 12-16 May 2014 (Lakatos et al.,
Eds.), WMO Climate Data and Monitoring WCDMP-No. 84, pp. 98-106.

Azorin-Molina C, Guijarro JA, McVicar TR, Vicente-Serrano SM, Chen D, Jerez S,
Espírito-Santo F (2016): Trends of daily peak wind gusts in Spain and
Portugal, 1961-2014. J. Geophys. Res. Atmos., 121, doi:10.1002/2015JD024485,
20 pp.
}
\author{Jose A. Guijarro}
\seealso{\code{\link{dahstat}}, \code{\link{dahgrid}}, \code{\link{outrename}},
\code{\link{dd2m}}}
\examples{
\dontrun{
#First we must generate the input files from example data:
data(Ptest)
write(dat,'Ptest_1951-2010.dat',ncolumns=12)
write.table(est.c,'Ptest_1951-2010.est',row.names=FALSE,col.names=FALSE)
rm(dat,est.c) #remove loaded data from memory space
#Now we can apply this function to these monthly precipitation data:
homogen('Ptest',1951,2010,std=2,dz.max=7,snht1=15,gp=4)
#
#Now with daily mean temperatures. Preparation of the input files:
data(Ttest)
write(dat,'Ttest_1981-2000.dat')
write.table(est.c,'Ttest_1981-2000.est',row.names=FALSE,col.names=FALSE)
rm(dat,est.c)
#Daily series are quite noisy. Let's aggregate them into monthly values:
dd2m('Ttest',1981,2000)
#Homogenization of the monthly averages:
homogen('Ttest-m',1981,2000)
# If you have a history of the stations, you can edit the file
# Ttest_1981-2000_brk.csv to adjust the dates of the breaks to known relevant
# changes.
# Now we can obtain the homogenized daily series adjusting them to the breaks
# found at the monthly scale:
homogen('Ttest',1981,2000,dz.max=7,metad=TRUE)
}
}
\keyword{datagen}
\keyword{ts}
\keyword{manip}

