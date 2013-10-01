\name{homogen}
\alias{homogen}
\title{Automatic homogenization of climatological series}
\description{
  Automatic homogenization of climatological series, including missing
  data filling and detection and correction of outliers and shifts in the
  mean of the series.
}
\usage{
homogen(varcli, anyi, anyf, nm=12,  nref=10, dz.max=5,
wd=c(0,0,100), swa=60, snht1=25, tVt=snht1, snht2=25, snhtt=snht2,
tol=0.02, tVf=tol, mxdif=.05, force=FALSE, a=0, b=1, wz=.001,
deg=TRUE, rtrans=0, std=3, ndec=1, mndat=0, leer=TRUE, gp=3,
na.strings="NA", nclust=100, maxite=50, ini="", vmin=NA, vmax=NA,
verb=TRUE) }
\arguments{
  \item{varcli}{Acronym of the name of the studied climatic variable,
    as in the data file name.}
  \item{anyi}{Initial year of the data present in the file.}
  \item{anyf}{Final year of the data present in the file.}
  \item{nm}{Number of data per year in each station. (Defaults to 12,
    for monthly data).}
  \item{nref}{Maximum number of references for data estimation. (10 by
    default).}
  \item{dz.max}{Threshold of outlier tolerance, in standard deviations.
    (5 by default).}
  \item{wd}{Distance (in km) at which reference data will weigh half that
    of another located at the same site of the series been estimated.
    (Defaults to \code{c(0,0,100)}, meaning that no weighting will be applied
    in the first two stages, and 100 km in the third).}
  \item{swa}{Size of the step forward to be applied to the staggered window
    application of SNHT.}
  \item{snht1}{Threshold value of the stepped SNHT window test (25 by default).}
  \item{tVt}{(Former name of \code{snht1} maintained for backward
    compatibility).}
  \item{snht2}{Threshold value for the SNHT test when applied to the complete
  series (25 by default).}
  \item{snhtt}{(Former name of \code{snht2} maintained for backward
    compatibility).}
  \item{tol}{Tolerance factor to split several series at a time.}
  \item{tVf}{(Former name of \code{tVf} maintained for backward
    compatibility).}
  \item{mxdif}{Maximum difference of any data item in consecutive iterations.
    (Defaults to 0.05).}
  \item{force}{Force break even when only one reference is available.
    (\code{FALSE} by default).}
  \item{a, b}{Parameters of the optional transformation \code{a+b*dat} to be
    applied to data when read from the files.}
  \item{wz}{Scale parameter of the vertical coordinate \code{Z}.}
  \item{deg}{Set to \code{FALSE} if the input coordinates are in km (the
    distance unit used internally in the package), or left in its default
    \code{TRUE} value if they are in geographical degrees.}
  \item{rtrans}{Root transformation to apply to the data. (2 for square
  root, 3 for cubic root, etc. Fractional numbers are allowed).}
  \item{std}{Type of normalization:
    \describe{
      \item{1:}{deviations from the mean,}
      \item{2:}{proportions of the mean (only for means greater than 1),}
      \item{3:}{standardization (subtract the mean and divide by the
        standard deviation of the sample).}
    }
  }
  \item{ndec}{Number of decimal digits to which the homogenized data must
  be rounded.}
  \item{mndat}{Minimum number of data for a split fragment to become
  a new series.}
  \item{leer}{set to \code{FALSE} if you read your data with your own
  R routines.}
  \item{gp}{Graphic parameter:
    \describe{
      \item{0:}{no graphic output,}
      \item{1:}{only descriptive graphics of the input data,}
      \item{2:}{as with 1, plus diagnostic graphics of anomalies,}
      \item{3:}{as with 2, plus graphics of running annual means and
      applied corrections,}
      \item{4:}{as with 3, but running annual totals (instead of means)
      will be plotted. (Better when working with precipitation data).}
    }
  }
  \item{na.strings}{Character string to be treated as a missing value.
  (It can be a vector of strings, if more than one is needed).}
  \item{nclust}{Maximum number of stations for the cluster analysis. (If much
    greater than 100, the default value, the process may be too long and the
    graphic too dense).}
  \item{maxite}{Maximum number of iterations when computing the means of
  the series. (50 by default).}
  \item{ini}{Initial date. Void by default, if set (with format
  \code{'YYYY-MM-DD'}), it will be assumed that the series contain daily data.}
  \item{vmin}{Minimum possible value (lower limit) of the studied variable.
  Unset by default, but note that \code{vmin=0} will be applied if \code{std}
  is set to 2.}
  \item{vmax}{Maximum possible value (upper limit) of the studied variable.
  (E.g., for relative humidity or relative sunshine hours it is advisable to
  set \code{vmax=100}).}
  \item{verb}{Verbosity. Set to \code{FALSE} to avoid messages being output to
  the console. (They will be in the output log file anyway).}
}
\details{
  This function reads the data from two files, and writes its results in
  other four. They are distinguished by their extension, since all have the
  same base name, \file{VAR_FIRSTY-LASTY}, where \code{VAR}
  stands for any acronym of the involved climatological variable, and
  \code{FIRSTY} and \code{LASTY} are the initial and final years of the
  data. Data are stored in a text file with extension \code{dat}, station
  by station, and chronologically within each station block, in a free
  blank separated format.

  The root transformation of the input data may be very useful in
  normalizing highly biased L-shape distributed variables (e.g.
  precipitation, wind speed, ...)

  The other input file must have \code{est} as its extension, and is also
  a blank separated text file with coordinates, codes and names of the
  observing stations. Each line identifies a single station, with structure
  \code{X Y Z CODE NAME}, where the coordinates \code{X} and \code{Y} may
  be in km or in geographical degrees (longitude and latitude) in decimal
  form (not in degrees, minutes and seconds). The other parameters are the
  altitude \code{Z} in m, an identification \code{CODE}, and the full
  \code{NAME} of the station (quoted if it contains blanks).

  This function may be applied to daily data, by setting \code{nm=0} and
  using \code{ini} to indicate the initial date (that can be different than
  January the first). The names of the input files are almost the same:
  \file{VAR-d_FIRSTY-LASTY}. As you can see, the only difference is the
  \code{-d} attached to the variable acronym, that helps to distinguish
  the daily file from a monthly file possibly coexisting in the same
  directory. Anyway, the homogenization of daily data is very challenging,
  and only should be addressed if the shifts in the mean are quite big.
  But it may be very useful to fill missing daily data (set \code{tVt=0}
  to skip the detection/correction of breaks).

  The default values of \code{dz.max}, \code{tVt} and \code{snhtt} try to
  be fairly conservative with the original series, but appropriate values
  may vary with the type of climatic variable, its spatial variability, the
  density of the observing network, and the kind of data (annual, seasonal,
  monthly, daily, ...), and therefore it is advisable to adjust them
  empirically with the help of the graphic diagnostic output of a first
  exploratory application. For a more extensive explanation and discussion,
  see the accompanying user's guide, also available at
  \file{http://www.climatol.eu/climatol-guide.pdf}.
}
\value{
  This function does not return any value, but the objects it creates
  (either directly or by calling internal functions) remain in the memory
  space during the rest of the R session if not explicitly removed,
  allowing the user to analyze them by their own means. This objects are:
  \describe{
    \item{dat:}{observed input data,} 
    \item{dat.d:}{working copy of observed data,} 
    \item{dat.z:}{normalized observed data,} 
    \item{dat.e:}{normalized estimated data,}
    \item{dat.c:}{computed (estimated) data,}
    \item{dat.m and dat.m0:}{data means,}
    \item{dat.s and dat.s0:}{data standard deviations,}
    \item{dat.na:}{missing data,}
    \item{datmed:}{overall data means,}
    \item{refmed:}{overall reference mean,}
    \item{refstd:}{overall reference standard deviation,}
    \item{refhom:}{index of preset homogeneous reference series,}
    \item{dah:}{homogenized output data,}
    \item{anom:}{anomalies,}
    \item{sanom:}{standardized anomalies,}
    \item{outan:}{standardized anomalies of the outliers,}
    \item{est.c:}{station data (coordinates, codes and names),}
    \item{est.d:}{inter-station distances,}
    \item{est.i:}{initial station data,}
    \item{est.p:}{inter-station proximity ranks,}
    \item{est.w:}{inter-station weights,}
    \item{iest:}{index pointing at original series,}
    \item{oneref:}{flag of data with only one reference,}
    \item{mindist:}{minimum distance to the reference data,}
    \item{used:}{flag of used reference series,}
    \item{anyi:}{first year of data,}
    \item{anyf:}{last year of data,}
    \item{na:}{number of years of the study period,}
    \item{nd:}{number of data per station,}
    \item{ne:}{number of stations,}
    \item{nei:}{initial number of stations,}
    \item{tVx:}{maximum value of the SNHT applications,}
    \item{kpx:}{location of each \code{tVx} in each series,}
    \item{verde:}{Green color used in several places.}
  }

  After successfully running this function, four output files will be saved
  in your working directory, with the same base name as the input files,
  and extensions:
  \describe{
    \item{.txt}{A text file that logs all the processing output to the
    console,}
    \item{.pdf}{A PDF file with a collection of diagnostic graphics,}
    \item{.dah}{A text file containing the homogenized data, with missing
    data filled and the same structure as the \code{.dat} input data file,}
    \item{.esh}{A text file with the coordinates and names of the stations
    of the homogenized data file.}
  }
}
\references{
  Guijarro, J.A. (2008): Homogenization of a dense thermo-pluviometric
  monthly database in the Balearic Islands using the free contributed R
  package CLIMATOL.  WMO Fifth Seminar for Homogenization and Quality
  Control in Climatological Databases, Budapest, Hungary, May 2006.
  WCDMP-No. 68, WMO-TD No. 1434.

  Guijarro, J.A. (2008): Homogeneización de las series termopluviométricas
  de Baleares y su influencia en la evaluación de las tendencias.
  En SIGRÓ J. et al. (Eds.), Cambio climático regional y sus impactos,
  Asociación Española de Climatología, A-6:151-160.
}
\author{Jose A. Guijarro}
\seealso{\code{\link{dahstat}}, \code{\link{outrename}}, \code{\link{dd2m}}}
\examples{
#As this function only works with files, you must uncompress the
#example data in your working directory first. (The example files
#can be obtained from http://www.climatol.eu/climatol-dat.zip).
#Afterwards, you can run the example:
#
#   homogen("Tmin", 1956, 2005, deg=FALSE)
#
#See the resulting four output files "Tmin_1956-2005.dah",
#"Tmin_1956-2005.esh", "Tmin_1956-2005.txt" and "Tmin_1956-2005.pdf".
#
#Another example with daily data, but only for filling missing data:
#
#   homogen('Tmax', 2007, 2010, nm=0, tVt=0, snhtt=0, ini='2007-01-01')
#
#The results will have been saved to four Tmax-d_2007-2010.* files.
}
\keyword{datagen}
\keyword{ts}
\keyword{manip}
