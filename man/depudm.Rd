\name{depudm}
\alias{depudm}
\title{Monthly data homogenization}
\description{
  Climatological series homogeneity studies, with missing data estimation,
  inhomogeneities (point errors, mean shifts and trends) detection, and
  graphical displays.
}
\usage{
  depudm(varcli, anyi, anyf, nm = 12, wa = 100, dz.max = 2, difumb = 0.05,
    leer = TRUE, a = 0, b = 1, wz=0.001, sqrtrans = FALSE, ttip = 3,
    refglob = FALSE, ndec = 1, pval = 0.05, graf = FALSE, auto = FALSE,
    verb=TRUE)
}
\arguments{
  \item{varcli}{Acronym of the name of the studied climatic variable,
    as in the data file name.}
  \item{anyi}{Initial year of the data present in the file (four digits)}
  \item{anyf}{Final year of the data present in the file (four digits)}
  \item{nm}{Number of series in each station. (Default=12, for monthly data).}
  \item{wa}{Shape parameter of the weighting function \code{1/(1+d^2/wa)},
    where d stands for distance. Low values (1-10) increase the weight of
    nearby stations, while high values (>1000) take in account also the
    far stations influence. Do \code{wa=0} if you want an unweighted average of
    all the stations.}
  \item{dz.max}{Threshold to accept differences between observed and estimated
    data, in standard deviation units. (Only used when \code{auto=TRUE}).}
  \item{difumb}{Maximum acceptable difference in the series average values
    between missing data filling iterations. This process will stop when
    the maximum difference is lower or equal to \code{difumb}.}
  \item{leer}{if \code{TRUE}, climatological data will be read from files.
    (Can be set to \code{FALSE} if data have already been read in a previous
    call to \code{depudm}).}
  \item{a, b}{Parameters of the optional transformation \code{a+b*dat} to be
    applied to data when read from the files.}
  \item{wz}{Scale parameter of the vertical coordinate \code{Z}. The default
    value assumes that \code{X} and \code{Y} are expressed in km, while
    \code{Z} is expressed in m.  Can be used to change the \code{Z} weight in
    inter-station distance computations.}
  \item{sqrtrans}{if \code{TRUE}, a square root transformation will  be applied
    to all data greater than 1. (Useful with e.g. precipitation data, to
    approximate their distribution to a Gauss one).}
  \item{ttip}{Type of standardization:
    \describe{
      \item{0:}{none,}
      \item{1:}{deviations from the mean,}
      \item{2:}{proportions of the mean (only for means greater than 1),}
      \item{3:}{full standardization (subtract the mean and divide by the
        standard deviation).}
    }
  }
  \item{refglob}{if \code{TRUE}, use annual averages for normalization of all
    the series of each station. Incompatible with \code{graf=TRUE}. Can be
    useful in arid places with frequent null precipitation mixed with high
    precipitation values in other years, making quite unstable the computation
    of their averages.}
  \item{ndec}{Number of decimal places of the purged data, to be saved in
    the file \file{VAR\_AI-AF.dep}.}
  \item{pval}{If greater than 0, t-test of mean difference will be applied
    to running windows of 10 and 20 terms (between samples of 5 and 10 terms),
    and an overall trend test. If \code{graf=TRUE}, p-values of the t-test
    will be plotted and, provided that the trend is significant (lower than
    \code{pval}), the regression line with time will be plotted as well.}
  \item{graf}{If \code{TRUE}, interactive graphs will be displayed for each
    series.}
  \item{auto}{If \code{TRUE}, data whose differences to their estimated
    (normalized) values are greater than \code{dz.max} will be substituted by
    their estimates.} \item{verb}{If \code{TRUE}, progress indications will be
    shown in the terminal.}
}
\details{
  This is the main function for error correction and homogeneity testing of the
  climatological series, and makes frequent calls to other subordinated
  functions.  The climatological series are read from files named
  \file{VAR\_AI-AF.dat}, and the coordinates and names of the stations from
  \file{VAR\_AI-AF.est}, where \code{VAR} stands for any acronym of the
  involved climatological variable, and \code{AI} and \code{AF} are the two
  last digits of the initial and final year of the data. Data are stored
  station by station, and chronologically within each station block.  Missing
  values are specified as \code{NA} (the usual way in R). In the stations file
  there will be a record (line) for each station, with structure \code{X Y Z ID
  NAME}, where \code{X} and \code{Y} are the UTM coordinates in km, \code{Z}
  the altitude in m, \code{ID} a station identifier, followed by the full
  \code{NAME} of the station, that must be put between double quotes if it
  contains any space character (see the example files in the package's data
  directory).  \code{depudm} will save the purged, filled series in files
  \file{VAR\_AI-AF.dep}, eventually overwriting the pre-existing ones. (Rename
  them if results of different calls to \code{depudm} are to be saved). A log
  of each run will also be appended in a file named \file{climatol.log}.
}
\value{
  Objects created by this function (original data, \code{dat.d}; normalized
  data, \code{dat.z}; estimated data, \code{dat.e}; ...) will remain resident
  in the memory space during the rest of the R session while not explicitly
  removed, therefore been susceptible to apply on them all the extended
  variety of statistical and graphic functions available in R.
}
\references{
  Climatol: Software libre para la depuración y homogeneización de datos
  climatológicos. IV Congress of the Spanish Climatology Society (Santander,
  2-5 of November 2004).
}
\author{Jose A. Guijarro}
\note{
  As inhomogeneities in one series will affect other nearby series, making
  them suspects of inhomogeneity even if they are good, it is advisable to
  proceed step by step, beginning with a fairly high \code{wa} (>=1000) to
  only correct the most prominent errors, and successively repeat the
  process with decreasing values of \code{wa}.
}
\seealso{\code{\link{grafanom}}, \code{\link{grabeps}}, \code{\link{depstat}}}
\examples{
  #The two files PTOT_51-00.* of the package's data directory must be
  #  first copied to the working directory. Then run:
  \dontrun{depudm("PTOT",1951,2000,b=.1,ttip=2,sqrtrans=TRUE,graf=TRUE)}
}
\keyword{ts}
\keyword{misc}
\keyword{hplot}
