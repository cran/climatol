\name{db2dat}
\alias{db2dat}
\title{Get daily or monthly data from a database and build input files *.dat and *.est}
\description{
  This function facilitates the creation of the input files needed by this package by retrieving the data from a database through an RODBC connection.
}
\usage{
db2dat(varcli, anyi, anyf, minny=5, daily=TRUE, ch, dformat='\%Y-\%m-\%d',
vtable, vcode, vdate, vval, stable, scode, sname, sx, sy, sz)
}
\arguments{
  \item{varcli}{Acronym of the name of the studied climatic variable,
    as it will appear in all data file names.}
  \item{anyi}{Initial year of the data to be included in the file.}
  \item{anyf}{Final year of the data to be included in the file.}
  \item{minny}{Minimum number of years with data for a series to be included in
  the file.}
  \item{daily}{Logical flag indicating whether the data are daily (the default) or monthly (set \code{daily=FALSE} in this case).}
  \item{ch}{Already open ODBC connection to the climatic database.}
  \item{dformat}{Date format in the database.}
  \item{vtable}{Name of the table containing our climatic variable.}
  \item{vcode}{Name of the variable containing station codes in the \code{vtable} table.}
  \item{vdate}{Name of the variable containing dates in the \code{vtable} table.}
  \item{vval}{Name of the climatic variable in the \code{vtable} table.}
  \item{stable}{Name of the table containing station information (metadata).}
  \item{scode}{Name of the variable containing station codes in the table \code{stable}.}
  \item{sname}{Name of the variable containing station names in the \code{stable} table.}
  \item{sx}{Name of the variable containing longitudes (degrees with decimals!) in the \code{stable} table.}
  \item{sy}{Name of the variable containing latitudes (degrees with decimals!) in the \code{stable} table.}
  \item{sz}{Name of the variable containing elevations (meters) in the \code{stable} table.}
}
\details{
This function creates the two input files needed by the homogenization functions of this package, \file{VAR_YEAR-YEAR.dat} (holding the data) and \file{VAR_YEAR-YEAR.est} (holding station coordinates, codes and names).

The table in the accessed database must contain either daily or monthly data (set \code{daily=FALSE} in this case). Otherwise the number of data per series will not be match the expected value and the function will fail.

Moreover, every data item must be in a different record in the database, as in this example table of monthly data (different variables for the same time step are O.K.):

\code{Station   Date     T.max T.min Rel.Hum Precip  Wind.speed}\cr
\code{S032    1991-01-01  12.1  -2.1    59    128.2     5.4}\cr
\code{S032    1991-02-01  13.2  -2.5    62     78.4     6.2}\cr
\code{...}

But if the table in the database arranges all monthly values of one year (or all daily values of one month) in a single record, then this function will not retreive the data correctly.
}
\value{
  This function does not return any value.
}
\seealso{\code{\link{homogen}}, \code{\link{homogsplit}}}
\examples{
\dontrun{
# First we must access our climatic database through RODBC, wich requires to
# have this package installed. System programs that allow ODBC connections to
# databases must also be installed and properly configured.
#
# For this example we will assume that our database is named "climate" and we
# access it with user "USER" and password "PASS".  Then we open the connection
# with:
#
library(RODBC)
ch <- odbcConnect("climate",uid="USER",pwd="PASS")
#
# Now we want to use this function to gather all monthly relative humidity
# averages for the period 1961-2015, requiring a minimum of 10 years of data
# (not necessarily consecutive). We must use the corresponding names of tables
# and headers existing the the database, and putting the parameters in the
# required order we avoid the need to name them:
#
db2dat('HRel',1961,2015,10,FALSE,ch,'%Y-%m-%d','monthly_relhum',
'Station','Date','Value','stations','Station','Name','Longitude',
'Latitude','Elevation')
#
odbcClose(ch) #close the connection if you do not need it anymore
#
# Our data would now be ready to be homogenized with the homogen function:
#
homogen('HRel',1961,2015,vmin=0,vmax=100)
#
}
}
\keyword{datagen}
