\name{xls2csv}
\alias{xls2csv}
\title{Join all data in *.xls or *.xlsx files into a single CSV file}
\description{
  This function reads all *.xls or *.xlsx files contained in a directory and dumps their data into a single CSV file.
}
\usage{
  xls2csv(tmpdir, archdir, var, datcols=1:4, codesep='-', dec='.', sep=',') 
}
\arguments{
\item{tmpdir}{temporal directory containing the files to read.}
\item{archdir}{directory where to archive files after processing.}
\item{var}{destination name of the variable.}
\item{datcols}{data columns to be written to the output file.}
\item{codesep}{character string separating the code from the rest of the file name (\code{'-'} by default).}
\item{dec}{character to use as decimal point in the output file (\code{'.'} by default).}
\item{sep}{character separating data in the output file (\code{','} by default).}
}
\details{
  File names must begin with their station code, which may optionally be followed by a hyphen ('-') or other code separator character (specified with the parameter \code{codesep}) and the name of the station or other characters.

  File contents must have one header line at the top. If they contain more,
supplementary header lines should have at least one empty cell in the columns
of date and data to be read.

  After their data have been dumped into the output \code{xls_*_data.csv} file, original files are moved to the \code{archdir} directory.

Note that data are appended to the output CSV files every time you run this function putting new files in the \code{tmpdir} directory.

Code and station names (if included in the file names) are appended to \code{xls_*_stations.csv}.

\code{climatol} input files can then be obtained from both output
\code{xls_*.csv} files with the \code{csv2climatol} function.
}
\seealso{\code{\link{csv2climatol}}, \code{\link{homogen}}}
\examples{
## Set a temporal working directory:
wd <- tempdir()
wd0 <- setwd(wd)

## Create origin and destination directories and copy example input files:
dir.create('dir1'); dir.create('dir2')
file.copy(exampleFiles('p064.xlsx'),'dir1')
file.copy(exampleFiles('p082.xlsx'),'dir1')
file.copy(exampleFiles('p084.xlsx'),'dir1')

## Now run the example:
xls2csv('dir1','dir2','TN',datcols=c(1:3,6))

## Return to user's working directory:
setwd(wd0)

## Input and output files can be found in the directory:
print(wd)
}
\keyword{manip}
