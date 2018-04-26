\name{datsubset}
\alias{datsubset}
\title{Subset Climatol input data into new files}
\description{
  This function allows saving a subset of Climatol input data into new input
files, by selecting a subperiod and/or series with a minimum number of years with data.
}
\usage{
datsubset(varcli, anyi, anyf, anyis=anyi, anyfs=anyf, minny=NA)
}
\arguments{
  \item{varcli}{Acronym of the name of the studied climatic variable.}
  \item{anyi}{Initial year of the data present in the file.}
  \item{anyf}{Final year of the data present in the file.}
  \item{anyis}{First year of the output subperiod. (Defaults to \code{anyi}).}
  \item{anyfs}{Last year of the output subperiod. (Defaults to \code{anyf}).}
  \item{minny}{Minimum number of years with data to retain the series.}
}
\details{
  Climatol homogenization requires that no time step be totally void of data in all stations simultaneously. This function allows subsetting already existing Climatol input files by selecting a subperiod and/or stations with a minimum number of years with data (may contain gaps).
}
\value{
  This function does not return any value.
}
\examples{
#Set a temporal working directory and write input files:
wd <- tempdir()
wd0 <- setwd(wd)
data(Ptest)
write(dat,'Ptest_1951-2010.dat')
write.table(est.c, 'Ptest_1951-2010.est', row.names=FALSE, col.names=FALSE)
rm(dat,est.c) #remove loaded data from memory space
datsubset('Ptest',1951,2010,1971,2000,20)
#Return to user's working directory:
setwd(wd0)
#Input and output files can be found in directory:
print(wd)
}
\keyword{datagen}
