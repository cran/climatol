\name{rclimdex2climatol}
\alias{rclimdex2climatol}
\title{Convert RClimDex daily data files to Climatol input format}
\description{
  This function can be useful to prepare the Climatol input files when the user have their daily data in RClimDex format.
}
\usage{
rclimdex2climatol(stfile, kvar, varcli='', chrcod=c(6,10), anyi=NA, anyf=NA,
  mis=-99.9, mindat=365, names=FALSE)
}
\arguments{
  \item{stfile}{Name of the file with the list of data file names and station coordinates, codes and names (in HOMER format).}  
  \item{kvar}{RClimDex variable to extract: 1(RR), 2(TX), 3(TN).}
  \item{varcli}{Acronym of the name of the studied climatic variable.}
  \item{chrcod}{Initial and final characters of data file names to be used as station codes. (Between the sixth and the tenth by default.)}
  \item{anyi}{First year to study. (Defaults to the first year of available data.)}
  \item{anyf}{Last year to study. (Defaults to the last year of available data.)}
  \item{mis}{Missing data code. (Defaults to -99.9.)}
  \item{mindat}{Minimum required number of data per station. (Defaults to 365 daily data.)}
  \item{names}{Are there a 9th column with the station names? (\code{FALSE} by default.)}
}
\details{
Users of the RClimDex program can convert their daily data files to the Climatol format. All files listed in \code{stfile} will be read, and the selected variable (precipitation, maximum or minimum temperature) will be stored in a unique *.dat file, with its companion *.est station file. Therefore, if you want to convert all three variables, you must run this function three times.
}
\value{
  This function does not return any value.
}
\seealso{\code{\link{homogen}}}
\examples{
#Set a temporal working directory and write input files:
wd <- tempdir()
wd0 <- setwd(wd)
data(rclimdex)
write.table(rcxsta,'stations.txt',sep='\t',row.names=FALSE,col.names=FALSE,
  quote=FALSE)
write.table(df1,'raxxd000WY003d.txt',row.names=FALSE,col.names=FALSE)
write.table(df2,'raxxd000WY018d.txt',row.names=FALSE,col.names=FALSE)
write.table(df3,'raxxd000WY020d.txt',row.names=FALSE,col.names=FALSE)
#Now run the examples:
rclimdex2climatol('stations.txt',1,'RR',c(9,13))
rclimdex2climatol('stations.txt',2,'TX',c(9,13))
rclimdex2climatol('stations.txt',3,'TN',c(9,13))
#Return to user's working directory:
setwd(wd0)
#Input and output files can be found in directory:
print(wd)
}

