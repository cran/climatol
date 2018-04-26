\name{climatol2rclimdex}
\alias{climatol2rclimdex}
\title{Convert DAILY data from Climatol to RClimDex input format}
\description{
  This function can be useful to apply the RClimDex program to daily data homogenized with Climatol.}
\usage{
climatol2rclimdex(varRR, varTX, varTN, yiRR, yfRR, yiTX=yiRR, yfTX=yfRR,
  yiTN=yiRR, yfTN=yfRR, prefix='hoclm', dir=NA, na='-99.9', nm=NA, dah=NA,
  nei=NA, est.c=NA)
}
\arguments{
  \item{varRR, varTX, varTN}{Name of the variables in the Climatol destination files. If some variable is not available, name it as ''.}  
  \item{yiRR, yfRR}{Initial and final years for the RR variable.}
  \item{yiTX, yfTX, yiTN, yfTN}{Initial and final years for the TX and TN variables. (By default they are the same for the three variables; otherwise, the output file will contain data for the common period only.)}
  \item{prefix}{Prefix to apply to the names of the RClimDex destination files. ('hoclm' by default.)}
  \item{dir}{Directory where to save the RClimDex destination files. (If not set, they will be saved into the R working directory.)}
  \item{na}{Missing data code to use in the ouput files. ('-99.9' by default.)}
  \item{nm, dah, nei, est.c}{Objects to be read from the Climatol
homogenization output. (The user does not need to care about them.)}
}
\details{
After homogenizing daily series with Climatol, the user may be interested in
applying the RClimDex program to the homogenized series. This function
automatizes the conversion of file formats between both programs.
}
\value{
  This function does not return any value.
}
\seealso{\code{\link{homogen}}}
\examples{
#Set a temporal working directory and generate input files:
wd <- tempdir()
wd0 <- setwd(wd)
data(tmax)
write.table(est.c,'Tmax_2001-2003.est',row.names=FALSE,col.names=FALSE)
write(dat,'Tmax_2001-2003.dat')
homogen('Tmax',2001,2003,snht1=30)
#Now run the example:
climatol2rclimdex('','Tmax','',2001,2003)
#Return to user's working directory:
setwd(wd0)
#Input and output files can be found in directory:
print(wd)
#(RR and TX columns are filled with missing data codes, because we do not have
#provided Climatol homogenization results for those variables.)
}

