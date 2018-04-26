\name{outrename}
\alias{outrename}
\title{Rename homogen's output files}
\description{
  This function inserts a suffix to the output file names of \code{homogen},
  to prevent them from being rewritten by any further function run.
}
\usage{
outrename(varcli, anyi, anyf, suffix, restore=FALSE)
}
\arguments{
  \item{varcli}{Acronym of the name of the studied climatic variable,
    as in the data file name.}
  \item{anyi}{Initial year of the study period}
  \item{anyf}{Final year of the study period}
  \item{suffix}{Suffix to be inserted (or removed) in the output file names.}
  \item{restore}{Set this parameter to \code{TRUE} to remove the \code{suffix}
    previously inserted by this function.}
}
\details{
  The variable (or file base) name is separated from the appended suffix by a
  hyphen. The purpose of this function is to allow a new application of
  \code{homogen} to the same data with different parameters without overwriting
  the previous results.
}
\value{
  This function does not return any value.
}
\seealso{\code{\link{homogen}}}
\examples{
#Set a temporal working directory and write input files:
wd <- tempdir()
wd0 <- setwd(wd)
data(Ptest)
dim(dat) <- c(720,20)
dat[601:720,5] <- dat[601:720,5]*1.8
write(dat[481:720,1:5],'pcp_1991-2010.dat')
write.table(est.c[1:5,1:5],'pcp_1991-2010.est',row.names=FALSE,col.names=FALSE)
homogen('pcp',1991,2010,std=2) #homogenization
#Now run the example:
outrename('pcp',1991,2010,'ok')
#Return to user's working directory:
setwd(wd0)
#Input and output files can be found in directory:
print(wd)
}
\keyword{utilities}
