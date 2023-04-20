\name{exampleFiles}
\alias{exampleFiles}
\title{Get the path to some example files}
\description{
  This function provides the path to files needed to run examples of some functions of the \code{climatol} package. 
}
\usage{
  exampleFiles(file=NULL)
}
\arguments{
  \item{file}{Name of the needed file. If NULL, all example files will be listed.}
}
\details{
  This function is an adaptation of \code{readxl_example}, of the \code{readxl} package.
}
\examples{
exampleFiles()
exampleFiles('Temp_1991-2000.rda')
}
\keyword{misc}
