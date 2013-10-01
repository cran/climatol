\name{outrename}
\alias{outrename}
\title{Rename homogen's output files}
\description{
  This function inserts a suffix to the output file names of \code{homogen},
  to prevent them from being rewritten by any further function run.
}
\usage{
outrename(varcli, anyi, anyf, suffix, varsuf=TRUE, restore=FALSE)
}
\arguments{
  \item{varcli}{Acronym of the name of the studied climatic variable,
    as in the data file name.}
  \item{anyi}{Initial year of the study period}
  \item{anyf}{Final year of the study period}
  \item{suffix}{Suffix to be inserted in the output file names.}
  \item{varsuf}{Logical value indicating that the suffix will be appended to
    the variable name. Set it to \code{FALSE} if you prefere to append the
    suffix to the file base name.}
  \item{restore}{Set this parameter to \code{TRUE} to remove a suffix
    previously inserted by this function.}
}
\details{
  The variable (or file base) name is separated from the appended suffix by
  a hyfen. The purpose of this function is to facilitate the comparison of
  two or more applications of \code{homogen} to the same data with
  different parameters.
}
\value{
  This function does not return any value.
}
\seealso{\code{\link{homogen}}}
\examples{
#After having run the first example of the homogen function, you can try:
#
#   outrename("Tmin", 1956, 2005, "old")
#
#The previous four output files will display 'Tmin-old' as the new
#variable name.
}
\keyword{utilities}
