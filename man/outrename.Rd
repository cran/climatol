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
\dontrun{
#After having run the example of the homogen function, you can do:
outrename("Ptest", 1951, 2010, "old")
#
#The previous output files will be renamed to Ptest-old_1951-2010*
}
}
\keyword{utilities}
