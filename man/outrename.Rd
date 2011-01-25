\name{outrename}
\alias{outrename}
\title{Rename homogen's output files}
\description{
  This function appends a suffix to the output files of \code{homogen},
  to prevent them from being rewritten by another run of the function.
}
\usage{
  outrename(varcli, anyi, anyf, suffix)
}
\arguments{
  \item{varcli}{Acronym of the name of the studied climatic variable,
    as in the data file name.}
  \item{anyi}{Initial year of the study period}
  \item{anyf}{Final year of the study period}
  \item{suffix}{Suffix to be appended to the base names of the files.}
}
\details{
  The base name is separeted from the appended suffix by a hyfen.
  The purpose of this function is to facilitate the comparison of
  two or more applications of \code{homogen} to the same data with
  different parameters.
}
\value{
  This function does not return any value.
}
\seealso{\code{\link{homogen}}}
\examples{
  #After having run the first example of the homogen function, you can try:
  \dontrun{outrename("Tmin", 1956, 2005, "old")}
  #The previous four output files will have '-old' appended to their base name.
}
\keyword{utilities}
