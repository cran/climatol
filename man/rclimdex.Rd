\name{rclimdex}
\alias{rclimdex}
\alias{rcxsta}
\alias{df1}
\alias{df2}
\alias{df3}
\docType{data}
\title{Daily temperature data}
\description{
  Daily temperatures (1981-1990) for three stations, in separate files.
}
\usage{data(rclimdex)}
\format{This data set contains:
\describe{
\item{rcxsta}{Data frame with file names, degrees, minutes and seconds of
latitude, degrees, minutes and seconds of longitude, and elevations (invented)}
\item{df1, df2, df3}{Data frames with years, months, days, precipitations,
maximum and minimum temperatures (invented)}
}
}
\source{Synthetic invented data.}
\seealso{\code{\link{rclimdex2climatol}}}
\examples{
data(rclimdex)
print(rcxsta)
head(df1)
#clean memory space:
rm(rcxsta, df1, df2, df3)
}
\keyword{datasets}

