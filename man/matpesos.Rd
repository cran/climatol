\name{matpesos}
\alias{matpesos}
\title{Inter-station weight matrix computation}
\description{
  Computation of the inter-station weight matrix that will be used by
  \code{datest} to estimate each station series from all other stations.
}
\usage{
  matpesos(wa=5, wz=.001, verb=TRUE)
}
\arguments{
  \item{wa}{Shape parameter of the weighting function \code{1/(1+d^2/wa)},
    where d stands for distance. Low values (1-10) increase the weight of
    nearby stations, while high values (>1000) take in account also the
    far stations influence. Do \code{wa=0} if you want an unweighted average of
    all the stations.}
  \item{wz}{Scale parameter of the vertical coordinate \code{Z}. The default
    value assumes that \code{X} and \code{Y} are expressed in km, while
    \code{Z} is expressed in m.  Can be used to change the \code{Z} weight in
    inter-station distance computations.}
  \item{verb}{If \code{TRUE}, progress indications will be shown in the
    terminal.}
}
\value{
  Object \code{est.w} is created, containing the computed weight matrix.
}
\seealso{\code{\link{depudm}}}
\keyword{internal}
